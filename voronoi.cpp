#include "voronoi.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "quickhull/QuickHull.cpp"

namespace voronoi {
namespace {

using Clock = std::chrono::steady_clock;

struct VoronoiCell {
    std::vector<Vec3> vertices;
    quickhull::HalfEdgeMesh<FloatType, size_t> mesh;
};

static inline Vec3 cross_product(Vec3 a, Vec3 b) {
    return Vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

static inline FloatType dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline FloatType norm2(Vec3 a) {
    return dot(a, a);
}

static inline bool solve_dual_vertex(Vec3 u1, Vec3 u2, Vec3 u3, Vec3& out_v) {
    Vec3 cp = cross_product(u2, u3);
    FloatType det = dot(u1, cp);
    if (std::abs(det) < 1e-16) return false;
    Vec3 sum_cp = cross_product(u2, u3) + cross_product(u3, u1) + cross_product(u1, u2);
    out_v = sum_cp / det;
    return true;
}

struct Grid {
    std::vector<int> bin_start;
    std::vector<int> ix, iy, iz;
    std::vector<int> orig_index;
    int cells_per_axis = 1;
    FloatType cell_size = 1.0;
    FloatType inv_cell_size = 1.0;

    inline int get_bin_index(FloatType v) const {
        int idx = static_cast<int>(v * inv_cell_size);
        if (idx < 0) idx = 0;
        if (idx >= cells_per_axis) idx = cells_per_axis - 1;
        return idx;
    }
};

static Grid build_grid_and_reorder(const std::vector<Vec3>& seeds_in, int target_per_cell, std::vector<Vec3>& seeds_out) {
    const int n = static_cast<int>(seeds_in.size());
    const FloatType ratio = static_cast<FloatType>(n) / static_cast<FloatType>(target_per_cell);
    const int cells_per_axis = std::max(1, static_cast<int>(std::cbrt(ratio)));

    Grid grid;
    grid.cells_per_axis = cells_per_axis;
    grid.inv_cell_size = static_cast<FloatType>(cells_per_axis);
    grid.cell_size = 1.0 / grid.inv_cell_size;

    const int num_bins = cells_per_axis * cells_per_axis * cells_per_axis;
    grid.bin_start.assign(num_bins + 1, 0);
    grid.ix.resize(n);
    grid.iy.resize(n);
    grid.iz.resize(n);
    grid.orig_index.resize(n);

    std::vector<int> particle_bin(n);
    for (int i = 0; i < n; ++i) {
        const int ix = grid.get_bin_index(seeds_in[i].x);
        const int iy = grid.get_bin_index(seeds_in[i].y);
        const int iz = grid.get_bin_index(seeds_in[i].z);
        const int bin_idx = iz * cells_per_axis * cells_per_axis + iy * cells_per_axis + ix;
        particle_bin[i] = bin_idx;
        grid.bin_start[bin_idx]++;
    }

    int accum = 0;
    for (int i = 0; i < num_bins; ++i) {
        const int count = grid.bin_start[i];
        grid.bin_start[i] = accum;
        accum += count;
    }
    grid.bin_start[num_bins] = accum;

    seeds_out.resize(static_cast<size_t>(n));
    std::vector<int> new_ix(n), new_iy(n), new_iz(n);
    std::vector<int> current_pos = grid.bin_start;

    for (int old_i = 0; old_i < n; ++old_i) {
        const int bin = particle_bin[old_i];
        const int dest = current_pos[bin]++;
        seeds_out[static_cast<size_t>(dest)] = seeds_in[static_cast<size_t>(old_i)];
        grid.orig_index[static_cast<size_t>(dest)] = old_i;

        const int ix = grid.get_bin_index(seeds_out[static_cast<size_t>(dest)].x);
        const int iy = grid.get_bin_index(seeds_out[static_cast<size_t>(dest)].y);
        const int iz = grid.get_bin_index(seeds_out[static_cast<size_t>(dest)].z);
        new_ix[dest] = ix;
        new_iy[dest] = iy;
        new_iz[dest] = iz;
    }

    grid.ix = std::move(new_ix);
    grid.iy = std::move(new_iy);
    grid.iz = std::move(new_iz);
    return grid;
}

static void get_grid_neighbors(
    const Grid& grid,
    const std::vector<Vec3>& seeds,
    int seed_idx,
    FloatType min_dist,
    FloatType max_dist,
    std::vector<int>& out_neighbors) {
    const int cx = grid.ix[seed_idx];
    const int cy = grid.iy[seed_idx];
    const int cz = grid.iz[seed_idx];
    const int cpa = grid.cells_per_axis;

    const int search_r = static_cast<int>(std::ceil(max_dist * grid.inv_cell_size));
    const int start_z = std::max(0, cz - search_r);
    const int end_z = std::min(cpa - 1, cz + search_r);
    const int start_y = std::max(0, cy - search_r);
    const int end_y = std::min(cpa - 1, cy + search_r);
    const int start_x = std::max(0, cx - search_r);
    const int end_x = std::min(cpa - 1, cx + search_r);

    int exclude_r = -1;
    if (min_dist > 0) {
        exclude_r = static_cast<int>((min_dist * grid.inv_cell_size) * 0.57) - 1;
    }

    const int in_start_z = cz - exclude_r;
    const int in_end_z = cz + exclude_r;
    const int in_start_y = cy - exclude_r;
    const int in_end_y = cy + exclude_r;
    const int in_start_x = cx - exclude_r;
    const int in_end_x = cx + exclude_r;

    const FloatType min_dist_sq = min_dist * min_dist;
    const FloatType max_dist_sq = max_dist * max_dist;

    for (int z = start_z; z <= end_z; ++z) {
        const bool z_in = (exclude_r >= 0) && (z >= in_start_z && z <= in_end_z);
        for (int y = start_y; y <= end_y; ++y) {
            const bool y_in = z_in && (y >= in_start_y && y <= in_end_y);
            for (int x = start_x; x <= end_x; ++x) {
                if (y_in && (x >= in_start_x && x <= in_end_x)) continue;

                const int bin_idx = z * cpa * cpa + y * cpa + x;
                const int start = grid.bin_start[bin_idx];
                const int end = grid.bin_start[bin_idx + 1];

                for (int idx = start; idx < end; ++idx) {
                    const FloatType d2 = norm2(seeds[static_cast<size_t>(idx)] - seeds[static_cast<size_t>(seed_idx)]);
                    if (d2 >= min_dist_sq && d2 <= max_dist_sq) {
                        out_neighbors.push_back(idx);
                    }
                }
            }
        }
    }
}

static VoronoiCell compute_cell_vertices(
    Vec3 seed,
    const std::vector<int>& candidate_indices,
    const std::vector<Vec3>& all_seeds,
    std::vector<Vec3>& dual_points,
    quickhull::QuickHull<FloatType>& qh,
    VoronoiStats& stats) {
    auto t0 = Clock::now();
    constexpr FloatType eps = 1e-16;
    if (seed.x > eps) dual_points.push_back(Vec3(-1.0 / seed.x, 0, 0));
    if (seed.x < 1.0 - eps) dual_points.push_back(Vec3(1.0 / (1.0 - seed.x), 0, 0));
    if (seed.y > eps) dual_points.push_back(Vec3(0, -1.0 / seed.y, 0));
    if (seed.y < 1.0 - eps) dual_points.push_back(Vec3(0, 1.0 / (1.0 - seed.y), 0));
    if (seed.z > eps) dual_points.push_back(Vec3(0, 0, -1.0 / seed.z));
    if (seed.z < 1.0 - eps) dual_points.push_back(Vec3(0, 0, 1.0 / (1.0 - seed.z)));

    for (int idx : candidate_indices) {
        Vec3 r = all_seeds[static_cast<size_t>(idx)] - seed;
        const FloatType len_sq = norm2(r);
        if (len_sq > 1e-16) {
            dual_points.push_back(r * (2.0 / len_sq));
        }
    }

    auto t1 = Clock::now();
    stats.time_ms_dual_setup += std::chrono::duration<FloatType, std::milli>(t1 - t0).count();

    auto t_qh1 = Clock::now();
    stats.total_quickhull_calls++;
    stats.total_dual_points += static_cast<uint64_t>(dual_points.size());
    auto mesh = qh.getConvexHullAsMesh(
        reinterpret_cast<const FloatType*>(dual_points.data()),
        dual_points.size(),
        true,
        quickhull::defaultEps<FloatType>());
    auto t_qh2 = Clock::now();
    const FloatType qh_as_mesh_ms = std::chrono::duration<FloatType, std::milli>(t_qh2 - t_qh1).count();
    stats.time_ms_quickhull += qh_as_mesh_ms;

    const auto& diag = qh.getDiagnostics();
    const FloatType qh_internal_total_ms = static_cast<FloatType>(diag.m_time_total_micro / 1000.0);
    stats.time_ms_qh_total_internal += qh_internal_total_ms;
    stats.time_ms_qh_setup_internal += diag.m_time_setup_micro / 1000.0;
    stats.time_ms_qh_process_internal += diag.m_time_process_micro / 1000.0;
    stats.time_ms_qh_mesh_convert += std::max<FloatType>(0.0, qh_as_mesh_ms - qh_internal_total_ms);

    auto t2 = Clock::now();

    VoronoiCell cell;
    cell.mesh = std::move(mesh);

    if (cell.mesh.m_faces.empty()) {
        auto t3 = Clock::now();
        stats.time_ms_primal_recovery += std::chrono::duration<FloatType, std::milli>(t3 - t2).count();
        return cell;
    }

    cell.vertices.resize(cell.mesh.m_faces.size());
    for (size_t i = 0; i < cell.mesh.m_faces.size(); ++i) {
        const auto& face = cell.mesh.m_faces[i];
        const size_t he_idx = face.m_halfEdgeIndex;
        const auto& he1 = cell.mesh.m_halfEdges[he_idx];
        const auto& he2 = cell.mesh.m_halfEdges[he1.m_next];
        const auto& he3 = cell.mesh.m_halfEdges[he2.m_next];

        const Vec3 u1 = cell.mesh.m_vertices[he3.m_endVertex];
        const Vec3 u2 = cell.mesh.m_vertices[he1.m_endVertex];
        const Vec3 u3 = cell.mesh.m_vertices[he2.m_endVertex];

        Vec3 v;
        if (solve_dual_vertex(u1, u2, u3, v)) {
            cell.vertices[i] = v;
        } else {
            cell.vertices[i] = (u1 + u2 + u3) * (1.0 / 3.0);
        }
    }

    auto t3 = Clock::now();
    stats.time_ms_primal_recovery += std::chrono::duration<FloatType, std::milli>(t3 - t2).count();
    return cell;
}

static VoronoiCell build_cell(
    int i,
    const std::vector<Vec3>& seeds,
    const Grid& grid,
    quickhull::QuickHull<FloatType>& qh,
    VoronoiStats& stats,
    std::vector<int>& neighbor_indices,
    std::vector<Vec3>& dual_points) {
    FloatType current_dist = grid.cell_size * 1.5;
    FloatType prev_dist = 0.0;
    int step = 0;

    dual_points.clear();

    while (true) {
        auto t_start = Clock::now();
        neighbor_indices.clear();
        get_grid_neighbors(grid, seeds, i, prev_dist, current_dist, neighbor_indices);
        auto t_neigh = Clock::now();
        stats.time_ms_neighbors += std::chrono::duration<FloatType, std::milli>(t_neigh - t_start).count();

        if (step >= static_cast<int>(stats.neighbors_count_per_step.size())) {
            stats.neighbors_count_per_step.resize(static_cast<size_t>(step + 1), 0);
            stats.step_count.resize(static_cast<size_t>(step + 1), 0);
        }
        stats.neighbors_count_per_step[static_cast<size_t>(step)] += neighbor_indices.size();
        stats.step_count[static_cast<size_t>(step)]++;

        auto t_copy = Clock::now();
        stats.time_ms_overhead_other += std::chrono::duration<FloatType, std::milli>(t_copy - t_neigh).count();

        stats.neighbors_total += neighbor_indices.size();
        VoronoiCell cell = compute_cell_vertices(seeds[static_cast<size_t>(i)], neighbor_indices, seeds, dual_points, qh, stats);

        FloatType max_dist_sq = 0.0;
        for (const auto& v : cell.vertices) {
            max_dist_sq = std::max(max_dist_sq, norm2(v));
        }
        const FloatType max_dist = std::sqrt(max_dist_sq);
        const FloatType safe_dist = current_dist * 0.5;

        if (max_dist < safe_dist || current_dist > grid.cell_size * grid.cells_per_axis * 1.5) {
            return cell;
        }

        prev_dist = current_dist;
        const FloatType target_radius = max_dist * 2;
        const FloatType min_step = grid.cell_size * 0.1;
        const FloatType next_radius = std::max(target_radius, current_dist + min_step);
        current_dist = next_radius * 1.1;

        dual_points = std::move(cell.mesh.m_vertices);
        stats.total_expansions++;
        step++;
    }
}

struct PolyhedraContext {
    std::vector<size_t> v_in_he;
};

static void convert_to_polyhedron_inplace(const VoronoiCell& cell, Vec3 seed, PolyhedraContext& ctx, Polyhedron& poly) {
    poly.vertices.assign(cell.vertices.begin(), cell.vertices.end());
    for (auto& v : poly.vertices) {
        v = v + seed;
    }
    poly.face_indices.clear();
    poly.face_degree.clear();
    poly.edges.clear();

    const auto& mesh_verts = cell.mesh.m_vertices;
    const auto& mesh_he = cell.mesh.m_halfEdges;
    constexpr size_t invalid = std::numeric_limits<size_t>::max();

    if (ctx.v_in_he.size() < mesh_verts.size()) {
        ctx.v_in_he.resize(mesh_verts.size());
    }
    std::fill(ctx.v_in_he.begin(), ctx.v_in_he.begin() + mesh_verts.size(), invalid);

    for (size_t i = 0; i < mesh_he.size(); ++i) {
        size_t v = mesh_he[i].m_endVertex;
        ctx.v_in_he[v] = i;
    }

    const size_t num_faces = mesh_verts.size();
    poly.face_degree.reserve(num_faces);
    poly.face_indices.reserve(mesh_he.size());
    poly.edges.reserve(mesh_he.size() / 2);

    for (size_t i = 0; i < num_faces; ++i) {
        size_t start_he = ctx.v_in_he[i];
        if (start_he == invalid) continue;

        uint8_t degree = 0;
        size_t curr = start_he;
        do {
            size_t f_idx = mesh_he[curr].m_face;
            poly.face_indices.push_back(static_cast<int>(f_idx));
            degree++;

            const size_t next_he = mesh_he[curr].m_next;
            curr = mesh_he[next_he].m_opp;
        } while (curr != start_he && curr != invalid);

        if (degree > 0) {
            poly.face_degree.push_back(degree);
        }
    }

    for (size_t i = 0; i < mesh_he.size(); ++i) {
        const size_t opp = mesh_he[i].m_opp;
        if (i < opp) {
            const size_t f1 = mesh_he[i].m_face;
            const size_t f2 = mesh_he[opp].m_face;
            poly.edges.push_back({static_cast<int>(f1), static_cast<int>(f2)});
        }
    }
}

}  // namespace

std::vector<Vec3> generate_random_points_box(size_t n, uint64_t seed, FloatType box_len) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<FloatType> dist(0.0, box_len);
    std::vector<Vec3> pts;
    pts.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        pts.push_back(Vec3(dist(rng), dist(rng), dist(rng)));
    }
    return pts;
}

void for_each_polyhedron(
    const std::vector<Vec3>& seeds_in,
    int target_per_cell,
    VoronoiStats& stats,
    const std::function<void(size_t seed_index, const Polyhedron& poly)>& on_polyhedron) {
    const auto t_total_start = Clock::now();
    std::cout << "Running for_each_polyhedron with " << seeds_in.size() << " seeds and target_per_cell = " << target_per_cell << "\n";

    std::vector<Vec3> seeds;
    const auto t_build_start = Clock::now();
    Grid grid = build_grid_and_reorder(seeds_in, target_per_cell, seeds);
    const auto t_build_end = Clock::now();
    stats.time_ms_build = std::chrono::duration<FloatType, std::milli>(t_build_end - t_build_start).count();

    quickhull::QuickHull<FloatType> qh;
    std::vector<int> neighbor_indices;
    neighbor_indices.reserve(512);
    std::vector<Vec3> dual_points;
    dual_points.reserve(128);
    PolyhedraContext poly_ctx;
    Polyhedron poly_tmp;

    const int n = static_cast<int>(seeds.size());
    for (int i = 0; i < n; ++i) {
        const auto t_cells_start = Clock::now();
        VoronoiCell cell = build_cell(i, seeds, grid, qh, stats, neighbor_indices, dual_points);
        const auto t_cells_end = Clock::now();
        stats.time_ms_cells += std::chrono::duration<FloatType, std::milli>(t_cells_end - t_cells_start).count();

        const auto t_poly_start = Clock::now();
        convert_to_polyhedron_inplace(cell, seeds[static_cast<size_t>(i)], poly_ctx, poly_tmp);
        const auto t_poly_end = Clock::now();
        stats.time_ms_polyhedra += std::chrono::duration<FloatType, std::milli>(t_poly_end - t_poly_start).count();

        const size_t original_index = static_cast<size_t>(grid.orig_index[static_cast<size_t>(i)]);
        on_polyhedron(original_index, poly_tmp);
    }

    const auto t_total_end = Clock::now();
    stats.total_ms = std::chrono::duration<FloatType, std::milli>(t_total_end - t_total_start).count();
}

FloatType compute_polyhedron_volume(const Polyhedron& poly) {
    if (poly.vertices.size() < 4 || poly.face_degree.empty()) return 0.0;

    Vec3 center(0, 0, 0);
    for (const auto& v : poly.vertices) {
        center = center + v;
    }
    center = center * (1.0 / static_cast<FloatType>(poly.vertices.size()));

    long double accum = 0.0L;
    size_t offset = 0;
    for (uint8_t deg_u8 : poly.face_degree) {
        const int deg = static_cast<int>(deg_u8);
        if (deg < 3) {
            offset += static_cast<size_t>(deg);
            continue;
        }

        Vec3 face_center(0, 0, 0);
        for (int k = 0; k < deg; ++k) {
            int idx = poly.face_indices[offset + static_cast<size_t>(k)];
            face_center = face_center + poly.vertices[static_cast<size_t>(idx)];
        }
        face_center = face_center * (1.0 / static_cast<FloatType>(deg));

        const int i0 = poly.face_indices[offset];
        const int i1 = poly.face_indices[offset + 1];
        const int i2 = poly.face_indices[offset + 2];
        const Vec3 v0 = poly.vertices[static_cast<size_t>(i0)];
        const Vec3 v1 = poly.vertices[static_cast<size_t>(i1)];
        const Vec3 v2 = poly.vertices[static_cast<size_t>(i2)];
        const Vec3 normal = cross_product(v1 - v0, v2 - v0);

        const bool reverse = dot(normal, face_center - center) < 0.0;
        for (int k = 1; k < deg - 1; ++k) {
            const int ia = poly.face_indices[offset + static_cast<size_t>(k)];
            const int ib = poly.face_indices[offset + static_cast<size_t>(k + 1)];
            const Vec3 p0 = v0;
            const Vec3 p1 = poly.vertices[static_cast<size_t>(reverse ? ib : ia)];
            const Vec3 p2 = poly.vertices[static_cast<size_t>(reverse ? ia : ib)];
            accum += static_cast<long double>(dot(p0, cross_product(p1, p2)));
        }

        offset += static_cast<size_t>(deg);
    }

    const long double vol = std::abs(accum) / 6.0L;
    return static_cast<FloatType>(vol);
}

void compute_polyhedron_face_areas(const Polyhedron& poly, std::vector<double>& out_areas) {
    out_areas.clear();
    out_areas.reserve(poly.face_degree.size());
    size_t offset = 0;
    for (uint8_t deg_u8 : poly.face_degree) {
        const int deg = static_cast<int>(deg_u8);
        if (deg < 3) {
            out_areas.push_back(0.0);
            offset += static_cast<size_t>(deg);
            continue;
        }

        const int i0 = poly.face_indices[offset];
        const Vec3 v0 = poly.vertices[static_cast<size_t>(i0)];
        
        double face_area_accum = 0.0;
        
        for (int k = 1; k < deg - 1; ++k) {
            const int i1 = poly.face_indices[offset + static_cast<size_t>(k)];
            const int i2 = poly.face_indices[offset + static_cast<size_t>(k + 1)];
            const Vec3 v1 = poly.vertices[static_cast<size_t>(i1)];
            const Vec3 v2 = poly.vertices[static_cast<size_t>(i2)];
            
            Vec3 cross = cross_product(v1 - v0, v2 - v0);
            double tri_area = std::sqrt(norm2(cross)) * 0.5;
            face_area_accum += tri_area;
        }
        out_areas.push_back(face_area_accum);
        offset += static_cast<size_t>(deg);
    }
}

TopologyOutput write_voro_compatible_output(std::ostream& out, size_t seed_index, const Vec3& seed, const Polyhedron& poly, double vol, std::vector<double>& area_buffer) {
    TopologyOutput topo;
    topo.index = static_cast<int>(seed_index);
    topo.volume = vol;
    topo.num_faces = static_cast<int>(poly.face_degree.size());

    compute_polyhedron_face_areas(poly, area_buffer);
    topo.face_areas = area_buffer;

    topo.face_edge_counts.clear();
    topo.face_edge_counts.reserve(poly.face_degree.size());
    for (uint8_t d : poly.face_degree) {
        topo.face_edge_counts.push_back(static_cast<int>(d));
    }

    out << std::scientific << std::setprecision(5);
    out << topo.index << " " << topo.volume << " " << topo.num_faces;

    for (double a : topo.face_areas) {
        out << " " << a;
    }

    for (int c : topo.face_edge_counts) {
        out << " " << c;
    }
    out << "\n";

    return topo;
}

void write_polyhedra_vtk(const std::string& path, const std::vector<Polyhedron>& polys, const std::vector<int>& cell_id_per_poly) {
    size_t total_points = 0;
    size_t total_polygons = 0;
    size_t polygons_index_count = 0;
    for (const auto& p : polys) {
        total_points += p.vertices.size();
        total_polygons += p.face_degree.size();
        for (uint8_t deg_u8 : p.face_degree) {
            polygons_index_count += 1 + static_cast<size_t>(deg_u8);
        }
    }

    std::ofstream out(path);
    out << "# vtk DataFile Version 3.0\n";
    out << "Voronoi polyhedra\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";
    out << "POINTS " << total_points << " double\n";
    out << std::setprecision(17);

    for (const auto& p : polys) {
        for (const auto& v : p.vertices) {
            out << v.x << " " << v.y << " " << v.z << "\n";
        }
    }

    out << "POLYGONS " << total_polygons << " " << polygons_index_count << "\n";

    std::vector<int> cell_id;
    cell_id.reserve(total_polygons);

    size_t base = 0;
    for (size_t pi = 0; pi < polys.size(); ++pi) {
        const auto& p = polys[pi];
        size_t off = 0;
        for (uint8_t deg_u8 : p.face_degree) {
            const int deg = static_cast<int>(deg_u8);
            out << deg;
            for (int k = 0; k < deg; ++k) {
                const int local_idx = p.face_indices[off + static_cast<size_t>(k)];
                out << " " << (base + static_cast<size_t>(local_idx));
            }
            out << "\n";
            off += static_cast<size_t>(deg);
            cell_id.push_back(cell_id_per_poly[pi]);
        }
        base += p.vertices.size();
    }

    out << "CELL_DATA " << total_polygons << "\n";
    out << "SCALARS cell_id int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int v : cell_id) {
        out << v << "\n";
    }
}

}  // namespace voronoi
