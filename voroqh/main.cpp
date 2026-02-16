#include "voronoi.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using FloatType = voronoi::FloatType;
using Vec3 = voronoi::Vec3;

struct SimulationConfig {
    int n;
    uint64_t seed_val;
    int target_per_cell;
    std::string out_path;
    int volume_flag;
    int voro_flag;
};

static bool ends_with(const std::string& s, const std::string& suffix) {
    if (s.size() < suffix.size()) return false;
    return std::equal(suffix.rbegin(), suffix.rend(), s.rbegin());
}

static bool run_voro_benchmark(const std::vector<Vec3>& seeds, FloatType& out_ms) {
    int available = std::system("command -v voro++ > /dev/null 2>&1");
    if (available != 0) return false;

    std::string path = "voro_input.txt";
    {
        std::ofstream out(path);
        out << std::setprecision(17);
        for (size_t i = 0; i < seeds.size(); ++i) {
            out << i << " " << seeds[i].x << " " << seeds[i].y << " " << seeds[i].z << "\n";
        }
    }

    auto t0 = std::chrono::steady_clock::now();
    std::string cmd = "voro++ -c \"%i %q %v %f %v\" 0 1 0 1 0 1 " + path + " > /dev/null 2>&1";
    int code = std::system(cmd.c_str());
    auto t1 = std::chrono::steady_clock::now();

    std::remove(path.c_str());
    std::remove((path + ".vol").c_str());

    if (code != 0) return false;
    out_ms = std::chrono::duration<FloatType, std::milli>(t1 - t0).count();
    return true;
}

static void run_simulation(const SimulationConfig& config) {
    const FloatType sys_length = 1.0;
    std::vector<Vec3> seeds = voronoi::generate_random_points_box(static_cast<size_t>(config.n), config.seed_val, sys_length);

    voronoi::VoronoiStats stats;
    double total_volume = 0.0;
    double time_ms_volume = 0.0;
    uint64_t faces_total = 0;
    uint64_t faces_max = 0;

    const bool write_vtk = (!config.out_path.empty() && ends_with(config.out_path, ".vtk"));
    std::vector<voronoi::Polyhedron> polys_out;
    std::vector<int> cell_ids;
    if (write_vtk) {
        polys_out.reserve(static_cast<size_t>(config.n));
        cell_ids.reserve(static_cast<size_t>(config.n));
    }

    voronoi::for_each_polyhedron(
        seeds,
        config.target_per_cell,
        stats,
        [&](size_t seed_index, const voronoi::Polyhedron& poly) {
            const uint64_t faces = static_cast<uint64_t>(poly.face_degree.size());
            faces_total += faces;
            faces_max = std::max(faces_max, faces);

            if (config.volume_flag) {
                auto tv0 = std::chrono::steady_clock::now();
                total_volume += voronoi::compute_polyhedron_volume(poly);
                auto tv1 = std::chrono::steady_clock::now();
                time_ms_volume += std::chrono::duration<double, std::milli>(tv1 - tv0).count();
            }

            if (write_vtk) {
                polys_out.push_back(poly);
                cell_ids.push_back(static_cast<int>(seed_index));
            }
        });

    if (write_vtk) {
        voronoi::write_polyhedra_vtk(config.out_path, polys_out, cell_ids);
    }

    FloatType voro_ms = 0.0;
    bool voro_ok = false;
    if (config.voro_flag != 0) {
        voro_ok = run_voro_benchmark(seeds, voro_ms);
    }

    std::cout << std::setprecision(15);
    std::cout << "n=" << config.n << " seed=" << config.seed_val << " target_per_cell=" << config.target_per_cell << "\n";
    if (config.volume_flag) {
        std::cout << "sum_volume=" << total_volume << " error=" << std::abs(total_volume - sys_length * sys_length * sys_length) << "\n";
    }

    std::cout << "neighbors_total=" << stats.neighbors_total << " neighbors_avg=" << (static_cast<double>(stats.neighbors_total) / config.n) << "\n";
    std::cout << "faces_total=" << faces_total << " faces_mean=" << (static_cast<double>(faces_total) / config.n) << " faces_max=" << faces_max << "\n";
    std::cout << "expansions_total=" << stats.total_expansions << " expansions_avg=" << (static_cast<double>(stats.total_expansions) / config.n) << "\n";
    std::cout << "dual_points_total=" << stats.total_dual_points << " dual_points_avg=" << (static_cast<double>(stats.total_dual_points) / stats.total_quickhull_calls) << "\n";
    for (size_t i = 0; i < stats.step_count.size(); ++i) {
        if (stats.step_count[i] > 0) {
            std::cout << "step " << i << ": avg_neighbors=" << (static_cast<double>(stats.neighbors_count_per_step[i]) / stats.step_count[i])
                      << " (count=" << stats.step_count[i] << ")\n";
        }
    }
    std::cout << "quickhull_calls_total=" << stats.total_quickhull_calls << " quickhull_calls_avg=" << (static_cast<double>(stats.total_quickhull_calls) / config.n) << "\n";
    std::cout << "time_ms_total=" << stats.total_ms << " time_ms_build=" << stats.time_ms_build << " time_ms_cells=" << stats.time_ms_cells
              << " time_us_per_cell=" << (stats.time_ms_cells * 1000.0 / config.n) << "\n";
    std::cout << "time_ms_quickhull=" << stats.time_ms_quickhull << " (" << (stats.time_ms_quickhull / stats.time_ms_cells * 100.0) << "% of cells time)\n";
    std::cout << "  [Detailed Breakdown]\n";
    std::cout << "  neighbors_search=" << stats.time_ms_neighbors << " ms (" << (stats.time_ms_neighbors / stats.time_ms_cells * 100.0) << "%)\n";
    std::cout << "  dual_setup=" << stats.time_ms_dual_setup << " ms (" << (stats.time_ms_dual_setup / stats.time_ms_cells * 100.0) << "%)\n";
    std::cout << "  primal_recovery=" << stats.time_ms_primal_recovery << " ms (" << (stats.time_ms_primal_recovery / stats.time_ms_cells * 100.0) << "%)\n";
    std::cout << "  overhead_other=" << stats.time_ms_overhead_other << " ms (" << (stats.time_ms_overhead_other / stats.time_ms_cells * 100.0) << "%)\n";
    std::cout << "  [QuickHull Internal]\n";
    std::cout << "  qh_total=" << stats.time_ms_qh_total_internal << " ms\n";
    std::cout << "  qh_setup=" << stats.time_ms_qh_setup_internal << " ms\n";
    std::cout << "  qh_process=" << stats.time_ms_qh_process_internal << " ms\n";
    std::cout << "  qh_mesh_convert=" << stats.time_ms_qh_mesh_convert << " ms\n";
    std::cout << "  [Polyhedra Conversion]\n";
    std::cout << "  total_time=" << stats.time_ms_polyhedra << " ms (" << (stats.time_ms_polyhedra / stats.time_ms_cells * 100.0) << "% of cells time)\n";
    if (config.volume_flag) {
        std::cout << "  [Volume]\n";
        std::cout << "  total_time=" << time_ms_volume << " ms\n";
    }

    if (config.voro_flag == 0) {
        std::cout << "voro_ms=off\n";
    } else if (voro_ok) {
        std::cout << "voro_ms=" << voro_ms << "\n";
    } else {
        std::cout << "voro_ms=na\n";
    }
}

int main(int argc, char** argv) {
    SimulationConfig config;
    config.n = 200;
    config.seed_val = 1;
    config.target_per_cell = 4;
    config.out_path = "polyhedra.vtk";
    config.volume_flag = 0;
    config.voro_flag = 1;

    if (argc > 1) config.n = std::max(1, std::atoi(argv[1]));
    if (argc > 2) config.seed_val = static_cast<uint64_t>(std::strtoull(argv[2], nullptr, 10));
    if (argc > 3) config.target_per_cell = std::max(1, std::atoi(argv[3]));
    if (argc > 4) config.out_path = argv[4];
    if (argc > 5) config.volume_flag = std::atoi(argv[5]) != 0 ? 1 : 0;
    if (argc > 6) config.voro_flag = std::atoi(argv[6]) != 0 ? 1 : 0;

    run_simulation(config);
    return 0;
}

#if 0

struct Grid {
    std::vector<int> bin_start;
    std::vector<int> ix, iy, iz;
    int cells_per_axis;
    FloatType cell_size;
    FloatType inv_cell_size;

    inline int get_bin_index(FloatType v) const {
        int idx = static_cast<int>(v * inv_cell_size);
        if (idx < 0) idx = 0;
        if (idx >= cells_per_axis) idx = cells_per_axis - 1;
        return idx;
    }
};

static Grid build_grid_and_reorder(std::vector<Vec3>& seeds, int target_per_cell) {
    int n = static_cast<int>(seeds.size());
    FloatType ratio = static_cast<FloatType>(n) / target_per_cell;
    int cells_per_axis = std::max(1, static_cast<int>(std::cbrt(ratio)));
    
    Grid grid;
    grid.cells_per_axis = cells_per_axis;
    grid.inv_cell_size = static_cast<FloatType>(cells_per_axis);
    grid.cell_size = 1.0 / grid.inv_cell_size;
    
    int num_bins = cells_per_axis * cells_per_axis * cells_per_axis;
    grid.bin_start.assign(num_bins + 1, 0);
    grid.ix.resize(n);
    grid.iy.resize(n);
    grid.iz.resize(n);

    std::vector<int> particle_bin(n);

    for (int i = 0; i < n; ++i) {
        int ix = grid.get_bin_index(seeds[i].x);
        int iy = grid.get_bin_index(seeds[i].y);
        int iz = grid.get_bin_index(seeds[i].z);
        int bin_idx = iz * cells_per_axis * cells_per_axis + iy * cells_per_axis + ix;
        particle_bin[i] = bin_idx;
        grid.bin_start[bin_idx]++;
    }

    int accum = 0;
    for (int i = 0; i < num_bins; ++i) {
        int count = grid.bin_start[i];
        grid.bin_start[i] = accum;
        accum += count;
    }
    grid.bin_start[num_bins] = accum;

    std::vector<Vec3> new_seeds(n);
    std::vector<int> new_ix(n), new_iy(n), new_iz(n);
    std::vector<int> current_pos = grid.bin_start;

    for (int i = 0; i < n; ++i) {
        int bin = particle_bin[i];
        int dest = current_pos[bin]++;
        new_seeds[dest] = seeds[i];
        
        int ix = grid.get_bin_index(new_seeds[dest].x);
        int iy = grid.get_bin_index(new_seeds[dest].y);
        int iz = grid.get_bin_index(new_seeds[dest].z);
        new_ix[dest] = ix;
        new_iy[dest] = iy;
        new_iz[dest] = iz;
    }
    
    seeds = std::move(new_seeds);
    grid.ix = std::move(new_ix);
    grid.iy = std::move(new_iy);
    grid.iz = std::move(new_iz);
    
    return grid;
}

static void get_grid_neighbors(const Grid& grid, const std::vector<Vec3>& seeds, int seed_idx, FloatType min_dist, FloatType max_dist, std::vector<int>& out_neighbors) {
    int cx = grid.ix[seed_idx];
    int cy = grid.iy[seed_idx];
    int cz = grid.iz[seed_idx];
    int cpa = grid.cells_per_axis;
    
    int search_r = static_cast<int>(std::ceil(max_dist * grid.inv_cell_size));
    int start_z = std::max(0, cz - search_r);
    int end_z = std::min(cpa - 1, cz + search_r);
    int start_y = std::max(0, cy - search_r);
    int end_y = std::min(cpa - 1, cy + search_r);
    int start_x = std::max(0, cx - search_r);
    int end_x = std::min(cpa - 1, cx + search_r);

    // Safe exclusion radius: max distance in box[-r, r] must be < min_dist.
    // Max dist is approx sqrt(3)*(r+1)*cell_size.
    // So r < (min_dist/cell_size)/sqrt(3) - 1.
    int exclude_r = -1;
    if (min_dist > 0) {
        exclude_r = static_cast<int>((min_dist * grid.inv_cell_size) * 0.57) - 1;
    }

    int in_start_z = cz - exclude_r;
    int in_end_z = cz + exclude_r;
    int in_start_y = cy - exclude_r;
    int in_end_y = cy + exclude_r;
    int in_start_x = cx - exclude_r;
    int in_end_x = cx + exclude_r;
    
    FloatType min_dist_sq = min_dist * min_dist;
    FloatType max_dist_sq = max_dist * max_dist;
    
    for (int z = start_z; z <= end_z; ++z) {
        bool z_in = (exclude_r >= 0) && (z >= in_start_z && z <= in_end_z);
        for (int y = start_y; y <= end_y; ++y) {
            bool y_in = z_in && (y >= in_start_y && y <= in_end_y);
            for (int x = start_x; x <= end_x; ++x) {
                if (y_in && (x >= in_start_x && x <= in_end_x)) continue;
                
                int bin_idx = z * cpa * cpa + y * cpa + x;
                int start = grid.bin_start[bin_idx];
                int end = grid.bin_start[bin_idx + 1];
                
                for (int idx = start; idx < end; ++idx) {
                    FloatType d2 = norm2(seeds[idx] - seeds[seed_idx]);
                    if (d2 >= min_dist_sq && d2 <= max_dist_sq) {
                        out_neighbors.push_back(idx);
                    }
                }
            }
        }
    }
}

static VoronoiCell compute_cell_vertices(Vec3 seed, const std::vector<int>& candidate_indices, const std::vector<Vec3>& all_seeds,
                             std::vector<Vec3>& dual_points, quickhull::QuickHull<FloatType>& qh, 
                             SimulationStats& stats) {

    auto t0 = std::chrono::steady_clock::now();
    // std::cout << candidates.size() + dual_points.size() << '\n'; 
    constexpr FloatType eps = 1e-16;
    if (seed.x > eps) dual_points.push_back(Vec3(-1.0 / seed.x, 0, 0));
    if (seed.x < 1.0 - eps) dual_points.push_back(Vec3(1.0 / (1.0 - seed.x), 0, 0));
    if (seed.y > eps) dual_points.push_back(Vec3(0, -1.0 / seed.y, 0));
    if (seed.y < 1.0 - eps) dual_points.push_back(Vec3(0, 1.0 / (1.0 - seed.y), 0));
    if (seed.z > eps) dual_points.push_back(Vec3(0, 0, -1.0 / seed.z));
    if (seed.z < 1.0 - eps) dual_points.push_back(Vec3(0, 0, 1.0 / (1.0 - seed.z)));

    for (int idx : candidate_indices) {
        Vec3 r = all_seeds[idx] - seed;
        FloatType len_sq = norm2(r);
        if (len_sq > 1e-16) {
            auto dp = r * (2.0 / len_sq);
            dual_points.push_back(dp);
        }
    }
    // for (auto dp : dual_points) {
    //     std::cout << "dual point " << dp.x << ' ' << dp.y << ' ' << dp.z << '\n'; 
    // }
    auto t1 = std::chrono::steady_clock::now();
    stats.time_ms_dual_setup += std::chrono::duration<FloatType, std::milli>(t1 - t0).count();

    auto t_qh1 = std::chrono::steady_clock::now();
    stats.total_quickhull_calls++;
    // std::cout << dual_points.size() << '\n'; 
    stats.total_dual_points += dual_points.size();
    
    auto mesh = qh.getConvexHullAsMesh(reinterpret_cast<const FloatType*>(dual_points.data()), dual_points.size(), true, quickhull::defaultEps<FloatType>());
    // for (auto mesh_point : mesh.m_vertices) {
    //     // std::cout << "mesh point " << mesh_point.x << ' ' << mesh_point.y << ' ' << mesh_point.z << '\n'; 
    // }
    
    auto t_qh2 = std::chrono::steady_clock::now();
    stats.time_ms_quickhull += std::chrono::duration<FloatType, std::milli>(t_qh2 - t_qh1).count();
    
    const auto& diag = qh.getDiagnostics();
    stats.time_ms_qh_total_internal += diag.m_time_total_micro / 1000.0;
    stats.time_ms_qh_setup_internal += diag.m_time_setup_micro / 1000.0;
    stats.time_ms_qh_process_internal += diag.m_time_process_micro / 1000.0;

    auto t2 = std::chrono::steady_clock::now();

    VoronoiCell cell;
    cell.mesh = std::move(mesh);
    
    if (cell.mesh.m_faces.empty()) {
        auto t3 = std::chrono::steady_clock::now();
        stats.time_ms_primal_recovery += std::chrono::duration<FloatType, std::milli>(t3 - t2).count();
        return cell;
    }

    // Ensure 1:1 mapping between faces and vertices
    cell.vertices.resize(cell.mesh.m_faces.size());

    for (size_t i = 0; i < cell.mesh.m_faces.size(); ++i) {
        const auto& face = cell.mesh.m_faces[i];
        size_t he_idx = face.m_halfEdgeIndex;
        const auto& he1 = cell.mesh.m_halfEdges[he_idx];
        const auto& he2 = cell.mesh.m_halfEdges[he1.m_next];
        const auto& he3 = cell.mesh.m_halfEdges[he2.m_next];
        
        Vec3 u1 = cell.mesh.m_vertices[he3.m_endVertex];
        Vec3 u2 = cell.mesh.m_vertices[he1.m_endVertex];
        Vec3 u3 = cell.mesh.m_vertices[he2.m_endVertex];
        
        Vec3 v;
        if (solve_dual_vertex(u1, u2, u3, v)) {
            cell.vertices[i] = v;
        } else {
            // Degenerate face or numerical issue. 
            // Assign a dummy value or the centroid (though centroid of a line is not helpful).
            // For now, keep it as (0,0,0) or last valid?
            // Since we initialized with default constructor, it is (0,0,0).
            // This case should be extremely rare.
            cell.vertices[i] = (u1 + u2 + u3) * (1.0/3.0); // Centroid of dual face as fallback
        }
    }
    auto t3 = std::chrono::steady_clock::now();
    stats.time_ms_primal_recovery += std::chrono::duration<FloatType, std::milli>(t3 - t2).count();
    return cell;
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

    // Build adjacency: Vertex -> One incoming HalfEdge
    if (ctx.v_in_he.size() < mesh_verts.size()) {
        ctx.v_in_he.resize(mesh_verts.size());
    }
    std::fill(ctx.v_in_he.begin(), ctx.v_in_he.begin() + mesh_verts.size(), invalid);
    
    for (size_t i = 0; i < mesh_he.size(); ++i) {
        size_t v = mesh_he[i].m_endVertex;
        // In valid mesh, v < mesh_verts.size()
        ctx.v_in_he[v] = i;
    }

    // Faces: Dual Vertex -> Voronoi Face
    size_t num_faces = mesh_verts.size();
    poly.face_degree.reserve(num_faces);
    poly.face_indices.reserve(mesh_he.size());
    poly.edges.reserve(mesh_he.size() / 2);

    for (size_t i = 0; i < num_faces; ++i) {
        size_t start_he = ctx.v_in_he[i];
        if (start_he == invalid) continue;

        uint8_t degree = 0;
        size_t curr = start_he;
        do {
            // The face of the current half-edge corresponds to a Voronoi vertex
            size_t f_idx = mesh_he[curr].m_face;
            poly.face_indices.push_back(static_cast<int>(f_idx));
            degree++;
            
            // Move to next half-edge around vertex 'i'
            size_t next_he = mesh_he[curr].m_next;
            curr = mesh_he[next_he].m_opp;
            
        } while (curr != start_he && curr != invalid);

        if (degree > 0) {
            poly.face_degree.push_back(degree);
        }
    }

    // Edges: Dual Edge -> Voronoi Edge
    for (size_t i = 0; i < mesh_he.size(); ++i) {
        size_t opp = mesh_he[i].m_opp;
        if (i < opp) {
            size_t f1 = mesh_he[i].m_face;
            size_t f2 = mesh_he[opp].m_face;
            poly.edges.push_back({static_cast<int>(f1), static_cast<int>(f2)});
        }
    }
}

static Polyhedron convert_to_polyhedron(const VoronoiCell& cell, Vec3 seed, PolyhedraContext& ctx) {
    Polyhedron poly;
    convert_to_polyhedron_inplace(cell, seed, ctx, poly);
    return poly;
}

static FloatType compute_polyhedron_volume(const Polyhedron& poly) {
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
        Vec3 normal = cross_product(v1 - v0, v2 - v0);

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

static void write_polyhedra_vtk(const std::string& path, const std::vector<Polyhedron>& polys) {
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
            cell_id.push_back(static_cast<int>(pi));
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

static void write_cells(const std::string& path, const std::vector<Vec3>& seeds, const std::vector<VoronoiCell>& cells) {
    std::ofstream out(path);
    out << std::setprecision(17);
    out << seeds.size() << "\n";
    for (size_t i = 0; i < seeds.size(); ++i) {
        out << seeds[i].x << " " << seeds[i].y << " " << seeds[i].z << " " << cells[i].vertices.size() << "\n";
        for (const auto& v : cells[i].vertices) {
            out << v.x << " " << v.y << " " << v.z << "\n";
        }
    }
}


#include <vector>
#include <algorithm>
#include <cmath>

static VoronoiCell build_cell(int i, const std::vector<Vec3>& seeds, const Grid& grid, quickhull::QuickHull<FloatType>& qh, SimulationStats& stats, std::vector<int>& neighbor_indices, std::vector<Vec3>& dual_points) {
    FloatType current_dist = grid.cell_size * 1.5;
    FloatType prev_dist = 0.0;
    int step = 0;

    dual_points.clear();

    while (true) {
        auto t_start = std::chrono::steady_clock::now();
        neighbor_indices.clear();
        get_grid_neighbors(grid, seeds, i, prev_dist, current_dist, neighbor_indices);
        auto t_neigh = std::chrono::steady_clock::now();
        stats.time_ms_neighbors += std::chrono::duration<FloatType, std::milli>(t_neigh - t_start).count();
        
        if (step >= stats.neighbors_count_per_step.size()) {
            stats.neighbors_count_per_step.resize(step + 1, 0);
            stats.step_count.resize(step + 1, 0);
        }
        stats.neighbors_count_per_step[step] += neighbor_indices.size();
        stats.step_count[step]++;

        // current_neighbors removal
        auto t_copy = std::chrono::steady_clock::now();
        stats.time_ms_overhead_other += std::chrono::duration<FloatType, std::milli>(t_copy - t_neigh).count();
        
        stats.neighbors_total += neighbor_indices.size();
        // std::cout << "step " << step << " Neighbors size: " << current_neighbors.size() << '\n';
        VoronoiCell cell = compute_cell_vertices(seeds[i], neighbor_indices, seeds, dual_points, qh, stats);
        
        FloatType max_dist_sq = 0.0;
        for (const auto& v : cell.vertices) {
            max_dist_sq = std::max(max_dist_sq, norm2(v));
        }
        FloatType max_dist = std::sqrt(max_dist_sq);
        FloatType safe_dist = current_dist * 0.5;
        
        if (max_dist < safe_dist || current_dist > grid.cell_size * grid.cells_per_axis * 1.5) {
            return cell;
        }
        
        prev_dist = current_dist;

        FloatType target_radius = max_dist * 2;

        FloatType min_step = grid.cell_size * 0.1;
        FloatType next_radius = std::max(target_radius, current_dist + min_step);
        
        // TODO: tweak 1.1
        current_dist = next_radius * 1.1;

        dual_points = std::move(cell.mesh.m_vertices);
        stats.total_expansions++; 
        step++;
    }
}

static std::vector<VoronoiCell> generate_voronoi_diagram(std::vector<Vec3>& seeds, int target_per_cell, SimulationStats& stats) {
    auto t_build_start = std::chrono::steady_clock::now();
    Grid grid = build_grid_and_reorder(seeds, target_per_cell);
    auto t_build_end = std::chrono::steady_clock::now();
    stats.time_ms_build = std::chrono::duration<FloatType, std::milli>(t_build_end - t_build_start).count();


    std::vector<VoronoiCell> all_cells;
    all_cells.reserve(int(seeds.size()));
    
    quickhull::QuickHull<FloatType> qh;
    
    auto t_cells_start = std::chrono::steady_clock::now();

    std::vector<int> neighbor_indices;
    neighbor_indices.reserve(512);
    std::vector<Vec3> dual_points;
    dual_points.reserve(128);

    for (int i = 0; i < int(seeds.size()); ++i) {
        all_cells.push_back(build_cell(i, seeds, grid, qh, stats, neighbor_indices, dual_points));
    }
    
    auto t_cells_end = std::chrono::steady_clock::now();
    stats.time_ms_build = std::chrono::duration<FloatType, std::milli>(t_build_end - t_build_start).count();
    stats.time_ms_cells = std::chrono::duration<FloatType, std::milli>(t_cells_end - t_cells_start).count();
    stats.total_ms = std::chrono::duration<FloatType, std::milli>(t_cells_end - t_build_start).count();
    return all_cells;
}

static bool ends_with(const std::string& s, const std::string& suffix) {
    if (s.size() < suffix.size()) return false;
    return std::equal(suffix.rbegin(), suffix.rend(), s.rbegin());
}

static void run_simulation(const SimulationConfig& config) {
    const FloatType sys_length = 1.0; 
    std::mt19937_64 rng(config.seed_val);
    std::uniform_real_distribution<FloatType> dist(0.0, sys_length);
    std::vector<Vec3> seeds;
    seeds.reserve(config.n);
    for (int i = 0; i < config.n; ++i) {
        seeds.push_back(Vec3(dist(rng), dist(rng), dist(rng)));
    }

    SimulationStats stats;
    auto t_total_start = std::chrono::steady_clock::now();
    auto t_build_start = std::chrono::steady_clock::now();
    Grid grid = build_grid_and_reorder(seeds, config.target_per_cell);
    auto t_build_end = std::chrono::steady_clock::now();
    stats.time_ms_build = std::chrono::duration<FloatType, std::milli>(t_build_end - t_build_start).count();

    double total_volume = 0; 
    uint64_t faces_total = 0;
    uint64_t faces_max = 0;

    quickhull::QuickHull<FloatType> qh;
    std::vector<int> neighbor_indices;
    neighbor_indices.reserve(512);
    std::vector<Vec3> dual_points;
    dual_points.reserve(128);
    PolyhedraContext poly_ctx;
    Polyhedron poly_tmp;
    std::vector<Polyhedron> polys_out;
    const bool write_vtk = (!config.out_path.empty() && ends_with(config.out_path, ".vtk"));
    if (write_vtk) {
        polys_out.reserve(static_cast<size_t>(config.n));
    }

    for (int i = 0; i < config.n; ++i) {
        auto t_cells_start = std::chrono::steady_clock::now();
        VoronoiCell cell = build_cell(i, seeds, grid, qh, stats, neighbor_indices, dual_points);
        auto t_cells_end = std::chrono::steady_clock::now();
        stats.time_ms_cells += std::chrono::duration<FloatType, std::milli>(t_cells_end - t_cells_start).count();

        auto t_poly_start = std::chrono::steady_clock::now();
        convert_to_polyhedron_inplace(cell, seeds[static_cast<size_t>(i)], poly_ctx, poly_tmp);
        auto t_poly_end = std::chrono::steady_clock::now();
        stats.time_ms_polyhedra += std::chrono::duration<FloatType, std::milli>(t_poly_end - t_poly_start).count();

        uint64_t faces = static_cast<uint64_t>(poly_tmp.face_degree.size());
        faces_total += faces;
        faces_max = std::max(faces_max, faces);

        if (config.volume_flag) {
            auto t_vol_start = std::chrono::steady_clock::now();
            total_volume += compute_polyhedron_volume(poly_tmp);
            auto t_vol_end = std::chrono::steady_clock::now();
            stats.time_ms_volume += std::chrono::duration<FloatType, std::milli>(t_vol_end - t_vol_start).count();
        }

        if (write_vtk) {
            polys_out.push_back(poly_tmp);
        }
    }
    auto t_total_end = std::chrono::steady_clock::now();
    stats.total_ms = std::chrono::duration<FloatType, std::milli>(t_total_end - t_total_start).count();

    if (config.voro_flag != 0) {
        stats.voro_ok = run_voro_benchmark(seeds, stats.voro_ms);
    }
    if (write_vtk) {
        write_polyhedra_vtk(config.out_path, polys_out);
    }

    std::cout << std::setprecision(15);
    std::cout << "n=" << config.n << " seed=" << config.seed_val << " target_per_cell=" << config.target_per_cell << "\n";
    if (config.volume_flag) {
        std::cout << "sum_volume=" << total_volume << " error=" << std::abs(total_volume - sys_length*sys_length*sys_length) << "\n";
    }
    std::cout << "neighbors_total=" << stats.neighbors_total << " neighbors_avg=" << (static_cast<FloatType>(stats.neighbors_total) / config.n) << "\n";
    std::cout << "faces_total=" << faces_total
              << " faces_mean=" << (config.n <= 0 ? 0.0 : (static_cast<double>(faces_total) / config.n))
              << " faces_max=" << faces_max << "\n";
    std::cout << "expansions_total=" << stats.total_expansions << " expansions_avg=" << (static_cast<double>(stats.total_expansions) / config.n) << "\n";
    std::cout << "dual_points_total=" << stats.total_dual_points << " dual_points_avg=" << (static_cast<double>(stats.total_dual_points) / stats.total_quickhull_calls) << "\n";
    for (size_t i = 0; i < stats.step_count.size(); ++i) {
        if (stats.step_count[i] > 0) {
            std::cout << "step " << i << ": avg_neighbors=" 
                      << (static_cast<double>(stats.neighbors_count_per_step[i]) / stats.step_count[i]) 
                      << " (count=" << stats.step_count[i] << ")\n";
        }
    }
    std::cout << "quickhull_calls_total=" << stats.total_quickhull_calls << " quickhull_calls_avg=" << (static_cast<double>(stats.total_quickhull_calls) / config.n) << "\n";
    std::cout << "time_ms_total=" << stats.total_ms << " time_ms_build=" << stats.time_ms_build << " time_ms_cells=" << stats.time_ms_cells
              << " time_us_per_cell=" << (stats.time_ms_cells * 1000.0 / config.n) << "\n";
    std::cout << "time_ms_quickhull=" << stats.time_ms_quickhull << " (" << (stats.time_ms_quickhull / stats.time_ms_cells * 100.0) << "% of cells time)\n";
    std::cout << "  [Detailed Breakdown]\n";
    std::cout << "  neighbors_search=" << stats.time_ms_neighbors << " ms (" << (stats.time_ms_neighbors / stats.time_ms_cells * 100.0) << "%)\n";
    std::cout << "  dual_setup=" << stats.time_ms_dual_setup << " ms (" << (stats.time_ms_dual_setup / stats.time_ms_cells * 100.0) << "%)\n";
    std::cout << "  primal_recovery=" << stats.time_ms_primal_recovery << " ms (" << (stats.time_ms_primal_recovery / stats.time_ms_cells * 100.0) << "%)\n";
    std::cout << "  overhead_other=" << stats.time_ms_overhead_other << " ms (" << (stats.time_ms_overhead_other / stats.time_ms_cells * 100.0) << "%)\n";
    std::cout << "  [QuickHull Internal]\n";
    std::cout << "  qh_total=" << stats.time_ms_qh_total_internal << " ms\n";
    std::cout << "  qh_setup=" << stats.time_ms_qh_setup_internal << " ms\n";
    std::cout << "  qh_process=" << stats.time_ms_qh_process_internal << " ms\n";
    std::cout << "  qh_mesh_convert=" << stats.time_ms_qh_mesh_convert << " ms\n";
    std::cout << "  [Polyhedra Conversion]\n";
    std::cout << "  total_time=" << stats.time_ms_polyhedra << " ms (" << (stats.time_ms_polyhedra / stats.time_ms_cells * 100.0) << "% of cells time)\n";
    if (config.volume_flag) {
        std::cout << "  [Volume]\n";
        std::cout << "  total_time=" << stats.time_ms_volume << " ms\n";
    }

    if (config.voro_flag == 0) {
        std::cout << "voro_ms=off\n";
    } else if (stats.voro_ok) {
        std::cout << "voro_ms=" << stats.voro_ms << "\n";
    } else {
        std::cout << "voro_ms=na\n";
    }
}

int main(int argc, char** argv) {
    SimulationConfig config;
    config.n = 200;
    config.seed_val = 1;
    config.target_per_cell = 4;
    config.out_path = "polyhedra.vtk";
    config.volume_flag = 0;
    config.voro_flag = 1;
    
    if (argc > 1) config.n = std::max(1, std::atoi(argv[1]));
    if (argc > 2) config.seed_val = static_cast<uint64_t>(std::strtoull(argv[2], nullptr, 10));
    if (argc > 3) config.target_per_cell = std::max(1, std::atoi(argv[3]));
    if (argc > 4) config.out_path = argv[4];
    if (argc > 5) config.volume_flag = std::atoi(argv[5]) != 0 ? 1 : 0;
    if (argc > 6) config.voro_flag = std::atoi(argv[6]) != 0 ? 1 : 0;

    run_simulation(config);
    
    return 0;
}
#endif
