#include "voronoi.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

static bool approx(double a, double b, double tol) {
    return std::abs(a - b) <= tol;
}

static voronoi::Vec3 cross3(voronoi::Vec3 a, voronoi::Vec3 b) {
    return voronoi::Vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

static double poly_surface_area(const voronoi::Polyhedron& poly) {
    double area = 0.0;
    size_t off = 0;
    for (uint8_t deg_u8 : poly.face_degree) {
        const int deg = static_cast<int>(deg_u8);
        if (deg < 3) {
            off += static_cast<size_t>(deg);
            continue;
        }
        const int i0 = poly.face_indices[off];
        const auto p0 = poly.vertices[static_cast<size_t>(i0)];
        for (int k = 1; k < deg - 1; ++k) {
            const int ia = poly.face_indices[off + static_cast<size_t>(k)];
            const int ib = poly.face_indices[off + static_cast<size_t>(k + 1)];
            const auto p1 = poly.vertices[static_cast<size_t>(ia)];
            const auto p2 = poly.vertices[static_cast<size_t>(ib)];
            const auto cp = cross3(p1 - p0, p2 - p0);
            area += 0.5 * std::sqrt(cp.x * cp.x + cp.y * cp.y + cp.z * cp.z);
        }
        off += static_cast<size_t>(deg);
    }
    return area;
}

struct VoroCellMetrics {
    bool present = false;
    double volume = 0.0;
    int faces = 0;
    int vertices = 0;
    int edges = 0;
    double surface_area = 0.0;
};

static void test_polyhedra_topology_and_euler() {
    const size_t n = 4096;
    const uint64_t seed = 42;

    std::cout << "N=" << n
              << " seeds in unit box, seed=" << seed
              << ", check V,E,F and Euler characteristic V-E+F=2 with E derived from sum(face_degree)";
    std::cout.flush();

    auto seeds = voronoi::generate_random_points_box(n, seed, 1.0);
    voronoi::VoronoiStats stats;

    size_t cell_count = 0;

    voronoi::for_each_polyhedron(
        seeds,
        4,
        stats,
        [&](size_t, const voronoi::Polyhedron& poly) {
            const int V = static_cast<int>(poly.vertices.size());
            const int F = static_cast<int>(poly.face_degree.size());
            assert(V >= 4);
            assert(F >= 4);

            size_t sum_deg = 0;
            for (uint8_t d_u8 : poly.face_degree) {
                const int d = static_cast<int>(d_u8);
                assert(d >= 3);
                sum_deg += static_cast<size_t>(d);
            }

            assert((sum_deg % 2) == 0);

            const int E = static_cast<int>(sum_deg / 2);
            assert(E >= 6);

            assert(V - E + F == 2);

            cell_count++;
        });

    assert(cell_count == n);
}

static void test_polyhedra_basic_invariants_and_bbox() {
    const size_t n = 4096;
    const uint64_t seed = 1;
    const double eps = 1e-7;

    std::cout << "N=" << n
              << " seeds in unit box, seed=" << seed
              << ", check volume_sum≈1.0, vertex indices, and global bbox within [0,1]^3 plus eps=" << eps;
    std::cout.flush();

    auto seeds = voronoi::generate_random_points_box(n, seed, 1.0);
    voronoi::VoronoiStats stats;

    double volume_sum = 0.0;
    size_t poly_count = 0;

    double global_min = std::numeric_limits<double>::infinity();
    double global_max = -std::numeric_limits<double>::infinity();

    voronoi::for_each_polyhedron(
        seeds,
        4,
        stats,
        [&](size_t, const voronoi::Polyhedron& poly) {
            assert(!poly.vertices.empty());
            assert(!poly.face_degree.empty());
            size_t face_indices_count = 0;
            for (uint8_t deg_u8 : poly.face_degree) {
                const int deg = static_cast<int>(deg_u8);
                assert(deg >= 3);
                face_indices_count += static_cast<size_t>(deg);
            }
            assert(face_indices_count == poly.face_indices.size());
            for (int idx : poly.face_indices) {
                assert(idx >= 0);
                assert(static_cast<size_t>(idx) < poly.vertices.size());
            }

            for (const auto& v : poly.vertices) {
                global_min = std::min(global_min, std::min(v.x, std::min(v.y, v.z)));
                global_max = std::max(global_max, std::max(v.x, std::max(v.y, v.z)));
            }

            volume_sum += voronoi::compute_polyhedron_volume(poly);
            poly_count++;
        });

    assert(poly_count == n);
    assert(approx(volume_sum, 1.0, 1e-6));

    assert(global_min >= -eps);
    assert(global_max <= 1.0 + eps);

    assert(stats.time_ms_quickhull >= 0.0);
    assert(stats.time_ms_qh_total_internal >= 0.0);
    assert(stats.time_ms_qh_mesh_convert >= 0.0);
    assert(stats.time_ms_quickhull + 1e-9 >= stats.time_ms_qh_total_internal);
}

static void test_polyhedra_surface_area_consistency() {
    const size_t n = 2048;
    const uint64_t seed = 5;
    const double tol_rel = 1e-9;

    std::cout << "N=" << n
              << " seeds in unit box, seed=" << seed
              << ", compare compute_polyhedron_face_areas sum with independent triangulated surface area from poly_surface_area using relative tol=" << tol_rel;
    std::cout.flush();

    auto seeds = voronoi::generate_random_points_box(n, seed, 1.0);
    voronoi::VoronoiStats stats;

    size_t cell_count = 0;

    voronoi::for_each_polyhedron(
        seeds,
        4,
        stats,
        [&](size_t, const voronoi::Polyhedron& poly) {
            std::vector<double> areas;
            voronoi::compute_polyhedron_face_areas(poly, areas);

            double sum_areas = 0.0;
            for (double a : areas) {
                assert(a > 0.0);
                sum_areas += a;
            }

            const double triangulated_area = poly_surface_area(poly);
            const double scale = std::max(1.0, triangulated_area);
            assert(approx(sum_areas, triangulated_area, tol_rel * scale));

            cell_count++;
        });

    assert(cell_count == n);
}

static void test_volume_and_centroid_translation_invariance() {
    const size_t n = 32;
    const uint64_t seed = 7;
    const double dx = 0.1;
    const double dy = -0.07;
    const double dz = 0.05;
    const double tol_vol = 1e-12;
    const double tol_centroid = 1e-10;

    std::cout << "N=" << n
              << " seeds in unit box, seed=" << seed
              << ", test volume and centroid invariance under translation delta=("
              << dx << "," << dy << "," << dz
              << "), volume tol=" << tol_vol
              << ", centroid tol=" << tol_centroid;
    std::cout.flush();

    auto seeds = voronoi::generate_random_points_box(n, seed, 1.0);
    voronoi::VoronoiStats stats;

    bool tested = false;

    voronoi::for_each_polyhedron(
        seeds,
        4,
        stats,
        [&](size_t, const voronoi::Polyhedron& poly) {
            if (tested) {
                return;
            }

            const double vol1 = voronoi::compute_polyhedron_volume(poly);
            const auto c1 = voronoi::compute_polyhedron_centroid(poly);

            voronoi::Polyhedron shifted = poly;
            const voronoi::Vec3 delta(dx, dy, dz);
            for (auto& v : shifted.vertices) {
                v = v + delta;
            }

            const double vol2 = voronoi::compute_polyhedron_volume(shifted);
            const auto c2 = voronoi::compute_polyhedron_centroid(shifted);

            assert(approx(vol1, vol2, tol_vol));
            assert(approx(c1.x + delta.x, c2.x, tol_centroid));
            assert(approx(c1.y + delta.y, c2.y, tol_centroid));
            assert(approx(c1.z + delta.z, c2.z, tol_centroid));

            tested = true;
        });

    assert(tested);
}

static void test_two_seed_voronoi_split() {
    const double left_x = 0.25;
    const double right_x = 0.75;
    const double split_x = 0.5;
    const double tol_total = 1e-12;
    const double tol_vol = 1e-9;

    std::cout << "N=2 seeds at x=" << left_x << " and x=" << right_x
              << " in unit box, expect split plane x=" << split_x
              << ", each cell volume≈0.5 within tol=" << tol_vol
              << " and total volume≈1.0 within tol=" << tol_total;
    std::cout.flush();

    std::vector<voronoi::Vec3> seeds;
    seeds.emplace_back(left_x, 0.5, 0.5);
    seeds.emplace_back(right_x, 0.5, 0.5);

    voronoi::VoronoiStats stats;

    double volumes[2] = {0.0, 0.0};
    voronoi::Vec3 centroids[2];
    bool present[2] = {false, false};

    voronoi::for_each_polyhedron(
        seeds,
        2,
        stats,
        [&](size_t index, const voronoi::Polyhedron& poly) {
            assert(index < 2);
            const double vol = voronoi::compute_polyhedron_volume(poly);
            const auto c = voronoi::compute_polyhedron_centroid(poly);
            assert(vol > 0.0);
            volumes[index] = vol;
            centroids[index] = c;
            present[index] = true;
        });

    assert(present[0] && present[1]);

    const double total_vol = volumes[0] + volumes[1];
    assert(approx(total_vol, 1.0, tol_total));
    assert(approx(volumes[0], 0.5, tol_vol));
    assert(approx(volumes[1], 0.5, tol_vol));
    assert(centroids[0].x < split_x);
    assert(centroids[1].x > split_x);
}

static void test_clustered_seeds_positive_volume() {
    const size_t n = 5;
    const double cx = 0.5;
    const double cy = 0.5;
    const double cz = 0.5;
    const double offset1 = 1e-6;
    const double offset2 = 2e-6;
    const double min_vol = 1e-15;

    std::cout << "N=" << n
              << " nearly coincident seeds clustered around ("
              << cx << "," << cy << "," << cz
              << ") with offsets up to " << offset2
              << ", require all cell volumes>" << min_vol;
    std::cout.flush();

    std::vector<voronoi::Vec3> seeds;
    seeds.emplace_back(cx, cy, cz);
    seeds.emplace_back(cx + offset1, cy, cz);
    seeds.emplace_back(cx, cy + offset1, cz);
    seeds.emplace_back(cx, cy, cz + offset1);
    seeds.emplace_back(cx + offset2, cy + offset2, cz);

    voronoi::VoronoiStats stats;

    size_t count = 0;
    voronoi::for_each_polyhedron(
        seeds,
        4,
        stats,
        [&](size_t, const voronoi::Polyhedron& poly) {
            const double vol = voronoi::compute_polyhedron_volume(poly);
            assert(vol > min_vol);
            count++;
        });

    assert(count == seeds.size());
}

static void test_volume_summation_stability_and_statistics() {
    const size_t n = 8192;
    const uint64_t seed = 99;
    const double tol_sum = 1e-12;
    const double tol_vol = 1e-6;
    const int faces_min = 4;
    const int faces_max = 128;
    const int edges_min = 6;
    const int edges_max = 256;
    const double diameter_max = 4.0;

    std::cout << "N=" << n
              << " seeds in unit box, seed=" << seed
              << ", check volume_sum≈1.0 within tol=" << tol_vol
              << " for two summation orders (tolerance for sum difference=" << tol_sum
              << ") and statistics on faces, edges, and diameter with faces∈["
              << faces_min << "," << faces_max << "], edges∈["
              << edges_min << "," << edges_max << "], diameter∈(0," << diameter_max << "]";
    std::cout.flush();

    auto seeds = voronoi::generate_random_points_box(n, seed, 1.0);
    voronoi::VoronoiStats stats;

    std::vector<double> volumes;
    volumes.reserve(n);

    int min_faces = std::numeric_limits<int>::max();
    int max_faces = 0;
    int min_edges = std::numeric_limits<int>::max();
    int max_edges = 0;
    double min_diameter = std::numeric_limits<double>::infinity();
    double max_diameter = 0.0;

    voronoi::for_each_polyhedron(
        seeds,
        8,
        stats,
        [&](size_t seed_index, const voronoi::Polyhedron& poly) {
            const double vol = voronoi::compute_polyhedron_volume(poly);
            assert(vol > 0.0);
            volumes.push_back(vol);

            const int F = static_cast<int>(poly.face_degree.size());
            size_t sum_deg = 0;
            for (uint8_t d_u8 : poly.face_degree) {
                const int d = static_cast<int>(d_u8);
                assert(d >= 3);
                sum_deg += static_cast<size_t>(d);
            }
            assert((sum_deg % 2) == 0);
            const int E = static_cast<int>(sum_deg / 2);

            min_faces = std::min(min_faces, F);
            max_faces = std::max(max_faces, F);
            min_edges = std::min(min_edges, E);
            max_edges = std::max(max_edges, E);

            const double diameter = voronoi::compute_polyhedron_diameter(poly, seeds[seed_index]);
            assert(diameter > 0.0);
            assert(diameter <= diameter_max);
            min_diameter = std::min(min_diameter, diameter);
            max_diameter = std::max(max_diameter, diameter);
        });

    assert(volumes.size() == n);

    double sum1 = 0.0;
    for (double v : volumes) {
        sum1 += v;
    }

    std::vector<double> sorted = volumes;
    std::sort(sorted.begin(), sorted.end(), [](double a, double b) { return std::abs(a) < std::abs(b); });

    double sum2 = 0.0;
    for (double v : sorted) {
        sum2 += v;
    }

    assert(approx(sum1, sum2, tol_sum));
    assert(approx(sum1, 1.0, tol_vol));

    assert(min_faces >= faces_min);
    assert(max_faces <= faces_max);
    assert(min_edges >= edges_min);
    assert(max_edges <= edges_max);
    assert(min_diameter > 0.0);
    assert(max_diameter <= diameter_max);
}

static void test_vertices_near_side_of_cutting_planes() {
    const size_t n = 10000;
    const uint64_t seed = 1234;
    const double tol = 1e-15;

    std::cout << "N=" << n
              << " seeds in unit box, seed=" << seed
              << ", for each generator i, every Voronoi cell vertex, and every other generator j, verify vertex lies on near side of perpendicular bisector cutting plane between seeds i and j within tol=" << tol;
    std::cout.flush();

    auto seeds = voronoi::generate_random_points_box(n, seed, 1.0);
    voronoi::VoronoiStats stats;

    voronoi::for_each_polyhedron(
        seeds,
        8,
        stats,
        [&](size_t i, const voronoi::Polyhedron& poly) {
            const voronoi::Vec3 si = seeds[i];
            const double six = si.x;
            const double siy = si.y;
            const double siz = si.z;
            const double si2 = six * six + siy * siy + siz * siz;

            for (const auto& v : poly.vertices) {
                for (size_t j = 0; j < seeds.size(); ++j) {
                    if (j == i) {
                        continue;
                    }
                    const voronoi::Vec3 sj = seeds[j];
                    const double sjx = sj.x;
                    const double sjy = sj.y;
                    const double sjz = sj.z;
                    const double sj2 = sjx * sjx + sjy * sjy + sjz * sjz;

                    const double nx = sjx - six;
                    const double ny = sjy - siy;
                    const double nz = sjz - siz;

                    const double lhs = v.x * nx + v.y * ny + v.z * nz;
                    const double rhs = 0.5 * (sj2 - si2);

                    assert(lhs <= rhs + tol);
                }
            }
        });
}

static void test_vtk_writer_smoke() {
    const int num_cells = 2;
    const int id0 = 7;
    const int id1 = 9;

    std::cout << "write " << num_cells
              << " tetrahedral Polyhedron cells with ids " << id0 << " and " << id1
              << " to test_polyhedra.vtk and verify presence of VTK POLYDATA headers and cell_id scalar field, then delete file";
    std::cout.flush();

    voronoi::Polyhedron tetra;
    tetra.vertices = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
    };
    tetra.face_degree = {3, 3, 3, 3};
    tetra.face_indices = {
        0, 2, 1,
        0, 1, 3,
        0, 3, 2,
        1, 2, 3,
    };

    std::vector<voronoi::Polyhedron> polys = {tetra, tetra};
    std::vector<int> ids = {id0, id1};

    const std::string path = "test_polyhedra.vtk";
    voronoi::write_polyhedra_vtk(path, polys, ids);

    std::ifstream in(path);
    assert(in.good());
    std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    assert(content.find("DATASET POLYDATA") != std::string::npos);
    assert(content.find("POINTS ") != std::string::npos);
    assert(content.find("POLYGONS ") != std::string::npos);
    assert(content.find("CELL_DATA ") != std::string::npos);
    assert(content.find("SCALARS cell_id int 1") != std::string::npos);

    std::remove(path.c_str());
}

static void run_test(int index, void (*test)()) {
    std::cout << "Test #" << index << ": ";
    std::cout.flush();
    test();
    std::cout << " OK\n\n";
}

int main() {
    int index = 1;
    run_test(index++, test_polyhedra_topology_and_euler);
    run_test(index++, test_polyhedra_basic_invariants_and_bbox);
    run_test(index++, test_polyhedra_surface_area_consistency);
    run_test(index++, test_volume_and_centroid_translation_invariance);
    run_test(index++, test_two_seed_voronoi_split);
    run_test(index++, test_clustered_seeds_positive_volume);
    run_test(index++, test_volume_summation_stability_and_statistics);
    run_test(index++, test_vertices_near_side_of_cutting_planes);
    run_test(index++, test_vtk_writer_smoke);
    return 0;
}
