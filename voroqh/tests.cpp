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

static void test_random_points_deterministic() {
    std::cout << "Test random_points_deterministic: N=32 seeds in unit box, seed=123, box_len=1.0 ...";
    auto a = voronoi::generate_random_points_box(32, 123, 1.0);
    auto b = voronoi::generate_random_points_box(32, 123, 1.0);
    assert(a.size() == b.size());
    for (size_t i = 0; i < a.size(); ++i) {
        assert(a[i].x == b[i].x);
        assert(a[i].y == b[i].y);
        assert(a[i].z == b[i].z);
    }
    std::cout << " OK\n";
}

static void test_polyhedra_topology_and_euler() {
    const size_t n = 128;
    std::cout << "Test polyhedra_topology_and_euler: N=128 seeds in unit box, seed=42, target_per_cell=4, check V,E,F and Euler characteristic V-E+F=2 with E derived from sum(face_degree) ...";
    auto seeds = voronoi::generate_random_points_box(n, 42, 1.0);
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
    std::cout << " OK\n";
}

static void test_polyhedra_basic_invariants_and_bbox() {
    const size_t n = 200;
    std::cout << "Test polyhedra_basic_invariants_and_bbox: N=200 seeds in unit box, seed=1, target_per_cell=4, check volume_sum≈1.0, vertex indices, and global bbox within [0,1]^3 plus eps=1e-7 ...";
    auto seeds = voronoi::generate_random_points_box(n, 1, 1.0);
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

    const double eps = 1e-7;
    assert(global_min >= -eps);
    assert(global_max <= 1.0 + eps);

    assert(stats.time_ms_quickhull >= 0.0);
    assert(stats.time_ms_qh_total_internal >= 0.0);
    assert(stats.time_ms_qh_mesh_convert >= 0.0);
    assert(stats.time_ms_quickhull + 1e-9 >= stats.time_ms_qh_total_internal);
    std::cout << " OK\n";
}

static void test_polyhedra_surface_area_consistency() {
    const size_t n = 64;
    std::cout << "Test polyhedra_surface_area_consistency: N=64 seeds in unit box, seed=5, target_per_cell=4, compare compute_polyhedron_face_areas sum with independent triangulated surface area from poly_surface_area using relative tol=1e-9 ...";
    auto seeds = voronoi::generate_random_points_box(n, 5, 1.0);
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
            assert(approx(sum_areas, triangulated_area, 1e-9 * scale));

            cell_count++;
        });

    assert(cell_count == n);
    std::cout << " OK\n";
}

static void test_volume_and_centroid_translation_invariance() {
    const size_t n = 32;
    std::cout << "Test volume_and_centroid_translation_invariance: N=32 seeds in unit box, seed=7, target_per_cell=4, test volume and centroid invariance under translation delta=(0.1,-0.07,0.05) for first polyhedron ...";
    auto seeds = voronoi::generate_random_points_box(n, 7, 1.0);
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
            const voronoi::Vec3 delta(0.1, -0.07, 0.05);
            for (auto& v : shifted.vertices) {
                v = v + delta;
            }

            const double vol2 = voronoi::compute_polyhedron_volume(shifted);
            const auto c2 = voronoi::compute_polyhedron_centroid(shifted);

            assert(approx(vol1, vol2, 1e-12));
            assert(approx(c1.x + delta.x, c2.x, 1e-10));
            assert(approx(c1.y + delta.y, c2.y, 1e-10));
            assert(approx(c1.z + delta.z, c2.z, 1e-10));

            tested = true;
        });

    assert(tested);
    std::cout << " OK\n";
}

static void test_two_seed_voronoi_split() {
    std::cout << "Test two_seed_voronoi_split: N=2 seeds at x=0.25 and x=0.75 in unit box, target_per_cell=2, expect plane x=0.5, each cell volume≈0.5 and centroids on opposite sides of x=0.5 ...";
    std::vector<voronoi::Vec3> seeds;
    seeds.emplace_back(0.25, 0.5, 0.5);
    seeds.emplace_back(0.75, 0.5, 0.5);

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
    assert(approx(total_vol, 1.0, 1e-12));
    assert(approx(volumes[0], 0.5, 1e-9));
    assert(approx(volumes[1], 0.5, 1e-9));
    assert(centroids[0].x < 0.5);
    assert(centroids[1].x > 0.5);
    std::cout << " OK\n";
}

static void test_clustered_seeds_positive_volume() {
    std::cout << "Test clustered_seeds_positive_volume: N=5 nearly coincident seeds clustered around (0.5,0.5,0.5) with offsets up to 2e-6, target_per_cell=4, require all cell volumes>1e-15 ...";
    std::vector<voronoi::Vec3> seeds;
    seeds.emplace_back(0.5, 0.5, 0.5);
    seeds.emplace_back(0.500001, 0.5, 0.5);
    seeds.emplace_back(0.5, 0.500001, 0.5);
    seeds.emplace_back(0.5, 0.5, 0.500001);
    seeds.emplace_back(0.500002, 0.500002, 0.5);

    voronoi::VoronoiStats stats;

    size_t count = 0;
    voronoi::for_each_polyhedron(
        seeds,
        4,
        stats,
        [&](size_t, const voronoi::Polyhedron& poly) {
            const double vol = voronoi::compute_polyhedron_volume(poly);
            assert(vol > 1e-15);
            count++;
        });

    assert(count == seeds.size());
    std::cout << " OK\n";
}

static void test_volume_summation_stability_and_statistics() {
    const size_t n = 256;
    std::cout << "Test volume_summation_stability_and_statistics: N=256 seeds in unit box, seed=99, target_per_cell=8, check volume_sum≈1.0 for two summation orders, and statistics on faces, edges, and diameter within bounds (faces∈[4,128], edges∈[6,256], diameter∈(0,4]) ...";
    auto seeds = voronoi::generate_random_points_box(n, 99, 1.0);
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
            assert(diameter <= 4.0);
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

    assert(approx(sum1, sum2, 1e-12));
    assert(approx(sum1, 1.0, 1e-6));

    assert(min_faces >= 4);
    assert(max_faces <= 128);
    assert(min_edges >= 6);
    assert(max_edges <= 256);
    assert(min_diameter > 0.0);
    assert(max_diameter <= 4.0);
    std::cout << " OK\n";
}

static void test_vtk_writer_smoke() {
    std::cout << "Test vtk_writer_smoke: write two identical tetrahedral Polyhedron cells with ids 7 and 9 to test_polyhedra.vtk and verify presence of VTK POLYDATA headers and cell_id scalar field, then delete file ...";
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
    std::vector<int> ids = {7, 9};

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
    std::cout << " OK\n";
}

int main() {
    test_random_points_deterministic();
    test_polyhedra_topology_and_euler();
    test_polyhedra_basic_invariants_and_bbox();
    test_polyhedra_surface_area_consistency();
    test_volume_and_centroid_translation_invariance();
    test_two_seed_voronoi_split();
    test_clustered_seeds_positive_volume();
    test_volume_summation_stability_and_statistics();
    test_vtk_writer_smoke();
    return 0;
}
