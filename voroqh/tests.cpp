#include "voronoi.hpp"

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
    auto a = voronoi::generate_random_points_box(32, 123, 1.0);
    auto b = voronoi::generate_random_points_box(32, 123, 1.0);
    assert(a.size() == b.size());
    for (size_t i = 0; i < a.size(); ++i) {
        assert(a[i].x == b[i].x);
        assert(a[i].y == b[i].y);
        assert(a[i].z == b[i].z);
    }
}

static void test_polyhedra_basic_invariants_and_bbox() {
    const size_t n = 200;
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
}

static void test_vtk_writer_smoke() {
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
}

int main() {
    test_random_points_deterministic();
    test_polyhedra_basic_invariants_and_bbox();
    test_vtk_writer_smoke();
    return 0;
}
