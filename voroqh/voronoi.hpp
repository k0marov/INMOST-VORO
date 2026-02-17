#pragma once

#include <cstdint>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "quickhull/Structs/Vector3.hpp"

namespace voronoi {

using FloatType = double;
using Vec3 = quickhull::Vector3<FloatType>;

struct Polyhedron {
    std::vector<Vec3> vertices;
    std::vector<int> face_indices;
    std::vector<uint8_t> face_degree;
    std::vector<std::pair<int, int>> edges;
};

struct VoronoiStats {
    uint64_t neighbors_total = 0;
    uint64_t total_expansions = 0;
    uint64_t total_quickhull_calls = 0;
    FloatType time_ms_quickhull = 0.0;
    FloatType time_ms_build = 0.0;
    FloatType time_ms_cells = 0.0;
    FloatType total_ms = 0.0;
    std::vector<uint64_t> neighbors_count_per_step;
    std::vector<uint64_t> step_count;
    uint64_t total_dual_points = 0;

    FloatType time_ms_neighbors = 0.0;
    FloatType time_ms_dual_setup = 0.0;
    FloatType time_ms_primal_recovery = 0.0;
    FloatType time_ms_overhead_other = 0.0;

    FloatType time_ms_qh_total_internal = 0.0;
    FloatType time_ms_qh_setup_internal = 0.0;
    FloatType time_ms_qh_process_internal = 0.0;
    FloatType time_ms_qh_mesh_convert = 0.0;

    FloatType time_ms_polyhedra = 0.0;
};

std::vector<Vec3> generate_random_points_box(size_t n, uint64_t seed, FloatType box_len = 1.0);

void for_each_polyhedron(
    const std::vector<Vec3>& seeds,
    int target_per_cell,
    VoronoiStats& stats,
    const std::function<void(size_t seed_index, const Polyhedron& poly)>& on_polyhedron);

FloatType compute_polyhedron_volume(const Polyhedron& poly);

void compute_polyhedron_face_areas(const Polyhedron& poly, std::vector<double>& out_areas);

void write_voro_compatible_output(std::ostream& out, size_t seed_index, const Vec3& seed, const Polyhedron& poly, double vol, std::vector<double>& area_buffer);

void write_polyhedra_vtk(
    const std::string& path,
    const std::vector<Polyhedron>& polys,
    const std::vector<int>& cell_id_per_poly);

}  // namespace voronoi
