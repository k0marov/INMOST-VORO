#include "voronoi_builder.h"
#include "voronoi.hpp"
#include <cmath>
#include <iostream>

VoronoiBuilder::VoronoiBuilder(const std::vector<std::tuple<double, double, double>>& seeds, int target_per_cell)
    : seeds(seeds), target_per_cell(target_per_cell) {
}

void add_polyhedron_to_mesh(Mesh& mesh, const voronoi::Polyhedron& poly) {
    ElementArray<Node> nodes(&mesh);
    nodes.reserve(poly.vertices.size());
    for (const auto& v : poly.vertices) {
        Storage::real coords[3] = {v.x, v.y, v.z};
        nodes.push_back(mesh.CreateNode(coords));
    }

    ElementArray<Face> faces(&mesh);
    size_t offset = 0;
    for (uint8_t deg : poly.face_degree) {
        ElementArray<Node> face_nodes(&mesh);
        face_nodes.reserve(deg);
        for (int k = 0; k < deg; ++k) {
            int node_idx = poly.face_indices[offset + k];
            face_nodes.push_back(nodes[node_idx]);
        }
        faces.push_back(mesh.CreateFace(face_nodes).first);
        offset += deg;
    }

    mesh.CreateCell(faces);
}

Mesh VoronoiBuilder::build(voronoi::VoronoiStats* stats_out) {
    Mesh global_mesh;
    global_mesh.RemTopologyCheck(DEFAULT_CHECK); 
    global_mesh.SetDimensions(3);

    // Convert seeds to voroqh format
    std::vector<voronoi::Vec3> voro_seeds;
    voro_seeds.reserve(seeds.size());
    for (const auto& s : seeds) {
        voro_seeds.emplace_back(std::get<0>(s), std::get<1>(s), std::get<2>(s));
    }

    voronoi::VoronoiStats stats_tmp;

    auto polys_out = voronoi::generate_voronoi_diagram(
        voro_seeds,
        target_per_cell,
        stats_tmp
    );

    const auto t_inmost_start = std::chrono::steady_clock::now();
    for (size_t i = 0; i < seeds.size(); ++i) {
        add_polyhedron_to_mesh(global_mesh, polys_out[i]);
    }
    const auto t_inmost_end = std::chrono::steady_clock::now();
    stats_tmp.time_ms_callback = std::chrono::duration<double, std::milli>(t_inmost_end - t_inmost_start).count();
    if (stats_out != nullptr) {
        *stats_out = stats_tmp;
    }

    return global_mesh;
}
