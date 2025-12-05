#include "voronoi_builder.h"
#include "planar_cutter.h" // Now includes the new Polyhedron structs
#include <cmath>
#include <iostream>
#include <map>
#include <tuple>

using std::get;

// Custom comparator for Point3D to use it in a map
//struct Point3DComparator {
//    bool operator()(const Point3D& a, const Point3D& b) const {
//        if (std::abs(a.x - b.x) > GEOMETRY_EPSILON) return a.x < b.x;
//        if (std::abs(a.y - b.y) > GEOMETRY_EPSILON) return a.y < b.y;
//        if (std::abs(a.z - b.z) > GEOMETRY_EPSILON) return a.z < b.z;
//        return false;
//    }
//};
//

//Plane bisector_plane(const Vec3 &seed, const Vec3 &neighbor) {
//  Vec3 midpoint = (seed + neighbor) * 0.5;
//  Vec3 normal = normalize(neighbor - seed); // points from seed toward neighbor
//  double d = -dot(normal, midpoint);         // plane passes through midpoint
//  // To keep the seed point "inside", flip normal if necessary
//  if(dot(normal, seed) + d > 0) normal = normal * -1;
//  return Plane(normal, d);
//}
//
VoronoiBuilder::VoronoiBuilder(const std::vector<std::tuple<double, double, double>>& seeds, SystemSize system_size)
    : seeds(seeds), system_size(system_size) {

    double min_x = std::get<0>(system_size).first;
    double max_x = std::get<0>(system_size).second;
    double min_y = std::get<1>(system_size).first;
    double max_y = std::get<1>(system_size).second;
    double min_z = std::get<2>(system_size).first;
    double max_z = std::get<2>(system_size).second;

    double lx = max_x - min_x;
    double ly = max_y - min_y;
    double lz = max_z - min_z;

    if (lx <= 0 || ly <= 0 || lz <= 0 || seeds.empty()) {
        nx = 1; ny = 1; nz = 1;
    } else {
        double system_volume = lx * ly * lz;
        double seed_density = seeds.size() / system_volume;
        const double seeds_per_container = 8.0;
        double ideal_container_volume = seeds_per_container / seed_density;
        double ideal_container_side = std::cbrt(ideal_container_volume);
        nx = static_cast<int>(std::round(lx / ideal_container_side));
        ny = static_cast<int>(std::round(ly / ideal_container_side));
        nz = static_cast<int>(std::round(lz / ideal_container_side));
    }

    if (nx <= 0) nx = 1;
    if (ny <= 0) ny = 1;
    if (nz <= 0) nz = 1;

    container_size_x = lx / nx;
    container_size_y = ly / ny;
    container_size_z = lz / nz;

    containers.resize(nx, std::vector<std::vector<std::vector<int>>>(ny, std::vector<std::vector<int>>(nz)));

    for (size_t i = 0; i < seeds.size(); ++i) {
        double x = std::get<0>(seeds[i]);
        double y = std::get<1>(seeds[i]);
        double z = std::get<2>(seeds[i]);
        int i_x = static_cast<int>((x - min_x) / container_size_x);
        int i_y = static_cast<int>((y - min_y) / container_size_y);
        int i_z = static_cast<int>((z - min_z) / container_size_z);
        if (i_x >= 0 && i_x < nx && i_y >= 0 && i_y < ny && i_z >= 0 && i_z < nz) {
            containers[i_x][i_y][i_z].push_back(i);
        }
    }
}

Mesh VoronoiBuilder::build() {
    Mesh global_mesh;
    global_mesh.SetDimensions(3);
    global_mesh.CreateTag("SEED_ID", DATA_INTEGER, CELL, NONE, 1);

    double min_x_sys = std::get<0>(system_size).first;
    double max_x_sys = std::get<0>(system_size).second;
    double min_y_sys = std::get<1>(system_size).first;
    double max_y_sys = std::get<1>(system_size).second;
    double min_z_sys = std::get<2>(system_size).first;
    double max_z_sys = std::get<2>(system_size).second;

    for (size_t seed_index = 0; seed_index < seeds.size(); ++seed_index) {
        const auto& current_seed = seeds[seed_index];
        double x = std::get<0>(current_seed);
        double y = std::get<1>(current_seed);
        double z = std::get<2>(current_seed);

        int ix = static_cast<int>((x - min_x_sys) / container_size_x);
        int iy = static_cast<int>((y - min_y_sys) / container_size_y);
        int iz = static_cast<int>((z - min_z_sys) / container_size_z);

        if (ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz) continue;

        // 1. Create initial lightweight Polyhedron
//        PolyhedronRaw current_poly = create_cube(min_x_sys, max_x_sys, min_y_sys, max_y_sys, min_z_sys, max_z_sys);

        // 2. Collect and sort neighbors
        std::vector<std::pair<double, int>> neighbors;
        for (int i = -1; i <= 1; ++i) {
            for (int j = -1; j <= 1; ++j) {
                for (int k = -1; k <= 1; ++k) {
                    int nix = ix + i, niy = iy + j, niz = iz + k;
                    if (nix >= 0 && nix < nx && niy >= 0 && niy < ny && niz >= 0 && niz < nz) {
                        for (int neighbor_seed_index : containers[nix][niy][niz]) {
                            if (neighbor_seed_index == seed_index) continue;
                            const auto& s = seeds[neighbor_seed_index];
                            double dist_sq = (std::get<0>(s) - x) * (std::get<0>(s) - x) +
                                             (std::get<1>(s) - y) * (std::get<1>(s) - y) +
                                             (std::get<2>(s) - z) * (std::get<2>(s) - z);
                            neighbors.push_back({dist_sq, neighbor_seed_index});
                        }
                    }
                }
            }
        }
        std::sort(neighbors.begin(), neighbors.end());
        std::vector<Vec3> npoints;
        for (auto [_, n] : neighbors) {
          npoints.push_back(Vec3(get<0>(seeds[n]), get<1>(seeds[n]), get<2>(seeds[n])));
        }

        GCell voro_cell = createVoronoiCell({x, y, z}, npoints);

        // 4. If the final polyhedron is valid, add it to the INMOST mesh
        if (!voro_cell.faces.empty()) {
//            std::cout << "--- Polyhedron for seed " << seed_index << " before converting to INMOST ---" << std::endl;
            add_polyhedron_to_inmost(&global_mesh, voro_cell);
//            add_polyhedron_to_mesh(&global_mesh, current_poly, seed_index);
        } else {
          std::cout << "Voronoi cell for seed " << seed_index << " was entirely cut away.\n";
        }
    }

    return global_mesh;
}