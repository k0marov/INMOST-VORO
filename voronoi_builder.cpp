#include "voronoi_builder.h"
#include "planar_cutter.h"
#include "helpers.h"
#include <cmath>
#include <iostream>

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

    double min_x_sys = std::get<0>(system_size).first;
    double min_y_sys = std::get<1>(system_size).first;
    double min_z_sys = std::get<2>(system_size).first;

    PlanarCutter cutter;

    for (size_t seed_index = 0; seed_index < seeds.size(); ++seed_index) {
        double x = std::get<0>(seeds[seed_index]);
        double y = std::get<1>(seeds[seed_index]);
        double z = std::get<2>(seeds[seed_index]);

        std::cout << "seed: " << x << ' ' << y << ' ' << z <<'\n';

        int ix = static_cast<int>((x - min_x_sys) / container_size_x);
        int iy = static_cast<int>((y - min_y_sys) / container_size_y);
        int iz = static_cast<int>((z - min_z_sys) / container_size_z);

        if (ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz) continue;

//        double center_x = min_x_sys + (ix + 0.5) * container_size_x;
//        double center_y = min_y_sys + (iy + 0.5) * container_size_y;
//        double center_z = min_z_sys + (iz + 0.5) * container_size_z;

        Cell current_voronoi_cell = CreateCubeAtPoint(&global_mesh, x, y, z, 2*container_size_x, 2*container_size_y, 2*container_size_z);
        bool cell_was_destroyed = false;

        // TODO: it's  a lot better to just construct the initial cell such that it does not overflow bounds
        current_voronoi_cell = cutter.Cut(current_voronoi_cell, 1, 0, 0, 1, true);
        current_voronoi_cell = cutter.Cut(current_voronoi_cell, -1, 0, 0, 0, true);

        current_voronoi_cell = cutter.Cut(current_voronoi_cell, 0, 1, 0, 1, true);
        current_voronoi_cell = cutter.Cut(current_voronoi_cell, 0, -1, 0, 0, true);

        current_voronoi_cell = cutter.Cut(current_voronoi_cell, 0, 0, 1, 1, true);
        current_voronoi_cell = cutter.Cut(current_voronoi_cell, 0, 0, -1, 0, true);

        for (int i = -1; i <= 1; ++i) {
            for (int j = -1; j <= 1; ++j) {
                for (int k = -1; k <= 1; ++k) {
                    int nix = ix + i;
                    int niy = iy + j;
                    int niz = iz + k;

                    if (nix >= 0 && nix < nx && niy >= 0 && niy < ny && niz >= 0 && niz < nz) {
                        for (int neighbor_seed_index : containers[nix][niy][niz]) {
                            if (neighbor_seed_index == seed_index) continue;

                            double nx_coord = std::get<0>(seeds[neighbor_seed_index]);
                            double ny_coord = std::get<1>(seeds[neighbor_seed_index]);
                            double nz_coord = std::get<2>(seeds[neighbor_seed_index]);

//                          std::cout << "neighbor: " << nx_coord << ' ' << ny_coord << ' ' << nz_coord <<'\n';

                            double a = 2 * (nx_coord - x);
                            double b = 2 * (ny_coord - y);
                            double c = 2 * (nz_coord - z);
                            double d = (nx_coord * nx_coord - x * x) +
                                       (ny_coord * ny_coord - y * y) +
                                       (nz_coord * nz_coord - z * z);
                            bool cut_positive =  a*x + b*y + c*z - d < 0;
                            current_voronoi_cell = cutter.Cut(current_voronoi_cell, a, b, c, d, cut_positive);

                            if (!current_voronoi_cell.isValid()) {
                                std::cerr << "Warning: Voronoi cell for seed " << seed_index << " was entirely cut away." << std::endl;
                                cell_was_destroyed = true;
                                break;
                            }
                        }
                    }
                    if(cell_was_destroyed) break;
                }
                if(cell_was_destroyed) break;
            }
            if(cell_was_destroyed) break;
        }
    }

    return global_mesh;
}
