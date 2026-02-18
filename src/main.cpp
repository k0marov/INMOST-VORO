#include "inmost.h"
#include "voronoi_builder.h"
#include "voroqh/voronoi.hpp"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <tuple>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <string>


const int target_per_cell = 5;
void direct_vtk_pipeline(std::ifstream& infile, std::string output_filename) {
    std::vector<voronoi::Vec3> seeds; 
    std::string line;
    double x, y, z;
    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        if (ss >> x >> y >> z) {
        seeds.emplace_back(voronoi::Vec3{x, y, z});
        }
    }
    
    std::cout << "Read " << seeds.size() << " seeds" << std::endl;
    int n = seeds.size();

    voronoi::VoronoiStats build_stats;
    auto polys_out = voronoi::generate_voronoi_diagram(
        seeds,
        target_per_cell,
        build_stats
    );
    std::vector<int> cell_ids(n, 0);
    for (int i = 0; i < n; ++i) {
        cell_ids[i] = i; 
    }
    
    const auto t_build_start = std::chrono::steady_clock::now();
    voronoi::write_polyhedra_vtk(output_filename, polys_out, cell_ids);
    const auto t_build_end = std::chrono::steady_clock::now();
    auto vtk_write_time = std::chrono::duration<double, std::milli>(t_build_end - t_build_start).count();
    std::cout << "Voronoi diagram generation via voroqh took " << build_stats.total_ms << " ms" << std::endl;
    std::cout << "Writing to VTK file took " << vtk_write_time << " ms" << std::endl;
    std::cout << "Saved Voronoi tessellation to " << output_filename << std::endl;
}

void inmost_pipeline(std::ifstream& infile, std::string output_filename) {
    std::vector<std::tuple<double, double, double>> seeds; 
    std::string line;
    double x, y, z;
    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        if (ss >> x >> y >> z) {
        seeds.emplace_back(x, y, z);
        }
    }
    
    std::cout << "Read " << seeds.size() << " seeds" << std::endl;
    int n = seeds.size();

    VoronoiBuilder builder(seeds, target_per_cell);
    voronoi::VoronoiStats build_stats;
    INMOST::Mesh voronoi_mesh = builder.build(&build_stats);
    std::cout << "Voronoi diagram generation via voroqh took " << build_stats.total_ms << " ms" << std::endl;
    std::cout << "INMOST::CreateCell() took " << build_stats.time_ms_callback << " ms" << std::endl;
    const auto t_build_start = std::chrono::steady_clock::now();
    voronoi_mesh.Save(output_filename);
    const auto t_build_end = std::chrono::steady_clock::now();
    auto vtk_write_time = std::chrono::duration<double, std::milli>(t_build_end - t_build_start).count();
    std::cout << "Saving to VTK file took " << vtk_write_time + build_stats.time_ms_callback << " ms" << std::endl;
    std::cout << "Saved Voronoi tessellation to " << output_filename << std::endl;
}



int main(int argc, char ** argv)
{
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <path_to_xyz_file> [--direct-vtk]" << std::endl;
        return 1;
    }
    bool direct_vtk = false; 
    if (argc == 3 && std::string(argv[2]) == "--direct-vtk") {
        direct_vtk = true;
    }

    std::string input_path = argv[1];
    std::string output_filename = input_path + ".vtk"; 
    std::vector<std::tuple<double, double, double>> seeds;
    
    std::ifstream infile(input_path);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open file " << input_path << std::endl;
        return 1;
    }

    // voroqh_pipeline(infile, output_filename);
    if (direct_vtk) {
        direct_vtk_pipeline(infile, output_filename);
    } else {
        inmost_pipeline(infile, output_filename);
    }
    return 0;
}
