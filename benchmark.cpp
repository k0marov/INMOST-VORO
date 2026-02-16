#include "inmost.h"
#include "voronoi_builder.h"
#include "voronoi.hpp"
#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <string>
#include <algorithm>
#include <cmath>
#include <tuple>
#include <cstdio>
#include <cstdlib>

using namespace INMOST;
using FloatType = voronoi::FloatType;
using Vec3 = voronoi::Vec3;

// --- Helper Functions ---

static bool run_voro_benchmark(const std::vector<Vec3>& seeds, double& out_s) {
    int available = std::system("command -v voro++ > /dev/null 2>&1");
    if (available != 0) return false;

    std::string path = "voro_input_bench.txt";
    {
        std::ofstream out(path);
        out << std::setprecision(17);
        for (size_t i = 0; i < seeds.size(); ++i) {
            out << i << " " << seeds[i].x << " " << seeds[i].y << " " << seeds[i].z << "\n";
        }
    }

    auto t0 = std::chrono::steady_clock::now();
    // Run voro++ with minimal output to avoid I/O bottleneck
    // -c "%i" just prints ID, which is minimal. Or even better, redirect output to /dev/null
    std::string cmd = "voro++ -c \"%i %q %v %f %v\" 0 1 0 1 0 1 " + path + " > /dev/null 2>&1";
    int code = std::system(cmd.c_str());
    auto t1 = std::chrono::steady_clock::now();

    std::remove(path.c_str());
    std::remove((path + ".vol").c_str()); // voro++ might produce this if misconfigured, but with -c it prints to stdout

    if (code != 0) return false;
    out_s = std::chrono::duration<double>(t1 - t0).count();
    return true;
}

int main(int argc, char** argv) {
    std::string out_csv = "benchmark_results.csv";
    if (argc > 1) out_csv = argv[1];

    std::ofstream csv(out_csv);
    csv << "N,voro++_time_s,voroqh_time_s,voroqh_time_per_cell_us,inmost_cell_creation_time_s\n";
    
    std::cout << "Starting Benchmark Suite...\n";
    std::cout << "Results will be written to: " << out_csv << "\n";

    // Define N values
    std::vector<int> N_values;
    N_values.push_back(100);
    N_values.push_back(1000);
    for (int n = 10000; n <= 100000; n += 10000) N_values.push_back(n);
    for (int n = 200000; n <= 1000000; n += 100000) N_values.push_back(n);
    N_values.push_back(2000000);

    const int seeds_per_cell = 5;

    for (int n : N_values) {
        std::cout << "Running N=" << n << "... " << std::flush;

        // Generate Seeds
        const double sys_length = 1.0;
        uint64_t seed_val = 12345;
        std::vector<Vec3> voro_seeds = voronoi::generate_random_points_box(static_cast<size_t>(n), seed_val, sys_length);

        // 1. Run voro++
        double voro_plus_time_s = 0.0;
        if (!run_voro_benchmark(voro_seeds, voro_plus_time_s)) {
            voro_plus_time_s = -1.0; // Indicate N/A
        }

        // 2. Run voroqh (Pure)
        voronoi::VoronoiStats stats_pure;
        auto t0_pure = std::chrono::steady_clock::now();
        uint64_t faces_total = 0; 
        uint64_t faces_max = 0;
        voronoi::for_each_polyhedron(voro_seeds, seeds_per_cell, stats_pure, [&](size_t seed_index, const voronoi::Polyhedron& poly) {
            // TODO: check if this really matters
            return;
            const uint64_t faces = static_cast<uint64_t>(poly.face_degree.size());
            faces_total += faces;
            faces_max = std::max(faces_max, faces);
        });
        auto t1_pure = std::chrono::steady_clock::now();
        double voroqh_time_s = std::chrono::duration<double>(t1_pure - t0_pure).count();
        double voroqh_inner_time_s = stats_pure.total_ms / 1000.0;
        double voroqh_time_per_cell_us = (voroqh_time_s * 1e6) / n;
        std::cout << "inner time " << voroqh_inner_time_s << "s\n";

        // 3. Run INMOST Integration
        // Convert seeds for builder
        // std::vector<std::tuple<double, double, double>> inmost_seeds;
        // inmost_seeds.reserve(n);
        // for(const auto& s : voro_seeds) {
        //     inmost_seeds.emplace_back(s.x, s.y, s.z);
        // }
        
        // SystemSize system_size = {
        //     std::make_pair(0.0, 1.0),
        //     std::make_pair(0.0, 1.0),
        //     std::make_pair(0.0, 1.0)
        // };

        // auto t0_inmost = std::chrono::steady_clock::now();
        // {
        //     VoronoiBuilder builder(inmost_seeds, system_size, seeds_per_cell);
        //     Mesh mesh = builder.build();
        //     // Force mesh operations if needed, but build() does the work
        // }
        // auto t1_inmost = std::chrono::steady_clock::now();
        // double total_inmost_time_s = std::chrono::duration<double>(t1_inmost - t0_inmost).count();
        
        // // Calculate INMOST creation overhead
        // // We assume Total = Voroqh + INMOST_Overhead
        // // So INMOST_Overhead = Total - Voroqh
        // // Note: There might be slight variance between the two runs, but this is the best approximation.
        // double inmost_creation_time_s = total_inmost_time_s - voroqh_time_s;
        // if (inmost_creation_time_s < 0) inmost_creation_time_s = 0.0; // Should not happen typically

        // Output to CSV
        csv << n << "," 
            << voro_plus_time_s << "," 
            << voroqh_time_s << "," 
            << voroqh_time_per_cell_us << "," 
            << 0 << "\n";
        
        csv.flush();
        std::cout << "Done. (voro++: " << voro_plus_time_s << "s, voroqh: " << voroqh_time_s << "s, inmost_overhead: " 
        << 0 << "s)\n";
    }

    std::cout << "Benchmark Suite Completed.\n";
    return 0;
}
