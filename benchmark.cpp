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
#include <sstream>
#include <voro++.hh>

using namespace INMOST;
using FloatType = voronoi::FloatType;
using Vec3 = voronoi::Vec3;

static double run_voro_library(const std::vector<Vec3>& seeds, int target_per_cell) {
    auto t0 = std::chrono::steady_clock::now();
    int n = static_cast<int>(seeds.size());
    if (n == 0) return 0.0;

    double ratio = static_cast<double>(n) / static_cast<double>(target_per_cell);
    int cells_per_axis = std::cbrt(n/10); 

    voro::container con(
        0.0, 1.0,
        0.0, 1.0,
        0.0, 1.0,
        cells_per_axis, cells_per_axis, cells_per_axis,
        false, false, false,
        8);

    for (int i = 0; i < n; ++i) {
        con.put(i, seeds[static_cast<size_t>(i)].x, seeds[static_cast<size_t>(i)].y, seeds[static_cast<size_t>(i)].z);
    }

    voro::voronoicell c;
    voro::c_loop_all cl(con);

    double total_volume = 0.0;
    if (cl.start()) do {
        if (con.compute_cell(c, cl)) {
            // total_volume += c.volume();
        }
    } while (cl.inc());
    auto t1 = std::chrono::steady_clock::now();

    (void)total_volume;
    return std::chrono::duration<double>(t1 - t0).count();
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
    // N_values.push_back(2000000);

    const int target_per_cell = 5;

    for (int n : N_values) {
        std::cout << "Running N=" << n << "... " << std::flush;

        const double sys_length = 1.0;
        uint64_t seed_val = 1234;
        std::vector<Vec3> seeds = voronoi::generate_random_points_box(static_cast<size_t>(n), seed_val, sys_length);


        voronoi::VoronoiStats stats;
        auto t_compute_start = std::chrono::steady_clock::now();
        double voroqh_total_volume = 0.0;
        voronoi::for_each_polyhedron(seeds, target_per_cell, stats, [&](size_t idx, const voronoi::Polyhedron& poly) {
            (void)idx;
            double vol = 0; // voronoi::compute_polyhedron_volume(poly);
            voroqh_total_volume += vol;
        });
        auto t_compute_end = std::chrono::steady_clock::now();
        double voroqh_time_s = std::chrono::duration<double>(t_compute_end - t_compute_start).count();
        double voroqh_time_per_cell_us = (voroqh_time_s / n) * 1e6;

        double inmost_time_s = 0.0;

        double voro_time_s = run_voro_library(seeds, target_per_cell);
        csv << n << "," << voro_time_s << "," << voroqh_time_s << "," << voroqh_time_per_cell_us << "," << inmost_time_s << "\n";
        std::cout << "Done. voro++: " << voro_time_s << "s, voroqh: " << voroqh_time_s << "s\n";
    }
    
    std::cout << "Benchmark complete.\n";
    return 0;
}
