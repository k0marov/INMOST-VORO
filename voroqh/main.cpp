#include "voronoi.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using FloatType = voronoi::FloatType;
using Vec3 = voronoi::Vec3;

struct SimulationConfig {
    int n;
    uint64_t seed_val;
    int target_per_cell;
    std::string out_path;
    int volume_flag;
    int voro_flag;
};

static bool ends_with(const std::string& s, const std::string& suffix) {
    if (s.size() < suffix.size()) return false;
    return std::equal(suffix.rbegin(), suffix.rend(), s.rbegin());
}

static bool run_voro_benchmark(const std::vector<Vec3>& seeds, FloatType& out_ms) {
    int available = std::system("command -v voro++ > /dev/null 2>&1");
    if (available != 0) return false;

    std::string path = "voro_input.txt";
    {
        std::ofstream out(path);
        out << std::setprecision(17);
        for (size_t i = 0; i < seeds.size(); ++i) {
            out << i << " " << seeds[i].x << " " << seeds[i].y << " " << seeds[i].z << "\n";
        }
    }

    auto t0 = std::chrono::steady_clock::now();
    std::string cmd = "voro++ -c \"%i %q %v %f %v\" 0 1 0 1 0 1 " + path + " > /dev/null 2>&1";
    int code = std::system(cmd.c_str());
    auto t1 = std::chrono::steady_clock::now();

    std::remove(path.c_str());
    std::remove((path + ".vol").c_str());

    if (code != 0) return false;
    out_ms = std::chrono::duration<FloatType, std::milli>(t1 - t0).count();
    return true;
}

static void run_simulation(const SimulationConfig& config) {
    const FloatType sys_length = 1.0;
    std::vector<Vec3> seeds = voronoi::generate_random_points_box(static_cast<size_t>(config.n), config.seed_val, sys_length);

    voronoi::VoronoiStats stats;
    double total_volume = 0.0;
    double time_ms_volume = 0.0;
    uint64_t faces_total = 0;
    uint64_t faces_max = 0;

    const bool write_vtk = (!config.out_path.empty() && ends_with(config.out_path, ".vtk"));
    std::vector<voronoi::Polyhedron> polys_out;
    std::vector<int> cell_ids;
    if (write_vtk) {
        polys_out.reserve(static_cast<size_t>(config.n));
        cell_ids.reserve(static_cast<size_t>(config.n));
    }

    voronoi::for_each_polyhedron(
        seeds,
        config.target_per_cell,
        stats,
        [&](size_t seed_index, const voronoi::Polyhedron& poly) {
            const uint64_t faces = static_cast<uint64_t>(poly.face_degree.size());
            faces_total += faces;
            faces_max = std::max(faces_max, faces);

            if (config.volume_flag) {
                auto tv0 = std::chrono::steady_clock::now();
                total_volume += voronoi::compute_polyhedron_volume(poly);
                auto tv1 = std::chrono::steady_clock::now();
                time_ms_volume += std::chrono::duration<double, std::milli>(tv1 - tv0).count();
            }

            if (write_vtk) {
                polys_out.push_back(poly);
                cell_ids.push_back(static_cast<int>(seed_index));
            }
        });

    if (write_vtk) {
        voronoi::write_polyhedra_vtk(config.out_path, polys_out, cell_ids);
    }

    FloatType voro_ms = 0.0;
    bool voro_ok = false;
    if (config.voro_flag != 0) {
        voro_ok = run_voro_benchmark(seeds, voro_ms);
    }

    std::cout << std::setprecision(15);
    std::cout << "n=" << config.n << " seed=" << config.seed_val << " target_per_cell=" << config.target_per_cell << "\n";
    if (config.volume_flag) {
        std::cout << "sum_volume=" << total_volume << " error=" << std::abs(total_volume - sys_length * sys_length * sys_length) << "\n";
    }

    std::cout << "neighbors_total=" << stats.neighbors_total << " neighbors_avg=" << (static_cast<double>(stats.neighbors_total) / config.n) << "\n";
    std::cout << "faces_total=" << faces_total << " faces_mean=" << (static_cast<double>(faces_total) / config.n) << " faces_max=" << faces_max << "\n";
    std::cout << "expansions_total=" << stats.total_expansions << " expansions_avg=" << (static_cast<double>(stats.total_expansions) / config.n) << "\n";
    std::cout << "dual_points_total=" << stats.total_dual_points << " dual_points_avg=" << (static_cast<double>(stats.total_dual_points) / stats.total_quickhull_calls) << "\n";
    for (size_t i = 0; i < stats.step_count.size(); ++i) {
        if (stats.step_count[i] > 0) {
            std::cout << "step " << i << ": avg_neighbors=" << (static_cast<double>(stats.neighbors_count_per_step[i]) / stats.step_count[i])
                      << " (count=" << stats.step_count[i] << ")\n";
        }
    }
    std::cout << "quickhull_calls_total=" << stats.total_quickhull_calls << " quickhull_calls_avg=" << (static_cast<double>(stats.total_quickhull_calls) / config.n) << "\n";
    std::cout << "time_ms_total=" << stats.total_ms << " time_ms_build=" << stats.time_ms_build << " time_ms_cells=" << stats.time_ms_cells
              << " time_us_per_cell=" << (stats.time_ms_cells * 1000.0 / config.n) << "\n";
    std::cout << "time_ms_quickhull=" << stats.time_ms_quickhull << " (" << (stats.time_ms_quickhull / stats.time_ms_cells * 100.0) << "% of cells time)\n";
    std::cout << "  [Detailed Breakdown]\n";
    std::cout << "  neighbors_search=" << stats.time_ms_neighbors << " ms (" << (stats.time_ms_neighbors / stats.time_ms_cells * 100.0) << "%)\n";
    std::cout << "  dual_setup=" << stats.time_ms_dual_setup << " ms (" << (stats.time_ms_dual_setup / stats.time_ms_cells * 100.0) << "%)\n";
    std::cout << "  primal_recovery=" << stats.time_ms_primal_recovery << " ms (" << (stats.time_ms_primal_recovery / stats.time_ms_cells * 100.0) << "%)\n";
    std::cout << "  overhead_other=" << stats.time_ms_overhead_other << " ms (" << (stats.time_ms_overhead_other / stats.time_ms_cells * 100.0) << "%)\n";
    std::cout << "  [QuickHull Internal]\n";
    std::cout << "  qh_total=" << stats.time_ms_qh_total_internal << " ms\n";
    std::cout << "  qh_setup=" << stats.time_ms_qh_setup_internal << " ms\n";
    std::cout << "  qh_process=" << stats.time_ms_qh_process_internal << " ms\n";
    std::cout << "  qh_mesh_convert=" << stats.time_ms_qh_mesh_convert << " ms\n";
    std::cout << "  [Polyhedra Conversion]\n";
    std::cout << "  total_time=" << stats.time_ms_polyhedra << " ms (" << (stats.time_ms_polyhedra / stats.time_ms_cells * 100.0) << "% of cells time)\n";
    if (config.volume_flag) {
        std::cout << "  [Volume]\n";
        std::cout << "  total_time=" << time_ms_volume << " ms\n";
    }

    if (config.voro_flag == 0) {
        std::cout << "voro_ms=off\n";
    } else if (voro_ok) {
        std::cout << "voro_ms=" << voro_ms << "\n";
    } else {
        std::cout << "voro_ms=na\n";
    }
}

int main(int argc, char** argv) {
    SimulationConfig config;
    config.n = 200;
    config.seed_val = 1;
    config.target_per_cell = 4;
    config.out_path = "polyhedra.vtk";
    config.volume_flag = 0;
    config.voro_flag = 1;

    if (argc > 1) config.n = std::max(1, std::atoi(argv[1]));
    if (argc > 2) config.seed_val = static_cast<uint64_t>(std::strtoull(argv[2], nullptr, 10));
    if (argc > 3) config.target_per_cell = std::max(1, std::atoi(argv[3]));
    if (argc > 4) config.out_path = argv[4];
    if (argc > 5) config.volume_flag = std::atoi(argv[5]) != 0 ? 1 : 0;
    if (argc > 6) config.voro_flag = std::atoi(argv[6]) != 0 ? 1 : 0;

    run_simulation(config);
    return 0;
}