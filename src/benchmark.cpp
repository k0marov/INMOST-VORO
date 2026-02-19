#include "inmost.h"
#include "voronoi_builder.h"
#include "voroqh/voronoi.hpp"
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

static std::string format_topology_line(const voronoi::TopologyOutput& t) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(5);
    oss << t.index << " " << t.volume << " " << t.num_faces;
    std::vector<std::pair<double,int>> faces;
    faces.reserve(t.num_faces);
    for (int i = 0; i < t.num_faces; ++i) {
        faces.emplace_back(t.face_areas[static_cast<size_t>(i)],
                           t.face_edge_counts[static_cast<size_t>(i)]);
    }
    std::sort(faces.begin(), faces.end(),
              [](const std::pair<double,int>& a, const std::pair<double,int>& b) {
                  return a.first < b.first;
              });
    for (const auto& f : faces) oss << " " << f.first;
    for (const auto& f : faces) oss << " " << f.second;
    return oss.str();
}

static bool validate_diagram(const std::vector<voronoi::TopologyOutput>& voropp_out,
                             const std::vector<voronoi::TopologyOutput>& voroqh_out,
                             int n) {
    if (static_cast<int>(voropp_out.size()) != n) return false;
    if (static_cast<int>(voroqh_out.size()) != n) return false;
    std::vector<const voronoi::TopologyOutput*> a(n, nullptr);
    std::vector<const voronoi::TopologyOutput*> b(n, nullptr);
    for (const auto& t : voropp_out) {
        if (t.index < 0 || t.index >= n) return false;
        a[static_cast<size_t>(t.index)] = &t;
    }
    for (const auto& t : voroqh_out) {
        if (t.index < 0 || t.index >= n) return false;
        b[static_cast<size_t>(t.index)] = &t;
    }
    const double volume_eps = 1e-8; // absolute error
    const double area_eps = 1e-8;
    std::vector<double> vols_a;
    std::vector<double> vols_b;
    vols_a.reserve(static_cast<size_t>(n));
    vols_b.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        const auto& va = *a[static_cast<size_t>(i)];
        const auto& vb = *b[static_cast<size_t>(i)];
        vols_a.push_back(va.volume);
        vols_b.push_back(vb.volume);
        // double denom = std::max({std::abs(va.volume), std::abs(vb.volume), 1e-14});
        double rel_err = std::abs(va.volume - vb.volume);// * n; // average cell volume ~ 1/N
        if (rel_err > volume_eps) {
            std::cerr << "Volume mismatch at cell " << i << "\n";
            std::cerr << "Absolute error is " << rel_err << '\n';
            std::cerr << "  voro++: " << format_topology_line(va) << "\n";
            std::cerr << "  voroqh: " << format_topology_line(vb) << "\n";
            return false;
        }
        double dx = va.centroid_x - vb.centroid_x;
        double dy = va.centroid_y - vb.centroid_y;
        double dz = va.centroid_z - vb.centroid_z;
        double dist = std::abs(dx) + std::abs(dy) + std::abs(dz); 
        // const double centroid_eps = 5e-6;
        // double D = std::max(va.diameter, vb.diameter);
        double threshold = 1e-5;
        if (dist > threshold) {
            std::cerr << "Centroid mismatch at cell " << i << "\n";
            std::cerr << "  distance = " << dist << ", threshold = " << threshold << "\n";
            std::cerr << "  voro++ centroid: " << va.centroid_x << " " << va.centroid_y << " " << va.centroid_z << "\n";
            std::cerr << "  voroqh centroid: " << vb.centroid_x << " " << vb.centroid_y << " " << vb.centroid_z << "\n";
            std::cerr << "  voro++: " << format_topology_line(va) << "\n";
            std::cerr << "  voroqh: " << format_topology_line(vb) << "\n";
            return false;
        }
    }
    std::cout << "Consistency checks between voro++ and voroqh passed for N = " << n << '\n'; 
    auto sum_sorted = [](const std::vector<double>& v) {
        std::vector<double> tmp = v;
        std::sort(tmp.begin(), tmp.end(), [](double x, double y) {
            return std::abs(x) < std::abs(y);
        });
        long double s = 0.0L;
        for (double x : tmp) s += static_cast<long double>(x);
        return static_cast<double>(s);
    };
    double sum_vol_a = sum_sorted(vols_a);
    double sum_vol_b = sum_sorted(vols_b);
    double diff_a = sum_vol_a - 1.0;
    double diff_b = sum_vol_b - 1.0;
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(5);
    oss << "Total volume diff from expected volume for unit cube (1.0):\nvoro++ = " << diff_a
        << ", voroqh = " << diff_b << "\n";
    std::cout << oss.str();
    return true;
}

static std::vector<voronoi::TopologyOutput> run_voro_library(const std::vector<Vec3>& seeds,
                                                             double& compute_time_s) {
    auto t0 = std::chrono::steady_clock::now();
    int n = static_cast<int>(seeds.size());
    if (n == 0) {
        compute_time_s = 0.0;
        return {};
    }

    double ratio = static_cast<double>(n) / 10.0;
    int cells_per_axis = std::max(1, static_cast<int>(std::cbrt(ratio)));
    // cells_per_axis = 6; 

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

    if (cl.start()) do {
        (void)con.compute_cell(c, cl);
    } while (cl.inc());
    auto t1 = std::chrono::steady_clock::now();
    compute_time_s = std::chrono::duration<double>(t1 - t0).count();

    std::vector<voronoi::TopologyOutput> out;
    out.reserve(static_cast<size_t>(n));
    voro::c_loop_all cl2(con);
    if (cl2.start()) do {
        if (con.compute_cell(c, cl2)) {
            int pid = cl2.pid();
            voronoi::TopologyOutput topo;
            topo.index = pid;
            topo.volume = c.volume();
            topo.num_faces = c.number_of_faces();
            double cx=0, cy=0, cz=0;
            double px=0, py=0, pz=0;
            cl2.pos(px, py, pz);
            c.centroid(cx, cy, cz);
            topo.centroid_x = cx + px;
            topo.centroid_y = cy + py;
            topo.centroid_z = cz + pz;
            std::vector<double> verts;
            c.vertices(verts);
            double D = 0.0;
            int nv = verts.size() / 3;

            for (int i = 0; i < nv; ++i) {
                for (int j = i+1; j < nv; ++j) {
                    double dx = verts[3*i]   - verts[3*j];
                    double dy = verts[3*i+1] - verts[3*j+1];
                    double dz = verts[3*i+2] - verts[3*j+2];
                    double d2 = dx*dx + dy*dy + dz*dz;
                    if (d2 > D) D = d2;
                }
            }
            topo.diameter = std::sqrt(D);
            std::vector<int> fo;
            std::vector<double> fa;
            c.face_orders(fo);
            c.face_areas(fa);
            topo.face_edge_counts = std::move(fo);
            topo.face_areas = std::move(fa);
            out.push_back(std::move(topo));
        }
    } while (cl2.inc());
    return out;
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
    for (int n = 1000000; n <= 3000000; n += 1000000) N_values.push_back(n);

    const int target_per_cell = 5;

    for (int n : N_values) {
        std::cout << "Running N=" << n << "...\n" << std::flush;

        const double sys_length = 1.0;
        uint64_t seed_val = 1234;
        std::vector<Vec3> seeds = voronoi::generate_random_points_box(static_cast<size_t>(n), seed_val, sys_length);

        voronoi::VoronoiStats stats;
        std::vector<voronoi::TopologyOutput> voroqh_topos;
        voroqh_topos.reserve(static_cast<size_t>(n));
        std::vector<double> area_buffer;
        double voroqh_topology_time_s = 0.0;
        auto t_compute_start = std::chrono::steady_clock::now();
        voronoi::for_each_polyhedron(seeds, target_per_cell, stats, [&](size_t idx, const voronoi::Polyhedron& poly) {
            auto topo_t0 = std::chrono::steady_clock::now();
            double vol = voronoi::compute_polyhedron_volume(poly);
            std::ostringstream null_out;
            auto topo = voronoi::write_voro_compatible_output(null_out, idx, seeds[idx], poly, vol, area_buffer);
            auto topo_t1 = std::chrono::steady_clock::now();
            voroqh_topology_time_s += std::chrono::duration<double>(topo_t1 - topo_t0).count();
            voroqh_topos.push_back(std::move(topo));
        });
        auto t_compute_end = std::chrono::steady_clock::now();
        double voroqh_total_time_s = std::chrono::duration<double>(t_compute_end - t_compute_start).count();
        double voroqh_time_s = voroqh_total_time_s - voroqh_topology_time_s;
        double voroqh_time_per_cell_us = (voroqh_time_s / n) * 1e6;

        double inmost_time_s = 0.0;

        double voro_time_s = 0.0;
        auto voro_topos = run_voro_library(seeds, voro_time_s);

        if (!validate_diagram(voro_topos, voroqh_topos, n)) {
            std::cerr << "Mismatch detected for N=" << n << "\n";
            return 1;
        }

        csv << n << "," << voro_time_s << "," << voroqh_time_s << "," << voroqh_time_per_cell_us << "," << inmost_time_s << "\n" << std::flush;
        std::cout << "Done. voro++: " << voro_time_s << "s, voroqh: " << voroqh_time_s << "s\n\n";
    }
    
    std::cout << "Benchmark complete.\n";
    return 0;
}
