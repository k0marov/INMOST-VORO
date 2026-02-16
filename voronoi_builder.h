#ifndef CODE_VORONOI_BUILDER_H
#define CODE_VORONOI_BUILDER_H

#include <vector>
#include <tuple>
#include "inmost.h"

using namespace INMOST;

// Represents the system size as a tuple of pairs for x, y, and z dimensions.
using SystemSize = std::tuple<std::pair<double, double>, std::pair<double, double>, std::pair<double, double>>;

class VoronoiBuilder {
public:
    // Constructor to initialize the Voronoi builder with seed points and system size.
    VoronoiBuilder(const std::vector<std::tuple<double, double, double>>& seeds, SystemSize system_size, int target_per_cell = 25);

    // Builds the Voronoi tessellation and returns it as an INMOST mesh.
    Mesh build();

private:
    // Computes a single Voronoi cell for a given seed.
    void compute_voronoi_cell(int seed_index, Mesh& mesh);

    // Vector of seed points.
    std::vector<std::tuple<double, double, double>> seeds;
    // System dimensions.
    SystemSize system_size;
    int target_per_cell;
};

#endif // CODE_VORONOI_BUILDER_H