#ifndef CODE_VORONOI_BUILDER_H
#define CODE_VORONOI_BUILDER_H

#include <vector>
#include <tuple>
#include "inmost.h"

using namespace INMOST;

using SystemSize = std::tuple<std::pair<double, double>, std::pair<double, double>, std::pair<double, double>>;

class VoronoiBuilder {
public:
    VoronoiBuilder(const std::vector<std::tuple<double, double, double>>& seeds, int target_per_cell = 5);

    Mesh build();

private:
    void compute_voronoi_cell(int seed_index, Mesh& mesh);

    std::vector<std::tuple<double, double, double>> seeds;
    int target_per_cell;
};

#endif 