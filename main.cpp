#include "inmost.h"
#include "voronoi_builder.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <tuple>

using std::get;

using namespace INMOST;

int main(int argc, char ** argv)
{
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <path_to_xyz_file>" << std::endl;
        return 1;
    }

    std::string filepath = argv[1];
    std::vector<std::tuple<double, double, double>> seeds;
    
    std::ifstream infile(filepath);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open file " << filepath << std::endl;
        return 1;
    }

    std::string line;
    double x, y, z;
    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        if (ss >> x >> y >> z) {
            seeds.emplace_back(x, y, z);
        }
    }
    
    std::cout << "Read " << seeds.size() << " seeds from " << filepath << std::endl;

    // Define a 1x1x1 system size
    SystemSize system_size = {
        std::make_pair(0.0, 1.0),
        std::make_pair(0.0, 1.0),
        std::make_pair(0.0, 1.0)
    };

    // Build the Voronoi tessellation
    VoronoiBuilder builder(seeds, system_size);
    Mesh voronoi_mesh = builder.build();

//    for (auto s : seeds) {
//      CreateCubeAtPoint(&voronoi_mesh, get<0>(s), get<1>(s), get<2>(s), 0.1, 0.1, 0.1);
////      voronoi_mesh.CreateNode({get<0>(s), get<1>(s), get<2>(s)})
//    }


    // Save the result to a VTK file
    std::string output_filename = "voronoi_output.vtk";
    voronoi_mesh.Save(output_filename);
    std::cout << "Saved Voronoi tessellation to " << output_filename << std::endl;

    return 0;
}