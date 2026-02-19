# voroqh

voroqh is a small parallelizable 3D Voronoi diagram generator designed as an alternative to voro++ for INMOST-VORO project. It builds Voronoi cells in a unit box using a custom pipeline based on a modified QuickHull implementation for convex hulls.

The focus is on:
- Being O(N) for uniform seed distributions 
- Producing correct results for non-uniform seed distributions 
- Computing each Voronoi Cell independently, so that it can be parallelized
- Detailed timing and statistics for algorithm profiling
- Simple, dependency-free C++17 code that can be dropped into other projects

The `quickhull` subdirectory contains a specific QuickHull variant tuned for this project.

## QuickHull variant

The QuickHull code under `quickhull/` is an adapted fork of an existing public-domain implementation from [https://github.com/akuukka/quickhull](https://github.com/akuukka/quickhull). In this project it is used in a very specific way:
- Microsecond optimized for tiny inputs, typically number of generators per cell neighborhood `N < 100` (because Voronoi cells in 3D have ~16 faces and thanks to an algorithm of dynamic neighbor selection on average we only need to build convex hulls for ~80 neighbors)
- Tight integration with Voronoi-specific statistics and timing
- Implementation of tracking of point IDs to efficiently determine Voronoi cell neighbors without extra computation. (TODO) 

The original QuickHull implementation itself is generic and can handle larger inputs, but the way voroqh drives it is deliberately tailored to small local hulls. But this modified version still works correctly on big inputs. 

## Layout

- `voronoi.hpp`, `voronoi.cpp`: core Voronoi types and algorithms
- `main.cpp`: command-line entry point for running a Voronoi simulation and exporting VTK
- `tests.cpp`: standalone test runner for Voronoi correctness and robustness
- `quickhull/`: embedded QuickHull library and its own tests

## Building the main binary

From the `voroqh` directory:

```bash
cd /Users/roshi/Desktop/inmost-voronoi/code/INMOST-VORO/voroqh

g++ -std=c++17 -O3 \
    main.cpp \
    voronoi.cpp \
    -o voroqh
```

This produces a `voroqh` executable in the same directory.

### Running the main binary

The main program simulates a Voronoi tessellation in the unit box and optionally writes VTK output (that can be viewed in programs like Paraview).

```bash
./voroqh [n] [seed] [target_per_cell] [out_path] [volume_flag]
```

Where:
- `n`: number of seeds (default `200`)
- `seed`: RNG seed (default `1`)
- `target_per_cell`: hyperparameter for the spacial grid that may influence performance: approximate target number of neighboring seeds per cell expansion (default `5`)
- `out_path`: output path for VTK (`polyhedra.vtk` by default). If empty or not ending in `.vtk`, no VTK file is written.
- `volume_flag`: `0` or `1` to disable/enable volume computations and reporting (default `0`)

Example:

```bash
./voroqh 200 1 5 polyhedra.vtk 1
```

This runs a 200-seed simulation, writes `polyhedra.vtk`, and prints detailed statistics and timings to stdout.

## Building and running Voronoi tests

The Voronoi tests live in `tests.cpp` and depend only on `voronoi.cpp` and the standard library.

From the `voroqh` directory:

```bash
cd /Users/roshi/Desktop/inmost-voronoi/code/INMOST-VORO/voroqh

g++ -std=c++17 -O3 \
    tests.cpp \
    voronoi.cpp \
    -o voroqh_tests

./voroqh_tests
```

The test binary prints a series of numbered tests, for example:

- Topology and Euler characteristic checks on many cells
- Bounding-box and basic invariants
- Surface area consistency checks
- Translation invariance of volume and centroid
- Two-seed split sanity check
- Degenerate clustered seeds with positive volume
- Volume summation stability and statistics
- Half-space checks for Voronoi cell vertices against cutting planes
- VTK writer smoke test

All tests are assertion-based; if any assertion fails, the process terminates.

## Building and running QuickHull tests

The embedded QuickHull implementation has its own test suite under `quickhull/TestsSrc`, with a CMake configuration.

From the `quickhull/TestsSrc` directory:

```bash
cd /Users/roshi/Desktop/inmost-voronoi/code/INMOST-VORO/voroqh/quickhull/TestsSrc

cmake -S . -B build
cmake --build build

./build/QuickHullTests
```