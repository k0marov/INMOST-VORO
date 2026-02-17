# **INMOST-VORO**

> [!WARNING]
> THIS IS WORK IN PROGRESS


![visualization.png](docs/imgs/visualization.png)

**A PDF overview of the project can be found [here](https://skomarov.com/static/distributed_voronoi_remeshing.pdf).**


* **VoroLib** — a library with Voronoi-related utilities
* **Voro** — a command-line executable using the library

---

## **Requirements**

* **C++ compiler** (GCC, Clang, or MSVC)
* **CMake ≥ 3.5**
* **INMOST library** installed on your system
  (and discoverable via `find_package(inmost REQUIRED)`)

---

## **Build Instructions**

### **0. Install INMOST**

Compile INMOST by using the guide from here: [https://github.com/INMOST-DEV/INMOST/wiki/0100-Compilation](https://github.com/INMOST-DEV/INMOST/wiki/0100-Compilation)

### **1. Clone the project**

```bash
git clone https://github.com/yourusername/InmostVoro.git
cd InmostVoro
```

### **2. Configure with CMake**

```bash
cmake -B build
```

If INMOST is installed in a custom location, set:

```bash
cmake -B build -DINMOST_DIR=/path/to/inmost
```

### **3. Build**

```bash
cmake --build build
```

This produces:

* `build/Voro` (the executable)
* `build/Benchmark` (the benchmarking suite)
* `build/libVoroLib.a` or `libVoroLib.so` (depending on platform)

---

## **Installation (Optional)**

To install the executable, library, and headers:

```bash
cmake --install build --prefix /your/install/prefix
```

Default install locations:

```
bin/          → Voro executable
lib/          → VoroLib library
include/      → public headers (voronoi_builder.h)
lib/cmake/InmostVoro/ → CMake package config (exported targets)
```

---

## **Example usage**

1. Use a helper script to generate a cloud of random seed points in with coords between 0 and 1:
```bash
python scripts/generate_random_points.py 1000 random_points_1000.txt
```
2. Run Voro to get a Voronoi tessellation with these points 
```bash
./Voro random_points_1000.txt
```
3.93. Open output file `voronoi_output.vtk` in ParaView to view the tessellation and use it

---

## **Benchmarks**

A dedicated benchmark suite is included to compare the performance of:
1. **voro++**: The external command-line tool.
2. **voroqh**: The internal engine used by this project.
3. **INMOST Integration**: The full pipeline including INMOST cell creation.

### **Running the Benchmark**

After building, run:
```bash
./build/Benchmark
```

This will:
- Iterate through various input sizes (N=100 to N=2,000,000).
- Measure execution times for each component.
- Output results to `benchmark_results.csv`.

**Note**: The benchmark simulates real-world conditions by including file I/O overhead for all methods to ensure a fair comparison.

---

## **Migration to Voroqh**

The project has been migrated from the original geometrical engine to **voroqh**. Key changes include:
*   **Engine Replacement**: Replaced the legacy Voronoi engine with `voroqh` for efficient parallel generation.
*   **Optimization**: Compiled with `-O3 -march=native` for maximum performance.
*   **Simplified Pipeline**: Removed `PlanarCutter`. Boundary clipping is now handled natively by `voroqh`, significantly reducing overhead (~2x faster conversion).
*   **INMOST Conversion**: Implemented direct conversion from `voroqh::Polyhedron` to `INMOST::Cell`.
