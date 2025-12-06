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
3. Open output file `voronoi_output.vtk` in ParaView to view the tessellation and use it
