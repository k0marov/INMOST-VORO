#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <functional>

const double GEOMETRY_EPSILON = 1e-9;

// --- Basic Geometric Data Structures ---

struct Point3D {
    double x = 0.0, y = 0.0, z = 0.0;
};

// A polygon is represented by an ordered list of vertices
struct Polygon {
    std::vector<Point3D> vertices;
};

// A polyhedron is a collection of faces (polygons)
struct Polyhedron {
    std::vector<Polygon> faces;
};

// A plane is defined by the equation: n.x * p.x + n.y * p.y + n.z * p.z - d = 0
struct Plane {
    Point3D n; // Normal vector
    double d = 0.0; // Distance from origin
};

// --- Function Declarations ---

Polyhedron create_cube(double min_x, double max_x, double min_y, double max_y, double min_z, double max_z);

Polyhedron clip_polyhedron(const Polyhedron& subject, const Plane& plane, bool keep_positive_side = false);

void print_polyhedron(const Polyhedron& poly);

// --- Inline Helpers ---

inline double dot(const Point3D& a, const Point3D& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

#endif // GEOMETRY_H