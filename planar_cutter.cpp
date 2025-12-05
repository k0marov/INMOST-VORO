#include "planar_cutter.h"
#include <iostream>

void print_polyhedron(const Polyhedron& poly) {
  std::cout << "Polyhedron has " << poly.faces.size() << " faces:" << std::endl;
  int face_num = 0;
  for (const auto& face : poly.faces) {
    std::cout << "  Face " << face_num++ << " has " << face.vertices.size() << " vertices:" << std::endl;
    for (const auto& v : face.vertices) {
      std::cout << "    (" << v.x << ", " << v.y << ", " << v.z << ")" << std::endl;
    }
  }
}

// --- Vector Math Helpers ---

static Point3D subtract(const Point3D& a, const Point3D& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

static Point3D add(const Point3D& a, const Point3D& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

static Point3D scale(const Point3D& a, double s) {
    return {a.x * s, a.y * s, a.z * s};
}

static Point3D cross(const Point3D& a, const Point3D& b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

// --- Clipping Algorithm ---

/**
 * @brief Calculates the signed distance of a point from a plane.
 */
static double classify_vertex(const Point3D& v, const Plane& plane) {
    return dot(plane.n, v) - plane.d;
}

/**
 * @brief Calculates the intersection point of a line segment with a plane.
 */
static Point3D intersect_edge(const Point3D& p1, const Point3D& p2, double d1, double d2) {
    Point3D line_dir = subtract(p2, p1);
    double t = d1 / (d1 - d2);
    return {p1.x + t * line_dir.x, p1.y + t * line_dir.y, p1.z + t * line_dir.z};
}


Polyhedron clip_polyhedron(const Polyhedron& subject, const Plane& plane, bool keep_positive_side) {
    Polyhedron result;
    std::vector<Point3D> cap_face_vertices;

    for (const auto& face : subject.faces) {
        if (face.vertices.empty()) continue;

        Polygon clipped_face;
        Point3D S = face.vertices.back();
        double dist_S = classify_vertex(S, plane);

        for (const Point3D& E : face.vertices) {
            double dist_E = classify_vertex(E, plane);

            bool s_is_kept, e_is_kept;
            if (keep_positive_side) {
                s_is_kept = dist_S > -GEOMETRY_EPSILON;
                e_is_kept = dist_E > -GEOMETRY_EPSILON;
            } else {
                s_is_kept = dist_S < GEOMETRY_EPSILON;
                e_is_kept = dist_E < GEOMETRY_EPSILON;
            }

            if (s_is_kept != e_is_kept) {
                Point3D I = intersect_edge(S, E, dist_S, dist_E);
                clipped_face.vertices.push_back(I);
                cap_face_vertices.push_back(I);
            }

            if (e_is_kept) {
                clipped_face.vertices.push_back(E);
            }

            S = E;
            dist_S = dist_E;
        }

        if (clipped_face.vertices.size() >= 3) {
            result.faces.push_back(clipped_face);
        }
    }

    if (cap_face_vertices.size() >= 3) {
        // First, remove duplicate vertices that may have been added from shared edges
        auto comparator = [](const Point3D& a, const Point3D& b) {
            if (std::abs(a.x - b.x) > GEOMETRY_EPSILON) return a.x < b.x;
            if (std::abs(a.y - b.y) > GEOMETRY_EPSILON) return a.y < b.y;
            return a.z < b.z;
        };
        std::sort(cap_face_vertices.begin(), cap_face_vertices.end(), comparator);
        cap_face_vertices.erase(std::unique(cap_face_vertices.begin(), cap_face_vertices.end(),
            [](const Point3D& a, const Point3D& b) {
                return std::abs(a.x - b.x) < GEOMETRY_EPSILON &&
                       std::abs(a.y - b.y) < GEOMETRY_EPSILON &&
                       std::abs(a.z - b.z) < GEOMETRY_EPSILON;
            }), cap_face_vertices.end());

        if (cap_face_vertices.size() >= 3) {
            Point3D centroid = {0, 0, 0};
            for (const auto& v : cap_face_vertices) {
                centroid = add(centroid, v);
            }
            centroid.x /= cap_face_vertices.size();
            centroid.y /= cap_face_vertices.size();
            centroid.z /= cap_face_vertices.size();
            
            Point3D sort_axis = keep_positive_side ? scale(plane.n, -1.0) : plane.n;

            std::sort(cap_face_vertices.begin(), cap_face_vertices.end(),
                      [&](const Point3D& a, const Point3D& b) {
                        Point3D vec_a = subtract(a, centroid);
                        Point3D vec_b = subtract(b, centroid);
                        double sign = dot(sort_axis, cross(vec_a, vec_b));
                        return sign > 0;
                      });
            
            result.faces.push_back({cap_face_vertices});
        }
    }

    return result;
}

Polyhedron create_cube(double min_x, double max_x, double min_y, double max_y, double min_z, double max_z) {
    Polyhedron cube;
    // Note: Winding order matters for face normals. Assuming outward-facing normals.
    cube.faces.push_back({{{min_x, min_y, min_z}, {max_x, min_y, min_z}, {max_x, max_y, min_z}, {min_x, max_y, min_z}}}); // Bottom (-z)
    cube.faces.push_back({{{min_x, min_y, max_z}, {min_x, max_y, max_z}, {max_x, max_y, max_z}, {max_x, min_y, max_z}}}); // Top (+z)
    cube.faces.push_back({{{min_x, min_y, min_z}, {min_x, min_y, max_z}, {max_x, min_y, max_z}, {max_x, min_y, min_z}}}); // Front (-y)
    cube.faces.push_back({{{min_x, max_y, min_z}, {max_x, max_y, min_z}, {max_x, max_y, max_z}, {min_x, max_y, max_z}}}); // Back (+y)
    cube.faces.push_back({{{min_x, min_y, min_z}, {min_x, max_y, min_z}, {min_x, max_y, max_z}, {min_x, min_y, max_z}}}); // Left (-x)
    cube.faces.push_back({{{max_x, min_y, min_z}, {max_x, min_y, max_z}, {max_x, max_y, max_z}, {max_x, max_y, min_z}}}); // Right (+x)
    return cube;
}
