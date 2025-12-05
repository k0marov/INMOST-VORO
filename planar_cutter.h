#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include "inmost.h" // Make sure INMOST include path is set
using namespace INMOST;

const double EPS = 1e-9;

// ===================== Vector and Plane =====================
struct Vec3 {
  double x, y, z;
  Vec3() : x(0), y(0), z(0) {}
  Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
  Vec3 operator+(const Vec3 &b) const { return Vec3(x+b.x, y+b.y, z+b.z); }
  Vec3 operator-(const Vec3 &b) const { return Vec3(x-b.x, y-b.y, z-b.z); }
  Vec3 operator*(double s) const { return Vec3(x*s, y*s, z*s); }
};

static Vec3 cross(const Vec3 &a, const Vec3 &b) {
  return Vec3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}
static double dot(const Vec3 &a,const Vec3 &b){ return a.x*b.x + a.y*b.y + a.z*b.z; }
static double norm(const Vec3 &a){ return sqrt(dot(a,a)); }
static Vec3 normalize(const Vec3 &a){ double n=norm(a); return n<1e-15 ? a : a*(1.0/n); }

struct Plane {
  Vec3 n; // normal
  double d;
  Plane() : n(), d(0) {}
  Plane(const Vec3 &N, double D) : n(N), d(D) {}
  double distance(const Vec3 &p) const { return dot(n,p)+d; }
};

// ===================== Polyhedron Representation =====================
using FaceRaw = std::vector<Vec3>;
using PolyhedronRaw = std::vector<FaceRaw>; // each face: list of points

struct Polyhedron {
  std::vector<Vec3> vertices;
  std::vector<std::vector<int>> faces; // faces use vertex indices
};

// ===================== Helper functions =====================
Vec3 intersect_segment_plane(const Vec3 &p, const Vec3 &q, const Plane &pl) {
  double dp = pl.distance(p);
  double dq = pl.distance(q);
  double t = dp - dq;
  if (fabs(t)<1e-15) return p;
  double alpha = dp / (dp - dq);
  if (alpha < 0) alpha = 0;
  if (alpha > 1) alpha = 1;
  return p + (q - p) * alpha;
}

FaceRaw clip_polygon_by_plane(const FaceRaw &poly, const Plane &pl, std::vector<Vec3> &intersections) {
  FaceRaw out;
  if (poly.empty()) return out;
  int m = (int)poly.size();
  for (int i=0;i<m;i++){
    Vec3 A=poly[i]; Vec3 B=poly[(i+1)%m];
    double da=pl.distance(A); double db=pl.distance(B);
    bool ina=(da<=EPS), inb=(db<=EPS);
    if (ina && inb) { out.push_back(B); }
    else if (ina && !inb) { Vec3 I=intersect_segment_plane(A,B,pl); out.push_back(I); intersections.push_back(I); }
    else if (!ina && inb) { Vec3 I=intersect_segment_plane(A,B,pl); out.push_back(I); out.push_back(B); intersections.push_back(I); }
  }
  // remove near-duplicate consecutive vertices
  FaceRaw cleaned;
  for (auto &v : out) {
    if (cleaned.empty()) cleaned.push_back(v);
    else if (norm(v-cleaned.back())>1e-12) cleaned.push_back(v);
  }
  if (cleaned.size()>=2 && norm(cleaned.front()-cleaned.back())<1e-12) cleaned.pop_back();
  return cleaned;
}

std::vector<Vec3> unique_points(const std::vector<Vec3> &pts, double tol=1e-9) {
  std::vector<Vec3> res;
  for (auto &p : pts) {
    bool found=false;
    for (auto &q : res) if (norm(p-q)<=tol){found=true; break;}
    if (!found) res.push_back(p);
  }
  return res;
}

FaceRaw sort_points_on_plane(const std::vector<Vec3> &pts, const Plane &pl) {
  FaceRaw ordered;
  if (pts.size()<3) return ordered;
  Vec3 centroid(0,0,0);
  for (auto &p: pts) centroid = centroid + p;
  centroid = centroid * (1.0/pts.size());
  Vec3 n = normalize(pl.n);
  Vec3 arbitrary = fabs(n.x)<0.9 ? Vec3(1,0,0) : Vec3(0,1,0);
  Vec3 u = normalize(cross(n, arbitrary));
  Vec3 v = cross(n,u);
  std::vector<std::pair<double, Vec3>> ang_pts;
  for (auto &p : pts) {
    Vec3 rel = p - centroid;
    double angle = atan2(dot(rel,v), dot(rel,u));
    ang_pts.emplace_back(angle,p);
  }
  std::sort(ang_pts.begin(), ang_pts.end(), [](auto &a, auto &b){return a.first<b.first;});
  for (auto &ap : ang_pts) ordered.push_back(ap.second);
  return ordered;
}

// ===================== Polyhedron clipping =====================
PolyhedronRaw clip_polyhedron(const PolyhedronRaw &poly, const Plane &pl) {
  PolyhedronRaw newpoly;
  std::vector<Vec3> all_intersections;
  for (auto &face : poly) {
    FaceRaw clipped = clip_polygon_by_plane(face, pl, all_intersections);
    if (clipped.size()>=3) newpoly.push_back(clipped);
  }
  std::vector<Vec3> uniq = unique_points(all_intersections,1e-9);
  if (uniq.size()>=3) {
    FaceRaw cap = sort_points_on_plane(uniq,pl);
    if (cap.size()>=3) {
      Vec3 e1=cap[1]-cap[0], e2=cap[2]-cap[0];
      Vec3 capn=normalize(cross(e1,e2));
      if (dot(capn,pl.n)>0) std::reverse(cap.begin(),cap.end());
      newpoly.push_back(cap);
    }
  }
  return newpoly;
}

// ===================== Convert to shared-vertex Polyhedron =====================
Polyhedron unify_vertices(const PolyhedronRaw &raw, double tol=1e-12) {
  Polyhedron out;
  auto find_or_add = [&](const Vec3 &v){
    for (int i=0;i<(int)out.vertices.size();++i)
      if (norm(out.vertices[i]-v)<tol) return i;
    out.vertices.push_back(v);
    return (int)out.vertices.size()-1;
  };
  for (auto &face_raw : raw) {
    std::vector<int> f;
    f.reserve(face_raw.size());
    for (auto &v: face_raw) f.push_back(find_or_add(v));
    out.faces.push_back(std::move(f));
  }
  return out;
}

// ===================== INMOST conversion =====================
Cell add_polyhedron_to_inmost(Mesh *mesh, const Polyhedron &poly) {
  mesh->BeginModification();  // Begin modification epoch

  // 1. Create nodes using ElementArray and const real* coordinates
  ElementArray<Node> nodes(mesh);
  nodes.reserve(poly.vertices.size());
  for (const auto &v : poly.vertices) {
    const Storage::real coords[3] = {v.x, v.y, v.z};
    auto n = mesh->CreateNode(coords);  // assuming pair return
    nodes.push_back(n);
    // optionally: check 'created'
  }

  // 2. Create faces using ElementArray
  ElementArray<Face> faces(mesh);
  faces.reserve(poly.faces.size());
  for (const auto &fidx : poly.faces) {
    ElementArray<Node> fnodes(mesh);
    fnodes.reserve(fidx.size());
    for (int vid : fidx) fnodes.push_back(nodes[vid]);
    auto [f, created] = mesh->CreateFace(fnodes);  // assuming pair return
    faces.push_back(f);
  }

  // 3. Create cell from faces
  auto [c, created] = mesh->CreateCell(faces);  // assuming pair return

  mesh->EndModification();  // Finish modification epoch

  return c;
}
// ===================== Cube helper =====================
// Create a rectangular cuboid/polyhedron with system bounds
PolyhedronRaw create_cube(double min_x, double max_x,
                          double min_y, double max_y,
                          double min_z, double max_z)
{
  PolyhedronRaw poly;

  // 8 vertices of the cuboid
  std::vector<Vec3> V = {
      Vec3(min_x, min_y, min_z),
      Vec3(max_x, min_y, min_z),
      Vec3(max_x, max_y, min_z),
      Vec3(min_x, max_y, min_z),
      Vec3(min_x, min_y, max_z),
      Vec3(max_x, min_y, max_z),
      Vec3(max_x, max_y, max_z),
      Vec3(min_x, max_y, max_z)
  };

  // 6 faces of the cuboid (each face is a vector of vertices in CCW order)
  poly.push_back({V[0], V[1], V[2], V[3]}); // bottom face (z = min_z)
  poly.push_back({V[4], V[5], V[6], V[7]}); // top face (z = max_z)
  poly.push_back({V[0], V[1], V[5], V[4]}); // front face (y = min_y)
  poly.push_back({V[1], V[2], V[6], V[5]}); // right face (x = max_x)
  poly.push_back({V[2], V[3], V[7], V[6]}); // back face (y = max_y)
  poly.push_back({V[3], V[0], V[4], V[7]}); // left face (x = min_x)

  return poly;
}

#endif // GEOMETRY_H