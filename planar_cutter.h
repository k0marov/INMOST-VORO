#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <map>

#include "inmost.h" // Make sure INMOST include path is set
using namespace INMOST;

const double EPS = 1e-12;

// ===================== Vector =====================
struct Vec3 {
  double x,y,z;
  Vec3() : x(0), y(0), z(0) {}
  Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
  Vec3 operator-(const Vec3& b) const { return {x-b.x, y-b.y, z-b.z}; }
  Vec3 operator+(const Vec3& b) const { return {x+b.x, y+b.y, z+b.z}; }
  Vec3 operator*(double s) const { return {x*s, y*s, z*s}; }
};

inline double dot(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}
inline double length(const Vec3& v) { return std::sqrt(dot(v,v)); }
inline Vec3 normalize(const Vec3& v) {
  double l = length(v);
  return l<EPS ? v : v*(1.0/l);
}

// ===================== GFace & GCell =====================
struct GFace {
  Vec3 normal;
  std::vector<Vec3> vertices;
};

struct GCell {
  std::vector<GFace> faces;
};

// ===================== Plane intersection =====================
inline bool segmentPlaneIntersection(const Vec3& p0, const Vec3& p1,
                                     const Vec3& planePoint, const Vec3& planeNormal,
                                     Vec3& outIntersect)
{
  Vec3 dir = p1 - p0;
  double denom = dot(planeNormal, dir);
  if (std::abs(denom) < EPS) return false;
  double t = dot(planePoint - p0, planeNormal) / denom;
  if (t < -EPS || t > 1.0+EPS) return false;
  outIntersect = p0 + dir*t;
  return true;
}

inline Vec3 centroid(const std::vector<Vec3>& pts) {
  Vec3 c{0,0,0};
  for(auto& p : pts) c = c + p;
  return c * (1.0/pts.size());
}

// ===================== Clip GCell by Plane =====================
inline void clipGCellByPlane(GCell& cell, const Vec3& planePoint, const Vec3& planeNormal) {
  std::vector<GFace> newGFaces;
  std::vector<Vec3> planeIntersections;

  for(auto& face : cell.faces) {
    std::vector<Vec3> newVerts;
    int N = face.vertices.size();
    for(int i=0;i<N;i++) {
      Vec3 curr = face.vertices[i];
      Vec3 next = face.vertices[(i+1)%N];
      double dCurr = dot(curr - planePoint, planeNormal);
      double dNext = dot(next - planePoint, planeNormal);

      if(dCurr >= -EPS) newVerts.push_back(curr);

      Vec3 inter;
      if((dCurr > EPS && dNext < -EPS) || (dCurr < -EPS && dNext > EPS)) {
        if(segmentPlaneIntersection(curr, next, planePoint, planeNormal, inter)) {
          newVerts.push_back(inter);
          planeIntersections.push_back(inter);
        }
      }
    }
    if(newVerts.size() >= 3) {
      face.vertices = newVerts;
      newGFaces.push_back(face);
    }
  }

  if(planeIntersections.size() >= 3) {
    Vec3 c = centroid(planeIntersections);
    Vec3 n = planeNormal;

    auto sortedVerts = planeIntersections;
    std::sort(sortedVerts.begin(), sortedVerts.end(),
              [&c,&n](const Vec3& a, const Vec3& b){
                Vec3 va = a-c, vb = b-c;
                double angle = atan2(dot(cross(va,vb), n), dot(va,vb));
                return angle > 0;
              });

    GFace newGFace;
    newGFace.normal = n;
    newGFace.vertices = sortedVerts;
    newGFaces.push_back(newGFace);
  }

  cell.faces = newGFaces;
}

// ===================== Bounding Cube =====================
inline GCell createBoundingCube(const Vec3& seed, const std::vector<Vec3>& neighbors, double margin=1.0) {
  double minX = seed.x, maxX = seed.x;
  double minY = seed.y, maxY = seed.y;
  double minZ = seed.z, maxZ = seed.z;

  for(auto& n : neighbors) {
    if(n.x<minX) minX=n.x; if(n.x>maxX) maxX=n.x;
    if(n.y<minY) minY=n.y; if(n.y>maxY) maxY=n.y;
    if(n.z<minZ) minZ=n.z; if(n.z>maxZ) maxZ=n.z;
  }

  minX-=margin; maxX+=margin;
  minY-=margin; maxY+=margin;
  minZ-=margin; maxZ+=margin;

  std::vector<Vec3> v = {
      {minX,minY,minZ},{maxX,minY,minZ},{maxX,maxY,minZ},{minX,maxY,minZ},
      {minX,minY,maxZ},{maxX,minY,maxZ},{maxX,maxY,maxZ},{minX,maxY,maxZ}
  };

  GCell cell;
  cell.faces = {
      {{0,0,-1},{v[0],v[1],v[2],v[3]}},
      {{0,0,1},{v[4],v[5],v[6],v[7]}},
      {{-1,0,0},{v[0],v[3],v[7],v[4]}},
      {{1,0,0},{v[1],v[2],v[6],v[5]}},
      {{0,-1,0},{v[0],v[1],v[5],v[4]}},
      {{0,1,0},{v[3],v[2],v[6],v[7]}}
  };
  return cell;
}

// ===================== Voronoi GCell =====================
inline GCell createVoronoiCell(const Vec3& seed, const std::vector<Vec3>& neighbors){
  GCell cell = createBoundingCube(seed, neighbors, 1.0);

  for(auto& n : neighbors) {
    Vec3 dir = seed - n;
    Vec3 planeNormal = normalize(dir);
    Vec3 planePoint = seed + (n-seed)*0.5;
    clipGCellByPlane(cell, planePoint, planeNormal);
  }

  return cell;
}

// ===================== INMOST conversion =====================
inline GCell add_polyhedron_to_inmost(Mesh* mesh, const GCell& cell) {
  mesh->BeginModification();

  ElementArray<Node> nodes(mesh);
  nodes.reserve(1000); // heuristic

  ElementArray<Face> faces(mesh);
  faces.reserve(cell.faces.size());

  for(const auto& f : cell.faces) {
    ElementArray<Node> fnodes(mesh);
    fnodes.reserve(f.vertices.size());
    for(const auto& v : f.vertices) {
      const Storage::real coords[3] = {v.x, v.y, v.z};
      fnodes.push_back(mesh->CreateNode(coords));
    }
    auto [faceObj, createdGFace] = mesh->CreateFace(fnodes);
    faces.push_back(faceObj);
  }

  auto [cellObj, createdGCell] = mesh->CreateCell(faces);

  mesh->EndModification();

  return cell;
}

#endif // GEOMETRY_H