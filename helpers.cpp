#include "helpers.h"
#include <iostream>
#include <vector>

using namespace INMOST;

// This function creates a single cube cell in the mesh and returns it.
Cell CreateCubeAtPoint(Mesh *m, double center_x, double center_y, double center_z, double size_x, double size_y, double size_z)
{
  double half_size_x = size_x / 2.0;
  double half_size_y = size_y / 2.0;
  double half_size_z = size_z / 2.0;

  // 1. Create 8 nodes for the cube corners relative to the center point
  ElementArray<Node> nodes(m);
  Storage::real xyz[3];

  double corners[8][3] = {
      {center_x - half_size_x, center_y - half_size_y, center_z - half_size_z}, // 0
      {center_x + half_size_x, center_y - half_size_y, center_z - half_size_z}, // 1
      {center_x + half_size_x, center_y + half_size_y, center_z - half_size_z}, // 2
      {center_x - half_size_x, center_y + half_size_y, center_z - half_size_z}, // 3
      {center_x - half_size_x, center_y - half_size_y, center_z + half_size_z}, // 4
      {center_x + half_size_x, center_y - half_size_y, center_z + half_size_z}, // 5
      {center_x + half_size_x, center_y + half_size_y, center_z + half_size_z}, // 6
      {center_x - half_size_x, center_y + half_size_y, center_z + half_size_z}  // 7
  };

  for(int i = 0; i < 8; ++i)
  {
    xyz[0] = corners[i][0];
    xyz[1] = corners[i][1];
    xyz[2] = corners[i][2];
    nodes.push_back(m->CreateNode(xyz));
  }

  // 2. Create 6 faces for the cube
  ElementArray<Face> faces(m);
  ElementArray<Node> face_nodes(m, 4);

  // Order of nodes matters for face orientation.
  // This order should result in outward-pointing normals for a right-handed system.
  int face_indices[6][4] = {
      {0, 3, 2, 1}, // Bottom face (z-negative)
      {4, 5, 6, 7}, // Top face (z-positive)
      {0, 1, 5, 4}, // Front face (y-negative)
      {2, 3, 7, 6}, // Back face (y-positive)
      {0, 4, 7, 3}, // Left face (x-negative)
      {1, 2, 6, 5}  // Right face (x-positive)
  };

  for(int i = 0; i < 6; ++i)
  {
    face_nodes[0] = nodes[face_indices[i][0]];
    face_nodes[1] = nodes[face_indices[i][1]];
    face_nodes[2] = nodes[face_indices[i][2]];
    face_nodes[3] = nodes[face_indices[i][3]];
    faces.push_back(m->CreateFace(face_nodes).first);
  }

  // 3. Create the cell from the faces and return it
  return m->CreateCell(faces).first;
}
