#ifndef HELPERS_H
#define HELPERS_H

#include "inmost.h"

using namespace INMOST;

// This function creates a single cube cell in the mesh and returns it.
Cell CreateCubeAtPoint(Mesh *m, double center_x, double center_y, double center_z, double size_x, double size_y, double size_z);

#endif // HELPERS_H
