#ifndef PLANAR_CUTTER_H
#define PLANAR_CUTTER_H

#include "inmost.h"

using namespace INMOST;

/// A self-contained class for cutting a single INMOST Cell with a plane.
class PlanarCutter {
public:
    PlanarCutter();

    /// Cuts a single cell by a plane and returns the resulting cell.
    /// The original cell is destroyed in the process.
    /// @param cell_to_cut The cell to be cut.
    /// @param a, b, c, d The coefficients of the plane equation ax + by + cz = d.
    /// @param cut_positive If true, the part of the cell where ax+by+cz-d > 0 is cut away.
    /// @return The new cell that results from the cut. Returns an invalid cell if the original was entirely cut away or not cut at all.
    Cell Cut(Cell& cell_to_cut, double a, double b, double c, double d, bool cut_positive);
};

#endif // PLANAR_CUTTER_H
