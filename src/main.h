
#ifndef MAIN_H
#define MAIN_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "CommonTypes.h"
#include "GroFileTools.h"
#include "SimulationBox.h"

// Vector norm
double norm(const Real3& x);
double norm(const std::vector<double>& x);

// Use the dot product to get the angle between the vectors
double angleBetweenVectors(const Real3& a, const Real3& b);

// Get linear cell indx from a triple of grid indices
int getLinearCellIndex(const Int3& grid_indices, const Int3& grid_dimensions) {
	return grid_indices[X_DIM]*grid_dimensions[Y_DIM]*grid_dimensions[Z_DIM] +
	       grid_indices[Y_DIM]*grid_dimensions[Z_DIM] +
	       grid_indices[Z_DIM];
}

#endif // MAIN_H
