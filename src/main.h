
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


// Contains the indices of two atoms which constitute a pair
typedef int Pair[2];

static const int X_DIM = 0;
static const int Y_DIM = 1;
static const int Z_DIM = 2;
static const int DIM_  = 3;

// Position vector
typedef double Rvec[DIM_];

using Int3  = std::array<int, DIM_>;
using Real3 = std::array<double, DIM_>;

// Apply the minimum image convention
void minImage(const Rvec x1, const Rvec x2, const Rvec boxL, Rvec x12, double& distSq);
void minImage(
	const std::vector<double>& x1, 
	const std::vector<double>& x2, 
	const std::vector<double>& boxL, 
	// Output
	std::vector<double>& x12, 
	double& distSq
);

// Keep the atom in the simulaton box
void keepInBox(const Rvec boxL, Rvec x);

void keepInBox(const std::vector<double>& boxL, std::vector<double>& x);

// Checks whether each oxygen has 2 hydrogens
bool areHydrogensCorrectlyPlaced(std::vector<int> hydrogenCounts);

// Vector norm
double norm(const Rvec x);
double norm(const std::vector<double>& x);

// Use the dot product to get the angle between the vectors
double angleBetweenVectors(const Rvec a, const Rvec b);

// Get linear cell indx from a triple of grid indices
int getLinearCellIndex(const Int3& grid_indices, const Int3& grid_dimensions) {
	return grid_indices[X_DIM]*grid_dimensions[Y_DIM]*grid_dimensions[Z_DIM] +
	       grid_indices[Y_DIM]*grid_dimensions[Z_DIM] +
	       grid_indices[Z_DIM];
}

#endif // MAIN_H
