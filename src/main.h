// MakeKaoliniteSlab
//
// - Creates a slab of kaolinite and a GROMACS topology file for it
//   according to the CLAY-FF force field
//   - Also makes a slab which is reflected across the xy-plane
// - Notes
//   - A "layer" of kaolinite is created by translating the unit cell in the xy-plane
//   - The hexagonal rings of Si and O are parallel to the face of the layer
//   - A "slab" of kaolinite is created by stacking slabs on top of each other

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

#include "ClayFF.h"
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


void writePositionRestraints(
	const std::string& posre_file,
	const std::string& header,
	const std::vector<int>& atom_indices,
	const std::string& parameter_comment,  // Comment describing each entry
	const std::string& parameters          // Contains all parameters for the restraints (except the atom serial)
);

void writeIndexFile(
	const std::string& ndx_file,
	const std::string& header,
	const std::string& group_name,
	const std::vector<int>& atom_indices,
	const int num_atoms_per_line = 20
);

#endif // MAIN_H
