
// gro file format: "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
//                       ^left-aligned
// - Number of characters per field: 5, 5, 5, 5, 8, 8, 8, 8, 8, 8

#ifndef GRO_FILE_TOOLS_H
#define GRO_FILE_TOOLS_H

#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "CommonTypes.h"
#include "StringTools.h"

class GroFileTools 
{
 public:
	static void readGroFile(
		const std::string& gro_file,
		// Output
		std::string&              title,
		std::vector<int>&         residue_numbers, 
		std::vector<std::string>& residue_names, 
		std::vector<std::string>& atom_names,
		std::vector<Real3>&       positions, 
		Matrix&                   box_matrix
	);

	static void writeGroFile(
		const std::string&              gro_file,
		const std::string&              title,
		const std::vector<int>&         residue_numbers, 
		const std::vector<std::string>& residue_names, 
		const std::vector<std::string>& atom_names,
		const std::vector<Real3>&       positions, 
		const Matrix&                   box_matrix
	);

	static void writeGroFile(
		const std::string&              gro_file,
		const std::string&              title,
		const std::vector<int>&         residue_numbers, 
		const std::vector<std::string>& residue_names, 
		const std::vector<std::string>& atom_names,
		const std::vector<Real3>&       positions, 
		const Real3&                    box_lengths
	) {
		Box box_matrix;
		for ( int a=0; a<DIM_; ++a ) {
			for ( int b=0; b<DIM_; ++b ) {
				if ( a == b ) {
					box_matrix[a][a] = box_lengths[a];
				}
				else {
					box_matrix[a][b] = 0;
				}
			}
		}

		writeGroFile(gro_file, title, residue_numbers, residue_names, atom_names, positions, box_matrix);
	}

	static void shiftPositions(
		std::vector<Real3>& positions,
		const Real3& shift
	);
};

#endif /* GRO_FILE_TOOLS_H */
