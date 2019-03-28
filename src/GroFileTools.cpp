#include "GroFileTools.h"


void GroFileTools::readGroFile(
	const std::string& gro_file, std::string& title,
	std::vector<int>& residue_numbers, std::vector<std::string>& residue_names, 
	std::vector<std::string>& atom_names, std::vector<Real3>& positions, Matrix& box_matrix
) 
{
	// Try to open the file
	std::ifstream ifs;
	ifs.open(gro_file);
	if ( ! ifs.is_open() ) {
		std::stringstream err_ss;
		err_ss << "GroFileTools::readGroFile - Unable to open .gro file.\n"
		       << "  Input: \"" + gro_file + "\"\n";
		throw std::runtime_error( err_ss.str() );
	}

	// Working variables
	std::string line;
	std::istringstream ss;

	// Read header
	int num_atoms_from_gro_file = -1; // sentinel value
	getline(ifs, title); // First line is the title
	getline(ifs, line);  // Second line contains the number of atoms 
	ss.str(line);
	ss >> num_atoms_from_gro_file;

	// Checks
	if ( num_atoms_from_gro_file <= 0 ) {
		std::stringstream err_ss;
		err_ss << "  GroFileTools::readGroFile - Invalid number of atoms ("
           << num_atoms_from_gro_file << ") reported by .gro file." << "\n"
		       << "  Input: " << gro_file << "\n";
		throw std::runtime_error( err_ss.str() );
	}

	// Allocate memory
	residue_numbers.resize(num_atoms_from_gro_file);
	residue_names.resize(num_atoms_from_gro_file);
	atom_names.resize(num_atoms_from_gro_file);
	positions.resize(num_atoms_from_gro_file);

	// Read information about atoms
	int atom_counter = 0, atom_serial;
	std::stringstream parsing_buffer;
	while ( (atom_counter<num_atoms_from_gro_file) and getline(ifs, line) ) {
		// Parse line
		ss.str(line);

		// GROMACS FORMAT: "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
		// - Number of characters per field: 5, 5, 5, 5, 8, 8, 8, 8, 8, 8
		// - Use a stringstream to neatly handle the different variable types

		// Clear the parsing buffer
		parsing_buffer.str( std::string() );
		parsing_buffer.clear();

		// Field 1: residue index (0-4, 5 char)
		parsing_buffer << line.substr(0, 5) << "\t";

		// Field 2: residue name (5-9, 5 char)
		parsing_buffer << line.substr(5, 5) << "\t";

		// Field 3: atom name (10-14, 5 char)
		parsing_buffer << line.substr(10, 5) << "\t";

		// Field 4: atom serial (15-19, 5 char)
		parsing_buffer << line.substr(15, 5) << "\t";

		// Field 5: x-position (20-27, 8 char)
		parsing_buffer << line.substr(20, 8) << "\t";

		// Field 6: y-position (28-35, 8 char)
		parsing_buffer << line.substr(28, 8) << "\t";

		// Field 7: z-position (36-43, 8 char)
		parsing_buffer << line.substr(36, 8) << "\t";

		// Store variables as needed
		parsing_buffer >> residue_numbers[atom_counter];;
		parsing_buffer >> residue_names[atom_counter];
		parsing_buffer >> atom_names[atom_counter];
		parsing_buffer >> atom_serial;
		for ( int d=0; d<DIM_; ++d ) {
			parsing_buffer >> positions[atom_counter][d];
		}

		atom_counter++;
	}

	if ( atom_counter != num_atoms_from_gro_file ) {
		ifs.close();
		std::stringstream ss;
		ss << "Error in GroFileTools::readGroFile: expected " << num_atoms_from_gro_file
		   << " atoms but found " << atom_counter;
		throw std::runtime_error( ss.str() );
	}

	// Last line has the box lengths
	// (FIXME: box matrix, if triclinic? how to check if extra, unexpected lines?)
	getline(ifs, line);
	ss.clear();   // Clear flags before using with a new line
	ss.str(line);
	for ( int d=0; d<3; ++d ) {
		ss >> box_matrix[d][d];
	}

	// Wrap up
	ifs.close();

	// Checks
	for ( int d=0; d<DIM_; ++d ) {
		if ( box_matrix[d][d] <= 0.0 ) {
			std::stringstream err_ss;
			err_ss << "GroFileTools::readGroFile - Invalid box length for dimension d=" 
						 << d+1 << " of 3 (reference_box_matrix[d][d] = " << box_matrix[d][d]
						 << ")." << "\n";
			throw std::runtime_error( err_ss.str() );
		}
	}
}


void GroFileTools::writeGroFile(
	const std::string& gro_file, const std::string& title,
	const std::vector<int>& residue_numbers, const std::vector<std::string>& residue_names, 
	const std::vector<std::string>& atom_names, const std::vector<Real3>& positions, 
	const Matrix& box_matrix
)
{
	// Try to open the file
	std::ofstream ofs(gro_file);
	if ( ofs.fail() ) {
		std::stringstream err_ss;
		err_ss << "GroFileTools::writeGroFile - Encountered an error while opening .gro file.\n"
		       << "  Input: \"" + gro_file + "\"\n";
		throw std::runtime_error( err_ss.str() );
	}

	// Checks (TODO)
	int num_atoms = positions.size();

	// Title
	int pos = title.find('\n');
	if ( pos != std::string::npos ) {
		ofs << title << "\n";
	}
	else {
		// Ignore everything from the newline character onwards
		ofs << title.substr(0, pos) << "\n";
	}

	// Number of atoms gets its own line
	ofs << num_atoms << "\n";

	// Print information for each atom
	ofs << std::fixed << std::right;
	for ( int i=0; i<num_atoms; ++i ) {
		// Field 1: residue number (0-4, 5 char)
		StringTools::printFixedWidthString(ofs, residue_numbers[i], 5);

		// Field 2: residue name (5-9, 5 char) 
		ofs << std::left;
		StringTools::printFixedWidthString(ofs, residue_names[i], 5);
		ofs << std::right;

		// Field 3: atom name (10-14, 5 char)
		StringTools::printFixedWidthString(ofs, atom_names[i], 5);

		// Field 4: atom serial (15-19, 5 char)
		StringTools::printFixedWidthString(ofs, i+1, 5);

		// Fields 5-7: position (20-27,28-35,36-43, 8 char each)
		// - TODO: Files with a x,y,z >999 nm can't be handled by the code below
		//   - GROMACS has started to allow .gro files with extended position fields,
		//     but support is not universal throughout the code
		ofs << std::setprecision(3);
		for ( int d=0; d<DIM_; ++d ) {
			ofs << std::setw(8) << positions[i][d];

			// This other approach works with gmx editconf (v2016.3), suggesting it's generally valid,
			// but produces files that look really ugly and are hard to read
			// - It pads leading zeros for some reason, rather than whitespace
			// - Also: what about handling the sign when some positions are negative?!
			//StringTools::printFixedWidthString(ofs, positions[i][d], 8);
		}

		// TODO velocity
		//ofs << std::setprecision(4);  // extra place

		ofs << "\n";
	}

	// Last line has the box 
	// (TODO assumes orthorhombic box matrix)
	ofs << std::setprecision(5);
	for ( int d=0; d<3; ++d ) {
		ofs << "  " << box_matrix[d][d];
	}
	ofs << "\n";

	// Wrap up
	ofs.close();
}


void GroFileTools::shiftPositions(std::vector<Real3>& positions, const Real3& shift)
{
	int num_atoms = positions.size();
	for ( int i=0; i<num_atoms; ++i ) {
		for ( int d=0; d<DIM_; ++d ) {
			positions[i][d] += shift[d];
		}
	}
}
