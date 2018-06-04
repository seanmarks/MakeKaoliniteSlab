// GenKaoliniteSlab.exe
//
// Updated Jun. 1, 2018
//
// - A "slab" of kaolinite is created by translating the unit cell in the xy-plane
//   - The hexagonal rings of Si and O are parallel to the face of the slab
// - A "layer" of kaolinite is created by stacking slabs on top of each other


// Project headers
#include "main.h"

int main(int argc, char* argv[])
{
	//----- Input -----//

	if ( argc < DIM_ + 1 ) {
		std::cout << "Usage:  ./GenKaoliniteSlab <nx> <ny> <nz>\n";
		return 0;
	}

	// Number of unit cells along x, y, and z
	std::array<int,DIM_> unit_cell_grid;
	for ( int d=0; d<DIM_; ++d ) {
		unit_cell_grid[d] = std::stod( argv[d+1] );
	}
	//unit_cell_grid = { 11, 7, 2 }; // Full grid
	//unit_cell_grid = { 1, 1, 1 }; // single unit cell
	//unit_cell_grid = { 4, 3, 2 }; // small grid


	//----- Constants -----//

	const double PI = 3.14159265358979323846;
	const double DEGREES_TO_RADIANS = PI/180.0;

	// Kaolinite unit cell
	double alpha = 91.926*DEGREES_TO_RADIANS;  // [rad]
	double beta  = 105.046*DEGREES_TO_RADIANS; // [rad]
	double gamma = 89.797*DEGREES_TO_RADIANS;  // [rad]

	double a = 0.51535; // [nm]
	double b = 0.89419; // [nm]
	double c = 0.73906; // [nm]
	std::vector<double> latticeConstants = {{ a, b, b }};

	//----- Fractional unit cell -----//

	// Atoms are stored in the following order:
	// 2x Al, 2x Si, 5x O (standalone), 4x O (hydroxyl), 4x H (hydroxyl)
	int num_atoms_per_primitive_cell = 17;
	int num_atoms_per_unit_cell = 2*num_atoms_per_primitive_cell;
	std::vector<std::vector<double>> fractional_unit_cell(num_atoms_per_unit_cell);

	fractional_unit_cell[0]  = {{ 0.28900, 0.49660, 0.46600 }}; // Al1 (octahedral)
	fractional_unit_cell[1]  = {{ 0.79300, 0.32880, 0.46500 }}; // Al2 (octahedral)
	fractional_unit_cell[2]  = {{ 0.98900, 0.33950, 0.09060 }}; // Si1 (tetrahedral)
	fractional_unit_cell[3]  = {{ 0.50700, 0.16650, 0.09380 }}; // Si2 (tetrahedral)

	fractional_unit_cell[4]  = {{ 0.04900, 0.34820, 0.31680 }}; // O1
	fractional_unit_cell[5]  = {{ 0.11300, 0.65990, 0.31880 }}; // O2
	fractional_unit_cell[6]  = {{ 0.00000, 0.50000, 0.00000 }}; // O3
	fractional_unit_cell[7]  = {{ 0.20400, 0.22910, 0.03000 }}; // O4
	fractional_unit_cell[8]  = {{ 0.19700, 0.76410, 0.00100 }}; // O5

	fractional_unit_cell[9]  = {{ 0.05000, 0.97100, 0.32500 }}; // O-hydroxyl-1
	fractional_unit_cell[10] = {{ 0.96000, 0.16580, 0.60700 }}; // O-hydroxyl-2
	fractional_unit_cell[11] = {{ 0.03700, 0.47260, 0.60460 }}; // O-hydroxyl-3
	fractional_unit_cell[12] = {{ 0.03800, 0.85820, 0.60900 }}; // O-hydroxyl-4

	fractional_unit_cell[13] = {{ 0.14500, 0.06510, 0.32600 }}; // H-hydroxyl-1
	fractional_unit_cell[14] = {{ 0.06300, 0.16380, 0.73900 }}; // H-hydroxyl-2
	fractional_unit_cell[15] = {{ 0.03600, 0.50570, 0.73200 }}; // H-hydroxyl-3
	fractional_unit_cell[16] = {{ 0.53400, 0.31540, 0.72800 }}; // H-hydroxyl-4

	//----- Atom type parameters for [ moleculetype ] ----//

	// Names of atoms for GROMACS
	// - Key (CLAYFF nomenclature):
	//     AlO  = aluminum, octahedrally coordinated
	//     SiT  = silicon, tetrahedrally coordinated
	//     OBTS = bridging oxygen bonded to a tetrahedral substitution
	//     OBSS = bridging oxygen bonded to two oct. subs OR double-subs.
	//     OH   = hydroxyl oxygen
	//     OHS  = hydroxyl oxygen bonded to a substituted site
	//     HH   = hydroxyl hydrogen
	// - Note: "substituted" --> non-oxygen
	// - Suffix "K" for kaolinite

	// Start with Al and Si
	std::vector<std::string> unit_cell_atom_names = { "AlOK", "AlOK", "SiTK", "SiTK" };
	std::vector<double> unit_cell_atom_masses = { 26.982, 26.982, 28.086, 28.086 };
	std::vector<double> unit_cell_atom_charges = { 1.575, 1.575, 2.1, 2.1 };

	// OBK: bridging O
	// - First 2 are "interior" (apical?)
	// - Other three are basal (form hexagonal rings with Si atoms)
	for ( int i=0; i<5; ++i ) {
		unit_cell_atom_names.push_back( std::string("OBK") );
		unit_cell_atom_masses.push_back( 16.000 );
		unit_cell_atom_charges.push_back( -1.05 );
	}

	// OHK: hydroxyl O coordinated with Al atoms
	for ( int i=0; i<4; ++i ) {
		unit_cell_atom_names.push_back( std::string("OHK") );
		unit_cell_atom_masses.push_back( 16.000 );
		unit_cell_atom_charges.push_back( -0.95 );
	}

	// HHK: hydroxyl hydrogens
	for ( int i=0; i<4; ++i ) {
		unit_cell_atom_names.push_back( std::string("HHK") );
		unit_cell_atom_masses.push_back( 1.008 );
		unit_cell_atom_charges.push_back( 0.425 );
	}

	//----- Use symmetry to construct a "full" unit cell -----//

	// Make a (+1/2*x, +1/2*y, z) cell
	std::vector<double> fractional_shift = { 0.5, 0.5, 0 };
	for ( int i=0; i<num_atoms_per_primitive_cell; ++i ) {
		// Allocate memory for new atom's coordinates
		int index = i + num_atoms_per_primitive_cell;
		fractional_unit_cell[index].resize(DIM_);

		for ( int j=0; j<DIM_; ++j ) {
			fractional_unit_cell[index][j] = fractional_unit_cell[i][j] + fractional_shift[j];
		}

		// Extend parameter arrays
		unit_cell_atom_names.push_back( unit_cell_atom_names[i] );
		unit_cell_atom_masses.push_back( unit_cell_atom_masses[i] );
		unit_cell_atom_charges.push_back( unit_cell_atom_charges[i] );
	}

	// Impose PBCs on the atoms outsize the (x,y,z) central unit cell
	std::vector<double> fractional_box_lengths = { 1.0, 1.0, 1.0 };
	for ( int i=0; i<num_atoms_per_unit_cell; ++i ) {
		for ( int j=0; j<3; ++j ) {
			keepInBox(fractional_box_lengths, fractional_unit_cell[i]);
		}
	}

	//----- Convert to orthonormal unit cell in Cartestian coordinates-----//

	// Volume of the unit cell parallelepiped [nm^3]
	double v_unit_cell = a*b*c*sqrt(1.0 - cos(alpha)*cos(alpha) 
	                                   - cos(beta)*cos(beta)
	                                   - cos(gamma)*cos(gamma) 
	                                   + 2.0*cos(alpha)*cos(beta)*cos(gamma));

	// Transformation matrix: [M^{-1}]^T (transpose of inverse of M)
	std::vector<std::vector<double>> inv_M(DIM_);
	inv_M[0] = { a,     b*cos(gamma),  c*cos(beta)                                      };
	inv_M[1] = { 0.0,   b*sin(gamma),  c*(cos(alpha) - cos(beta)*cos(gamma))/sin(gamma) };
	inv_M[2] = { 0.0,   0.0,           v_unit_cell/(a*b*sin(gamma))                      };

	// Parallelepiped unit cell in Cartesian coordinates
	std::vector<std::vector<double>> unit_cell(num_atoms_per_unit_cell);
	for ( int i=0; i<num_atoms_per_unit_cell; ++i ) {
		unit_cell[i].assign(DIM_, 0.0);

		for ( int j=0; j<DIM_; ++j ) { // x,y,z coordinates of each atom (rows of M^{-1})
			for ( int k=0; k<DIM_; ++k ) {
				unit_cell[i][j] += inv_M[j][k]*fractional_unit_cell[i][k];
			}
		}
	}

	// Box lengths of the orthogonal unit cell are along the diagonal of inv_M
	std::vector<double> unit_cell_box_lengths(DIM_); // [nm]
	for ( int i=0; i<DIM_; ++i ) {
		unit_cell_box_lengths[i] = inv_M[i][i];
	}

	// Adjust the atom positions so that they lie within the orthogonal unit cell
	for ( int i=0; i<num_atoms_per_unit_cell; ++i ) {
		keepInBox(unit_cell_box_lengths, unit_cell[i]);
	}

	// Shift the whole box by (-x_H, -y_H, 0.0) (H --> i=13)
	// so that O-H bonds don't cross the boundary (purely for convenience)
	std::vector<double> x_shift(DIM_);
	int shiftAtom = 13;
	double one_plus_eps = 1.01;
	x_shift[X_DIM] = -1.0*one_plus_eps*unit_cell[shiftAtom][X_DIM];
	x_shift[Y_DIM] = -1.0*one_plus_eps*unit_cell[shiftAtom][Y_DIM];
	x_shift[Z_DIM] = 0.0;   // no need to shift along z

	for ( int i=0; i<num_atoms_per_unit_cell; ++i ) {
		for ( int j=0; j<DIM_; ++j ) {
			unit_cell[i][j] += x_shift[j];
		}

		// Ensure everything stays in the simulation box
		keepInBox(unit_cell_box_lengths, unit_cell[i]);
	}

	//----- Dipole of the unit cell -----//

	std::vector<double> unit_cell_dipole(DIM_, 0.0);

	for ( int i=0; i<num_atoms_per_unit_cell; ++i ) {
		for ( int j=0; j<DIM_; ++j ) {
			unit_cell_dipole[j] += unit_cell_atom_charges[i]*unit_cell[i][j];
		}
	}

	// Convert to Debyes
	const double E_NM_PER_DEBYE = 0.0208194; // [(e*nm)/Debye]
	for ( int j=0; j<DIM_; ++j ) { 
		unit_cell_dipole[j] /= E_NM_PER_DEBYE; 
	}

	//----- Make a layer of kaolinite  -----//

	// Number of unit cells along x, y, and z
	//std::vector<int> unit_cell_grid = { 11, 7, 3 };
	//unit_cell_grid = { 1, 1, 1 };
	//unit_cell_grid = { 4, 3, 2 };

	// Slab dimensions [nm]
	std::vector<double> slab_lengths(DIM_);
	for ( int j=0; j<DIM_; ++j ) {
		slab_lengths[j] = unit_cell_grid[j]*unit_cell_box_lengths[j];
	}

	// Make the "normal" slab
	int num_unit_cells_per_slab = 1; 
	for ( int d=0; d<DIM_; ++d ) {
		num_unit_cells_per_slab *= unit_cell_grid[d];
	}
	int num_atoms_per_slab = num_atoms_per_unit_cell*num_unit_cells_per_slab;

	std::vector<std::vector<double>> slab(num_atoms_per_slab, std::vector<double>(DIM_));
	std::vector<double> offset(DIM_); // Position of bottom-left corner of cell (l,m,n)

	int atom_counter = 0;
	for ( int l=0; l<unit_cell_grid[X_DIM]; ++l ) { // cells along x
		offset[X_DIM] = l*unit_cell_box_lengths[X_DIM];
		for ( int m=0; m<unit_cell_grid[Y_DIM]; ++m ) { // cells along y
			offset[Y_DIM] = m*unit_cell_box_lengths[Y_DIM];
			for ( int n=0; n<unit_cell_grid[Z_DIM]; ++n ) { // cells along z
				offset[Z_DIM] = n*unit_cell_box_lengths[Z_DIM];

				// Place atoms for this unit cell
				for ( int i=0; i<num_atoms_per_unit_cell; ++i ) {
					for ( int j=0; j<DIM_; ++j ) {
						slab[atom_counter][j] = unit_cell[i][j] + offset[j];
					}
					++atom_counter;
				}
			}
		}
	}

	//----- Make a slab reflected across the xy-plane -----//

	std::vector<std::vector<double>> reflected_slab(num_atoms_per_slab, std::vector<double>(DIM_));

	for ( int i=0; i<num_atoms_per_slab; ++i ) {
		// x and y are unchanged
		for ( int j=0; j<Z_DIM; ++j ) {
			reflected_slab[i][j] = slab[i][j];
		}

		// Reflect across xy-plane, and offset z so that all 
		// coordinates are positive
		reflected_slab[i][Z_DIM] = slab_lengths[Z_DIM] - slab[i][Z_DIM];
	}

	//----- Write a GRO file for the normal slab -----//
	
	std::string gro_file("kaolinite_slab.gro");
	std::ofstream ofs(gro_file);

	ofs << "Kaolinite slab" << "\n";
	ofs << "    " << num_atoms_per_slab << "\n";

	atom_counter = 0;
	int molecule_counter = 0;
	for ( int i=0; i<num_unit_cells_per_slab; ++i ) {
		// Each unit cell:
		// 2x Al, 2x Si, 5x O (standalone), 4x O (hydroxyl), 4x H (hydroxyl)

		for ( int n=0; n<num_atoms_per_unit_cell; ++n ) {
			// Residue/Molecule number
			ofs << std::right << std::setw(5) << molecule_counter + 1;

			// Residue/Molecule name
			ofs << std::left << std::setw(5) << "KAO";

			// Atom name
			ofs << std::right << std::setw(5) << unit_cell_atom_names[n];

			// Atom number
			ofs << std::right << std::setw(5) << atom_counter + 1;

			// (x, y, z) in nm
			for ( int j=0; j<3; ++j ) {
				// Right-justified fixed-point number with 3 decimal places (total width: 8)
				ofs << std::right << std::setw(8) << std::fixed << std::setprecision(3)
				    << slab[atom_counter][j];
			}

			// Done line
			ofs << "\n";

			++atom_counter;
		}
	}

	// Box lengths
	for ( int j=0; j<3; ++j ) {
		ofs << "   " << slab_lengths[j];
	}
	ofs << "\n";

	// Done
	ofs.close();


	//----- Write a GRO file for the reflected slab -----//
	
	ofs.open("kaolinite_reflected_slab.gro");

	ofs << "Kaolinite slab reflected across xy-plane" << "\n";
	ofs << "    " << num_atoms_per_slab << "\n";

	atom_counter = 0;
	molecule_counter = 0;
	for ( int i=0; i<num_unit_cells_per_slab; ++i ) {
		// Each unit cell:
		// 2x Al, 2x Si, 5x O (standalone), 4x O (hydroxyl), 4x H (hydroxyl)

		for ( int n=0; n<num_atoms_per_unit_cell; ++n ) {
			// Residue/Molecule number
			ofs << std::right << std::setw(5) << molecule_counter + 1;

			// Residue/Molecule name
			ofs << std::left << std::setw(5) << "KAO";

			// Atom name
			ofs << std::right << std::setw(5) << unit_cell_atom_names[n];

			// Atom number
			ofs << std::right << std::setw(5) << atom_counter + 1;

			// (x, y, z) in nm
			for ( int j=0; j<DIM_; ++j ) {
				// Right-justified fixed-point number with 3 decimal places (total width: 8)
				ofs << std::right << std::setw(8) << std::fixed << std::setprecision(3)
				    << reflected_slab[atom_counter][j];
			}

			// Done line
			ofs << "\n";

			++atom_counter;
		}
	}

	// Box lengths
	for ( int j=0; j<DIM_; ++j ) {
		ofs << "   " << slab_lengths[j];
	}
	ofs << "\n";

	// Done
	ofs.close();


	//----- Write a .itp file with the "molecule" information -----//

	std::string itp_file("kaolinite_slab.itp");
	ofs.open( itp_file );

	ofs << "; Topology file fragment for KAO molecule" << "\n";
	ofs << "; " << "\n";
	ofs << "; ifdef Flags" << "\n";
	ofs << ";     AL_SURFACE_FLEXIBLE_OH: O-H bonds at the Al surface are flexible." << "\n";
	ofs << ";     AL_SURFACE_FREE_H: free hydroxyl hydrogens at Al surface" << "\n";
	ofs << ";     AL_SURFACE_FREE_OH: free hydroxyl groups at Al surface" << "\n";
	ofs << "; " << "\n";
	ofs << "; NOTE: For fully rigid kaolinite, use a freeze group" << "\n";
	ofs << ";" << "\n";
	ofs << "; Dipole of unit cell: " << "\n"; 

	ofs << ";     mu = {";
	for ( int j=0; j<DIM_; ++j ) {
		ofs << std::setw(8) << std::fixed << std::setprecision(3) 
		    << unit_cell_dipole[j];

		if ( j != Z_DIM ) { ofs << ","; }
	}
	ofs << "} [Debye]" << "\n";

	ofs << ";    |mu| = " << norm(unit_cell_dipole) << " [Debye]" << "\n";
	ofs << ";" << "\n";
	ofs << "\n";

	// "Molecule" type declaration: kaolinite (KAO)
	ofs << "[ moleculetype ]" << "\n";
	ofs << "; molname      nrexcl" << "\n";
	ofs << "  KAO          3" << "\n";
	ofs << "\n";

	// Atoms in molecule
	ofs << "[ atoms ]" << "\n";
	ofs << ";  nr   type  resnr   residue  name    cgnr    charge      mass" << "\n";

	atom_counter = 0;
	int charge_group_index = 0;
	for ( int i=0; i<num_unit_cells_per_slab; ++i ) {
		// Each primitive unit cell: 17 atoms
		//   2x Al, 2x Si, 5x O (bridging), 4x O (hydroxyl), 4x H (hydroxyl)
		// Each unit cell: 2x as many atoms --> 34 total

		for ( int n=0; n<num_atoms_per_unit_cell; ++n ) {
			// Use charge groups of about 4 atoms
			// - Unit cell: 34 atoms
			// - Max. charge group size: 32 atoms
			if ( ((n % 4) == 0) && (atom_counter > 0) ) {
				++charge_group_index;
			}

			// Indent for readability (2 spaces)
			ofs << "  ";

			// Atom index
			ofs << std::setw(6) << (atom_counter + 1) << " ";

			// Atom type
			ofs << std::left << std::setw(6) << unit_cell_atom_names[n] << " ";

			// Residue/Molecule number
			ofs << std::setw(6) << (molecule_counter + 1) << " ";

			// Residue/Molecule name
			ofs << std::left << std::setw(5) << "KAO" << " ";

			// Atom name 
			ofs << std::left << std::setw(7) << unit_cell_atom_names[n] << " ";

			// Charge group index
			ofs << std::setw(6) << (atom_counter + 1) << " ";

			// Charge
			ofs << std::right << std::setw(8) << std::fixed << std::setprecision(5)
			    << unit_cell_atom_charges[n] << " ";

			// Mass
			ofs << std::right << std::setw(8) << std::fixed << std::setprecision(3)
			    << unit_cell_atom_masses[n] << " ";

			// Done line
			ofs << "\n";

			++atom_counter;
		}
	}
	ofs << "\n";

	// Bond potentials

	// - Parameters
	double k_bond = 231850;  // [kJ/mol*nm^2]
	double r_0    = 0.1;     // [nm]

	// - Functional forms
	int    bond_function       = 1; // Functional form for harmonic bond
	int    constraint_function = 1; // Constraints will generate nonbonded exclusions

	// - Relative indices of O-H bond pairs
	std::vector<int> rel_indices_O = {  9, 10, 11, 12, 26, 27, 28, 29 };
	std::vector<int> rel_indices_H = { 13, 14, 15, 33, 30, 31, 32, 16 };
	int num_OH_bonds_per_unit_cell = rel_indices_O.size();

	// - Bond section header for flexible O-H bonds
	ofs << "#ifdef AL_SURFACE_FLEXIBLE_OH" << "\n";
	ofs << "[ bonds ]" << "\n";
	ofs << "; ai    aj    funct     b0=r0(nm)    k(kJ/mol/nm^2)" << "\n";

	for ( int i=0; i<num_unit_cells_per_slab; ++i ) {
		int atom_index_offset = num_atoms_per_unit_cell*i;
		int index_O, index_H; 

		for ( int j=0; j<num_OH_bonds_per_unit_cell; ++j ) {
			index_O = atom_index_offset + rel_indices_O[j];
			index_H = atom_index_offset + rel_indices_H[j];

			ofs << " " << std::setw(5) << (index_O + 1)
			    << " " << std::setw(5) << (index_H + 1)
			    << " " << std::setw(3) << bond_function
			    << " " << std::setw(8) << r_0
			    << " " << std::setw(8) << k_bond
			    << "\n";
		}
	}
	ofs << "\n";

	// - Consraints for rigid O-H bonds
	ofs << "; Else O-H bonds are rigid" << "\n";
	ofs << "#else" << "\n";
	ofs << "[ constraints ]" << "\n";
	ofs << "; ai    aj    funct     b0=r0" << "\n";
	for ( int i=0; i<num_unit_cells_per_slab; ++i ) {
		int atom_index_offset = num_atoms_per_unit_cell*i;
		int index_O, index_H; 

		for ( int j=0; j<num_OH_bonds_per_unit_cell; ++j ) {
			index_O = atom_index_offset + rel_indices_O[j];
			index_H = atom_index_offset + rel_indices_H[j];

			ofs << " " << std::setw(5) << (index_O + 1)
			    << " " << std::setw(5) << (index_H + 1)
			    << " " << std::setw(3) << constraint_function
			    << " " << std::setw(8) << r_0
			    << "\n";
		}
	}

	ofs << "\n";
	ofs << "; End #ifdef AL_SURFACE_FLEXIBLE_OH" << "\n";
	ofs << "#endif" << "\n";
	ofs << "\n";

	// Bond angle potentials

	ofs << "[ angles ]" << "\n";
	ofs << "; ai   aj   ak   funct   theta0    k_bend" << "\n";

	// - Parameters
	double k_bend        = 251.04; // [kJ/mol*rad^2]
	double theta_0       = 109.47; // [degrees]
	int    angleFunction = 1;      // Functional form for harmonic bend 

	Int3 grid_indices, grid_indices_minus_x, grid_indices_plus_y;
	int  cell_index, cell_minus_x, cell_plus_y;
	int  atom_index_offset;

	for ( int l=0; l<unit_cell_grid[X_DIM]; ++l ) {
		for ( int m=0; m<unit_cell_grid[Y_DIM]; ++m ) {
			for ( int n=0; n<unit_cell_grid[Z_DIM]; ++n ) {
				// Current unit cell and index of its first atom
				grid_indices = {{ l, m, n }};
				cell_index = getLinearCellIndex(grid_indices, unit_cell_grid);
				atom_index_offset = num_atoms_per_unit_cell*cell_index;

				// Index of the "-x" unit cell
				grid_indices_minus_x = grid_indices;
				if ( grid_indices[X_DIM] > 0 ) {
					grid_indices_minus_x[X_DIM] -= 1;
				}
				else {
					// crossing the x=0 boundary of the box
					grid_indices_minus_x[X_DIM] = unit_cell_grid[X_DIM] - 1;
				}
				cell_minus_x = getLinearCellIndex(grid_indices_minus_x, unit_cell_grid);

				// Index of the "+y" unit cell
				grid_indices_plus_y = grid_indices;
				if ( grid_indices[Y_DIM] < unit_cell_grid[Y_DIM] - 1 ) {
					grid_indices_plus_y[Y_DIM] += 1;
				}
				else {
					// crossing the upper y-boundary of the box
					grid_indices_plus_y[Y_DIM] = 0;
				}
				cell_plus_y = getLinearCellIndex(grid_indices_plus_y, unit_cell_grid);

				// Note: Atom indexing in parantheses is derived from examining the 
				// first unit cell
				int index_M, index_O, index_H;

				//----- Al(0) -----//

				std::vector<std::array<int,4>> bond_angles = {{
					// Key:  Al_atom_offset,  cell_with_OH,  O_atom_offset,  H_atom_offset
					// Al(0)
					{{  0, cell_minus_x, 11, 15 }}, // OH(11), HH(15) [-x cell]
					{{  0, cell_index,   26, 30 }}, // OH(26), HH(30)
					{{  0, cell_index,   27, 31 }}, // OH(27), HH(31)
					{{  0, cell_index,   29, 16 }}, // OH(29), HH(16)
					// Al(1)
					{{  1, cell_index,   10, 14 }}, // OH(10), HH(14)
					{{  1, cell_index,   11, 15 }}, // OH(11), HH(15)
					{{  1, cell_index,   26, 30 }}, // OH(26), HH(30)
					{{  1, cell_index,   29, 16 }}, // OH(29), HH(16)
					// Al(17)
					{{ 17, cell_index,    9, 13 }}, // OH(9),  HH(13)
					{{ 17, cell_plus_y,  10, 14 }}, // OH(10), HH(14) [+y cell]
					{{ 17, cell_index,   12, 33 }}, // OH(12), HH(33)
					{{ 17, cell_index,   28, 32 }}, // OH(28), HH(32)
					// Al(18)
					{{ 18, cell_minus_x,  9, 13 }}, // OH(9),  HH(13) [-x cell]
					{{ 18, cell_minus_x, 12, 33 }}, // OH(12), HH(33) [-x cell]
					{{ 18, cell_index,   27, 31 }}, // OH(27), HH(31)
					{{ 18, cell_index,   28, 32 }}  // OH(28), HH(32)
				}};

				int num_bond_angles = bond_angles.size();
				for ( int k=0; k<num_bond_angles; ++k ) {
					// Metal atom (M)
					index_M = atom_index_offset + bond_angles[k][0];

					// Hydroxyl group (not necessarily in the same cell!)
					index_O = num_atoms_per_unit_cell*bond_angles[k][1] + bond_angles[k][2];
					index_H = num_atoms_per_unit_cell*bond_angles[k][1] + bond_angles[k][3];

					ofs << " " << std::setw(5) << (index_M + 1)
							<< " " << std::setw(5) << (index_O + 1)
					    << " " << std::setw(5) << (index_H + 1)
							<< " " << std::setw(3) << angleFunction
							<< " " << std::setw(8) << theta_0
					    << " " << std::setw(8) << k_bend
							<< "\n";
				}
			}
		}
	}
	ofs << "\n";

	ofs.close();


	//----- Index file for partially-frozen surfaces -----//

	// Indices (in the unit cell) of the Al surface's O and H atoms
	// - Note: O(9)--H(13) and O(26)--H(30) are surface hydroxyl groups
	std::vector<int> unit_cell_surface_O = { 10, 11, 12, 27, 28, 29 };
	std::vector<int> unit_cell_surface_H = { 14, 15, 33, 31, 32, 16 };   // corresponding H
	std::sort( unit_cell_surface_H.begin(), unit_cell_surface_H.end() ); // sort for convenience

	// Arrays of frozen atom indices
	std::vector<int> frozen_atoms_free_H;
	std::vector<int> frozen_atoms_free_OH;

	frozen_atoms_free_H.reserve(num_atoms_per_slab);
	frozen_atoms_free_OH.reserve(num_atoms_per_slab);

	// Number of surface hydroxyl groups per unit cell
	int num_surface_hydroxyl_groups = unit_cell_surface_O.size();

	// Index of topmost layer (along z)
	int top_layer = unit_cell_grid[2] - 1;

	// "Free H" surface: only hydroxyl H at the Al surface can move
	bool is_frozen;
	for ( int l=0; l<unit_cell_grid[0]; ++l ) {
		for ( int m=0; m<unit_cell_grid[1]; ++m ) {
			for ( int n=0; n<unit_cell_grid[2]; ++n ) {
				// Current unit cell and index of its first atom
				int cell_index = getLinearCellIndex( Int3{{l, m, n}}, unit_cell_grid);
				int atom_index_offset = num_atoms_per_unit_cell*cell_index;
				
				for ( int i=0; i<num_atoms_per_unit_cell; ++i ) {
					// Reset the flag
					is_frozen = true;

					// If in top layer, might be a free H
					if ( n == top_layer ) {
						for ( int j=0; j<num_surface_hydroxyl_groups; ++j ) {
							if ( i == unit_cell_surface_H[j] ) {
								// Atom is a surface H: don't freeze it
								is_frozen = false;
							}
						}
					}

					if ( is_frozen ) {
						// Store index
						int atom_index = atom_index_offset + i;
						frozen_atoms_free_H.push_back( atom_index );
					}
				}
			}
		}
	}

	// Free OH surface: only hydroxyl groups at exposed Al surface free 
	for ( int l=0; l<unit_cell_grid[0]; ++l ) {
		for ( int m=0; m<unit_cell_grid[1]; ++m ) {
			for ( int n=0; n<unit_cell_grid[2]; ++n ) {
				// Current unit cell and index of its first atom
				int cell_index = getLinearCellIndex( Int3{{l, m, n}}, unit_cell_grid);
				int atom_index_offset = num_atoms_per_unit_cell*cell_index;
				
				for ( int i=0; i<num_atoms_per_unit_cell; ++i ) {
					// Reset the flag
					is_frozen = true;

					// If in top layer, might be an H
					if ( n == top_layer ) {
						for ( int j=0; j<num_surface_hydroxyl_groups; ++j ) {
							if ( i == unit_cell_surface_H[j] ) {
								// Atom is a surface H: don't freeze it
								is_frozen = false;
							}
							else if ( i == unit_cell_surface_O[j] ) {
								// Atom is a surface O: don't freeze it
								is_frozen = false;
							}
						}
					}

					if ( is_frozen ) {
						// Store index
						int atom_index = atom_index_offset + i;
						frozen_atoms_free_OH.push_back( atom_index );
					}
				}
			}
		}
	}

	//----- Write a freeze group file for a single kaolinite slab -----//

	ofs.open("freeze_groups_monolayer.ndx");

	// Multiple indices per line for compact file
	int line_counter = 0;
	int num_atoms_per_line = 5;

	// Free H surface
	int num_frozen_atoms_free_H = frozen_atoms_free_H.size();

	ofs << "[ KAO_FROZEN_ATOMS_FREE_H ]" << "\n";

	for ( int i=0; i<num_frozen_atoms_free_H; ++i ) {
		ofs << std::setw(8) << frozen_atoms_free_H[i] << "\t";
		if ( line_counter < num_atoms_per_line ) {
			++line_counter;
		}
		else {
			// Done line
			ofs << "\n";
			line_counter = 0;
		}
	}
	ofs << "\n";

	// Free OH surface
	int num_frozen_atoms_free_OH = frozen_atoms_free_OH.size();

	ofs << "[ KAO_FROZEN_ATOMS_FREE_OH ]" << "\n";

	line_counter = 0;

	for ( int i=0; i<num_frozen_atoms_free_OH; ++i ) {
		ofs << std::setw(8) << frozen_atoms_free_OH[i] << "\t";
		if ( line_counter < num_atoms_per_line ) {
			++line_counter;
		}
		else {
			// Done line
			ofs << "\n";
			line_counter = 0;
		}
	}
	ofs << "\n";

	ofs.close();

	//----- Write a second freeze group file for a kaolinite bilayer -----//
	
	// Assumes layers (slabs) are stored consecutively
	std::vector<int> layer_atom_offsets = { 0, num_atoms_per_slab };
	int num_layers = layer_atom_offsets.size();

	ofs.open("freeze_groups_bilayer.ndx");

	// Free H surface
	ofs << "[ KAO_FROZEN_ATOMS_FREE_H ]" << "\n";

	for ( int l=0; l<num_layers; ++l ) {
		for ( int i=0; i<num_frozen_atoms_free_H; ++i ) {
			ofs << std::setw(8) 
			    << (frozen_atoms_free_H[i] + layer_atom_offsets[l] + 1) 
			    << "\t";

			if ( line_counter < num_atoms_per_line ) {
				++line_counter;
			}
			else {
				// Done line
				ofs << "\n";
				line_counter = 0;
			}
		}
	}
	ofs << "\n";

	// Free OH surface

	ofs << "[ KAO_FROZEN_ATOMS_FREE_OH ]" << "\n";

	line_counter = 0;

	for ( int l=0; l<num_layers; ++l ) {
		for ( int i=0; i<num_frozen_atoms_free_OH; ++i ) {
			ofs << std::setw(8) 
			    << (frozen_atoms_free_OH[i] + layer_atom_offsets[l] + 1) 
			    << "\t";

			if ( line_counter < num_atoms_per_line ) {
				++line_counter;
			}
			else {
				// Done line
				ofs << "\n";
				line_counter = 0;
			}
		}
	}
	ofs << "\n";

	ofs.close();

}



// Applies the minimum image convention
void minImage(const Rvec x1, const Rvec x2, const Rvec boxL, Rvec x12, double& dist_sq)
{
	dist_sq = 0.0;

	for ( int d=0; d<3; d++ ) {
		x12[d] = x2[d] - x1[d];
		// Apply minimum image convention
		if      ( x12[d] >  0.5*boxL[d] ) { x12[d] -= boxL[d]; }
		else if ( x12[d] < -0.5*boxL[d] ) { x12[d] += boxL[d]; }

		dist_sq += x12[d]*x12[d];
	}

	return;
}



void minImage(
	const std::vector<double>& x1, const std::vector<double>& x2, 
	const std::vector<double>& boxL,
	// Output
	std::vector<double>& x12, double& dist_sq)
{
	dist_sq = 0.0;

	for ( int d=0; d<3; d++ ) {
		x12[d] = x2[d] - x1[d];
		// Apply minimum image convention
		if      ( x12[d] >  0.5*boxL[d] ) { x12[d] -= boxL[d]; }
		else if ( x12[d] < -0.5*boxL[d] ) { x12[d] += boxL[d]; }

		dist_sq += x12[d]*x12[d];
	}

	return;
}



// Keeps the atom in the simulaton box by applying PBCs
void keepInBox(const Rvec boxL, Rvec x)
{
	for ( int d=0; d<3; d++ ) {
		// Apply minimum image convention
		if      ( x[d] > boxL[d] ) { x[d] -= boxL[d]; }
		else if ( x[d] < 0.0 )     { x[d] += boxL[d]; }
	}

	return;
}



void keepInBox(const std::vector<double>& boxL, std::vector<double>& x)
{
	for ( int d=0; d<3; d++ ) {
		// Apply minimum image convention
		if      ( x[d] > boxL[d] ) { x[d] -= boxL[d]; }
		else if ( x[d] < 0.0 )     { x[d] += boxL[d]; }
	}

	return;
}



// Checks whether each oxygen has 2 hydrogens
bool areHydrogensCorrectlyPlaced(std::vector<int> hydrogenCounts)
{
	int numWaters = hydrogenCounts.size();
	for ( int i=0; i<numWaters; i++ ) {
		if ( hydrogenCounts[i] != 2 ) {
			return false;
		}
	}
	return true;
}



// Vector norm
double norm(const Rvec x) {
	return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}



double norm(const std::vector<double>& x) {
	return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}



// Use the dot product to get the angle between the vectors
double angleBetweenVectors(const Rvec a, const Rvec b)
{
	
	double norm_a = 0.0, norm_b = 0.0, a_dot_b = 0.0;
	for ( int i=0; i<3; ++i ) {
		norm_a  += a[i]*a[i];
		norm_b  += b[i]*b[i];
		a_dot_b += a[i]*b[i];
	}
	norm_a = sqrt(norm_a);
	norm_b = sqrt(norm_b);

	double theta = acos(a_dot_b/(norm_a*norm_b));

	return theta;
}

