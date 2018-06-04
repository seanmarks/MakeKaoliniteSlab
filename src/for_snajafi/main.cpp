// GenKaoliniteSlab.exe
//
// Updated Feb. 5, 2018
//

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

// Project headers
#include "main.h"

int main(int argc, char* argv[])
{
	//----- Input -----//

	// Number of unit cells along x, y, and z
	std::vector<int> unitCellGrid(3); 
	unitCellGrid = { 11, 7, 2 }; // Full grid
	//unitCellGrid = { 1, 1, 1 }; // single unit cell
	//unitCellGrid = { 4, 3, 2 }; // small grid


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
	std::vector<double> latticeConstants = { a, b, b };

	//----- Fractional unit cell -----//

	// Atoms are stored in the following order:
	// 2x Al, 2x Si, 5x O (standalone), 4x O (hydroxyl), 4x H (hydroxyl)
	int numAtomsPerPrimitiveCell = 17;
	int numAtomsPerUnitCell = 2*numAtomsPerPrimitiveCell;
	std::vector<std::vector<double>> fractionalUnitCell(numAtomsPerUnitCell);

	fractionalUnitCell[0]  = { 0.28900, 0.49660, 0.46600 }; // Al1 (octahedral)
	fractionalUnitCell[1]  = { 0.79300, 0.32880, 0.46500 }; // Al2 (octahedral)
	fractionalUnitCell[2]  = { 0.98900, 0.33950, 0.09060 }; // Si1 (tetrahedral)
	fractionalUnitCell[3]  = { 0.50700, 0.16650, 0.09380 }; // Si2 (tetrahedral)

	fractionalUnitCell[4]  = { 0.04900, 0.34820, 0.31680 }; // O1
	fractionalUnitCell[5]  = { 0.11300, 0.65990, 0.31880 }; // O2
	fractionalUnitCell[6]  = { 0.00000, 0.50000, 0.00000 }; // O3
	fractionalUnitCell[7]  = { 0.20400, 0.22910, 0.03000 }; // O4
	fractionalUnitCell[8]  = { 0.19700, 0.76410, 0.00100 }; // O5

	fractionalUnitCell[9]  = { 0.05000, 0.97100, 0.32500 }; // O-hydroxyl-1
	fractionalUnitCell[10] = { 0.96000, 0.16580, 0.60700 }; // O-hydroxyl-2
	fractionalUnitCell[11] = { 0.03700, 0.47260, 0.60460 }; // O-hydroxyl-3
	fractionalUnitCell[12] = { 0.03800, 0.85820, 0.60900 }; // O-hydroxyl-4

	fractionalUnitCell[13] = { 0.14500, 0.06510, 0.32600 }; // H-hydroxyl-1
	fractionalUnitCell[14] = { 0.06300, 0.16380, 0.73900 }; // H-hydroxyl-2
	fractionalUnitCell[15] = { 0.03600, 0.50570, 0.73200 }; // H-hydroxyl-3
	fractionalUnitCell[16] = { 0.53400, 0.31540, 0.72800 }; // H-hydroxyl-4

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
	std::vector<std::string> unitCellAtomNames = { "AlOK", "AlOK", "SiTK", "SiTK" };
	std::vector<double> unitCellAtomMasses = { 26.982, 26.982, 28.086, 28.086 };
	std::vector<double> unitCellAtomCharges = { 1.575, 1.575, 2.1, 2.1 };

	// OBK: bridging O
	// - First 2 are "interior" (apical?)
	// - Other three are basal (form hexagonal rings with Si atoms)
	for ( int i=0; i<5; ++i )
	{
		unitCellAtomNames.push_back( std::string("OBK") );
		unitCellAtomMasses.push_back( 16.000 );
		unitCellAtomCharges.push_back( -1.05 );
	}

	// OHK: hydroxyl O coordinate to Al atoms
	for ( int i=0; i<4; ++i )
	{
		unitCellAtomNames.push_back( std::string("OHK") );
		unitCellAtomMasses.push_back( 16.000 );
		unitCellAtomCharges.push_back( -0.95 );
	}

	// HHK: hydroxyl hydrogens
	for ( int i=0; i<4; ++i )
	{
		unitCellAtomNames.push_back( std::string("HHK") );
		unitCellAtomMasses.push_back( 1.008 );
		unitCellAtomCharges.push_back( 0.425 );
	}

	//----- Use symmetry to construct a "full" unit cell -----//

	// Make a (+1/2*x, +1/2*y, z) cell
	std::vector<double> fractionalShift = { 0.5, 0.5, 0 };
	for ( int i=0; i<numAtomsPerPrimitiveCell; ++i )
	{
		// Allocate memory for new atom's coordinates
		int index = i + numAtomsPerPrimitiveCell;
		fractionalUnitCell[index].resize(3);

		for ( int j=0; j<3; ++j )
		{
			fractionalUnitCell[index][j] = fractionalUnitCell[i][j] + fractionalShift[j];
		}

		// Extend parameter arrays
		unitCellAtomNames.push_back( unitCellAtomNames[i] );
		unitCellAtomMasses.push_back( unitCellAtomMasses[i] );
		unitCellAtomCharges.push_back( unitCellAtomCharges[i] );
	}

	// Impose PBCs on the atoms outsize the (x,y,z) central unit cell
	std::vector<double> fractionalBoxL = { 1.0, 1.0, 1.0 };
	for ( int i=0; i<numAtomsPerUnitCell; ++i )
	{
		for ( int j=0; j<3; ++j )
		{
			keepInBox(fractionalBoxL, fractionalUnitCell[i]);
		}
	}

	//----- Convert to orthonormal unit cell in Cartestian coordinates-----//

	// Volume of the unit cell parallelepiped [nm^3]
	double v_unitCell = a*b*c*sqrt(1.0 - cos(alpha)*cos(alpha) 
	                                   - cos(beta)*cos(beta)
	                                   - cos(gamma)*cos(gamma) 
	                                   + 2.0*cos(alpha)*cos(beta)*cos(gamma));

	// Transformation matrix: [M^{-1}]^T (transpose of inverse of M)
	std::vector<std::vector<double>> invM(3);
	invM[0] = { a,     b*cos(gamma),  c*cos(beta)                                      };
	invM[1] = { 0.0,   b*sin(gamma),  c*(cos(alpha) - cos(beta)*cos(gamma))/sin(gamma) };
	invM[2] = { 0.0,   0.0,           v_unitCell/(a*b*sin(gamma))                      };

	// Parallelepiped unit cell in Cartesian coordinates
	std::vector<std::vector<double>> unitCell(numAtomsPerUnitCell);
	for ( int i=0; i<numAtomsPerUnitCell; ++i )
	{
		unitCell[i].assign(3, 0.0);

		for ( int j=0; j<3; ++j ) // x,y,z coordinates of each atom (rows of M^{-1})
		{
			for ( int k=0; k<3; ++k )
			{
				unitCell[i][j] += invM[j][k]*fractionalUnitCell[i][k];
			}
		}
	}

	// Box lengths of the orthogonal unit cell are along the diagonal of invM
	std::vector<double> unitCellBoxL(3); // [nm]
	for ( int i=0; i<3; ++i )
	{
		unitCellBoxL[i] = invM[i][i];
	}

	// Adjust the atom positions so that they lie within the orthogonal unit cell
	for ( int i=0; i<numAtomsPerUnitCell; ++i )
	{
		keepInBox(unitCellBoxL, unitCell[i]);
	}

	// Shift the whole box by (-x_H, -y_H, 0.0) (H --> i=13)
	// so that O-H bonds don't cross the boundary (purely for convenience)
	std::vector<double> x_shift(3);
	int shiftAtom = 13;
	double one_plus_eps = 1.01;
	x_shift[0] = -1.0*one_plus_eps*unitCell[shiftAtom][0];
	x_shift[1] = -1.0*one_plus_eps*unitCell[shiftAtom][1];
	x_shift[2] = 0.0;   // no need to shift along z

	for ( int i=0; i<numAtomsPerUnitCell; ++i )
	{
		for ( int j=0; j<3; ++j )
		{
			unitCell[i][j] += x_shift[j];
		}

		// Ensure everything stays in the simulation box
		keepInBox(unitCellBoxL, unitCell[i]);
	}

	//----- Dipole of the unit cell -----//

	std::vector<double> unitCellDipole(3,0.0);

	for ( int i=0; i<numAtomsPerUnitCell; ++i )
	{
		for ( int j=0; j<3; ++j )
		{
			unitCellDipole[j] += unitCellAtomCharges[i]*unitCell[i][j];
		}
	}

	// Convert to Debyes
	const double E_NM_PER_DEBYE = 0.0208194; // [(e*nm)/Debye]

	for ( int j=0; j<3; ++j ) { unitCellDipole[j] /= E_NM_PER_DEBYE; }

	//----- Make a layer of kaolinite  -----//

	// Number of unit cells along x, y, and z
	//std::vector<int> unitCellGrid = { 11, 7, 3 };
	//unitCellGrid = { 1, 1, 1 };
	//unitCellGrid = { 4, 3, 2 };

	// Slab dimensions [nm]
	std::vector<double> slabL(3);
	for ( int j=0; j<3; ++j )
	{
		slabL[j] = unitCellGrid[j]*unitCellBoxL[j];
	}

	// Make the "normal" slab
	int numUnitCellsPerSlab = unitCellGrid[0]*unitCellGrid[1]*unitCellGrid[2];
	int numAtomsPerSlab     = numAtomsPerUnitCell*numUnitCellsPerSlab;

	std::vector<std::vector<double>> slab(numAtomsPerSlab, std::vector<double>(3));
	std::vector<double>              offset(3); // Position of bottom-left corner of cell (l,m,n)

	int atomCounter = 0;
	for ( int l=0; l<unitCellGrid[0]; ++l ) // cells along x
	{
		offset[0] = l*unitCellBoxL[0];
		for ( int m=0; m<unitCellGrid[1]; ++m ) // cells along y
		{
			offset[1] = m*unitCellBoxL[1];
			for ( int n=0; n<unitCellGrid[2]; ++n ) // cells along z
			{
				offset[2] = n*unitCellBoxL[2];

				// Place atoms for this unit cell
				for ( int i=0; i<numAtomsPerUnitCell; ++i )
				{
					for ( int j=0; j<3; ++j )
					{
						slab[atomCounter][j] = unitCell[i][j] + offset[j];
					}
					++atomCounter;
				}
			}
		}
	}

	//----- Make a slab reflected across the xy-plane -----//

	std::vector<std::vector<double>> reflectedSlab(
			numAtomsPerSlab, std::vector<double>(3));

	for ( int i=0; i<numAtomsPerSlab; ++i )
	{
		// x and y are unchanged
		for ( int j=0; j<2; ++j )
		{
			reflectedSlab[i][j] = slab[i][j];
		}

		// Reflect across xy-plane, and offset z so that all 
		// coordinates are positive
		reflectedSlab[i][2] = slabL[2] - slab[i][2];
	}

	//----- Write a GRO file for the normal slab -----//
	
	std::string xyzFileName("kaolinite_slab.gro");
	std::ofstream ofs(xyzFileName);

	ofs << "Kaolinite slab" << "\n";
	ofs << "    " << numAtomsPerSlab << "\n";

	atomCounter = 0;
	int moleculeCounter = 0;
	for ( int i=0; i<numUnitCellsPerSlab; ++i )
	{
		// Each unit cell:
		// 2x Al, 2x Si, 5x O (standalone), 4x O (hydroxyl), 4x H (hydroxyl)

		for ( int n=0; n<numAtomsPerUnitCell; ++n )
		{
			// Residue/Molecule number
			ofs << std::right << std::setw(5) << moleculeCounter + 1;

			// Residue/Molecule name
			ofs << std::left << std::setw(5) << "KAO";

			// Atom name
			ofs << std::right << std::setw(5) << unitCellAtomNames[n];

			// Atom number
			ofs << std::right << std::setw(5) << atomCounter + 1;

			// (x, y, z) in nm
			for ( int j=0; j<3; ++j )
			{
				// Right-justified fixed-point number with 3 decimal places (total width: 8)
				ofs << std::right << std::setw(8) << std::fixed << std::setprecision(3)
				    << slab[atomCounter][j];
			}

			// Done line
			ofs << "\n";

			++atomCounter;
		}
	}

	// Box lengths
	for ( int j=0; j<3; ++j )
	{
		ofs << "   " << slabL[j];
	}
	ofs << "\n";

	// Done
	ofs.close();

	//----- Write a GRO file for the reflected slab -----//
	
	ofs.open("kaolinite_reflected_slab.gro");

	ofs << "Kaolinite slab reflected across xy-plane" << "\n";
	ofs << "    " << numAtomsPerSlab << "\n";

	atomCounter = 0;
	moleculeCounter = 0;
	for ( int i=0; i<numUnitCellsPerSlab; ++i )
	{
		// Each unit cell:
		// 2x Al, 2x Si, 5x O (standalone), 4x O (hydroxyl), 4x H (hydroxyl)

		for ( int n=0; n<numAtomsPerUnitCell; ++n )
		{
			// Residue/Molecule number
			ofs << std::right << std::setw(5) << moleculeCounter + 1;

			// Residue/Molecule name
			ofs << std::left << std::setw(5) << "KAO";

			// Atom name
			ofs << std::right << std::setw(5) << unitCellAtomNames[n];

			// Atom number
			ofs << std::right << std::setw(5) << atomCounter + 1;

			// (x, y, z) in nm
			for ( int j=0; j<3; ++j )
			{
				// Right-justified fixed-point number with 3 decimal places (total width: 8)
				ofs << std::right << std::setw(8) << std::fixed << std::setprecision(3)
				    << reflectedSlab[atomCounter][j];
			}

			// Done line
			ofs << "\n";

			++atomCounter;
		}
	}

	// Box lengths
	for ( int j=0; j<3; ++j )
	{
		ofs << "   " << slabL[j];
	}
	ofs << "\n";

	// Done
	ofs.close();

	//----- Write a .itp file with the "molecule" information -----//

	std::string itpFile("kaolinite_slab.itp");
	ofs.open( itpFile );

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
	for ( int j=0; j<3; ++j )
	{
		ofs << std::setw(8) << std::fixed << std::setprecision(3) 
		    << unitCellDipole[j];

		if ( j != 2 ) { ofs << ","; }
	}
	ofs << "} [Debye]" << "\n";

	ofs << ";    |mu| = " << norm(unitCellDipole) << " [Debye]" << "\n";
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

	atomCounter = 0;
	int chargeGroupIndex = 0;
	for ( int i=0; i<numUnitCellsPerSlab; ++i )
	{
		// Each primitive unit cell: 17 atoms
		//   2x Al, 2x Si, 5x O (bridging), 4x O (hydroxyl), 4x H (hydroxyl)
		// Each unit cell: 2x as many atoms --> 34 total

		for ( int n=0; n<numAtomsPerUnitCell; ++n )
		{
			// Use charge groups of about 4 atoms
			// - Unit cell: 34 atoms
			// - Max. charge group size: 32 atoms
			if ( ((n % 4) == 0) && (atomCounter > 0) )
			{
				++chargeGroupIndex;
			}

			// Indent for readability (2 spaces)
			ofs << "  ";

			// Atom index
			ofs << std::setw(6) << (atomCounter + 1) << " ";

			// Atom type
			ofs << std::left << std::setw(6) << unitCellAtomNames[n] << " ";

			// Residue/Molecule number
			ofs << std::setw(6) << (moleculeCounter + 1) << " ";

			// Residue/Molecule name
			ofs << std::left << std::setw(5) << "KAO" << " ";

			// Atom name 
			ofs << std::left << std::setw(7) << unitCellAtomNames[n] << " ";

			// Charge group index
			ofs << std::setw(6) << (atomCounter + 1) << " ";

			// Charge
			ofs << std::right << std::setw(8) << std::fixed << std::setprecision(5)
			    << unitCellAtomCharges[n] << " ";

			// Mass
			ofs << std::right << std::setw(8) << std::fixed << std::setprecision(3)
			    << unitCellAtomMasses[n] << " ";

			// Done line
			ofs << "\n";

			++atomCounter;
		}
	}
	ofs << "\n";

	// Bond potentials

	// - Parameters
	double k_bond = 231850;  // [kJ/mol*nm^2]
	double r_0    = 0.1;     // [nm]

	// - Functional forms
	int    bondFunction       = 1; // Functional form for harmonic bond
	int    constraintFunction = 1; // Constraints will generate nonbonded exclusions

	// - Relative indices of O-H bond pairs
	std::vector<int> rel_indices_O = {  9, 10, 11, 12, 26, 27, 28, 29 };
	std::vector<int> rel_indices_H = { 13, 14, 15, 33, 30, 31, 32, 16 };
	int num_OH_bonds_per_unit_cell = rel_indices_O.size();

	// - Bond section header for flexible O-H bonds
	ofs << "#ifdef AL_SURFACE_FLEXIBLE_OH" << "\n";
	ofs << "[ bonds ]" << "\n";
	ofs << "; ai    aj    funct     b0=r0      k" << "\n";

	for ( int i=0; i<numUnitCellsPerSlab; ++i )
	{
		int atomIndexOffset = numAtomsPerUnitCell*i;
		int index_O, index_H; 

		for ( int j=0; j<num_OH_bonds_per_unit_cell; ++j )
		{
			index_O = atomIndexOffset + rel_indices_O[j];
			index_H = atomIndexOffset + rel_indices_H[j];

			ofs << " " << std::setw(5) << (index_O + 1)
                << " " << std::setw(5) << (index_H + 1)
				<< " " << std::setw(3) << bondFunction
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
	for ( int i=0; i<numUnitCellsPerSlab; ++i )
	{
		int atomIndexOffset = numAtomsPerUnitCell*i;
		int index_O, index_H; 

		for ( int j=0; j<num_OH_bonds_per_unit_cell; ++j )
		{
			index_O = atomIndexOffset + rel_indices_O[j];
			index_H = atomIndexOffset + rel_indices_H[j];

			ofs << " " << std::setw(5) << (index_O + 1)
                << " " << std::setw(5) << (index_H + 1)
				<< " " << std::setw(3) << constraintFunction
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

	//int numUnitCellsPer_xy_Layer = unitCellGrid[0]*unitCellGrid[1]
	int numUnitCellsPer_yz_Layer = unitCellGrid[1]*unitCellGrid[2];

	for ( int l=0; l<unitCellGrid[0]; ++l )
	{
		for ( int m=0; m<unitCellGrid[1]; ++m )
		{
			for ( int n=0; n<unitCellGrid[2]; ++n )
			{
				// Current unit cell and index of its first atom
				int cellIndex       = l*unitCellGrid[1]*unitCellGrid[2] 
				                      + m*unitCellGrid[2] + n;
				int atomIndexOffset = numAtomsPerUnitCell*cellIndex;

				// Index of the "-x" unit cell
				int cell_minus_x = cellIndex - numUnitCellsPer_yz_Layer;
				if ( cell_minus_x < 0 ) // crossing the x=0 boundary of the box
				{ 
					cell_minus_x += numUnitCellsPerSlab;
				}

				// Index of the "+y" unit cell
				int cell_plus_y = cellIndex + unitCellGrid[2];
				if ( m == (unitCellGrid[1] - 1) ) // crossing +y boundary of box
				{
					cell_plus_y -= numUnitCellsPer_yz_Layer;
				}

				// Note: Atom indexing in parantheses is derived from examining the 
				// first unit cell
				int index_M, index_O, index_H;

				//----- Al(0) -----//

				index_M = atomIndexOffset + 0;

				// OH(11), HH(15) [-x cell]
				index_O = numAtomsPerUnitCell*cell_minus_x + 11;
				index_H = index_O + 4;
				ofs << " " << std::setw(5) << (index_M + 1)
				    << " " << std::setw(5) << (index_O + 1)
                    << " " << std::setw(5) << (index_H + 1)
				    << " " << std::setw(3) << angleFunction
				    << " " << std::setw(8) << theta_0
					<< " " << std::setw(8) << k_bend
				    << "\n";

				// OH(26), HH(30)
				index_O = atomIndexOffset + 26;
				index_H = index_O + 4;
				ofs << " " << std::setw(5) << (index_M + 1)
				    << " " << std::setw(5) << (index_O + 1)
                    << " " << std::setw(5) << (index_H + 1)
				    << " " << std::setw(3) << angleFunction 
				    << " " << std::setw(8) << theta_0
					<< " " << std::setw(8) << k_bend
				    << "\n";
				
				// OH(27), HH(31)
				index_O = atomIndexOffset + 27;
				index_H = index_O + 4;
				ofs << " " << std::setw(5) << (index_M + 1)
				    << " " << std::setw(5) << (index_O + 1)
                    << " " << std::setw(5) << (index_H + 1)
				    << " " << std::setw(3) << angleFunction
				    << " " << std::setw(8) << theta_0
					<< " " << std::setw(8) << k_bend
				    << "\n";

				// OH(29), HH(16)
				index_O = atomIndexOffset + 29;
				index_H = index_O - 13;
				ofs << " " << std::setw(5) << (index_M + 1)
				    << " " << std::setw(5) << (index_O + 1)
                    << " " << std::setw(5) << (index_H + 1)
				    << " " << std::setw(3) << angleFunction
				    << " " << std::setw(8) << theta_0
					<< " " << std::setw(8) << k_bend
				    << "\n";

				//----- Al(1) -----//

				index_M = atomIndexOffset + 1;

				// OH(10), HH(14)
				index_O = atomIndexOffset + 10;
				index_H = index_O + 4;
				ofs << " " << std::setw(5) << (index_M + 1)
				    << " " << std::setw(5) << (index_O + 1)
                    << " " << std::setw(5) << (index_H + 1)
				    << " " << std::setw(3) << angleFunction
				    << " " << std::setw(8) << theta_0
					<< " " << std::setw(8) << k_bend
				    << "\n";

				// OH(11), HH(15)
				index_O = atomIndexOffset + 11;
				index_H = index_O + 4;
				ofs << " " << std::setw(5) << (index_M + 1)
				    << " " << std::setw(5) << (index_O + 1)
                    << " " << std::setw(5) << (index_H + 1)
				    << " " << std::setw(3) << angleFunction
				    << " " << std::setw(8) << theta_0
					<< " " << std::setw(8) << k_bend
				    << "\n";

				// OH(26), HH(30)
				index_O = atomIndexOffset + 26;
				index_H = index_O + 4;
				ofs << " " << std::setw(5) << (index_M + 1)
				    << " " << std::setw(5) << (index_O + 1)
                    << " " << std::setw(5) << (index_H + 1)
				    << " " << std::setw(3) << angleFunction
				    << " " << std::setw(8) << theta_0
					<< " " << std::setw(8) << k_bend
				    << "\n";

				// OH(29), HH(16)
				index_O = atomIndexOffset + 29;
				index_H = index_O - 13;
				ofs << " " << std::setw(5) << (index_M + 1)
				    << " " << std::setw(5) << (index_O + 1)
                    << " " << std::setw(5) << (index_H + 1)
				    << " " << std::setw(3) << angleFunction
				    << " " << std::setw(8) << theta_0
					<< " " << std::setw(8) << k_bend
				    << "\n";

				//----- Al(17) -----//

				index_M = atomIndexOffset + 17;

				// OH(9), HH(13)
				index_O = atomIndexOffset + 9;
				index_H = index_O + 4;
				ofs << " " << std::setw(5) << (index_M + 1)
				    << " " << std::setw(5) << (index_O + 1)
                    << " " << std::setw(5) << (index_H + 1)
				    << " " << std::setw(3) << angleFunction
				    << " " << std::setw(8) << theta_0
					<< " " << std::setw(8) << k_bend
				    << "\n";

				// OH(10), HH(14) [+y cell]
				index_O = numAtomsPerUnitCell*cell_plus_y + 10;
				index_H = index_O + 4;
				ofs << " " << std::setw(5) << (index_M + 1)
				    << " " << std::setw(5) << (index_O + 1)
                    << " " << std::setw(5) << (index_H + 1)
				    << " " << std::setw(3) << angleFunction
				    << " " << std::setw(8) << theta_0
					<< " " << std::setw(8) << k_bend
				    << "\n";
				
				// OH(12), HH(33)
				index_O = atomIndexOffset + 12;
				index_H = index_O + 21;
				ofs << " " << std::setw(5) << (index_M + 1)
				    << " " << std::setw(5) << (index_O + 1)
                    << " " << std::setw(5) << (index_H + 1)
				    << " " << std::setw(3) << angleFunction
				    << " " << std::setw(8) << theta_0
					<< " " << std::setw(8) << k_bend
				    << "\n";

				// OH(28), HH(32)
				index_O = atomIndexOffset + 28;
				index_H = index_O + 4;
				ofs << " " << std::setw(5) << (index_M + 1)
				    << " " << std::setw(5) << (index_O + 1)
                    << " " << std::setw(5) << (index_H + 1)
				    << " " << std::setw(3) << angleFunction
				    << " " << std::setw(8) << theta_0
					<< " " << std::setw(8) << k_bend
				    << "\n";

				//----- Al(18) -----//

				index_M = atomIndexOffset + 18;

				// OH(9), HH(13) [-x cell]
				index_O = numAtomsPerUnitCell*cell_minus_x + 9;
				index_H = index_O + 4;
				ofs << " " << std::setw(5) << (index_M + 1)
				    << " " << std::setw(5) << (index_O + 1)
                    << " " << std::setw(5) << (index_H + 1)
				    << " " << std::setw(3) << angleFunction
				    << " " << std::setw(8) << theta_0
					<< " " << std::setw(8) << k_bend
				    << "\n";

				// OH(12), HH(33) [-x cell]
				index_O = numAtomsPerUnitCell*cell_minus_x + 12;
				index_H = index_O + 21;
				ofs << " " << std::setw(5) << (index_M + 1)
				    << " " << std::setw(5) << (index_O + 1)
                    << " " << std::setw(5) << (index_H + 1)
				    << " " << std::setw(3) << angleFunction
				    << " " << std::setw(8) << theta_0
					<< " " << std::setw(8) << k_bend
				    << "\n";
				
				// OH(27), HH(31)
				index_O = atomIndexOffset + 27;
				index_H = index_O + 4;
				ofs << " " << std::setw(5) << (index_M + 1)
				    << " " << std::setw(5) << (index_O + 1)
                    << " " << std::setw(5) << (index_H + 1)
				    << " " << std::setw(3) << angleFunction
				    << " " << std::setw(8) << theta_0
					<< " " << std::setw(8) << k_bend
				    << "\n";

				// OH(28), HH(32)
				index_O = atomIndexOffset + 28;
				index_H = index_O + 4;
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
	ofs << "\n";

	ofs.close();

	//----- Index file for partially frozen surfaces -----//

	// Indices (in the unit cell) of the Al surface's O and H atoms
	// - Note: O(9)--H(13) and O(26)--H(30) are surface hydroxyl groups
	std::vector<int> unitCellSurfaceO = { 10, 11, 12, 27, 28, 29 };
	//        Corresponding H:          { 14, 15, 33, 31, 32, 16 };

	// - Ordered H atoms:
	std::vector<int> unitCellSurfaceH = { 14, 15, 16, 31, 32, 33 };

	// Arrays of frozen atom indices
	std::vector<int> frozenAtomsFreeH;
	std::vector<int> frozenAtomsFreeOH;

	frozenAtomsFreeH.reserve(numAtomsPerSlab);
	frozenAtomsFreeOH.reserve(numAtomsPerSlab);

	// Number of surface hydroxyl groups per unit cell
	int numSurfaceHydroxylGroups = unitCellSurfaceO.size();

	// Index of topmost layer (along z)
	int topLayer = unitCellGrid[2] - 1;

	// "Free H" surface: only hydroxyl H at the Al surface can move

	bool isFrozen;

	for ( int l=0; l<unitCellGrid[0]; ++l )
	{
		for ( int m=0; m<unitCellGrid[1]; ++m )
		{
			for ( int n=0; n<unitCellGrid[2]; ++n )
			{
				// Current unit cell and index of its first atom
				int cellIndex       = l*unitCellGrid[1]*unitCellGrid[2] 
				                      + m*unitCellGrid[2] + n;
				int atomIndexOffset = numAtomsPerUnitCell*cellIndex;
				
				for ( int i=0; i<numAtomsPerUnitCell; ++i )
				{
					// Reset the flag
					isFrozen = true;

					// If in top layer, might be a free H
					if ( n == topLayer )
					{
						for ( int j=0; j<numSurfaceHydroxylGroups; ++j )
						{
							if ( i == unitCellSurfaceH[j] )
							{
								// Atom is a surface H: don't freeze it
								isFrozen = false;
							}
						}
					}

					if ( isFrozen )
					{
						// Store index
						int atomIndex = atomIndexOffset + i;
						frozenAtomsFreeH.push_back( atomIndex );
					}
				}
			}
		}
	}

	// Free OH surface: only hydroxyl groups at exposed Al surface free 

	for ( int l=0; l<unitCellGrid[0]; ++l )
	{
		for ( int m=0; m<unitCellGrid[1]; ++m )
		{
			for ( int n=0; n<unitCellGrid[2]; ++n )
			{
				// Current unit cell and index of its first atom
				int cellIndex       = l*unitCellGrid[1]*unitCellGrid[2] 
				                      + m*unitCellGrid[2] + n;
				int atomIndexOffset = numAtomsPerUnitCell*cellIndex;
				
				for ( int i=0; i<numAtomsPerUnitCell; ++i )
				{
					// Reset the flag
					isFrozen = true;

					// If in top layer, might be an H
					if ( n == topLayer )
					{
						for ( int j=0; j<numSurfaceHydroxylGroups; ++j )
						{
							if ( i == unitCellSurfaceH[j] )
							{
								// Atom is a surface H: don't freeze it
								isFrozen = false;
							}
							else if ( i == unitCellSurfaceO[j] )
							{
								// Atom is a surface O: don't freeze it
								isFrozen = false;
							}
						}
					}

					if ( isFrozen )
					{
						// Store index
						int atomIndex = atomIndexOffset + i;
						frozenAtomsFreeOH.push_back( atomIndex );
					}
				}
			}
		}
	}

	//----- Write a freeze group file for a single kaolinite slab -----//

	ofs.open("freeze_groups_monolayer.ndx");

	// Multiple indices per line for compact file
	int lineCounter = 0;
	int numAtomsPerLine = 5;

	// Free H surface
	int numFrozenAtomsFreeH = frozenAtomsFreeH.size();

	ofs << "[ KAO_FROZEN_ATOMS_FREE_H ]" << "\n";

	for ( int i=0; i<numFrozenAtomsFreeH; ++i )
	{
		ofs << std::setw(8) << frozenAtomsFreeH[i] << "\t";
		if ( lineCounter < numAtomsPerLine )
		{
			++lineCounter;
		}
		else
		{
			// Done line
			ofs << "\n";
			lineCounter = 0;
		}
	}
	ofs << "\n";

	// Free OH surface
	int numFrozenAtomsFreeOH = frozenAtomsFreeOH.size();

	ofs << "[ KAO_FROZEN_ATOMS_FREE_OH ]" << "\n";

	lineCounter = 0;

	for ( int i=0; i<numFrozenAtomsFreeOH; ++i )
	{
		ofs << std::setw(8) << frozenAtomsFreeOH[i] << "\t";
		if ( lineCounter < numAtomsPerLine )
		{
			++lineCounter;
		}
		else
		{
			// Done line
			ofs << "\n";
			lineCounter = 0;
		}
	}
	ofs << "\n";

	ofs.close();

	//----- Write a second freeze group file for a kaolinite bilayer -----//
	
	// Assumes layers (slabs) are stored consecutively
	std::vector<int> layerAtomOffsets = { 0, numAtomsPerSlab };
	int numLayers = layerAtomOffsets.size();

	ofs.open("freeze_groups_bilayer.ndx");

	// Free H surface
	ofs << "[ KAO_FROZEN_ATOMS_FREE_H ]" << "\n";

	for ( int l=0; l<numLayers; ++l )
	{
		for ( int i=0; i<numFrozenAtomsFreeH; ++i )
		{
			ofs << std::setw(8) 
			    << (frozenAtomsFreeH[i] + layerAtomOffsets[l] + 1) 
			    << "\t";

			if ( lineCounter < numAtomsPerLine )
			{
				++lineCounter;
			}
			else
			{
				// Done line
				ofs << "\n";
				lineCounter = 0;
			}
		}
	}
	ofs << "\n";

	// Free OH surface

	ofs << "[ KAO_FROZEN_ATOMS_FREE_OH ]" << "\n";

	lineCounter = 0;

	for ( int l=0; l<numLayers; ++l )
	{
		for ( int i=0; i<numFrozenAtomsFreeOH; ++i )
		{
			ofs << std::setw(8) 
			    << (frozenAtomsFreeOH[i] + layerAtomOffsets[l] + 1) 
			    << "\t";

			if ( lineCounter < numAtomsPerLine )
			{
				++lineCounter;
			}
			else
			{
				// Done line
				ofs << "\n";
				lineCounter = 0;
			}
		}
	}
	ofs << "\n";

	ofs.close();

}



// Applies the minimum image convention
void minImage(const Rvec x1, const Rvec x2, const Rvec boxL, Rvec x12, double& distSq)
{
	distSq = 0.0;

	for ( int d=0; d<3; d++ )
	{
		x12[d] = x2[d] - x1[d];
		// Apply minimum image convention
		if      ( x12[d] >  0.5*boxL[d] ) { x12[d] -= boxL[d]; }
		else if ( x12[d] < -0.5*boxL[d] ) { x12[d] += boxL[d]; }

		distSq += x12[d]*x12[d];
	}

	return;
}



void minImage(
	const std::vector<double>& x1, const std::vector<double>& x2, 
	const std::vector<double>& boxL,
	// Output
	std::vector<double>& x12, double& distSq)
{
	distSq = 0.0;

	for ( int d=0; d<3; d++ )
	{
		x12[d] = x2[d] - x1[d];
		// Apply minimum image convention
		if      ( x12[d] >  0.5*boxL[d] ) { x12[d] -= boxL[d]; }
		else if ( x12[d] < -0.5*boxL[d] ) { x12[d] += boxL[d]; }

		distSq += x12[d]*x12[d];
	}

	return;
}



// Keeps the atom in the simulaton box by applying PBCs
void keepInBox(const Rvec boxL, Rvec x)
{
	for ( int d=0; d<3; d++ )
	{
		// Apply minimum image convention
		if      ( x[d] > boxL[d] ) { x[d] -= boxL[d]; }
		else if ( x[d] < 0.0 )     { x[d] += boxL[d]; }
	}

	return;
}



void keepInBox(const std::vector<double>& boxL, std::vector<double>& x)
{
	for ( int d=0; d<3; d++ )
	{
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
	for ( int i=0; i<numWaters; i++ )
	{
		if ( hydrogenCounts[i] != 2 )
		{
			return false;
		}
	}
	return true;
}



// Vector norm
double norm(const Rvec x)
{
	return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}



double norm(const std::vector<double>& x)
{
	return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}



// Use the dot product to get the angle between the vectors
double angleBetweenVectors(const Rvec a, const Rvec b)
{
	
	double norm_a = 0.0, norm_b = 0.0, a_dot_b = 0.0;
	for ( int i=0; i<3; ++i )
	{
		norm_a  += a[i]*a[i];
		norm_b  += b[i]*b[i];
		a_dot_b += a[i]*b[i];
	}
	norm_a = sqrt(norm_a);
	norm_b = sqrt(norm_b);

	double theta = acos(a_dot_b/(norm_a*norm_b));

	return theta;
}

