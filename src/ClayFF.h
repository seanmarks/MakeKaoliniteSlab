// ClayFF
// - Organizes parameters for ClayFF force field
//   - See Cygan, Liang, & Kalinichev (J. Phys. Chem. B 2004, 108, 1255-1266)
//   - Converts parameters from the functional forms and units used in the 
//     original article, to the format used by GROMACS
//
// Names of atoms for GROMACS
// - Key (ClayFF nomenclature):
//     AlO  = aluminum, octahedrally coordinated
//     SiT  = silicon, tetrahedrally coordinated
//     OBTS = bridging oxygen bonded to a tetrahedral substitution
//     OBSS = bridging oxygen bonded to two oct. subs OR double-subs.
//     OH   = hydroxyl oxygen
//     OHS  = hydroxyl oxygen bonded to a substituted site
//     HH   = hydroxyl hydrogen
// - Note: "substituted" --> non-oxygen

#ifndef CLAYFF_H
#define CLAYFF_H

#include <array>
#include <cmath>
#include <exception>
#include <iomanip>
#include <fstream>
#include <map>
#include <vector>
#include <stdexcept>
#include <string>

#include "CommonTypes.h"

struct AtomType
{
	std::string type;
	std::string name;
	double mass;     // [amu]
	double charge;   // [e]
	std::string ptype = "A";  // indicates the "type" of atom (e.g. real/virtual site)
	// Lennard-Jones parameters (GROMACS style)
	double sigma;    // [nm]
	double epsilon;  // [kJ/mol]
};

class ClayFF
{
 public:
	ClayFF();

	// Throws an exception if type cannot be found
	const AtomType& getAtomType(const std::string& type) const;

	const std::map<std::string, AtomType>& get_atom_type_map() const { return atom_type_map_; }

	// Print atomtype parameters to file in GROMACS format
	void printAtomTypes(const std::string& atom_types_file) const;

 private:
	std::map<std::string, AtomType> atom_type_map_;

	// Converts from potential form and units used in ClayFF paper to GROMACS style,
	// and registers the atomtype
	void registerAtomType(
		const std::string& type,
		const std::string& name,
		const double mass,          // [amu]
		const double charge,        // [e]
		const double energy_param,  // [kcal/mol]
		const double size_param,    // location of LJ minimum [Angstroms]
		const std::string& ptype = "A"
	);

	// Unit conversions
	static constexpr double kJ_per_kcal_     = 4.184;
	static constexpr double nm_per_Angstrom_ = 0.1;
};

#endif /* CLAYFF_H */
