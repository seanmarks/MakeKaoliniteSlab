#ifndef CLAYFF_H
#define CLAYFF_H

#include <array>
#include <cmath>
#include <exception>
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
		const double size_param     // location of LJ minimum [Angstroms]
	);

	// Unit conversions
	static constexpr double kJ_per_kcal_     = 4.184;
	static constexpr double nm_per_Angstrom_ = 0.1;
};

#endif /* CLAYFF_H */
