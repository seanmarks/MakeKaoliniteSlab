#include "ClayFF.h"

ClayFF::ClayFF()
{
	// Taken from NIST
	// - If a range is given for the "Standard Atomic Weight", the midpoint of that 
	//   range is used
	const double mass_Al = 26.982;
	const double mass_Si = 0.5*(28.086 + 28.084);
	const double mass_O  = 16.000;
	const double mass_H  = 1.008;
	/* Comment out for now to silence 'unused variable' warnings
	const double mass_Mg = 0.5*(24.307 + 24.304);
	const double mass_Ca = 40.078;
	const double mass_Fe = 55.845;
	const double mass_Li = 0.5*(6.997 + 6.938);
	*/

	// TODO Include remaining types
	//               type          name   mass     charge   D0         R0
	registerAtomType("HH_ClayFF",  "HH",  mass_H,  +0.4250, 0.0000,    0.0000);  // H, hydroxyl
	registerAtomType("OH_ClayFF",  "OH",  mass_O,  -0.9500, 0.1554,    3.5532);  // O, hydroxyl 
	registerAtomType("OB_ClayFF",  "OB",  mass_O,  -1.0500, 0.1554,    3.5532);  // O, "bridging"
	registerAtomType("AlO_ClayFF", "AlO", mass_Al, +1.5750, 1.3298e-6, 4.7943);  // Al, octahedral
	registerAtomType("SiT_ClayFF", "SiT", mass_Si, +2.1000, 1.8405e-6, 3.7064);  // Si, tetrahedral
}


void ClayFF::registerAtomType(
	const std::string& type, const std::string& name, const double mass, const double charge,
	const double energy_param, const double size_param, const std::string& ptype)
{
	// Organize data and perform conversions
	AtomType new_type;
	new_type.type    = type;
	new_type.name    = name;
	new_type.mass    = mass;
	new_type.charge  = charge;
	new_type.sigma   = (size_param*nm_per_Angstrom_)/pow(2.0, 1.0/6.0);
	new_type.epsilon = energy_param*kJ_per_kcal_;
	new_type.ptype   = ptype;

	// Add to map
	auto insert_pair = atom_type_map_.insert( std::make_pair(new_type.type, new_type) );
	if ( not insert_pair.second ) {
		throw std::runtime_error("Error in ClayFF: Atom type \"" + type + "\" already exists");
	}
}


const AtomType& ClayFF::getAtomType(const std::string& type) const
{
	const auto find_it = atom_type_map_.find(type);
	if ( find_it != atom_type_map_.end() ) {
		return find_it->second;
	}
	else {
		throw std::runtime_error("Error in ClayFF: Atom type \"" + type + "\" does not exist");
	}
}

void ClayFF::printAtomTypes(const std::string& atom_types_file) const
{
	std::ofstream ofs(atom_types_file);
	ofs << "; ClayFF atom types\n"
	    << "[ atomtypes ]\n"
	    << "; type        mass       charge      ptype     sigma[nm]         epsilon[kJ/mol]\n";

	for ( auto it = atom_type_map_.begin(); it != atom_type_map_.end(); ++it ) {
		const AtomType& atom_type = it->second;
		ofs << "  " << atom_type.type 
		    << "   " << std::fixed << std::setprecision(3) << atom_type.mass 
		    << "   " << std::showpos << std::fixed << std::setprecision(3) << atom_type.charge << std::noshowpos
		    << "   " << atom_type.ptype 
		    << "   " << std::scientific << std::setprecision(5) << atom_type.sigma 
		    << "   " << std::scientific << std::setprecision(5) << atom_type.epsilon 
		    << "\n";
	}
	ofs << "\n";

	ofs.close();
}
