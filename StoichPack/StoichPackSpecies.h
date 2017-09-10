#ifndef __H_STOICHPACK_SPECIES__
#define __H_STOICHPACK_SPECIES__

#include "StoichPackDefines.h"
#include <cassert>

namespace StoichPack{

 enum species_type {mobile, immobile, unknown};

 class Species{
 private:
	std::string name; // name of the species
	species_type type;
	size_t pos_global; // global position of species
	size_t pos_type; // position within species with the same type
	
	Species(); //FORBID

 public:
	explicit Species(const std::string& n) : name(n), type(species_type::unknown), pos_global(-1), pos_type(-1) {
		if(n=="") throw StoichPackException("Name must not be empty!");
	}
	Species(const std::string& n, species_type t) : name(n), type(t), pos_global(-1), pos_type(-1){
		if(n=="") throw StoichPackException("Name must not be empty!");
		if(t==species_type::unknown) throw StoichPackException("Type must not be unknown!");
	}
	
	bool Initialized() const {
		return pos_global!=size_t(-1);
	}

	const std::string& Name() const { return name; }
	species_type Type() const { assert(type!=species_type::unknown); return type; }
	size_t GlobalPos() const { assert(Initialized()); return pos_global; }
	size_t TypePos() const { assert(Initialized()); return pos_type; }
	
	void Update(const Species& other) {
		if(Initialized()) throw StoichPackException("Initialized species must not be changed!");
		if(name!=other.name) throw StoichPackException("Cannot use operator= for species "+name+" and "+other.name);

		type=other.type;
		pos_global=other.pos_global;
		pos_type=other.pos_type;
	}

	bool operator == (const Species& other) const { return name==other.name; }

	void Initialize(size_t pos, size_t n_mobile, size_t n_immobile){
		if(Initialized()) throw StoichPackException("Species "+name+" already initialized!");
		if(type == species_type::unknown) throw StoichPackException("Cannot initialize species "+name+" with unknown type!");

		if(type == species_type::mobile){
			if(pos>=n_mobile) throw StoichPackException("Invalid type for desired position!");
			pos_type=pos;
		} else if (type == species_type::immobile){
			if(pos<n_mobile || pos>=n_mobile+n_immobile) throw StoichPackException("Invalid type for desired position!");
			pos_type=pos-n_mobile;
		/*} else if (type == species_type::constant){
			if(pos<n_mobile+n_immobile) throw StoichPackException("Invalid type for desired position!");
			pos_type=pos-n_mobile-n_immobile;*/
		} else { //??
			assert(false);
		}

		pos_global=pos;
	}

	template<typename ITV>
	sp_scalar GetValue(ITV v) const {
		return *(v+pos_global);
	}

	template<typename ITM, typename ITIMM>
	sp_scalar GetValue(ITM m, ITIMM imm) const{
		if(type==species_type::immobile) return *(imm+pos_type);
		return *(m+pos_type);
	}
 };

} //namespace

#endif
