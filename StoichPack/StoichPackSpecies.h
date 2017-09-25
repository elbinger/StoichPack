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
	Species(const Species& other) : name(other.name), type(other.type), pos_global(other.pos_global), pos_type(other.pos_type) {}
	Species(const std::string& n, species_type t=species_type::unknown) {
		if(n=="") throw StoichPackException("Name must not be empty!");
		name=n;
		type=t;
	}

	void Update(const Species& other){
		if(Name()!=other.Name()) throw StoichPackException("Cannot update from different species!");
		if(Type()!=species_type::unknown) throw StoichPackException("Species already updated!");
		if(other.Type()==species_type::unknown) throw StoichPackException("Cannot update from species with unknown type!");

		type=other.Type();
		pos_global=other.pos_global;
		pos_type=other.pos_type;
	}
	
	const std::string& Name() const { return name; }
	species_type Type() const { return type; }
	size_t GlobalPos() const { return pos_global; }
	size_t TypePos() const { return pos_type; }
	
	bool operator == (const Species& other) const { return name==other.name; }

	void Initialize(size_t pos, size_t n_mobile, size_t n_immobile){
		if(Type()==species_type::unknown) throw StoichPackException("Cannot initialize species with unknown type!");

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
		assert(Type()!=species_type::unknown);
		return *(v+pos_global);
	}

	template<typename ITM, typename ITIMM>
	sp_scalar GetValue(ITM m, ITIMM imm) const{
		assert(Type()!=species_type::unknown);
		if(type==species_type::immobile) return *(imm+pos_type);
		return *(m+pos_type);
	}
 };

} //namespace

#endif
