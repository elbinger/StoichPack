/*File: StoichPackSpecies.h
 * Purpose: Provide data structures for storing biochemical species.
 * Author: Tobias Elbinger <elbinger@math.fau.de> */

#ifndef __H_STOICHPACK_SPECIES__
#define __H_STOICHPACK_SPECIES__

#include "StoichPackDefines.h"
#include <cassert>
#include <vector>
#include <algorithm>
#include <sstream>

namespace StoichPack{
 /* Biochemical species may have 4 different types:
  * * mobile: mobile species can be transported
  * * immobile: immobile species cannot be transported
  * * microscopic: not supported in this version.
  * * constant: not supported in this version. */

 enum species_type {mobile, immobile, microscopic, constant};
#define N_SPECIES_TYPES 4

 /* class Basic Species:
  * * store the basic informations (name and type)
  * * provide functions for updating species
  * * provide functions for sorting species */
 class BasicSpecies{
 private:
	std::string name;
	species_type type;

	/* Helper function to check input data. */
	void checktype() const {
		if(Type()==constant) throw StoichPackException("Constant species are not supported in this version!");
		if(Type()==microscopic) throw StoichPackException("Microscopic species are not supported in this version!");

		if(Type()==mobile || Type()==immobile) return;
		throw StoichPackException("unkown type!");
	}

	//FORBID the following operators and constructors
	BasicSpecies();
	bool operator==(const BasicSpecies&) const;
	bool operator!=(const BasicSpecies&) const;

 public:
	/* Create a species by defining its name in n and its type in t */
	BasicSpecies(const std::string& n, species_type t) : name(n), type(t){
		if(n=="") throw StoichPackException("species must have a name!");
		checktype();
	}

	/* Can the type be changed? */
	static bool FixedType(species_type t) { return t!=constant; }

	/* getters */
	const std::string& Name() const { return name; }
	species_type Type() const { return type; }
	std::string TypeName() const { // convert species_type to string
		switch(Type()){
			case mobile: return "mobile";
			case immobile: return "immobile";
			case microscopic: return "microscopic";
			case constant: return "constant";
		}
		return "unknown type";
	}

	/* Write Information about the species into a string. */
	std::string Info() const {
		std::stringstream info;
		info<<Name()<<" (type="<<TypeName()<<")";
		return info.str();
	}

	/* Update from another species, i.e. change the type if it can be changed,
 	 * throw an exception if types differ and cannot be changed.
	 * Remark: This is particullary useful if 2 biochemical systems are combined. Imagine one system that models the consumtion of
	 * a certain species but not its production. In this scenario one might consider a constant value for that species. If this
	 * model is combined with a model for the production of this species, this update function becomes very useful. */
	void Update(const BasicSpecies& other){
		if(Name()!=other.Name()) throw StoichPackException("Cannot update from different species!");

		if(FixedType(Type())) { //Can my type be changed?
			//No. --> Other species must have same type or a type that can be changed
			if(FixedType(other.Type()) && Type()!=other.Type()) throw StoichPackException("Type mismatch!");
		} else type=other.Type(); //update type from other species
	}

	/* Specifies if this species has a lower index than the species in other. */
	bool operator<(const BasicSpecies& other) const {
		if(Type()==other.Type()) // same type --> make sure species are alphabetically ordered
			return Name().compare(other.Name())<0;
		else return Type()<other.Type(); //sort according to type
	}

	bool operator==(const std::string& n) const { return Name()==n; }
 };

 //Container for biochemical species
 class BasicSpeciesArray{
  private:
   std::vector<BasicSpecies> data;
  public:
	/* Add functions: Try to add a new species to the container:
	 * If a species with the same value is already contained, update this species if possible or throw an excption.
	 * Otherwise add to the container. */
	void Add(const std::string& n, species_type t){
		BasicSpecies S(n,t);
		Add(S);
	}

	void Add(const BasicSpecies& S){
		std::vector<BasicSpecies>::iterator it = std::find(data.begin(),data.end(),S.Name());
		if(it==data.end()) data.push_back(S);
		else it->Update(S);
	}

	/* getters */
	size_t size() const { return data.size(); }

	const BasicSpecies& operator[](size_t i) const {
		assert(i<size());
		return data[i];
	}

	/* Return sorted species. This makes sure that the positions of the species are independant of the order of they have been added. */
	std::vector<BasicSpecies> SortedData() const {
		std::vector<BasicSpecies> result(data);
		std::sort(result.begin(),result.end());
		return result;
	}
 };

 /* class InitializedSpecies:
  * * assigns 2 integer valued indeces to each species:
  * * * a global position.
  * * * a position within all species with the same type
  * * provides functions for reading the concentration of the species from a given vector of concentrations. */
 class InitializedSpecies{
 private:
	//FORBID
	InitializedSpecies();
	bool operator==(const InitializedSpecies&);
	bool operator!=(const InitializedSpecies&);

 public:
	const BasicSpecies Base;
	const size_t GlobalPos; // global position of species
	const size_t TypePos; // position within species with the same type

	/* Construct new InitializedSpecies, given a BasicSpecies and informations about already initialized species and
	 * update those informations.
	 * Parameters:
	 * s: a the according BasicSpecies object.
	 * pos: the number of already initialized species
	 * typepos: the number of already initialized mobile, immobile, microscopic and constant species */
	InitializedSpecies(const BasicSpecies& s, size_t& pos, std::vector<size_t>& typepos)
	                   : Base(s), GlobalPos(pos++), TypePos(typepos[s.Type()]++){}

	/* getters */
	const std::string& Name() const { return Base.Name(); }
	species_type Type() const { return Base.Type(); }
	const std::string TypeName() const { return Base.TypeName(); }

	/* Write information about the species into a string. */
	std::string Info() const {
		std::stringstream info;
		info<<Base.Name()<<" (type="<<TypeName()<<", globalpos="<<GlobalPos<<", typepos="<<TypePos<<")";
		return info.str();
	}
	
	bool operator == (const std::string& other) const { return Base==other; }

	/* Read concentration of the species out of a container.
	 * Parameters:
	 * v: an iterator to the beginning of the container. */
	template<typename ITV>
	sp_scalar GetValue(ITV v) const {
		return *(v+GlobalPos);
	}

	/* Read concentration of the species out of a container for mobile species and a container for immobile species.
	 * Parameters:
	 * m: a container with the concentrations of all mobile species
	 * imm: a container with the concentration of all immobile species */
	template<typename ITM, typename ITIMM>
	sp_scalar GetValue(ITM m, ITIMM imm) const{
		if(Type()==species_type::immobile) return *(imm+TypePos);
		return *(m+TypePos);
	}
 };

} //namespace

#endif
