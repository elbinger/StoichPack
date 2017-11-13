/* File: StoichPackBiochemicalSystem.h
 * Purpose: Describe biochemical systems.
 * Author: Tobias Elbinger <elbinger@math.fau.de> */

#ifndef __H_STOICPACK_BIOCHEMICAL_SYSTEM__
#define __H_STOICPACK_BIOCHEMICAL_SYSTEM__

#include "StoichPackIKineticReaction.h"
#include <iostream>

namespace StoichPack{

/* class BiochemicalSystem:
 * Store all reactions and biochemical species of a biochemical system.
 * Template parameters:
 *  * KineticType: The storage type for kinetic reactions (default = IKineticReaction). */
template<class KineticType = IKineticReaction>
class BiochemicalSystem{
	typedef std::shared_ptr<KineticType> KPTR;
private:
	std::vector<KPTR> kinetic_reactions; //all participating kinetic reactions
	BasicSpeciesArray species; //all participating species

public:
	//AddSpecies: add one or more species to the system
	void AddSpecies(const std::string& name, species_type type){ species.Add(name,type); }
	void AddSpecies(const std::vector<std::string>& names, species_type type){
		for(auto n : names) AddSpecies(n,type);
	}
	void AddSpecies(const std::initializer_list<std::string>& names, species_type type) {
		for(auto n : names) AddSpecies(n,type);
	}

	//Add a kinetic reaction to the system
	void AddKineticReaction(KineticType* r){ kinetic_reactions.push_back(KPTR(r)); }

	//getters
	const BasicSpeciesArray& Participants() const { return species; }
	const std::vector<KPTR>& KineticReactions() const { return kinetic_reactions; }
 };

/* class InitializedBiochemicalSystem:
 * Store an aequivalent to a BiochemicalSystem object with InitializedReactions (cf. StoichPackIReaction.h/StoichPackIKineticReaction.h)
 * and InitializedSpecies (cf. StoichPackSpecies.h).
 * Template parameters:
 *  * KineticType: storage type for kinetic reactions (default = IKineticReaction).*/
template<class KineticType = IKineticReaction>
class InitializedBiochemicalSystem{
private:
	std::vector<InitializedKineticReaction<KineticType> > kinetic_reactions; //participating kinetic reactions
	std::vector<InitializedSpecies> species; //participating species
	std::vector<size_t> counter; // the number of all species with a certain type (e.g. species_type::mobile)
	InitializedBiochemicalSystem(); //FORBID
public:
	typedef InitializedKineticReaction<KineticType> KineticReactionType;

	//Create InitializedBiochemicalSystem from a BiochemicalSystem
	explicit InitializedBiochemicalSystem(BiochemicalSystem<KineticType> S) : counter(N_SPECIES_TYPES,0) {
		std::vector<BasicSpecies> tmp = S.Participants().SortedData(); //sorted vector of species
		size_t pos=0;
		for(auto p : tmp) // for each participant p
			species.push_back(InitializedSpecies(p,pos,counter)); //add and update counters/positions
		for(auto kr : S.KineticReactions()) // for each kinetic reaction kr
			kinetic_reactions.push_back(KineticReactionType(kr,species)); //initialize and add
	}

	/* getters */
	size_t Count(species_type t) const { return counter[t]; }
	const std::vector<KineticReactionType>& KineticReactions() const { return kinetic_reactions; }
	const std::vector<InitializedSpecies>& Participants() const { return species; }
};

//print InitializedBiochemicalSystem
template<typename ReactionType>
std::ostream& operator<<(std::ostream& os, const InitializedBiochemicalSystem<ReactionType>& system){
	os<<"Species:"<<std::endl;
	for(auto x : system.Participants()) os<<"| "<<x.Info()<<std::endl;
	os<<"Kinetic Reactions:"<<std::endl;
	for(auto x : system.KineticReactions()) os<<"| "<<x.Info()<<std::endl;
	return os;
}

} //namespace

#endif

