/* File: StoichPackBiochemicalSystem.h
 * Purpose: Describe biochemical systems.
 * Author: Tobias Elbinger <elbinger@math.fau.de> */

#ifndef __H_STOICPACK_BIOCHEMICAL_SYSTEM__
#define __H_STOICPACK_BIOCHEMICAL_SYSTEM__

#include "StoichPackIKineticReaction.h"
#include <memory>
#include <iostream>

namespace StoichPack{

/* class BiochemicalSystem:
 * Store all reactions and biochemical species of a biochemical system. */
template<class KineticType = IKineticReaction>
class BiochemicalSystem{
	typedef KineticType* KPTR;
	typedef std::shared_ptr<KineticType> SHARED_KPTR;
	typedef const KineticType& KREF;
private:
	std::vector<SHARED_KPTR> kinetic_reactions;
	BasicSpeciesArray species;

public:
	void AddSpecies(const std::string& name, species_type type){ species.Add(name,type); }
	void AddSpecies(const std::vector<std::string>& names, species_type type){
		for(auto n : names) AddSpecies(n,type);
	}
	void AddSpecies(const std::initializer_list<std::string>& names, species_type type) {
		for(auto n : names) AddSpecies(n,type);
	}

	void AddKineticReaction(KPTR r){ kinetic_reactions.push_back(SHARED_KPTR(r)); }

	const BasicSpeciesArray& Participants() const { return species; }
	const std::vector<SHARED_KPTR>& KineticReactions() const { return kinetic_reactions; }
 };

/* class InitializedBiochemicalSystem:
 * Store an aequivalent to a BiochemicalSystem object with InitializedReactions (cf. StoichPackIReaction.h/StoichPackIKineticReaction.h)
 * and InitializedSpecies (cf. StoichPackSpecies.h). */
template<class KineticType = IKineticReaction>
class InitializedBiochemicalSystem{
private:
	std::vector<InitializedKineticReaction<KineticType> > kinetic_reactions;
	std::vector<InitializedSpecies> species;
	std::vector<size_t> counter; // the number of all species with a certain type (e.g. species_type::mobile)
	InitializedBiochemicalSystem(); //FORBID
public:
	typedef InitializedKineticReaction<KineticType> KineticReactionType;

	explicit InitializedBiochemicalSystem(BiochemicalSystem<KineticType> S) : counter(N_SPECIES_TYPES,0) {
		std::vector<BasicSpecies> tmp = S.Participants().SortedData(); //sorted vector of species
		size_t pos=0;
		for(auto p : tmp) // for each participant p
			species.push_back(InitializedSpecies(p,pos,counter));
		for(auto kr : S.KineticReactions()) // for each kinetic reaction kr
			kinetic_reactions.push_back(KineticReactionType(kr,species));
	}

	/* getters */
	size_t Count(species_type t) const { return counter[t]; }
	const std::vector<KineticReactionType>& KineticReactions() const { return kinetic_reactions; }
	const std::vector<InitializedSpecies>& Participants() const { return species; }
};

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

