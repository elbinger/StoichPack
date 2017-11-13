/* File: StoichPackIReaction.h
 * Purpose: Define interfaces for reactions.
 * Author: Tobias Elbinger <elbinger@math.fau.de> */

#ifndef __H_STOICHPACK_REACTION__
#define __H_STOICHPACK_REACTION__

#include "StoichPackSpecies.h"
#include <memory>

namespace StoichPack{
 /* class IReaction:
  * Provide basic information about a reactions: participants and their Stoichiometric coefficients. */
 class IReaction{
 private:
	std::vector<std::string> participants;
	std::vector<sp_scalar> coefficients;

	IReaction& operator=(const IReaction&); // FORBID

 protected:
	/* AddSpecies: add a new participant.
	 * The first added species can later be accessed with index 0, the second with index 1 ... */
	void AddSpecies(const std::string& name, sp_scalar weight) {
		participants.push_back(name);
		coefficients.push_back(weight);
	}
	void AddSpecies(const std::vector<std::string>& names, const std::vector<sp_scalar>& weights){
		if(names.size()!=weights.size()) throw StoichPackException("names and weights have different size!");
		for(size_t i=0;i<names.size();++i) AddSpecies(names[i],weights[i]);
	}
 public:
	/* getters */
	const std::vector<std::string>& Participants() const { return participants; }
	const std::vector<sp_scalar>& Coefficients() const { return coefficients; }

	/* Write information about a reaction into a string. */
	std::string Info() const {
		std::stringstream info;
		for(size_t i=0;i<participants.size();++i){
			if(coefficients[i]>=0) info<<"+";
			info<<coefficients[i]<<" "<<participants[i]<<" ";
		}
		info<<"= 0 ("<<Name()<<")";
		return info.str();
	}

	virtual std::string Name() const { return "<Reaction>"; } // provide a name if you want to
	virtual ~IReaction() { }
 };

 /* class InitializedReaction:
  * * Store a reaction and the InitializedSpecies (cf. StoichPackSpecies.h) of its participants.
  * * Provide functions to read the concentrations of the participants from containers.
  * * Provide functions to write values to the correct position in a container. */
 template<class ReactionType = IReaction>
 class InitializedReaction{
  private:
	typedef ReactionType* RPTR;
	typedef const ReactionType* CRPTR;

	const std::vector<InitializedSpecies> participants;
	const std::shared_ptr<ReactionType> reaction;
	//FORBID
	InitializedReaction();
	InitializedReaction(CRPTR, const std::vector<InitializedSpecies>&);
	InitializedReaction(RPTR, const std::vector<InitializedSpecies>&);
  public:
	/* Find the participants of an reaction in a vector of InitializedSpecies. */
	static std::vector<InitializedSpecies> FindParticipants(const ReactionType& r, const std::vector<InitializedSpecies>& s){
		std::vector<InitializedSpecies> result;
		for(auto p : r.Participants()){ //for each participant p
			//find p in s
			std::vector<InitializedSpecies>::const_iterator it = std::find(s.begin(),s.end(),p);
			if(it==s.end()) throw StoichPackException("Could not find species!"); //not found
			else result.push_back(*it); // store it
		}
		return result;
	}

	InitializedReaction(const std::shared_ptr<ReactionType>& r, const std::vector<InitializedSpecies>& s)
	                    : participants(FindParticipants(*r,s)), reaction(r) {}

	/* getters */
	std::string Name() const { return reaction->Name(); }
	const std::vector<InitializedSpecies>& Participants() const { return participants; }
	const std::vector<sp_scalar>& Coefficients() const { return reaction->Coefficients(); }
	const ReactionType& Reaction() const { return *reaction; }

	std::string Info() const { return Reaction().Info(); }

	/* Get concentrations of the participants out of a container with concentrations.
	 * Parameters:
	 * v: an iterator to the beginning of the container. */
	template<typename ITV>
	std::vector<sp_scalar> GetValues(ITV v) const {
		const size_t s=participants.size();
		std::vector<sp_scalar> ret(s);
		for(size_t i=0;i<s;++i) ret[i]=participants[i].GetValue(v);
		return ret;
	}

	/* Get concentrations of the participants out of a container for mobile concentrations and a container for immobile concentrations.
	 * Parameters:
	 * itm: an iterator to the beginning of the container with the mobile concentrations.
	 * itimm: an iterator the the beginning of the container with the immobile concentrations. */
	template<typename ITM, typename ITIMM>
	std::vector<sp_scalar> GetValues(ITM m, ITIMM imm) const {
		const size_t s=participants.size();
		std::vector<sp_scalar> ret(s);
		for(size_t i=0;i<s;++i) ret[i]=participants[i].GetValue(m,imm);
		return ret;
	}

	/* Write (more preciesly additive write) values for each participant to a container for mobile species values and a container
	 * for immobile species values.
	 * Parameters:
	 * itm: an iterator to the beginning of a container for mobile species values.
	 * itimm: an iterator to the beginning of a container for immobile species values.
	 * c: a vector containing one value for each participant. */
	template<typename ITM, typename ITIMM>
	void Add(ITM itm, ITIMM itimm, const std::vector<sp_scalar>& c) const {
		for(size_t i=0;i<participants.size();++i){
			if(participants[i].Type()==species_type::mobile) *(itm+participants[i].TypePos)+=c[i];
			else if(participants[i].Type()==species_type::immobile) *(itimm+participants[i].TypePos)+=c[i];
		}
	}

	/* Write (more preciesly additive write) values for each participant to a container for species values.
	 * Parameters:
	 * it: an iterator to the beginning of a container for species values. */
	template<typename IT>
	void Add(IT it, const std::vector<sp_scalar>& c) const {
		for(size_t i=0;i<participants.size();++i){
			*(it+participants[i].GlobalPos)+=c[i];
		}
	}

	virtual ~InitializedReaction() {}
};

} //namespace

#endif
