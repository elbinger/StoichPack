/* File: StoichPackIKineticReaction.h
 * Purpose: Provide interfaces for kinetic reactions.
 * Author: Tobias Elbinger <elbinger@math.fau.de> */

#ifndef __H_STOICHPACK_IKINETICREACTION__
#define __H_STOICHPACK_IKINETICREACTION__

#include "StoichPackIReaction.h"
#include <cassert>

namespace StoichPack{

 /* class IKineticReaction: Interface for kinetic reactions.
  * Provide basic functionality for kinetic reactions: the rate law, its derivative, dependencies. */
 class IKineticReaction : public IReaction {
 private:
	//FORBID
	IKineticReaction& operator=(const IKineticReaction&);

 protected:
	using IReaction::AddSpecies;

 public:
	//pure virtual functions

	/* RateImpl: the reaction rate given the concentrations of all participants. c contains only the concentrations of the participants.
	 * The values of c are ordered according to the order defined by AddSpecies. */
	virtual sp_scalar RateImpl(const std::vector<sp_scalar>& c) const =0;
	/* The derivative of the reaction rate with respect to each participant. */
	virtual std::vector<sp_scalar> DiffRateImpl(const std::vector<sp_scalar>& c) const =0;
	/* Dependencies()[i]=false means the i-th participant does not influence the reaction rate. Set all values true if you are not sure. */
	virtual std::vector<bool> Dependencies() const =0;

	//dtor
	virtual ~IKineticReaction() { }

	virtual std::string Name() const { return "<KineticReaction>"; }
 };

/* class InitializedKineticReaction:
 * Extend the InitializedReaction class for treating kinetic reactions, i.e. add function Rate, DiffRate and Dependencies. */
template<typename ReactionType = IKineticReaction>
class InitializedKineticReaction : public InitializedReaction<ReactionType> {
private:
	typedef ReactionType* RPTR;
	typedef const ReactionType* CRPTR;

public:
	InitializedKineticReaction(const std::shared_ptr<ReactionType>& r, const std::vector<InitializedSpecies>& s)
	                           : InitializedReaction<ReactionType>(r,s) {}

	template<typename ITV>
	sp_scalar Rate(ITV v) const {
		//1) Get values of participants
		//2) Call RateImpl of the underlying kinetic reaction.
		return this->Reaction().RateImpl(this->GetValues(v));
	}
	template<typename ITM, typename ITIMM>
	sp_scalar Rate(ITM m, ITIMM imm) const { return this->Reaction().RateImpl(this->GetValues(m,imm)); }

	template<typename ITV>
	std::vector<sp_scalar> DiffRate(ITV v) const { return this->Reaction().DiffRateImpl(this->GetValues(v)); }
	template<typename ITM, typename ITIMM>
	std::vector<sp_scalar> DiffRate(ITM m, ITIMM imm) const { return this->Reaction().DiffRateImpl(this->GetValues(m,imm)); }

	/* Return the result of the Dependencies function of the underlying kinetic reaction as vector of scalars.
	 * (1 if true, 0 if false) */
	std::vector<sp_scalar> Dependencies() const {
		const std::vector<bool> tmp=this->Reaction().Dependencies();
		std::vector<sp_scalar> result(tmp.size());
		for(size_t i=0;i<tmp.size();++i) result[i]=tmp[i] ? 1. : 0.;
		return result;
	}
};
} //namespace

#endif

