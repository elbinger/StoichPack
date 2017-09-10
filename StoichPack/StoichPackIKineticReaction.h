#ifndef __H_STOICHPACK_IKINETICREACTION__
#define __H_STOICHPACK_IKINETICREACTION__

#include "StoichPackIReaction.h"
#include <cassert>

namespace StoichPack{

 /* Interface for kinetic reactions */
 class IKineticReaction : public IReaction {
 private:
	//FORBID
	IKineticReaction& operator=(const IKineticReaction&);

 protected:
	using IReaction::AddSpecies;

 public:
	template<typename ITV>
	sp_scalar Rate(ITV v) const { return RateImpl(IReaction::GetValues(v)); }
	template<typename ITM, typename ITIMM>
	sp_scalar Rate(ITM m, ITIMM imm) const { return RateImpl(IReaction::GetValues(m,imm)); }

	template<typename ITV>
	std::vector<sp_scalar> DiffRate(ITV v) const { return DiffRateImpl(IReaction::GetValues(v)); }
	template<typename ITM, typename ITIMM>
	std::vector<sp_scalar> DiffRate(ITM m, ITIMM imm) const { return DiffRateImpl(IReaction::GetValues(m,imm)); }

	//pure virtual functions
	virtual sp_scalar RateImpl(const std::vector<sp_scalar>& c) const =0;
	virtual std::vector<sp_scalar> DiffRateImpl(const std::vector<sp_scalar>& c) const =0;

	//dtor
	virtual ~IKineticReaction() { }

	virtual std::string Name() const { return "<KineticReaction>"; }
	virtual std::vector<sp_scalar> Dependencies() const =0;
 };
} //namespace

#endif

