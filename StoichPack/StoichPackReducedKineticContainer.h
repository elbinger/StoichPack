#ifndef __H_STOICHPACK_REDUCED_KINETIC_CONTAINER__
#define __H_STOICHPACK_REDUCED_KINETIC_CONTAINER__

#include "StoichPackIHierarchicalLinearKineticContainer.h"

namespace StoichPack{
template<typename EXT, class BT = IKineticContainer<EXT> >
class ReducedKineticContainer : public IHierarchicalLinearKineticContainer<EXT,BT>{
 public:
	typedef typename EXT::VectorType VectorType;
	typedef typename EXT::VectorPairType VectorPairType;
	typedef typename EXT::VectorArrayType VectorArrayType;
	typedef typename EXT::VectorArrayPairType VectorArrayPairType;
	typedef typename EXT::MatrixType MatrixType;
	typedef typename EXT::MatrixPairType MatrixPairType;
	typedef typename EXT::OrthogonalDecompositionType OrthType;

	ReducedKineticContainer(const BT& base) : IHierarchicalLinearKineticContainer<EXT,BT>(base) {
	 for(size_t s=0;s<base.Stages();++s){
		const OrthType decomp1(Base().MobileStoichiometricMatrices()[s]);
		const OrthType decomp2(Base().ImmobileStoichiometricMatrices()[s]);
		const size_t coupled_mobile = EXT::rows(decomp1.R());
		const size_t coupled_immobile = EXT::rows(decomp2.R());
		for(size_t i=coupled_immobile;i<Base().ImmobileSpecies(s);++i){
		 const MatrixType Q = SubCols<EXT>(decomp2.Q(),std::vector<size_t>(1,i));
		 this->AddStage(EXT::CreateMatrix(1,0),s,false,EXT::CreateMatrix(Base().MobileSpecies(s),0),Q,
		                EXT::CreateMatrix(0,Base().MobileSpecies(s)), EXT::Transposed(Q));
	 	}
	 	for(size_t i=coupled_mobile;i<Base().MobileSpecies(s);++i){
		 const MatrixType Q = SubCols<EXT>(decomp1.Q(),std::vector<size_t>(1,i));
		 this->AddStage(EXT::CreateMatrix(1,0),s,false,Q,EXT::CreateMatrix(Base().ImmobileSpecies(s),0),
		                EXT::Transposed(Q), EXT::CreateMatrix(0,Base().ImmobileSpecies(s)));
	 	}

	 	const MatrixType Q1 = DivideCols<EXT>(decomp1.Q(),coupled_mobile).Mobile();
	 	const MatrixType Q2 = DivideCols<EXT>(decomp2.Q(),coupled_immobile).Mobile();

	 	this->AddStage(CombineRows<EXT>(decomp1.R(),decomp2.R()),s,true,Q1,Q2,EXT::Transposed(Q1),EXT::Transposed(Q2));
	 }
	 this->Finish();
	}

	VectorType ReactionRatesImpl1(const VectorType& all, size_t stage) const { 
		if(!this->CorrectionAllowed(stage)) return EXT::CreateVector(0);
		else return Base().ReactionRatesImpl1(all,this->SubStage(stage));
	}
	VectorType ReactionRatesImpl2(const VectorPairType& all, size_t stage) const {
		if(!this->CorrectionAllowed(stage)) return EXT::CreateVector(0);
		else return Base().ReactionRatesImpl2(all,this->SubStage(stage));
	}
	VectorType SubReactionRatesImpl1(const VectorType& all, const std::vector<size_t>& I, size_t stage) const { 
		if(!this->CorrectionAllowed(stage)) return EXT::CreateVector(0);
		else return Base().SubReactionRatesImpl1(all,I,this->SubStage(stage));
	}
	VectorType SubReactionRatesImpl2(const VectorPairType& all, const std::vector<size_t>& I, size_t stage) const {
		if(!this->CorrectionAllowed(stage)) return EXT::CreateVector(0);
		else return Base().SubReactionRatesImpl2(all,I,this->SubStage(stage));
	}

	MatrixPairType EmptyPair(size_t stage) const {
		return MatrixPairType(EXT::CreateMatrix(0,this->MobileSpecies(stage)),EXT::CreateMatrix(0,this->ImmobileSpecies(stage)));
	}

	MatrixType DiffReactionRatesImpl1(const VectorType& all, size_t stage) const {
		if(!this->CorrectionAllowed(stage)) return EXT::CreateMatrix(0,1);
		else return Base().DiffReactionRatesImpl1(all,this->SubStage(stage))*this->toBase(stage);
	}
	MatrixPairType DiffReactionRatesImpl2(const VectorPairType& all, size_t stage) const {
		if(!this->CorrectionAllowed(stage)) return EmptyPair(stage);
		const MatrixPairType tmp = Base().DiffReactionRatesImpl2(all,this->SubStage(stage));
		return MatrixPairType(tmp.Mobile()*this->toBaseMobile(stage),tmp.Immobile()*this->toBaseImmobile(stage));
	}

	MatrixType DiffSubReactionRatesImpl1(const VectorType& all, const std::vector<size_t>& I, size_t stage) const {
		if(!this->CorrectionAllowed(stage)) return EXT::CreateMatrix(0,1);
		else return Base().DiffSubReactionRatesImpl1(all,I,this->SubStage(stage))*this->toBase(stage);
	}
	MatrixPairType DiffSubReactionRatesImpl2(const VectorPairType& all, const std::vector<size_t>& I, size_t stage) const {
		if(!this->CorrectionAllowed(stage)) return EmptyPair(stage);
		const MatrixPairType tmp = Base().DiffSubReactionRatesImpl2(all,I,this->SubStage(stage));
		return MatrixPairType(tmp.Mobile()*this->toBaseMobile(stage),tmp.Immobile()*this->toBaseImmobile(stage));
	}

	VectorType ConstSpeciesRatesImpl1(const VectorType& all, size_t stage) const {
		return this->fromBase(stage)*Base().ConstSpeciesRatesImpl1(all,this->SubStage(stage));
	}

	VectorPairType ConstSpeciesRatesImpl2(const VectorPairType& all, size_t stage) const {
		const VectorPairType tmp = Base().ConstSpeciesRatesImpl2(all,this->SubStage(stage));
		return VectorPairType(this->fromBaseMobile(stage)*tmp.Mobile(),this->fromBaseImmobile(stage)*tmp.Immobile());
	}

	VectorType ConstMobileSpeciesRatesImpl(const VectorPairType& all, size_t stage) const {
		return this->fromBaseMobile(stage)*Base().ConstMobileSpeciesRatesImpl(all,this->SubStage(stage));
	}
	VectorType ConstImmobileSpeciesRatesImpl(const VectorPairType& all, size_t stage) const {
		return this->fromBaseImmobile(stage)*Base().ConstImmobileSpeciesRatesImpl(all,this->SubStage(stage));
	}

	bool ApplyMobileCorrectionImpl(VectorArrayType& mobile, const VectorArrayType& allowed) const {
		return true;		
	}
	bool ApplyImmobileCorrectionImpl(VectorArrayType& immobile, const VectorArrayType& allowed) const {
		return true;
	}
		
	bool ApplyCorrectionImpl(VectorArrayType& all, const VectorArrayType& allowed) const {
		return true;
	}

	MatrixType SpeciesStructure(size_t stage) const { return EXT::CreateMatrix(1,1); }

	IKineticContainer<EXT>* copy() const { return new ReducedKineticContainer(*this); }

	private:
	using IHierarchicalLinearKineticContainer<EXT,BT>::Base;
};
}

#endif
