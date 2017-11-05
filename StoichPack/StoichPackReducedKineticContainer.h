#ifndef __H_STOICHPACK_REDUCED_KINETIC_CONTAINER__
#define __H_STOICHPACK_REDUCED_KINETIC_CONTAINER__

#include "StoichPackIHierarchicalKineticContainer.h"

namespace StoichPack{
template<typename EXT>
class ReducedStageInfo : public HierarchicalStageInfo<EXT>{
public:
	typedef typename EXT::MatrixType MatrixType;
	ReducedStageInfo(const IKineticContainer<EXT>& Base) : base(Base) {
		for(size_t i=0;i<base.Stages();++i) ProcessStage(i);
	}

	size_t Stages() const { return toBase.size(); }
	size_t SubStage(size_t stage) const { return substage[stage]; }
	std::vector<size_t> SubReactions(size_t stage) const {
		const size_t reactions = EXT::cols(stoich[stage]);
		std::vector<size_t> result(reactions);
		for(size_t i=0;i<reactions;++i) result[i]=i;
		return result;
	}
	MatrixType StoichiometricMatrix(size_t stage) const { return stoich[stage]; }	
	MatrixType ToBase(size_t stage) const { return toBase[stage]; }
	MatrixType FromBaseGlobal(size_t stage) const { return fromBase[stage]; }
	size_t MobileSpecies(size_t stage) const { return mobile[stage]; }
	bool ForceNoCorrection(size_t stage) const { return linear[stage]; }
	const IKineticContainer<EXT>& Base() const { return base; }
private:
	const IKineticContainer<EXT>& base;
	std::vector<MatrixType> toBase, fromBase, stoich;
	std::vector<size_t> substage, mobile;
	std::vector<bool> linear;

	typedef typename EXT::OrthogonalDecompositionType OrthType;

	void ProcessStage(size_t stage){
		const OrthType decomp1(Base().MobileStoichiometricMatrix(stage));
		const OrthType decomp2(Base().ImmobileStoichiometricMatrix(stage));
		const size_t coupled_mobile = EXT::rows(decomp1.R());
		const size_t coupled_immobile = EXT::rows(decomp2.R());
		
		const MatrixType mobile_add = EXT::CreateZeroMatrix(EXT::rows(decomp2.Q()),1);
		const MatrixType immobile_add = EXT::CreateZeroMatrix(EXT::rows(decomp1.Q()),1);

		for(size_t i=coupled_immobile;i<Base().ImmobileSpecies(stage);++i){
		 const MatrixType Q = SubCols<EXT>(decomp2.Q(),std::vector<size_t>(1,i));
		 toBase.push_back(CombineRows<EXT>(immobile_add,Q));
		 fromBase.push_back(EXT::Transposed(toBase.back()));
		 stoich.push_back(EXT::CreateMatrix(1,0));
		 substage.push_back(stage);
		 mobile.push_back(0);
		 linear.push_back(true);
		}
	 	for(size_t i=coupled_mobile;i<Base().MobileSpecies(stage);++i){
		 const MatrixType Q = SubCols<EXT>(decomp1.Q(),std::vector<size_t>(1,i));
		 toBase.push_back(CombineRows<EXT>(Q,mobile_add));
		 fromBase.push_back(EXT::Transposed(toBase.back()));
		 stoich.push_back(EXT::CreateMatrix(1,0));
		 substage.push_back(stage);
		 mobile.push_back(1);
		 linear.push_back(true);
	 	}

	 	const MatrixType Q1 = DivideCols<EXT>(decomp1.Q(),coupled_mobile).Mobile();
	 	const MatrixType Q2 = DivideCols<EXT>(decomp2.Q(),coupled_immobile).Mobile();

		const MatrixType Q = CombineBlocks<EXT>(Q1,Q2);
		toBase.push_back(Q);
		fromBase.push_back(EXT::Transposed(Q));
		stoich.push_back(CombineRows<EXT>(decomp1.R(),decomp2.R()));
		substage.push_back(stage);
		mobile.push_back(coupled_mobile);
		linear.push_back(false);
	}

	ReducedStageInfo(); //FORBID
};

template<typename EXT, class BT = IKineticContainer<EXT> >
class ReducedKineticContainer : public IHierarchicalKineticContainer<EXT,BT>{
 public:
	typedef typename EXT::VectorType VectorType;
	typedef typename EXT::VectorArrayType VectorArrayType;
	typedef typename EXT::MatrixType MatrixType;

	ReducedKineticContainer(const BT& base)
		: IHierarchicalKineticContainer<EXT,BT>(ReducedStageInfo<EXT>(base),base){}

	VectorType ConstSpeciesRatesImpl(const VectorType& all, size_t stage) const {
		return FromBaseGlobal()[stage]*Base().ConstSpeciesRatesImpl(all,SubStages()[stage]);
	}
	VectorPair<EXT> ConstSpeciesRatesImpl(const VectorType& mobile, const VectorType& immobile, size_t stage) const {
		const VectorPair<EXT> tmp = Base().ConstSpeciesRatesImpl(mobile,immobile,SubStages()[stage]);
		return VectorPair<EXT>(FromBaseMobile()[stage]*tmp.Mobile(),FromBaseImmobile()[stage]*tmp.Immobile());
	}

	MatrixType RateStructure(size_t stage) const {
		ReducedStageInfo<EXT> I(Base());
		if(EXT::cols(I.StoichiometricMatrix(stage))==0) return EXT::CreateMatrix(0,1);
		else return Base().RateStructure(stage)*BooleanMatrix<EXT>(I.ToBase(stage));
	}

	IKineticContainer<EXT>* copy() const { return new ReducedKineticContainer(*this); }

	private:
	using IHierarchicalKineticContainer<EXT,BT>::Base;
	using IHierarchicalKineticContainer<EXT,BT>::FromBaseGlobal;
	using IHierarchicalKineticContainer<EXT,BT>::FromBaseImmobile;
	using IHierarchicalKineticContainer<EXT,BT>::FromBaseMobile;
	using IHierarchicalKineticContainer<EXT,BT>::SubStages;
};
}

#endif
