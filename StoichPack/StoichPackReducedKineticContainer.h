/* File: StoichPackReducedKineticStoichiometry.h
 * Author: Tobias Elbinger <elbinger@math.fau.de>
 * Purpose: Provide a preprocessing scheme using linear combinations of species. */

/* This file provides a reaction container using preprocessing as described in [1]. A short summary of [1]:
 *
 * We want to solve the system
 * 	dt c_m + L c_m = S_m*r(c_m,c_i) + F_m,
 * 	dt c_i = S_i*r(c_m,c_i) + F_i,                                                (eq. 1)
 * 	+ boundary condition and inital values.
 * In this case c_m is a vector of mobile concentrations, c_i is a vector of immobile concentrations, S_m and S_r
 * are stoichiometric matrices, r is a vector valued function representing the reaction rates, L is differential
 * operator that permutes with all linear transformations, i.e. A*L c_m = L A*c_m for all matrices A with fitting
 * dimension and F_m and F_i are source terms. We define matrices Q_m, R_m, Q_i, R_i, Z_m and Z_i such that
 * *	Q_m*R_m = S_m and Q_i*R_i=S_i,               (this yields image(Q_x)=image(S_x))
 * *	rank(Q_m)=rank(S_m) and rank(Q_i)=rank(S_i),
 * *	[Q_m Z_m] and [Q_i Z_i] are orthonormal.
 * (In fact, the conditions above only meet a special kind of condition that is needed in [1]. The choice of the
 *  conditions mentioned above yields benefits for implementation.)
 *
 * Let X' denote the transposed of X. We define:
 * * xi_m = Q_m'*c_m and xi_i = Q_i'*c_i,
 * * eta_m = Z_m'*c_m and eta_i = Z_i'*c_i.
 * In this case, the conditions 
 * * c_m = Q_m*xi_m + Z_m*eta_m and
 * * c_i = Q_i*xi_i + Z_i*eta_i
 * hold true. Consequently, solving (eq. 1) is equivalent to solving
 * 	dt eta_m + L eta_m = Z_m'*F_m,
 *                                                                                    (eq. 2a)
 * 	dt eta_ i = Z_i'*F_i,
 * and
 * 	dt xi_m + L xi_m = R_m*r(Q_m*xi_m+Z_m*eta_m,Q_i*xi_i,Z_i*eta_i)+Q_m'*F_m,
 *                                                                                    (eq. 2b)
 *	dt xi_i = R_i*r(Q_m*xi_m+Z_m*eta_m,Q_i*xi_i,Z_i*eta_i)+Q_i'*F_i.
 * Obviously, (eq. 2a) can be solved without knowing xi_m and xi_i. Therefore, we can first solve the (linear)
 * problem (eq. 2a) and then solve the (nonlinear) problem (eq. 2a).
 *
 *
 * In terms of StoichPack this means dividing the problem in (eq. 1) into dim(eta_m)+dim(eta_i)+1 stages:
 * * one stage for every entry of eta_m.
 * * one stage for every entry of eta_i.
 * * one stage that involves xi_m and xi_i.
 *
 *
 * [1]: S. Kräutle, P. Knabner: A new numerical reduction scheme for fully coupled multicomponent
 *      transport-reaction problems in porous media. WRR 2005, DOI: 10.1029/2004WR003624 */

#ifndef __H_STOICHPACK_REDUCED_KINETIC_CONTAINER__
#define __H_STOICHPACK_REDUCED_KINETIC_CONTAINER__

#include "StoichPackIHierarchicalKineticContainer.h"

namespace StoichPack{

//basic information about the stages of a ReducedKineticContainer
template<typename EXT>
class ReducedStageInfo : public IHierarchicalStageInfo<EXT>{
public:
	//construct in dependency of the original container
	typedef typename EXT::MatrixType MatrixType;
	ReducedStageInfo(const IKineticContainer<EXT>& Base) : base(Base) {
		//preprocess every stage of the original container
		for(size_t i=0;i<base.Stages();++i) ProcessStage(i);
	}

	/* functions inherited from interface */

	size_t Stages() const { return toBase.size(); }
	size_t SubStage(size_t stage) const { return substage[stage]; }
	std::vector<size_t> SubReactions(size_t stage) const {
		const size_t reactions = EXT::cols(stoich[stage]);
		std::vector<size_t> result(reactions);
		for(size_t i=0;i<reactions;++i) result[i]=i;

		// case 1: we are on a linear stage / a eta stage.
		//         the stoichiometric matrix will have 0 columns, so result will be empty
		// case 2: we are on nonlinear stage / a xi stage.
		//         result will contain all reactions

		return result;
	}

	MatrixType StoichiometricMatrix(size_t stage) const { return stoich[stage]; }	
	MatrixType ToBase(size_t stage) const { return toBase[stage]; }
	MatrixType FromBaseGlobal(size_t stage) const { return fromBase[stage]; }
	size_t MobileSpecies(size_t stage) const { return mobile[stage]; }
	//bool ForceNoCorrection(size_t stage) const { return linear[stage]; } //not needed in this version
	const IKineticContainer<EXT>& Base() const { return base; }

private:
	const IKineticContainer<EXT>& base; // the original container
	std::vector<MatrixType> toBase; //contributions of this stage to the values of the corresponding stage of the original container
	std::vector<MatrixType> fromBase; //"inverse" of toBase
	std::vector<MatrixType> stoich; //the stoichiometric matrices for each stage
	std::vector<size_t> substage; //corresponding stage of the original container for each stage
	std::vector<size_t> mobile; //mobile species involved in each stage
	//std::vector<bool> linear; //not needed in this version

	typedef typename EXT::OrthogonalDecompositionType OrthType;

	// apply preprocessing to one stage of the original container
	void ProcessStage(size_t stage){
		/* Orthogonal decompositions for ReducedKineticContainer need the following functions:
		 * * Constructor( (const) MatrixType (&) M): construct orthogonal decomposition for M.
		 * * Q(): returns a (const) matrix (reference) with orthonormal columns such that Q()*[R()]=M.
		 *                                                                                    [ 0 ]
		 * * R(): returns a (const) matrix (reference) with rank(M) rows,
		 *        fullfilling the formula mentioned for Q(). */

		// decompositions of stoichiometric matrices
		const OrthType decomp1(Base().MobileStoichiometricMatrix(stage));
		const OrthType decomp2(Base().ImmobileStoichiometricMatrix(stage));

		//ranks
		const size_t coupled_mobile = EXT::rows(decomp1.R());
		const size_t coupled_immobile = EXT::rows(decomp2.R());

		//zero matrices filling blocks		
		const MatrixType mobile_add = EXT::CreateMatrix(EXT::rows(decomp2.Q()),1,0);
		const MatrixType immobile_add = EXT::CreateMatrix(EXT::rows(decomp1.Q()),1,0);

		// add one stage for each entry of eta_i
		for(size_t i=coupled_immobile;i<Base().ImmobileSpecies(stage);++i){
		 const MatrixType Q = SubCols<EXT>(decomp2.Q(),std::vector<size_t>(1,i)); //corresponding col of Q_i
		 toBase.push_back(CombineRows<EXT>(immobile_add,Q)); //corresponding col of [Q_m \n Q_i]
		 fromBase.push_back(EXT::Transposed(toBase.back())); //"inverse"
		 stoich.push_back(EXT::CreateMatrix(1,0)); //one species, no reaction
		 substage.push_back(stage); //corresponding stage of original container
		 mobile.push_back(0); // no mobile species (eta_I!!!!)
		 //linear.push_back(true); //not needed in this version
		}
		// add one stage for each entry of eta_m
	 	for(size_t i=coupled_mobile;i<Base().MobileSpecies(stage);++i){
		 const MatrixType Q = SubCols<EXT>(decomp1.Q(),std::vector<size_t>(1,i)); //corresponding col of Q_m
		 toBase.push_back(CombineRows<EXT>(Q,mobile_add)); //corresponding col of [Q_m \n Q_i]
		 fromBase.push_back(EXT::Transposed(toBase.back())); //"inverse"
		 stoich.push_back(EXT::CreateMatrix(1,0)); //one species, no reaction
		 substage.push_back(stage); //corresponding stage of original container
		 mobile.push_back(1); // 1 mobile species (eta_M!!!!)
		 //linear.push_back(true); //not needed in this version
	 	}

		/* finally add a stage for xi_m and xi_i */

	 	const MatrixType Q1 = DivideCols<EXT>(decomp1.Q(),coupled_mobile).Mobile(); // cols for xi_m
	 	const MatrixType Q2 = DivideCols<EXT>(decomp2.Q(),coupled_immobile).Mobile(); // cols for xi_i

		const MatrixType Q = CombineBlocks<EXT>(Q1,Q2); // [Q_m 0 \n 0 Q_i]
		toBase.push_back(Q);
		fromBase.push_back(EXT::Transposed(Q)); //"inverse"
		stoich.push_back(CombineRows<EXT>(decomp1.R(),decomp2.R())); //stoichiometric matrix (all species)
		substage.push_back(stage); // corresponding stage of original container
		mobile.push_back(coupled_mobile); //number of mobile species in last stage
		//linear.push_back(false); //not needed in this version
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

	/* functions inherited from interface */

	/* ConstSpeciesRates: Base().ConstSpeciesRates with a transfromation */

	VectorType ConstSpeciesRatesImpl(const VectorType& all, size_t stage) const {
		return FromBaseGlobal()[stage]*Base().ConstSpeciesRatesImpl(all,SubStages()[stage]);
	}
	VectorPair<EXT> ConstSpeciesRatesImpl(const VectorType& mobile, const VectorType& immobile, size_t stage) const {
		const VectorPair<EXT> tmp = Base().ConstSpeciesRatesImpl(mobile,immobile,SubStages()[stage]);
		return VectorPair<EXT>(FromBaseMobile()[stage]*tmp.Mobile(),FromBaseImmobile()[stage]*tmp.Immobile());
	}

	MatrixType RateStructure(size_t stage) const {
		//application of chain rule
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
} //namespace

#endif
