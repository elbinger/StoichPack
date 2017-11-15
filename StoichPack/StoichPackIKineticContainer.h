/* File: StoichPackIKineticContainer.h
 * Purpose: Define an interface for kinetic reaction containers.
 * Author: Tobias Elbinger <elbinger@math.fau.de> */

#ifndef __H_STOICHPACK_IKINETIC_CONTAINER__
#define __H_STOICHPACK_IKINETIC_CONTAINER__

#include "StoichPackStageInfo.h"
#include "StoichPackSpecies.h"

namespace StoichPack{
 //general container interface
 template<typename EXT >
 class IKineticContainer{
	typedef typename EXT::VectorType VectorType;
	typedef typename EXT::MatrixType MatrixType;
	typedef typename EXT::VectorArrayType VectorArrayType;

 private:
	StageInfoArray<EXT> stageinfo; //information about the stages

	//MatrixType toOriginal_global, toOriginal_mobile, toOriginal_immobile; //not needed in this version
	//MatrixType fromOriginal_global, fromOriginal_mobile, fromOriginal_immobile; //not needed in this version

	IKineticContainer& operator=(const IKineticContainer<EXT>&);	//FORBID
	IKineticContainer(); //FORBID

 protected:
	IKineticContainer(const StageInfoArray<EXT>& info) : stageinfo(info) {}

	/*//not needed in this version:
	IKineticContainer(const StageInfoArray<EXT>& info) : stageinfo(info), toOriginal_global(EXT::CreateMatrix(info.GlobalSpecies())),
		toOriginal_mobile(EXT::CreateMatrix(info.MobileSpecies())), toOriginal_immobile(EXT::CreateMatrix(info.ImmobileSpecies())),
		fromOriginal_global(EXT::CreateMatrix(info.GlobalSpecies())), fromOriginal_mobile(EXT::CreateMatrix(info.MobileSpecies())),
		fromOriginal_immobile(EXT::CreateMatrix(info.ImmobileSpecies()))
		{

		std::vector<MatrixType> to_global, to_mobile, to_immobile, from_global, from_mobile, from_immobile;

		for(auto x : info){
			to_global.push_back(x.Global().ToOriginal());
			to_mobile.push_back(x.Mobile().ToOriginal());
			to_immobile.push_back(x.Immobile().ToOriginal());
			from_global.push_back(x.Global().FromOriginal());
			from_mobile.push_back(x.Mobile().FromOriginal());
			from_immobile.push_back(x.Immobile().FromOriginal());
		}

		typedef typename std::vector<MatrixType>::iterator IT;

		toOriginal_global = ColCat<EXT,IT>(to_global.begin(),to_global.end());
		toOriginal_mobile = ColCat<EXT,IT>(to_mobile.begin(),to_mobile.end());
		toOriginal_immobile = ColCat<EXT,IT>(to_immobile.begin(),to_immobile.end());

		fromOriginal_global = RowCat<EXT,IT>(from_global.begin(),from_global.end());
		fromOriginal_mobile = RowCat<EXT,IT>(from_mobile.begin(),from_mobile.end());
		fromOriginal_immobile = RowCat<EXT,IT>(from_immobile.begin(),from_immobile.end());
	}*/

 public:
	/* getters */
	const StageInfoArray<EXT>& StageInformation() const { return stageinfo; }
	size_t Stages() const { return stageinfo.size(); }

	const BasicStageInfo<EXT>& Global(size_t stage) const { return stageinfo[stage].Global(); }
	const BasicStageInfo<EXT>& Mobile(size_t stage) const { return stageinfo[stage].Mobile(); }
	const BasicStageInfo<EXT>& Immobile(size_t stage) const { return stageinfo[stage].Immobile(); }

	const MatrixType& StoichiometricMatrix(size_t stage) const { return Global(stage).StoichiometricMatrix(); }
	const MatrixType& MobileStoichiometricMatrix(size_t stage) const { return Mobile(stage).StoichiometricMatrix(); }
	const MatrixType& ImmobileStoichiometricMatrix(size_t stage) const { return Immobile(stage).StoichiometricMatrix(); }

	size_t GlobalSpecies(size_t stage) const { return Global(stage).GlobalSpecies(); }
	size_t GlobalSpecies() const { return stageinfo.GlobalSpecies(); }
	size_t MobileSpecies(size_t stage) const { return Mobile(stage).GlobalSpecies(); }
	size_t MobileSpecies() const { return stageinfo.MobileSpecies(); }
	size_t ImmobileSpecies(size_t stage) const { return Immobile(stage).GlobalSpecies(); }
	size_t ImmobileSpecies() const { return stageinfo.ImmobileSpecies(); }

	size_t Reactions(size_t stage) const { return EXT::cols(StoichiometricMatrix(stage)); }

	//const MatrixType& ToOriginalGlobal() const { return toOriginal_global; } // not needed in this version
	VectorType ToOriginalGlobal(const VectorArrayType& all) const { return stageinfo.ToOriginalGlobal(all); }

	//const MatrixType& ToOriginalMobile() const { return toOriginal_mobile; } // not needed in this version
	VectorType ToOriginalMobile(const VectorArrayType& mobile) const { return stageinfo.ToOriginalMobile(mobile); }
	
	//const MatrixType& ToOriginalImmobile() const { return toOriginal_immobile(); } // not needed in this version
	VectorType ToOriginalImmobile(const VectorArrayType& immobile) const { return stageinfo.ToOriginalImmobile(immobile); }

	VectorPair<EXT> ToOriginal(const VectorArrayType& mobile, const VectorArrayType& immobile) const {
		return VectorPair<EXT>(ToOriginalMobile(mobile),ToOriginalImmobile(immobile));
	}
	VectorPair<EXT> ToOriginal(const VectorArrayPair<EXT>& all) const {
		return ToOriginal(all.Mobile(),all.Immobile());
	}

	//const MatrixType& FromOriginalGlobal() const { return fromOriginal_global; } // not needed in this version
	VectorArrayType FromOriginalGlobal(const VectorType& all) const { return stageinfo.FromOriginalGlobal(all); }

	//const MatrixType& FromOriginalMobile() const { return fromOriginal_mobile; } // not needed in this version
	VectorArrayType FromOriginalMobile(const VectorType& mobile) const { return stageinfo.FromOriginalMobile(mobile); }
	
	//const MatrixType& FromOriginalImmobile() const { return fromOriginal_immobile; }// not needed in this version
	VectorArrayType FromOriginalImmobile(const VectorType& immobile) const { return stageinfo.FromOriginalImmobile(immobile); }

	VectorArrayPair<EXT> FromOriginal(const VectorType& mobile, const VectorType& immobile) const {
		return VectorArrayPair<EXT>(FromOriginalMobile(mobile), FromOriginalImmobile(immobile));
	}

	VectorArrayPair<EXT> FromOriginal(const VectorPair<EXT>& all) const {
		return FromOriginal(all.Mobile(),all.Immobile());
	}

	/* pure virtual functions: */

	/* Provide a function to evaluate all reaction rates of a stage, given the values of all concentrations (original presentation). */
	virtual VectorType ReactionRatesImpl(const VectorType& all, size_t stage) const =0;
	/* Provide a funtion that only evaluates the reactions speciefied in I. */
	virtual VectorType SubReactionRatesImpl(const VectorType& all, const std::vector<size_t>& I, size_t stage) const =0;
	/* Provide a function to evaluate all reaction rates of a stage, given the values of all mobile concentrations
	 * and all immobile concentrations separately (original presentation). */
	virtual VectorType ReactionRatesImpl(const VectorType& mobile, const VectorType& immobile, size_t stage) const =0;
	/* Provide a funtion that only evaluates the reactions speciefied in I. */
	virtual VectorType SubReactionRatesImpl(const VectorType& mobile, const VectorType& immobile,
	                                        const std::vector<size_t>& I, size_t stage) const =0;

	/* Provide a function for the Jacobian of ReactionRatesImpl. */
	virtual MatrixType DiffReactionRatesImpl(const VectorType& all, size_t stage) const =0;
	/* Provide a function for the Jacobian of SubReactionRatesImpl. */
	virtual MatrixType DiffSubReactionRatesImpl(const VectorType& all, const std::vector<size_t>& I, size_t stage) const=0;
	/* Provide a function for the Jacobian of ReactionRatesImpl. */
	virtual MatrixPair<EXT> DiffReactionRatesImpl(const VectorType& mobile, const VectorType& immobile, size_t stage) const =0;
	/* Provide a function for the Jacobian of SubReactionRatesImpl. */
	virtual MatrixPair<EXT> DiffSubReactionRatesImpl(const VectorType& mobile, const VectorType& immobile,
	                                                 const std::vector<size_t>& I, size_t stage) const=0;

	/* ConstSpeciesRates: Provide functions for right-hand side contributions that do not depend on concentrations of this stage
	 * (e.g. source terms). */
	virtual VectorType ConstSpeciesRatesImpl(const VectorType& all, size_t stage) const =0;
	virtual VectorPair<EXT> ConstSpeciesRatesImpl(const VectorType& mobile, const VectorType& immobile, size_t stage) const =0;

	/* Provide a function that returns all biochemical species of the underlying biochemical problem. The order must not be changed! */
	virtual const std::vector<InitializedSpecies>& Participants() const =0;

	/* Provide a function to determine if a certain reaction rate depends on a certain species, i.e. the entry in the i-th row
	 * and j-th column of the return value is 0 if reaction i does not depend on the concentration of species j
	 * and >0 otherwise (preprocessed presentation).*/
	virtual MatrixType RateStructure(size_t stage) const = 0;

	/* Create a copy of *this. */
	virtual IKineticContainer<EXT>* copy() const =0;

	/* Reaction rates of all reactions in a stage, given preprocessed concentration values */
	VectorType ReactionRates(const VectorArrayType& all, size_t stage) const {
		return ReactionRatesImpl(ToOriginalGlobal(all),stage);
	}

	/* Reaction rates of the subreactions specified in I of a stage, given preprocessed concentration values */
	VectorType SubReactionRates(const VectorArrayType& all, const std::vector<size_t>& I, size_t stage) const {
		return SubReactionRatesImpl(ToOriginalGlobal(all),I,stage);
	}

	/* Reaction rates of all reactions in a stage, given preprocessed concentration values */
	VectorType ReactionRates(const VectorArrayType& mobile, const VectorArrayType& immobile, size_t stage) const {
		return ReactionRatesImpl(ToOriginalMobile(mobile),ToOriginalImmobile(immobile),stage);
	}

	/* Reaction rates of all reactions in a stage, given preprocessed concentration values */
	VectorType ReactionRates(const VectorArrayPair<EXT>& all, size_t stage) const {
		return ReactionRates(all.Mobile(),all.Immobile(),stage);
	}

	/* Reaction rates of the subreactions specified in I of a stage, given preprocessed concentration values */
	VectorType SubReactionRates(const VectorArrayType& mobile, const VectorArrayType& immobile,
	                                const std::vector<size_t>& I, size_t stage) const {
		return SubReactionRatesImpl(ToOriginalMobile(mobile),ToOriginalImmobile(immobile),I,stage);
	}

	/* Reaction rates of the subreactions specified in I of a stage, given preprocessed concentration values */
	VectorType SubReactionRates(const VectorArrayPair<EXT>& all, const std::vector<size_t>& I, size_t stage) const {
		return SubReactionRates(all.Mobile(),all.Immobile(),I,stage);
	}

	/* DiffReactionRates and DiffSubReactionRates: return the Jacobians of ReactionRates and SubReactionRates.
	 * We have to use the chain rule here. */

	MatrixType DiffReactionRates(const VectorArrayType& all, size_t stage) const {
		return DiffReactionRatesImpl(ToOriginalGlobal(all),stage)*Global(stage).ToOriginal();
	}

	MatrixType DiffSubReactionRates(const VectorArrayType& all, const std::vector<size_t>& I, size_t stage) const {
		return DiffSubReactionRatesImpl(ToOriginalGlobal(all),I,stage)*Global(stage).ToOriginal();
	}

	MatrixPair<EXT> DiffReactionRates(const VectorArrayType& mobile, const VectorArrayType& immobile, size_t stage) const {
		const MatrixPair<EXT> tmp = DiffReactionRatesImpl(ToOriginalMobile(mobile),ToOriginalImmobile(immobile),stage);
		return MatrixPair<EXT>(tmp.Mobile()*Mobile(stage).ToOriginal(),tmp.Immobile()*Immobile(stage).ToOriginal());
	}

	MatrixPair<EXT> DiffReactionRates(const VectorArrayPair<EXT>& all, size_t stage) const {
		return DiffReactionRates(all.Mobile(),all.Immobile(),stage);
	}

	MatrixPair<EXT> DiffSubReactionRates(const VectorArrayType& mobile, const VectorArrayType& immobile,
	                                            const std::vector<size_t>& I, size_t stage) const {
		const MatrixPair<EXT> tmp = DiffSubReactionRatesImpl(ToOriginalMobile(mobile),ToOriginalImmobile(immobile),I,stage);
		return MatrixPair<EXT>(tmp.Mobile()*Mobile(stage).ToOriginal(),tmp.Immobile()*Immobile(stage).ToOriginal());
	}

	MatrixPair<EXT> DiffSubReactionRates(const VectorArrayPair<EXT> all, const std::vector<size_t>& I, size_t stage) const {
		return DiffSubReactionRates(all.Mobile(),all.Immobile(),I,stage);
	}


	/* ConstSpeciesRates: Provide right hand side contributions that do not depend on concentration values of the respective stage
	 * (e.g. source terms).*/
	VectorType ConstSpeciesRates(const VectorArrayType& all, size_t stage) const {
		return ConstSpeciesRatesImpl(ToOriginalGlobal(all),stage);
	}

	VectorPair<EXT> ConstSpeciesRates(const VectorArrayType& mobile, const VectorArrayType& immobile, size_t stage) const {
		return ConstSpeciesRatesImpl(ToOriginalMobile(mobile),ToOriginalImmobile(immobile),stage);
	}

	VectorPair<EXT> ConstSpeciesRates(const VectorArrayPair<EXT>& all, size_t stage) const {
		return ConstSpeciesRates(all.Mobile(),all.Immobile(),stage);
	}

	/* SpeciesRates: Evaluate the right hand side for all species in a stage. This is equivalent to
		stoichiometric_matrix * rate_vector + constants */
	VectorType SpeciesRates(const VectorArrayType& all, size_t stage) const {
		return Global(stage).StoichiometricMatrix()*ReactionRates(all,stage)+ConstSpeciesRates(all,stage);
	}

	VectorPair<EXT> SpeciesRates(const VectorArrayType& mobile, const VectorArrayType& immobile, size_t stage) const {
		const VectorPair<EXT> crates = ConstSpeciesRates(mobile,immobile,stage);
		const VectorType rrates = ReactionRates(mobile,immobile,stage);
		return VectorPair<EXT>(Mobile(stage).StoichiometricMatrix()*rrates+crates.Mobile(),
		                       Immobile(stage)*rrates+crates.Immobile());
	}

	VectorPair<EXT> SpeciesRates(const VectorArrayPair<EXT>& all, size_t stage) const {
		return SpeciesRates(all.Mobile(),all.Immobile(),stage);
	}

	/* Mobile/ImmobileSpeciesRates: Evaluate the right hand side for mobile/immobile species in a stage. */

	VectorType MobileSpeciesRates(const VectorArrayType& mobile, const VectorArrayType& immobile, size_t stage) const {
		return Mobile(stage).StoichiometricMatrix()*ReactionRates(mobile,immobile,stage)
		       +ConstSpeciesRates(mobile,immobile,stage).Mobile();
	}

	VectorType MobileSpeciesRates(const VectorArrayPair<EXT>& all, size_t stage) const {
		return MobileSpeciesRates(all.Mobile(),all.Immobile(),stage);
	}

	VectorType ImmobileSpeciesRates(const VectorArrayType& mobile, const VectorArrayType& immobile, size_t stage) const {
		return Immobile(stage).StoichiometricMatrix()*ReactionRates(mobile,immobile,stage)
		       +ConstSpeciesRates(mobile,immobile,stage).Immobile();
	}

	VectorType ImmobileSpeciesRates(const VectorArrayPair<EXT>& all, size_t stage) const {
		return ImmobileSpeciesRates(all.Mobile(),all.Immobile(),stage);
	}

	MatrixType DiffSpeciesRates(const VectorArrayType& all, size_t stage) const {
		return Global(stage).StoichiometricMatrix()*DiffReactionRates(all,stage);
	}

	/* DiffSpeciesRates/DiffMobileSpeciesRates/DiffImmobileSpeciesRates:
	 * return the Jacobians of SpeciesRates/MobileSpeciesRates/ImmobileSpeciesRates.*/
	MatrixQuad<EXT> DiffSpeciesRates(const VectorArrayType& mobile, const VectorArrayType& immobile, size_t stage) const {
		const MatrixPair<EXT> tmp = DiffReactionRates(mobile,immobile,stage);
		const MatrixType& smob = Mobile(stage).StoichiometricMatrix();
		const MatrixType& simmob = Immobile(stage).StoichiometricMatrix();

		return MatrixQuad<EXT>(MatrixPair<EXT>(smob*tmp.Mobile()  , smob*tmp.Immobile()),
		                       MatrixPair<EXT>(simmob*tmp.Mobile(), simmob*tmp.Immobile()));
	}

	MatrixQuad<EXT> DiffSpeciesRates(const VectorArrayPair<EXT>& all, size_t stage) const {
		return DiffSpeciesRates(all.Mobile(),all.Immobile(),stage);
	}


	MatrixPair<EXT> DiffMobileSpeciesRates(const VectorArrayType& mobile, const VectorArrayType& immobile,size_t stage) const{
		const MatrixPair<EXT> tmp = DiffReactionRates(mobile,immobile,stage);
		const MatrixType& smob = Mobile(stage).StoichiometricMatrix();
		return MatrixPair<EXT>(smob*tmp.Mobile(),smob*tmp.Immobile());
	}

	MatrixPair<EXT> DiffMobileSpeciesRates(const VectorArrayPair<EXT>& all, size_t stage) const{
		return DiffMobileSpeciesRates(all.Mobile(),all.Immobile(),stage);
	}

	MatrixPair<EXT> DiffImmobileSpeciesRates(const VectorArrayType& mobile, const VectorArrayType& immobile,size_t stage) const{
		const MatrixPair<EXT> tmp = DiffReactionRates(mobile,immobile,stage);
		const MatrixType& simmob = Immobile(stage).StoichiometricMatrix();
		return MatrixPair<EXT>(simmob*tmp.Mobile(),simmob*tmp.Immobile());
	}

	MatrixPair<EXT> DiffImmobileSpeciesRates(const VectorArrayPair<EXT>& all, size_t stage) const{
		return DiffImmobileSpeciesRates(all.Mobile(),all.Immobile(),stage);
	}

	/* specify whether there is a direct infulence of species j on species i, i.e. the entry in the i-th row and j-th column is
	 * 0 if there is a reaction that influences species i and is influenced by species j (preprocessed representation). */
	MatrixType SpeciesStructure(size_t stage) const { 
		return BooleanMatrix<EXT>(Global(stage).StoichiometricMatrix())*RateStructure(stage);
	}

	//dtor
	virtual ~IKineticContainer(){}	
 };

 //ensure backward compatibility for older naming convention
 template<typename EXT>
 using IKineticStoichiometry = IKineticContainer<EXT>;
} //namespace StoichPack

#endif
