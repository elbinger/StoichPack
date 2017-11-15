/* File: StoichPackStageInfo.h
 * Purpose: Provide basic informations about stages.
 * Author: Tobias Elbinger <elbinger@math.fau.de> */

//READ README_CONCEPTS

#ifndef __H_STOICHPACK_STAGEINFO__
#define __H_STOICHPACK_STAGEINFO__

#include "StoichPackUtility.h"
/* class BasicStageInfo:
 * Provide information about one stage. This information concerns either only mobile species, only immobile species or all species. */
namespace StoichPack{
 template<typename EXT>
 class BasicStageInfo{
 public:
	typedef typename EXT::MatrixType MatrixType;

	/* Constructor:
	 * Parameters:
	 * * Stoichiometric: the stoichiometric matrix of the stage.
	 * * ToOriginal: a matrix describing the contribution of this stage to the original solution.
	 * * FromOriginal: a matrix describing the contribution of an original solution the values in this stage.
	 * * n_mobile: number of mobile species treated in this class.
	 * * allmobile: number of mobile species in this stage (this differs from n_mobile if the information concerns at least
	 *              1 immobile species).
	 * * bool Correct: is it allowed to correct the values in this stage between 2 iterations?
	 *                 (not needed in this version)*/
	BasicStageInfo(const MatrixType& Stoichiometric, const MatrixType& ToOriginal, const MatrixType& FromOriginal,
			size_t n_mobile, size_t allmobile/*, bool Correct*/)
		: stoichiometric(Stoichiometric), toOriginal(ToOriginal), fromOriginal(FromOriginal),
		  mobile(n_mobile), origmobile(allmobile)/*, correct(Correct)*/ {
		//check
		if(MobileSpecies()>GlobalSpecies()) throw StoichPackException("Inconsitent parameters");
		if(EXT::rows(fromOriginal)!=GlobalSpecies()) throw StoichPackException("Inconsitent parameters");
		if(EXT::cols(toOriginal)!=GlobalSpecies()) throw StoichPackException("Inconsitent parameters");
	}

	/* return a BasicStageInfo object representing only the information about the mobile species. The matrices in the return value
	 * are subblocks of the matrices in *this (see StoichPackSpecies.h: mobile species are stored in front of immobile species). */
	BasicStageInfo<EXT> MobilePart() const {
		const MatrixType stoich_mobile = DivideRows<EXT>(stoichiometric,mobile).Mobile();
		const MatrixType toOriginal_mobile = SplitBlocks<EXT>(toOriginal,origmobile,mobile).Mobile();
		const MatrixType fromOriginal_mobile = SplitBlocks<EXT>(fromOriginal,mobile,origmobile).Mobile();
		return BasicStageInfo<EXT>(stoich_mobile,toOriginal_mobile,fromOriginal_mobile,mobile,origmobile/*,correct*/);
	}

	/* return a BasicStageInfo object representing only the information about the immobile species. The matrices in the return value
	 * are subblocks of the matrices in *this (see StoichPackSpecies.h: mobile species are stored in front of immobile species). */
	BasicStageInfo<EXT> ImmobilePart() const {
		const MatrixType stoich_immobile = DivideRows<EXT>(stoichiometric,mobile).Immobile();
		const MatrixType toOriginal_immobile = SplitBlocks<EXT>(toOriginal,origmobile,mobile).Immobile();
		const MatrixType fromOriginal_immobile = SplitBlocks<EXT>(fromOriginal,mobile,origmobile).Immobile();
		return BasicStageInfo<EXT>(stoich_immobile,toOriginal_immobile,fromOriginal_immobile,0,origmobile/*,correct*/);
	}

	//getters
	const MatrixType& StoichiometricMatrix() const { return stoichiometric; }
	const MatrixType& ToOriginal() const { return toOriginal; }
	const MatrixType& FromOriginal() const { return fromOriginal; }

	size_t GlobalSpecies() const { return EXT::rows(stoichiometric); }
	size_t MobileSpecies() const { return mobile; }
	size_t ImmobileSpecies() const { return GlobalSpecies()-MobileSpecies(); }

	//bool Correct() const { return correct; } //not needed in this version
 private:
	MatrixType stoichiometric; //stoichiometric matrix of this stage
	MatrixType toOriginal; //a matrix describing the contribution of this stage to the original solution
	MatrixType fromOriginal; //a matrix describing the contribution of an original solution the values in this stage
	size_t mobile, origmobile; //number of mobile species treated in this class / in this stage
	                           //(this differs if the information concerns at least 1 immobile species).

	//bool correct; //is it allowed to correct the values in this stage between 2 iterations? // not needed in this version

	BasicStageInfo(); //FORBID
 };

 // a combination of 3 BasicStageInfo objects: 1 for all species, 1 for mobile species only, one for immobile species only.
 template<typename EXT>
 class StageInfo{
 public:
	explicit StageInfo(const BasicStageInfo<EXT>& Original)
	                   : global(Original), mobile(Original.MobilePart()), immobile(Original.ImmobilePart()) {}

	const BasicStageInfo<EXT>& Global() const { return global; }
	const BasicStageInfo<EXT>& Mobile() const { return mobile; }
	const BasicStageInfo<EXT>& Immobile() const { return immobile; }

	//bool Correct() const { return global.Correct(); } //not needed in this version

 private:
	BasicStageInfo<EXT> global, mobile, immobile;
	StageInfo(); //FORBID
 };

 //Information about all stages. Basically, this is a std::vector<StageInfo> with additional functions
 template<typename EXT>
 class StageInfoArray{
 public:
	typedef StageInfo<EXT> entry_t;
	typedef std::vector<entry_t> container_t;
	typedef typename EXT::VectorType VectorType;
	typedef typename EXT::VectorArrayType VectorArrayType;

	//copy from std::vector
	StageInfoArray(const std::vector<BasicStageInfo<EXT> >& Original) {
		if(Original.size()==0) throw StoichPackException("At least 1 stage!");
		entries.reserve(Original.size());
		for(auto x : Original) entries.push_back(entry_t(x));
	}

	//sum over all stages of all species (=number of species in system)
	size_t GlobalSpecies() const {
		size_t result = 0;
		for(const entry_t& x : entries) result+=x.Global().GlobalSpecies();
		return result;
	}
	//sum over all stages of mobile species (=number of mobile species in system)
	size_t MobileSpecies() const {
		size_t result = 0;
		for(const entry_t& x : entries) result+=x.Global().MobileSpecies();
		return result;
	}
	//sum over all stages of immobile species (=number of immobile species in system)
	size_t ImmobileSpecies() const {
		size_t result = 0;
		for(const entry_t& x : entries) result+=x.Global().ImmobileSpecies();
		return result;
	}

	//Calculate the corresponding original solution for a given preprocessed solution (all species)
	VectorType ToOriginalGlobal(const VectorArrayType& all) const {
		VectorType result = EXT::CreateVector(GlobalSpecies(),0);
		for(size_t i=0;i<Stages();++i) {
			result+=entries[i].Global().ToOriginal()*all[i];
		}
		return result;
	}

	//Calculate the corresponding original solution for a given preprocessed solution (mobile species)
	VectorType ToOriginalMobile(const VectorArrayType& mobile) const {
		VectorType result = EXT::CreateVector(MobileSpecies(),0);
		for(size_t i=0;i<Stages();++i) result+=entries[i].Mobile().ToOriginal()*mobile[i];
		return result;
	}

	//Calculate the corresponding original solution for a given preprocessed solution (immobile species)	
	VectorType ToOriginalImmobile(const VectorArrayType& immobile) const {
		VectorType result = EXT::CreateVector(ImmobileSpecies(),0);
		for(size_t i=0;i<Stages();++i) result+=entries[i].Immobile().ToOriginal()*immobile[i];
		return result;
	}

	//Calculate the corresponding preprocessed solution for a given original solution (all species)
	VectorArrayType FromOriginalGlobal(const VectorType& all) const {
		VectorArrayType result = EXT::ReserveVectorArray(Stages());
		for(const StageInfo<EXT>& x : entries) EXT::PushBack(result,x.Global().FromOriginal()*all);
		return result;
	}

	//Calculate the corresponding preprocessed solution for a given original solution (mobile species)
	VectorArrayType FromOriginalMobile(const VectorType& mobile) const {
		VectorArrayType result = EXT::ReserveVectorArray(Stages());
		for(const StageInfo<EXT>& x : entries) EXT::PushBack(result,x.Mobile().FromOriginal()*mobile);
		return result;
	}

	//Calculate the corresponding preprocessed solution for a given original solution (immobile species)
	VectorArrayType FromOriginalImmobile(const VectorType& immobile) const {
		VectorArrayType result = EXT::ReserveVectorArray(Stages());
		for(const StageInfo<EXT>& x : entries) EXT::PushBack(result,x.Immobile().FromOriginal()*immobile);
		return result;
	}

	//getters
	const entry_t& operator[](size_t i) const { return entries[i]; }
	size_t size() const { return entries.size(); }
	size_t Stages() const { return entries.size(); }
	typename container_t::const_iterator begin() const { return entries.begin(); }
	typename container_t::const_iterator end() const { return entries.end(); }

 private:
	container_t entries;
 };

}

#endif
