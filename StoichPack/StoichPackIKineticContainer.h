#ifndef __H_STOICHPACK_IKINETIC_CONTAINER__
#define __H_STOICHPACK_IKINETIC_CONTAINER__

#include "StoichPackUtility.h"
#include "StoichPackSpecies.h"

namespace StoichPack{

 template<typename EXT >
 class IKineticContainer{
	typedef typename EXT::VectorType VectorType;
	typedef typename EXT::VectorPairType VectorPairType;
	typedef typename EXT::VectorArrayType VectorArrayType;
	typedef typename EXT::VectorArrayPairType VectorArrayPairType;
	typedef typename EXT::MatrixType MatrixType;
	typedef typename EXT::MatrixPairType MatrixPairType;
	typedef typename EXT::MatrixQuadType MatrixQuadType;

 private:
	std::vector<MatrixType> stoich_mobile, stoich_immobile, stoich_all;
	size_t mobile_species, immobile_species;
	std::vector<bool> correction;

	IKineticContainer& operator=(const IKineticContainer<EXT>&);	//FORBID

 protected:
	void AddStage(const MatrixType& stoich, size_t n_mobile, bool correct){
		if(n_mobile>EXT::rows(stoich)) throw StoichPackException("Illegal values in IKineticContainer::AddStage!");

		if(EXT::rows(stoich)+EXT::cols(stoich)==0) throw StoichPackException("Cannot add empty stage");

		stoich_all.push_back(stoich);

		MatrixPairType tmp = EXT::DivideRows(stoich,n_mobile);
		stoich_mobile.push_back(tmp.Mobile());
		stoich_immobile.push_back(tmp.Immobile());

		mobile_species+=n_mobile;
		immobile_species+=ImmobileSpecies(Stages()-1);
		correction.push_back(correct);
	}

	void AddStage(const MatrixType& stoich_mob, const MatrixType& stoich_immob, bool correct){
		if(EXT::cols(stoich_mob)!=EXT::cols(stoich_immob)) throw StoichPackException("Illegal values in IKineticContainer::AddStage!");
		AddStage(EXT::CombineRows(stoich_mob,stoich_immob),EXT::rows(stoich_mob),correct);
	}

 public:
	IKineticContainer() : mobile_species(0), immobile_species(0) {}

	size_t Stages() const { return stoich_all.size(); }

	const std::vector<MatrixType>& StoichiometricMatrices() const { return stoich_all; }
	const std::vector<MatrixType>& MobileStoichiometricMatrices() const { return stoich_mobile; }
	const std::vector<MatrixType>& ImmobileStoichiometricMatrices() const { return stoich_immobile; }

	size_t MobileSpecies(size_t stage) const { assert(stage<Stages()); return EXT::rows(stoich_mobile[stage]); }
	size_t MobileSpecies() const { return mobile_species; }
	size_t ImmobileSpecies() const { return immobile_species; }
	size_t ImmobileSpecies(size_t stage) const { assert(stage<Stages()); return EXT::rows(stoich_immobile[stage]); }
	size_t AllSpecies(size_t stage) const { return MobileSpecies(stage)+ImmobileSpecies(stage); }
	size_t AllSpecies() const { return MobileSpecies()+ImmobileSpecies(); }

	size_t Reactions(size_t stage) const { assert(stage<Stages()); return EXT::cols(stoich_all[stage]); }

  //debug
	bool CheckSize(const VectorType& x) const {
		return EXT::size(x)==AllSpecies();
	}
	bool CheckSizeMobile(const VectorType& x) const {
		return EXT::size(x)==MobileSpecies();
	}
	bool CheckSizeImmobile(const VectorType& x) const {
		return EXT::size(x)==ImmobileSpecies();
	}


	bool CheckSize(const VectorType& x, const VectorType& y) const {
		return EXT::size(x)==MobileSpecies() && EXT::size(y)==ImmobileSpecies();
	}

	bool CheckSize(const VectorPairType& x) const {
		return CheckSize(x.Mobile(),x.Immobile());
	}

	bool CheckSize(const VectorArrayType& x) const {
		if(EXT::size(x)!=Stages()) return false;

		typename EXT::ConstVectorArrayIteratorType it = EXT::ConstBegin(x);
		for(size_t i=0; i<Stages(); ++it, ++i) {
			if(EXT::size(*it)!=AllSpecies(i)) return false;
		}
		return true;
	}

	bool CheckSize(const VectorArrayType& x, const VectorArrayType& y) const {
		if(EXT::size(x)!=Stages() || EXT::size(y)!=Stages()) return false;

		typename EXT::ConstVectorArrayIteratorType it1 = EXT::ConstBegin(x);
		typename EXT::ConstVectorArrayIteratorType it2 = EXT::ConstBegin(x);
		for(size_t i=0; i<Stages(); ++it1, ++it2, ++i) {
			if(EXT::size(*it1)!=MobileSpecies(i) || EXT::size(*it2)!=ImmobileSpecies(i)) return false;
		}
		return true;
	}

	bool CheckSize(const VectorArrayPairType& x) const {
		return CheckSize(x.Mobile(),x.Immobile());
	}

	virtual void FromOriginalImpl(const VectorType& all, VectorArrayType& ret) const =0;
	virtual void FromOriginalMobileImpl(const VectorType& mobile, VectorArrayType& ret) const =0;
	virtual void FromOriginalImmobileImpl(const VectorType& immobile, VectorArrayType& ret) const=0;

	VectorArrayType FromOriginal(const VectorType& all) const {
		VectorArrayType ret = EXT::ReserveVectorArray(Stages());
		FromOriginalImpl(all,ret);
		return ret;
	}
	VectorArrayType FromOriginalMobile(const VectorType& mobile) const {
		VectorArrayType ret = EXT::ReserveVectorArray(Stages());
		FromOriginalMobileImpl(mobile,ret);
		return ret;
	}
	VectorArrayType FromOriginalImmobile(const VectorType& immobile) const {
		VectorArrayType ret = EXT::ReserveVectorArray(Stages());
		FromOriginalImmobileImpl(immobile,ret);
		return ret;
	}
	VectorArrayPairType FromOriginal(const VectorType& mobile, const VectorType& immobile) const {
		const size_t s=Stages();
		VectorArrayPairType ret(EXT::ReserveVectorArray(s),EXT::ReserveVectorArray(s));
		FromOriginalMobileImpl(mobile,ret.Mobile());
		FromOriginalImmobileImpl(immobile,ret.Immobile());
		return ret;
	}
	VectorArrayPairType FromOriginal(const VectorPairType& all) const { return FromOriginal(all.Mobile(),all.Immobile()); }

	virtual void ToOriginalImpl(const VectorArrayType& all, VectorType& orig) const =0;
	virtual void ToOriginalMobileImpl(const VectorArrayType& mobile, VectorType& orig) const =0;
	virtual void ToOriginalImmobileImpl(const VectorArrayType& immobile, VectorType& orig) const =0;

	VectorType ToOriginal(const VectorArrayType& all) const {
		VectorType ret = EXT::CreateVector(AllSpecies());
		ToOriginalImpl(all,ret);
		return ret;
	}
	VectorType ToOriginalMobile(const VectorArrayType& mobile) const{
		VectorType ret = EXT::CreateVector(MobileSpecies());
		ToOriginalMobileImpl(mobile,ret);
		return ret;
	}
	VectorType ToOriginalImmobile(const VectorArrayType& immobile) const{
		VectorType ret = EXT::CreateVector(ImmobileSpecies());
		ToOriginalImmobileImpl(immobile,ret);
		return ret;
	}
	VectorPairType ToOriginal(const VectorArrayType& mobile, const VectorArrayType& immobile) const {
		return VectorPairType(ToOriginalMobile(mobile),ToOriginalImmobile(immobile));
	}
	VectorPairType ToOriginal(const VectorArrayPairType& all) const {
		return ToOriginal(all.Mobile(),all.Immobile());
	}

	virtual VectorType ReactionRatesImpl1(const VectorType& all, size_t stage) const =0;
	virtual VectorType ReactionRatesImpl2(const VectorPairType& all, size_t stage) const =0;
	virtual VectorType SubReactionRatesImpl1(const VectorType& all, const std::vector<size_t>& I, size_t stage) const =0;
	virtual VectorType SubReactionRatesImpl2(const VectorPairType& all, const std::vector<size_t>& I, size_t stage) const =0;

	virtual MatrixType DiffReactionRatesImpl1(const VectorType& all, size_t stage) const =0;
	virtual MatrixPairType DiffReactionRatesImpl2(const VectorPairType& all, size_t stage) const =0;
	virtual MatrixType DiffSubReactionRatesImpl1(const VectorType& all, const std::vector<size_t>& I, size_t stage) const=0;
	virtual MatrixPairType DiffSubReactionRatesImpl2(const VectorPairType& all, const std::vector<size_t>& I, size_t stage) const=0;

	virtual VectorType ConstSpeciesRatesImpl1(const VectorType& all, size_t stage) const =0;
	virtual VectorPairType ConstSpeciesRatesImpl2(const VectorPairType& all, size_t stage) const =0;
	virtual VectorType ConstMobileSpeciesRatesImpl(const VectorPairType& all, size_t stage) const =0;
	virtual VectorType ConstImmobileSpeciesRatesImpl(const VectorPairType& all, size_t stage) const =0;

	virtual bool ApplyMobileCorrectionImpl(VectorArrayType& mobile, const VectorArrayType& allowed) const =0;
	virtual bool ApplyImmobileCorrectionImpl(VectorArrayType& immobile, const VectorArrayType& allowed) const =0;
	virtual bool ApplyCorrectionImpl(VectorArrayType& all, const VectorArrayType& allowed) const =0;

	virtual const std::vector<Species>& Participants() const =0;

	//virtual MatrixType RateStructure(size_t stage, const std::vector<size_t>& I) const =0;
	//virtual MatrixType RateStructure(size_t stage) const =0;
	virtual MatrixType SpeciesStructure(size_t stage) const=0;

	virtual MatrixType MobileTransformation() const =0;
	virtual MatrixType ImmobileTransformation() const =0;
	virtual MatrixType Transformation() const =0;

	virtual IKineticContainer<EXT>* copy() const =0;

	VectorType ReactionRates(const VectorArrayType& all, size_t stage) const {
		return ReactionRatesImpl1(ToOriginal(all),stage);
	}

	VectorType SubReactionRates(const VectorArrayType& all, const std::vector<size_t>& I, size_t stage) const {
		return SubReactionRatesImpl1(ToOriginal(all),I,stage);
	}

	VectorType ReactionRates(const VectorArrayType& mobile, const VectorArrayType& immobile, size_t stage) const {
		return ReactionRatesImpl2(ToOriginal(mobile,immobile),stage);
	}

	VectorType SubReactionRates(const VectorArrayType& mobile, const VectorArrayType& immobile,
	                                const std::vector<size_t>& I, size_t stage) const {
		return SubReactionRatesImpl2(ToOriginal(mobile,immobile),I,stage);
	}

	MatrixType DiffReactionRates(const VectorArrayType& all, size_t stage) const {
		return DiffReactionRatesImpl1(ToOriginal(all),stage);
	}

	MatrixType DiffSubReactionRates(const VectorArrayType& all, const std::vector<size_t>& I, size_t stage) const {
		return DiffSubReactionRatesImpl1(ToOriginal(all),I,stage);
	}

	MatrixPairType DiffReactionRates(const VectorArrayType& mobile, const VectorArrayType& immobile, size_t stage) const {
		return DiffReactionRatesImpl2(ToOriginal(mobile,immobile),stage);
	}

	MatrixPairType DiffSubReactionRates(const VectorArrayType& mobile, const VectorArrayType& immobile,
	                                            const std::vector<size_t>& I, size_t stage) const {
		return DiffSubReactionRatesImpl2(ToOriginal(mobile,immobile),I,stage);
	}

	VectorType ConstSpeciesRates(const VectorArrayType& all, size_t stage) const {
		return ConstSpeciesRatesImpl1(ToOriginal(all),stage);
	}

	VectorPairType ConstSpeciesRates(const VectorArrayType& mobile, const VectorArrayType& immobile, size_t stage) const {
		return ConstSpeciesRatesImpl2(ToOriginal(mobile,immobile),stage);
	}

	VectorType ConstMobileSpeciesRates(const VectorArrayType& mobile, const VectorArrayType& immobile, size_t stage) const {
		return ConstMobileSpeciesRatesImpl(ToOriginal(mobile,immobile),stage);
	}

	VectorType ConstImmobileSpeciesRates(const VectorArrayType& mobile, const VectorArrayType& immobile, size_t stage) const {
		return ConstImmobileSpeciesRatesImpl(ToOriginal(mobile,immobile),stage);
	}
	
	VectorType SpeciesRates(const VectorArrayType& all, size_t stage) const {
		const VectorType orig = ToOriginal(all);
		return stoich_all[stage]*ReactionRatesImpl1(orig,stage)+ConstSpeciesRatesImpl1(orig,stage);
	}

	VectorPairType SpeciesRates(const VectorArrayType& mobile, const VectorArrayType& immobile, size_t stage) const {
		const VectorPairType orig = ToOriginal(mobile,immobile);
		const VectorPairType crates = ConstSpeciesRatesImpl2(orig,stage);
		const VectorType rrates = ReactionRatesImpl2(orig,stage);
		return VectorPairType(stoich_mobile[stage]*rrates+crates.Mobile(), stoich_immobile[stage]*rrates+crates.Immobile());
	}

	VectorType MobileSpeciesRates(const VectorArrayType& mobile, const VectorArrayType& immobile, size_t stage) const {
		const VectorPairType orig = ToOriginal(mobile,immobile);
		return stoich_mobile[stage]*ReactionRatesImpl2(orig,stage)+ConstMobileSpeciesRatesImpl(orig,stage);
	}

	VectorType ImmobileSpeciesRates(const VectorArrayType& mobile, const VectorArrayType& immobile, size_t stage) const {
		const VectorPairType orig = ToOriginal(mobile,immobile);
		return stoich_immobile[stage]*ReactionRatesImpl2(orig,stage)+ConstImmobileSpeciesRatesImpl(orig,stage);
	}

	MatrixType DiffSpeciesRates(const VectorArrayType& all, size_t stage) const {
		return stoich_all[stage]*DiffReactionRatesImpl1(ToOriginal(all),stage);
	}

	MatrixQuadType DiffSpeciesRates(const VectorArrayType& mobile, const VectorArrayType& immobile, size_t stage) const {
		MatrixPairType tmp = DiffReactionRatesImpl2(ToOriginal(mobile,immobile),stage);
		return MatrixQuadType(MatrixPairType(stoich_mobile[stage]*tmp.Mobile(),stoich_mobile[stage]*tmp.Immobile()),
		                      MatrixPairType(stoich_immobile[stage]*tmp.Mobile(),stoich_immobile[stage]*tmp.Immobile()));
	}

	MatrixPairType DiffMobileSpeciesRates(const VectorArrayType& mobile, const VectorArrayType& immobile,size_t stage) const{
		MatrixPairType tmp = DiffReactionRatesImpl2(ToOriginal(mobile,immobile),stage);
		return MatrixPairType(MatrixPairType(stoich_mobile[stage]*tmp.Mobile(),stoich_mobile[stage]*tmp.Immobile()));
	}

	MatrixPairType DiffImmobileSpeciesRates(const VectorArrayType& mobile, const VectorArrayType& immobile,size_t stage) const{
		MatrixPairType tmp = DiffReactionRatesImpl2(ToOriginal(mobile,immobile),stage);
		return MatrixPairType(MatrixPairType(stoich_immobile[stage]*tmp.Mobile(),stoich_immobile[stage]*tmp.Immobile()));
	}

	bool ApplyCorrection(VectorArrayType& all, size_t stage) const {
		const size_t s=Stages();
		assert(stage<s);

		if(!CorrectionAllowed(stage)) return true;

		VectorArrayType allowed = EXT::ReserveVectorArray(s);
		for(size_t i=0;i<stage;++i) EXT::PushBack(allowed,EXT::CreateZeroVector(AllSpecies(i)));
		for(size_t i=stage;i<s;++i) EXT::PushBack(allowed,EXT::CreateVector(AllSpecies(i),1));

		return ApplyCorrectionImpl(all,allowed);
	}

	bool ApplyMobileCorrection(VectorArrayType& mobile, size_t stage) const {
		const size_t s=Stages();
		assert(stage<s);

		if(!CorrectionAllowed(stage)) return true;

		VectorArrayType allowed = EXT::ReserveVectorArray(s);
		for(size_t i=0;i<stage;++i) EXT::PushBack(allowed,EXT::CreateZeroVector(MobileSpecies(i)));
		for(size_t i=stage;i<s;++i) EXT::PushBack(allowed,EXT::CreateVector(MobileSpecies(i),1));

		return ApplyMobileCorrectionImpl(mobile,allowed);
	}

	bool ApplyImmobileCorrection(VectorArrayType& immobile, size_t stage) const {
		const size_t s=Stages();
		assert(stage<s);

		if(!CorrectionAllowed(stage)) return true;

		VectorArrayType allowed = EXT::ReserveVectorArray(s);
		for(size_t i=0;i<stage;++i) EXT::PushBack(allowed,EXT::CreateZeroVector(ImmobileSpecies(i)));
		for(size_t i=stage;i<s;++i) EXT::PushBack(allowed,EXT::CreateVector(ImmobileSpecies(i),1));

		return ApplyImmobileCorrectionImpl(immobile,allowed);
	}

	bool CorrectionAllowed(size_t stage) const { return correction[stage]; }
	//dtor
	virtual ~IKineticContainer(){}	
 };

 template<typename EXT>
 using IKineticStoichiometry = IKineticContainer<EXT>;
} //namespace StoichPack

#endif
