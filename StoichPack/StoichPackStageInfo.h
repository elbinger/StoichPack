#ifndef __H_STOICHPACK_STAGEINFO__
#define __H_STOICHPACK_STAGEINFO__

#include "StoichPackUtility.h"

namespace StoichPack{
 template<typename EXT>
 class BasicStageInfo{
 public:
	typedef typename EXT::MatrixType MatrixType;

	BasicStageInfo(const MatrixType& Stoichiometric, const MatrixType& ToOriginal, const MatrixType& FromOriginal,
			size_t n_mobile, size_t allmobile, bool c)
		: stoichiometric(Stoichiometric), toOriginal(ToOriginal), fromOriginal(FromOriginal),
		  mobile(n_mobile), origmobile(allmobile), correct(c) {
		if(MobileSpecies()>GlobalSpecies()) throw StoichPackException("Inconsitent parameters");
		if(EXT::rows(fromOriginal)!=GlobalSpecies()) throw StoichPackException("Inconsitent parameters");
		if(EXT::cols(toOriginal)!=GlobalSpecies()) throw StoichPackException("Inconsitent parameters");
	}

	BasicStageInfo(const BasicStageInfo<EXT>& mob, const BasicStageInfo<EXT>& immob)
			: stoichiometric(CombineRows<EXT>(mob.stoichiometric,immob.stoichiometric)),
			  toOriginal(CombineBlocks<EXT>(mob.ToOriginal(),immob.ToOriginal())),
			  fromOriginal(CombineBlocks<EXT>(mob.FromOriginal(),immob.FromOriginal())),
			  mobile(mob.mobile), origmobile(mob.origmobile), correct(mob.correct)
		{ assert(immob.correct==correct); }

	const MatrixType& StoichiometricMatrix() const { return stoichiometric; }
	const MatrixType& ToOriginal() const { return toOriginal; }
	const MatrixType& FromOriginal() const { return fromOriginal; }

	size_t GlobalSpecies() const { return EXT::rows(stoichiometric); }
	size_t MobileSpecies() const { return mobile; }
	size_t ImmobileSpecies() const { return GlobalSpecies()-MobileSpecies(); }

	BasicStageInfo<EXT> MobilePart() const {
		const MatrixType stoich_mobile = DivideRows<EXT>(stoichiometric,mobile).Mobile();
		const MatrixType toOriginal_mobile = SplitBlocks<EXT>(toOriginal,origmobile,mobile).Mobile();
		const MatrixType fromOriginal_mobile = SplitBlocks<EXT>(fromOriginal,mobile,origmobile).Mobile();
		return BasicStageInfo<EXT>(stoich_mobile,toOriginal_mobile,fromOriginal_mobile,mobile,origmobile,correct);
	}

	BasicStageInfo<EXT> ImmobilePart() const {
		const MatrixType stoich_immobile = DivideRows<EXT>(stoichiometric,mobile).Immobile();
		const MatrixType toOriginal_immobile = SplitBlocks<EXT>(toOriginal,origmobile,mobile).Immobile();
		const MatrixType fromOriginal_immobile = SplitBlocks<EXT>(fromOriginal,mobile,origmobile).Immobile();
		return BasicStageInfo<EXT>(stoich_immobile,toOriginal_immobile,fromOriginal_immobile,0,origmobile,correct);
	}

	bool Correct() const { return correct; }
 private:
	MatrixType stoichiometric, toOriginal, fromOriginal;
	size_t mobile, origmobile;
	bool correct;
	BasicStageInfo(); //FORBID
 };

 template<typename EXT>
 class StageInfo{
 public:
	explicit StageInfo(const BasicStageInfo<EXT>& Original)
	                   : global(Original), mobile(Original.MobilePart()), immobile(Original.ImmobilePart()) {}

	StageInfo(const BasicStageInfo<EXT>& mob, const BasicStageInfo<EXT>& immob, bool coorection)
		: global(mob,immob), mobile(mob), immobile(immob) {}

	const BasicStageInfo<EXT>& Global() const { return global; }
	const BasicStageInfo<EXT>& Mobile() const { return mobile; }
	const BasicStageInfo<EXT>& Immobile() const { return immobile; }

	bool Correct() const { return global.Correct(); }

 private:
	BasicStageInfo<EXT> global, mobile, immobile;
	StageInfo(); //FORBID
 };

 template<typename EXT>
 class StageInfoArray{
 public:
	typedef StageInfo<EXT> entry_t;
	typedef std::vector<entry_t> container_t;
	typedef typename EXT::VectorType VectorType;
	typedef typename EXT::VectorArrayType VectorArrayType;

	StageInfoArray(const std::vector<BasicStageInfo<EXT> >& Original) {
		if(Original.size()==0) throw StoichPackException("At least 1 stage!");
		entries.reserve(Original.size());
		for(auto x : Original) entries.push_back(StageInfo<EXT>(x));
	}

	const entry_t& operator[](size_t i) const { return entries[i]; }
	size_t size() const { return entries.size(); }

	typename container_t::const_iterator begin() const { return entries.begin(); }
	typename container_t::const_iterator end() const { return entries.end(); }

	size_t GlobalSpecies() const {
		size_t result = 0;
		for(const entry_t& x : entries) result+=x.Global().GlobalSpecies();
		return result;
	}
	size_t MobileSpecies() const {
		size_t result = 0;
		for(const entry_t& x : entries) result+=x.Global().MobileSpecies();
		return result;
	}
	size_t ImmobileSpecies() const {
		size_t result = 0;
		for(const entry_t& x : entries) result+=x.Global().ImmobileSpecies();
		return result;
	}

	size_t Stages() const { return entries.size(); }

	VectorType ToOriginalGlobal(const VectorArrayType& all) const {
		VectorType result = EXT::CreateVector(GlobalSpecies(),0);
		for(size_t i=0;i<Stages();++i) {
			result+=entries[i].Global().ToOriginal()*all[i];
		}
		return result;
	}

	VectorType ToOriginalMobile(const VectorArrayType& mobile) const {
		VectorType result = EXT::CreateVector(MobileSpecies(),0);
		for(size_t i=0;i<Stages();++i) result+=entries[i].Mobile().ToOriginal()*mobile[i];
		return result;
	}
	
	VectorType ToOriginalImmobile(const VectorArrayType& immobile) const {
		VectorType result = EXT::CreateVector(ImmobileSpecies(),0);
		for(size_t i=0;i<Stages();++i) result+=entries[i].Immobile().ToOriginal()*immobile[i];
		return result;
	}

	VectorArrayType FromOriginalGlobal(const VectorType& all) const {
		VectorArrayType result = EXT::ReserveVectorArray(Stages());
		for(const StageInfo<EXT>& x : entries) EXT::PushBack(result,x.Global().FromOriginal()*all);
		return result;
	}

	VectorArrayType FromOriginalMobile(const VectorType& mobile) const {
		VectorArrayType result = EXT::ReserveVectorArray(Stages());
		for(const StageInfo<EXT>& x : entries) EXT::PushBack(result,x.Mobile().FromOriginal()*mobile);
		return result;
	}

	VectorArrayType FromOriginalImmobile(const VectorType& immobile) const {
		VectorArrayType result = EXT::ReserveVectorArray(Stages());
		for(const StageInfo<EXT>& x : entries) EXT::PushBack(result,x.Immobile().FromOriginal()*immobile);
		return result;
	}


 private:
	container_t entries;
 };

}

#endif
