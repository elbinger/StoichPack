#ifndef __H_STOICHPACK_HIERARCHICAL_STOICHIOMETRY__
#define __H_STOICHPACK_HIERARCHICAL_STOICHIOMETRY__

#include "StoichPackIKineticStoichiometry.h"
#include <memory>

namespace StoichPack{

template<typename EXT, typename BT>
class BaseStorage{
	private:
	const BT storage;
	public:
	BaseStorage(const BT& x) : storage(x) {}
	const BT& get() const { return storage; }
};

template<typename EXT>
class BaseStorage<EXT,IKineticStoichiometry<EXT> >{
	private:
	const std::shared_ptr<IKineticStoichiometry<EXT> > storage;
	public:
	BaseStorage(const IKineticStoichiometry<EXT>& x) : storage(x.copy()) {}
	const IKineticStoichiometry<EXT>& get() const { return *storage; }
};

template<typename EXT, typename BT = IKineticStoichiometry<EXT> >
class IHierarchicalStoichiometry : public IKineticStoichiometry<EXT>{
 	typedef typename EXT::VectorType VectorType;
	typedef typename EXT::VectorPairType VectorPairType;
	typedef typename EXT::VectorArrayType VectorArrayType;
	typedef typename EXT::VectorArrayPairType VectorArrayPairType;
	typedef typename EXT::MatrixType MatrixType;
private:
	std::vector<size_t> substages, substages_first;
	const BaseStorage<EXT,BT> base;

	void UpdateSubStages(size_t subs) {
		assert(subs<Base().Stages());
		if(substages.size()==0){
			assert(subs==0);
			substages_first.push_back(0);
		} else if(substages.back()!=subs){
			assert(substages.back()==subs-1);
			substages_first.push_back(Stages());
		}
		substages.push_back(subs);
	}

	//FORBID
	IHierarchicalStoichiometry();
protected:
	IHierarchicalStoichiometry(const BT& bt) : base(bt) {}

	void Finish() {
		assert(substages_first.size()==Base().Stages());
		substages_first.push_back(Stages());
	}

	void AddStage(const MatrixType& stoich, size_t n_mobile, size_t subs, bool correction){
		UpdateSubStages(subs);
		IKineticStoichiometry<EXT>::AddStage(stoich,n_mobile,correction);
	}

	void AddStage(const MatrixType& stoich_mob, const MatrixType& stoich_immob, size_t subs, bool correction){
		UpdateSubStages(subs);
		IKineticStoichiometry<EXT>::AddStage(stoich_mob,stoich_immob, correction);
	}

	size_t SubStage(size_t i) const { assert(i<Stages()); return substages[i]; }
	size_t SubStageFirst(size_t i) const { assert(i<substages_first.size()); return substages_first[i]; }
	const BT& Base() const { return base.get(); }

public:
	using IKineticStoichiometry<EXT>::Stages;
	using IKineticStoichiometry<EXT>::CheckSize;

	virtual void FromBaseImpl(const VectorArrayType& all, VectorArrayType& ret) const =0;
	virtual void FromBaseMobileImpl(const VectorArrayType& mobile, VectorArrayType& ret) const =0;
	virtual void FromBaseImmobileImpl(const VectorArrayType& immobile, VectorArrayType& ret) const=0;

	VectorArrayType FromBase(const VectorArrayType& all) const {
		VectorArrayType ret = EXT::ReserveVectorArray(Stages());
		FromBaseImpl(all,ret);
		return ret;
	}

	VectorArrayType FromBaseMobile(const VectorArrayType& mobile) const {
		VectorArrayType ret = EXT::ReserveVectorArray(Stages());
		FromBaseMobileImpl(mobile,ret);
		return ret;
	}
	VectorArrayType FromBaseImmobile(const VectorArrayType& immobile) const {
		VectorArrayType ret = EXT::ReserveVectorArray(Stages());
		FromBaseImmobileImpl(immobile,ret);
		return ret;
	}
	void FromOriginalImpl(const VectorType& all, VectorArrayType& ret) const {
		const VectorArrayType fbase=Base().FromOriginal(all);
		FromBaseImpl(fbase,ret);
	}
		
	void FromOriginalMobileImpl(const VectorType& mobile, VectorArrayType& ret) const {
		const VectorArrayType fbase=Base().FromOriginalMobile(mobile);
		FromBaseMobileImpl(fbase,ret);
	}

	void FromOriginalImmobileImpl(const VectorType& immobile, VectorArrayType& ret) const{
		const VectorArrayType fbase=Base().FromOriginalImmobile(immobile);
		FromBaseImmobileImpl(fbase,ret);
	}
	
	virtual void ToBaseImpl(const VectorArrayType& all, VectorArrayType& ret) const =0;
	virtual void ToBaseMobileImpl(const VectorArrayType& mobile, VectorArrayType& ret) const=0;
	virtual void ToBaseImmobileImpl(const VectorArrayType& immobile, VectorArrayType& ret) const=0;

	VectorArrayType ToBase(const VectorArrayType& all) const {
		VectorArrayType ret=EXT::ReserveVectorArray(Base().Stages());
		ToBaseImpl(all,ret);
		return ret;
	}
	VectorArrayType ToBaseMobile(const VectorArrayType& mobile) const {
		VectorArrayType ret=EXT::ReserveVectorArray(Base().Stages());
		ToBaseMobileImpl(mobile,ret);
		return ret;
	}
	VectorArrayType ToBaseImmobile(const VectorArrayType& immobile) const {
		VectorArrayType ret=EXT::ReserveVectorArray(Base().Stages());
		ToBaseImmobileImpl(immobile,ret);
		return ret;
	}

	VectorArrayPairType ToBase(const VectorArrayType& mobile, const VectorArrayType& immobile) const {
		const size_t s=Base().Stages();
		VectorArrayPairType ret(EXT::ReserveVectorArray(s),EXT::ReserveVectorArray(s));
		ToBaseMobileImpl(mobile,ret.Mobile());
		ToBaseImmobileImpl(immobile,ret.Immobile());
		return ret;
	}

	VectorArrayPairType ToBase(const VectorArrayPairType& all) const {
		return ToBase(all.Mobile(),all.Immobile());
	}

	void ToOriginalImpl(const VectorArrayType& all, VectorType& ret) const {
		Base().ToOriginalImpl(ToBase(all),ret);
	}
	void ToOriginalMobileImpl(const VectorArrayType& mobile, VectorType& ret) const {
		Base().ToOriginalMobileImpl(ToBaseMobile(mobile),ret);
	}
	void ToOriginalImmobileImpl(const VectorArrayType& immobile, VectorType& ret) const {
		Base().ToOriginalImmobileImpl(ToBaseImmobile(immobile),ret);
	}
	
	virtual MatrixType MobileBaseTransformation() const =0;
	virtual MatrixType ImmobileBaseTransformation() const =0;
	virtual MatrixType BaseTransformation() const =0;
	
	MatrixType MobileTransformation() const { return Base().MobileTransformation()*MobileBaseTransformation(); }
	MatrixType ImmobileTransformation() const { return Base().ImmobileTransformation()*ImmobileBaseTransformation(); }
	MatrixType Transformation() const { return Base().Transformation()*BaseTransformation(); }


	const std::vector<Species>& Participants() const { return Base().Participants(); }

	virtual ~IHierarchicalStoichiometry() {}
};
}

#endif
