#ifndef __H_STOICHPACK_HIERARCHICAL_LINEAR_STOICHIOMETRY__
#define __H_STOICHPACK_HIERARCHICAL_LINEAR_STOICHIOMETRY__

#include "StoichPackIHierarchicalStoichiometry.h"

namespace StoichPack {
template<typename EXT, typename BT>
class IHierarchicalLinearStoichiometry : public IHierarchicalStoichiometry<EXT,BT>{
 	typedef typename EXT::VectorType VectorType;
	typedef typename EXT::VectorPairType VectorPairType;
	typedef typename EXT::VectorArrayType VectorArrayType;
	typedef typename EXT::VectorArrayPairType VectorArrayPairType;
	typedef typename EXT::MatrixType MatrixType;
private:
	std::vector<MatrixType> toBase_global, toBase_mobile, toBase_immobile;
	std::vector<MatrixType> fromBase_global, fromBase_mobile, fromBase_immobile;

	void FromBaseImplImpl(const VectorArrayType& base, VectorArrayType& ret, const std::vector<MatrixType>& R) const{
		const size_t s=R.size();
		assert(s==this->Stages());
		typename EXT::ConstVectorArrayIteratorType it = EXT::ConstBegin(base);
		for(size_t i=0;i<s;++i) {
			const size_t ss=this->SubStage(i);
			EXT::PushBack(ret,R[i]**(it+ss));
		}
	}

	void ToBaseImplImpl(VectorArrayType& base, const VectorArrayType& in, const std::vector<MatrixType>& R) const{
		const size_t s=Base().Stages();
		assert(R.size()==this->Stages());
		typename EXT::ConstVectorArrayIteratorType it = EXT::ConstBegin(in);
		for(size_t i=0;i<s;++i) {
			const size_t start = SubStageFirst(i);
			const size_t stop = SubStageFirst(i+1);
			EXT::PushBack(base,R[start]**(it+start));
			VectorType& value = *(EXT::Begin(base)+i);
			for(size_t j=start+1;j<stop;++j) value+=R[j]**(it+j);
		}
	}

protected:
	using IHierarchicalStoichiometry<EXT,BT>::Finish;
	using IHierarchicalStoichiometry<EXT,BT>::SubStage;
	using IHierarchicalStoichiometry<EXT,BT>::SubStageFirst;
	using IHierarchicalStoichiometry<EXT,BT>::Base;

	IHierarchicalLinearStoichiometry(const BT* oth) : IHierarchicalStoichiometry<EXT,BT>(oth) {}
	void AddStage(const MatrixType& stoich, size_t substage, bool correction, const MatrixType& toBasemobile, const MatrixType& toBaseimmobile, const MatrixType& fromBasemobile, const MatrixType& fromBaseimmobile){
		const size_t basemobile=EXT::rows(toBasemobile);
		const size_t baseimmobile=EXT::rows(toBaseimmobile);
		const size_t mobile=EXT::rows(fromBasemobile);
		const size_t immobile=EXT::rows(fromBaseimmobile);

		assert(mobile==EXT::cols(toBasemobile) && immobile==EXT::cols(toBaseimmobile));
		assert(EXT::cols(fromBasemobile)==basemobile && EXT::cols(fromBaseimmobile)==baseimmobile);

		toBase_mobile.push_back(toBasemobile);
		toBase_immobile.push_back(toBaseimmobile);
		fromBase_mobile.push_back(fromBasemobile);
		fromBase_immobile.push_back(fromBaseimmobile);

		MatrixType a = CombineCols<EXT>(toBasemobile,EXT::CreateZeroMatrix(basemobile,immobile));
		MatrixType b = CombineCols<EXT>(EXT::CreateZeroMatrix(baseimmobile,mobile),toBaseimmobile);
		toBase_global.push_back(CombineRows<EXT>(a,b));

		MatrixType c = CombineCols<EXT>(fromBasemobile,EXT::CreateZeroMatrix(mobile,baseimmobile));
		MatrixType d = CombineCols<EXT>(EXT::CreateZeroMatrix(immobile,basemobile),fromBaseimmobile);
		fromBase_global.push_back(CombineRows<EXT>(c,d));
		
		IHierarchicalStoichiometry<EXT,BT>::AddStage(stoich,mobile,substage,correction);
	}

	const MatrixType& toBase(size_t stage) const { return toBase_global[stage]; }
	const MatrixType& toBaseMobile(size_t stage) const { return toBase_mobile[stage]; }
	const MatrixType& toBaseImmobile(size_t stage) const { return toBase_immobile[stage]; }
	const MatrixType& fromBase(size_t stage) const { return fromBase_global[stage]; }
	const MatrixType& fromBaseMobile(size_t stage) const { return fromBase_mobile[stage]; }
	const MatrixType& fromBaseImmobile(size_t stage) const { return fromBase_immobile[stage]; }

public:
	void FromBaseImpl(const VectorArrayType& all, VectorArrayType& ret) const {
		FromBaseImplImpl(all,ret,fromBase_global);
	}
	void FromBaseMobileImpl(const VectorArrayType& mobile, VectorArrayType& ret) const {
		FromBaseImplImpl(mobile,ret,fromBase_mobile);
	}
	void FromBaseImmobileImpl(const VectorArrayType& immobile, VectorArrayType& ret) const {
		FromBaseImplImpl(immobile,ret,fromBase_immobile);
	}

	void ToBaseImpl(const VectorArrayType& all, VectorArrayType& ret) const {
		ToBaseImplImpl(ret,all,toBase_global);
	}
	void ToBaseMobileImpl(const VectorArrayType& mobile, VectorArrayType& ret) const {
		ToBaseImplImpl(ret,mobile,toBase_mobile);
	}
	void ToBaseImmobileImpl(const VectorArrayType& immobile, VectorArrayType& ret) const {
		ToBaseImplImpl(ret,immobile,toBase_immobile);
	}
	
	MatrixType MobileBaseTransformation() const { return EXT::CreateMatrix(1,1); }
	MatrixType ImmobileBaseTransformation() const { return EXT::CreateMatrix(1,1); }
	MatrixType BaseTransformation() const { return EXT::CreateMatrix(1,1); }

};
	
}

#endif
