#ifndef __H_STOICHPACK_BASIC_KINETIC_CONTAINER__
#define __H_STOICHPACK_BASIC_KINETIC_CONTAINER__

#include "StoichPackIKineticContainer.h"
#include "StoichPackKineticReactionArray.h"

namespace StoichPack{
 template<typename EXT, typename ReactionType = IKineticReaction>
 class BasicKineticContainer : public IKineticContainer<EXT>{
	typedef std::shared_ptr<ReactionType> SHARED_RPTR;
	typedef const ReactionType& RREF;
 private:
	template<typename OP, typename INDEX>
	typename OP::ReturnType ReactionLoop(OP& op, const INDEX& index) const {
		const size_t s=index.size();
		const std::vector<SHARED_RPTR>& r=this->reactions.Reactions();
		for(size_t i=0;i<s;++i,++op){
			assert(index[i]<r.size());
			op.apply(*r[index[i]]);
		}
		return op.result;
	}

	class FullIndex{
	 public:
	 const size_t s;
	 FullIndex(size_t si) : s(si) {}
	 size_t size() const { return s; }
	 size_t operator[](size_t i) const { return i; }
	};
	class Impl1RateOp{
	 public:
	 typedef typename EXT::VectorType& ReturnType;
	 const typename EXT::ConstVectorIteratorType data;
	 ReturnType result;
	 typename EXT::VectorIteratorType it;
	 Impl1RateOp(const typename EXT::VectorType& all, ReturnType res)
	             : data(EXT::ConstBegin(all)), result(res), it(EXT::Begin(res)) {}
	 void apply(RREF r){ *it = r.Rate(data); }
	 Impl1RateOp& operator++() { ++it; return *this; }
	};
	class Impl2RateOp{
	 public:
	 typedef typename EXT::VectorType& ReturnType;
	 const typename EXT::ConstVectorIteratorType data1;
	 const typename EXT::ConstVectorIteratorType data2;
	 ReturnType result;
	 typename EXT::VectorIteratorType it;
	 Impl2RateOp(const typename EXT::VectorPairType& all, ReturnType res)
	             : data1(EXT::ConstBegin(all.Mobile())), data2(EXT::ConstBegin(all.Immobile())), result(res),
	               it (EXT::Begin(res)) {}
	 void apply(RREF r){ *it = r.Rate(data1,data2); }
	 Impl2RateOp& operator++() { ++it; return *this; }
	};
	class Impl1DiffRateOp{
	 public:
	 typedef typename EXT::MatrixType ReturnType;
	 const typename EXT::ConstVectorIteratorType data;
	 ReturnType result;
	 size_t row;
	 Impl1DiffRateOp(const typename EXT::VectorType& all, size_t rows, size_t cols)
	                 : data(EXT::ConstBegin(all)), result(EXT::CreateZeroMatrix(rows,cols)), row(0) {
	 }
	 void apply(RREF r){
		r.Add(EXT::BeginRowWise(result,row),r.DiffRate(data));
	 }
	 Impl1DiffRateOp& operator++() { ++row; return *this; }
	};
	class Impl2DiffRateOp{
	 public:
	 typedef typename EXT::MatrixPairType ReturnType;
	 const typename EXT::ConstVectorIteratorType data1;
	 const typename EXT::ConstVectorIteratorType data2;
	 ReturnType result;
	 size_t row;
	 Impl2DiffRateOp(const typename EXT::VectorPairType& all, size_t rows, size_t cols1, size_t cols2)
	                 : data1(EXT::ConstBegin(all.Mobile())), data2(EXT::ConstBegin(all.Immobile())),
	                   result(EXT::CreateZeroMatrix(rows,cols1),EXT::CreateZeroMatrix(rows,cols2)), row(0) {
	 }
	 void apply(RREF r){
		r.Add(EXT::BeginRowWise(result.Mobile(),row),EXT::BeginRowWise(result.Immobile(),row),r.DiffRate(data1,data2));
	 }
	 Impl2DiffRateOp& operator++() { ++row; return *this; }
	};
	class StructureOp{
	 public:
	 typedef typename EXT::MatrixType ReturnType;
	 ReturnType result;
	 size_t row;
	 StructureOp(size_t rows, size_t cols)
	                 : result(EXT::CreateZeroMatrix(rows,cols)), row(0) {
	 }
	 void apply(RREF r){
		r.Add(EXT::BeginRowWise(result,row),r.Dependencies());
	 }
	 StructureOp& operator++() { ++row; return *this; }
	};

 public:
 	typedef typename EXT::VectorType VectorType;
	typedef typename EXT::VectorPairType VectorPairType;
	typedef typename EXT::VectorArrayType VectorArrayType;
	typedef typename EXT::VectorArrayPairType VectorArrayPairType;
	typedef typename EXT::MatrixType MatrixType;
	typedef typename EXT::MatrixPairType MatrixPairType;
	typedef typename EXT::MatrixQuadType MatrixQuadType;

	KineticReactionArray<ReactionType> reactions;

	void Initialize(){
		if(!reactions.Initialized()) reactions.Initialize();
		if(this->Stages()==1) return;

		if(this->Stages()!=0) throw StoichPackException("Something is really messed up!");

		//assemble stoichiometric matrix
		const std::vector<SHARED_RPTR>& r=reactions.Reactions();
		const size_t n = reactions.MobileSpecies()+reactions.ImmobileSpecies();
		typename EXT::MatrixType stoich = EXT::CreateZeroMatrix(n,r.size());
		for(size_t i=0;i<r.size();++i){
			r[i]->Add(EXT::BeginColWise(stoich,i),r[i]->Coefficients());
		}
		IKineticContainer<EXT>::AddStage(stoich,reactions.MobileSpecies(),true);
	}

	//interface
	void FromOriginalImpl(const VectorType& all, VectorArrayType& ret) const {
		assert(this->CheckSize(all));
		EXT::PushBack(ret,all);
	}
	void FromOriginalMobileImpl(const VectorType& mobile, VectorArrayType& ret) const{
		assert(this->CheckSizeMobile(mobile));
		EXT::PushBack(ret,mobile);
	}
	void FromOriginalImmobileImpl(const VectorType& immobile, VectorArrayType& ret) const {
		assert(this->CheckSizeImmobile(immobile));
		EXT::PushBack(ret,immobile);
	}

	void ToOriginalImpl(const VectorArrayType& all, VectorType& ret) const {
		ret=all[0];
	}
	void ToOriginalMobileImpl(const VectorArrayType& mobile, VectorType& ret) const {
		ret=mobile[0];
	}
	void ToOriginalImmobileImpl(const VectorArrayType& immobile, VectorType& ret) const {
		ret=immobile[0];
	}
	VectorType ReactionRatesImpl1(const VectorType& all, size_t stage) const {
		VectorType ret =EXT::CreateVector(reactions.Reactions().size());
		Impl1RateOp tmp(all,ret);
		ReactionLoop(tmp,FullIndex(reactions.Reactions().size()));
		return ret;
		/*assert(this->CheckSize(all,stage));
		const std::vector<SHARED_RPTR>& r=reactions.Reactions();
		const size_t s= r.size();
		VectorType ret = EXT::CreateVector(s);
		typename EXT::VectorIteratorType it = EXT::Begin(ret);
		const typename EXT::ConstVectorIteratorType data = EXT::ConstBegin(all);
		for(size_t i=0;i<s;++i, ++it){
			*it = r[i]->Rate(data);
		}
		return ret;*/
	}

	VectorType ReactionRatesImpl2(const VectorPairType& all, size_t stage) const {
		VectorType ret = EXT::CreateVector(reactions.Reactions().size());
		Impl2RateOp tmp(all,ret);
		ReactionLoop(tmp,FullIndex(reactions.Reactions().size()));
		return ret;
		/*assert(this->CheckSize(all,stage));
		const std::vector<SHARED_RPTR>& r=reactions.Reactions();
		const size_t s= r.size();
		VectorType ret = EXT::CreateVector(s);
		typename EXT::VectorIteratorType it = EXT::Begin(ret);
		const typename EXT::ConstVectorIteratorType data1 = EXT::ConstBegin(all.Mobile());
		const typename EXT::ConstVectorIteratorType data2 = EXT::ConstBegin(all.Immobile());
		for(size_t i=0;i<s;++i, ++it){
			*it = r[i]->Rate(data1,data2);
		}
		return ret;*/
	}


	VectorType SubReactionRatesImpl1(const VectorType& all, const std::vector<size_t>& I, size_t stage) const {
		VectorType ret = EXT::CreateVector(I.size());
		Impl1RateOp tmp(all,ret);
		ReactionLoop(tmp,I);
		return ret;
		/*assert(this->CheckSize(all,stage));
		const std::vector<SHARED_RPTR>& r=reactions.Reactions();
		const size_t s= I.size();
		VectorType ret = EXT::CreateVector(s);
		typename EXT::VectorIteratorType it = EXT::Begin(ret);
		const typename EXT::ConstVectorIteratorType data = EXT::ConstBegin(all);
		for(size_t i=0;i<s;++i, ++it){
			assert(I[i]<r.size());
			*it = r[I[i]]->Rate(data);
		}
		return ret;*/
	}

	VectorType SubReactionRatesImpl2(const VectorPairType& all, const std::vector<size_t>& I, size_t stage) const {
		VectorType ret = EXT::CreateVector(I.size());
		Impl2RateOp tmp(all,ret);
		ReactionLoop(tmp,I);
		return ret;
		/*assert(this->CheckSize(all,stage));
		const std::vector<SHARED_RPTR>& r=reactions.Reactions();
		const size_t s= I.size();
		VectorType ret = EXT::CreateVector(s);
		typename EXT::VectorIteratorType it = EXT::Begin(ret);
		const typename EXT::ConstVectorIteratorType data1 = EXT::ConstBegin(all.Mobile());
		const typename EXT::ConstVectorIteratorType data2 = EXT::ConstBegin(all.Immobile());
		for(size_t i=0;i<s;++i, ++it){
			assert(I[i]<r.size());
			*it = r[I[i]]->Rate(data1,data2);
		}
		return ret;*/
	}
	
	MatrixType DiffReactionRatesImpl1(const VectorType& all, size_t stage) const {
		const size_t s = reactions.Reactions().size();
		Impl1DiffRateOp tmp(all,s,this->AllSpecies());
		return ReactionLoop(tmp,FullIndex(s));
		/*assert(this->CheckSize(all,stage));
		const std::vector<SHARED_RPTR>& r=reactions.Reactions();
		const size_t s= r.size();
		MatrixType ret = EXT::CreateZeroMatrix(s,this->AllSpecies());
		const typename EXT::ConstVectorIteratorType data = EXT::ConstBegin(all);
		for(size_t i=0;i<s;++i){
			r[i]->Add(EXT::BeginRowWise(ret,i),r[i]->DiffRate(data));
		}
		return ret;*/
	}

	MatrixPairType DiffReactionRatesImpl2(const VectorPairType& all, size_t stage) const {
		const size_t s = reactions.Reactions().size();
		Impl2DiffRateOp tmp(all,s,this->MobileSpecies(),this->ImmobileSpecies());
		return ReactionLoop(tmp,FullIndex(s));
		/*assert(this->CheckSize(all,stage));
		const std::vector<SHARED_RPTR>& r=reactions.Reactions();
		const size_t s= r.size();
		MatrixPairType ret(EXT::CreateZeroMatrix(s,this->MobileSpecies()),EXT::CreateZeroMatrix(s,this->ImmobileSpecies()));
		const typename EXT::ConstVectorIteratorType data1 = EXT::ConstBegin(all.Mobile());
		const typename EXT::ConstVectorIteratorType data2 = EXT::ConstBegin(all.Immobile());

		for(size_t i=0;i<s;++i){
			r[i]->Add(EXT::BeginRowWise(ret.Mobile(),i),EXT::BeginRowWise(ret.Immobile(),i),r[i]->DiffRate(data1,data2));
		}
		return ret;*/
	}

	MatrixType DiffSubReactionRatesImpl1(const VectorType& all, const std::vector<size_t>& I, size_t stage) const {
		Impl1DiffRateOp tmp(all,I.size(),this->AllSpecies());
		return ReactionLoop(tmp,I);
	}

	MatrixPairType DiffSubReactionRatesImpl2(const VectorPairType& all, const std::vector<size_t>& I, size_t stage) const {
		Impl2DiffRateOp tmp(all,I.size(),this->MobileSpecies(),this->ImmobileSpecies());
		return ReactionLoop(tmp,I);
	}

	VectorType ConstSpeciesRatesImpl1(const VectorType& all, size_t stage) const {
		return EXT::CreateZeroVector(this->AllSpecies(0));
	}

	VectorPairType ConstSpeciesRatesImpl2(const VectorPairType& all, size_t stage) const {
		return VectorPairType(EXT::CreateZeroVector(this->MobileSpecies()), EXT::CreateZeroVector(this->ImmobileSpecies()));
	}

	VectorType ConstMobileSpeciesRatesImpl(const VectorPairType& all, size_t stage) const {
		return EXT::CreateZeroVector(this->MobileSpecies());
	}

	VectorType ConstImmobileSpeciesRatesImpl(const VectorPairType& all, size_t stage) const {
		return EXT::CreateZeroVector(this->ImmobileSpecies());
	}

	bool ApplyCorrectionImpl(VectorArrayType& all, const VectorArrayType& allowed) const {
		assert(EXT::size(allowed)==1 && EXT::size(*EXT::ConstBegin(allowed))==this->AllSpecies());
		return Correct(all,allowed);
	}

	bool ApplyMobileCorrectionImpl(VectorArrayType& mobile, const VectorArrayType& allowed) const {
		assert(EXT::size(allowed)==1 && EXT::size(*EXT::ConstBegin(allowed))==this->MobileSpecies());
		return Correct(mobile,allowed);
	}

	bool ApplyImmobileCorrectionImpl(VectorArrayType& immobile, const VectorArrayType& allowed) const {
		assert(EXT::size(allowed)==1 && EXT::size(*EXT::ConstBegin(allowed))==this->ImmobileSpecies());
		return Correct(immobile,allowed);
	}

	bool Correct(VectorArrayType& X, const VectorArrayType& allowed) const {
		assert(EXT::size(X)==1);
		VectorType& x=*EXT::Begin(X);
		assert(EXT::size(x)==EXT::size(*EXT::ConstBegin(allowed)));
		typename EXT::ConstVectorIteratorType a=EXT::ConstBegin(*EXT::ConstBegin(allowed));
		size_t i=0;
		for(typename EXT::VectorIteratorType it = EXT::Begin(x); it!=EXT::End(x); ++it, ++i){
			if(*it<0){
				if(*(a+i)==0) return false;
				*it=0;
			}
		}
		return true;
	}

	const std::vector<Species>& Participants() const { return reactions.Participants(); }

	MatrixType RateStructure(size_t stage, const std::vector<size_t>& I) const {
		assert(stage==0);
		StructureOp tmp(reactions.Reactions().size(),Participants().size());
		return ReactionLoop(tmp,I);
	}
	MatrixType RateStructure(size_t stage) const {
		assert(stage==0);
		StructureOp tmp(reactions.Reactions().size(),Participants().size());
		return ReactionLoop(tmp,FullIndex(reactions.Reactions().size()));
	}		

	MatrixType SpeciesStructure(size_t stage) const {
		StructureOp tmp(reactions.Reactions().size(),Participants().size());
		return BooleanMatrix<EXT>(this->StoichiometricMatrices()[0])*RateStructure(0);
	}

	MatrixType MobileTransformation() const { return Diag<EXT>(std::vector<sp_scalar>(this->MobileSpecies(),1)); }
	MatrixType ImmobileTransformation() const { return Diag<EXT>(std::vector<sp_scalar>(this->ImmobileSpecies(),1)); }
	MatrixType Transformation() const { return Diag<EXT>(std::vector<sp_scalar>(this->AllSpecies(),1)); }

	IKineticContainer<EXT>* copy() const { return new BasicKineticContainer(*this); }
 };

 template<typename EXT>
 using BasicKineticStoichiometry = BasicKineticContainer<EXT>;
} // namespace StoichPack

#endif
