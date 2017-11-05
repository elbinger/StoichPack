/* File: StoichPackBasicKineticContainer.h
 * Purpose: Provide a reaction container for the original representation of the problem.
 * Author: Tobias Elbinger <elbinger@math.fau.de> */

#ifndef __H_STOICHPACK_BASIC_KINETIC_CONTAINER__
#define __H_STOICHPACK_BASIC_KINETIC_CONTAINER__

#include "StoichPackIKineticContainer.h"
#include "StoichPackBiochemicalSystem.h"

namespace StoichPack{
 template<typename EXT, typename KineticType = IKineticReaction>
 class BasicKineticContainer : public IKineticContainer<EXT>{
 typedef InitializedBiochemicalSystem<KineticType> BiochemistryType;
 typedef typename BiochemistryType::KineticReactionType KineticReactionType;
 typedef typename EXT::VectorType VectorType;
 typedef typename EXT::MatrixType MatrixType;

 private:
	const BiochemistryType problem;
	BasicKineticContainer();

	static StageInfoArray<EXT> Init(const BiochemicalSystem<KineticType>& system){
		BiochemistryType problem(system);
		//assemble stoichiometric matrix
		const std::vector<KineticReactionType>& r=problem.KineticReactions();
		const size_t n = problem.Count(species_type::mobile)+problem.Count(species_type::immobile);
		MatrixType stoich = EXT::CreateZeroMatrix(n,r.size());
		for(size_t i=0;i<r.size();++i){
			r[i].Add(EXT::BeginColWise(stoich,i),r[i].Coefficients());
		}
		//assemble toBase / fromBase (identity)
		MatrixType identity = Diag<EXT>(std::vector<sp_scalar>(n,1.));

		const size_t nmob = problem.Count(species_type::mobile);
		const std::vector<BasicStageInfo<EXT> > tmp(1, BasicStageInfo<EXT>(stoich,identity,identity,nmob,nmob,true));
		return StageInfoArray<EXT>(tmp);
	}

	template<typename OP, typename INDEX>
	typename OP::ReturnType ReactionLoop(OP& op, const INDEX& index) const {
		const size_t s=index.size();
		const std::vector<KineticReactionType>& r=this->problem.KineticReactions();
		for(size_t i=0;i<s;++i,++op){
			assert(index[i]<r.size());
			op.apply(r[index[i]]);
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
	 void apply(const KineticReactionType& r){ *it = r.Rate(data); }
	 Impl1RateOp& operator++() { ++it; return *this; }
	};
	class Impl2RateOp{
	 public:
	 typedef typename EXT::VectorType& ReturnType;
	 const typename EXT::ConstVectorIteratorType data1;
	 const typename EXT::ConstVectorIteratorType data2;
	 ReturnType result;
	 typename EXT::VectorIteratorType it;
	 Impl2RateOp(const VectorType& mobile, const VectorType& immobile, ReturnType res)
	             : data1(EXT::ConstBegin(mobile)), data2(EXT::ConstBegin(immobile)), result(res),
	               it (EXT::Begin(res)) {}
	 void apply(const KineticReactionType& r){ *it = r.Rate(data1,data2); }
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
	 void apply(const KineticReactionType& r){
		r.Add(EXT::BeginRowWise(result,row),r.DiffRate(data));
	 }
	 Impl1DiffRateOp& operator++() { ++row; return *this; }
	};
	class Impl2DiffRateOp{
	 public:
	 typedef MatrixPair<EXT> ReturnType;
	 const typename EXT::ConstVectorIteratorType data1;
	 const typename EXT::ConstVectorIteratorType data2;
	 ReturnType result;
	 size_t row;
	 Impl2DiffRateOp(const VectorType& mobile, const VectorType& immobile, size_t rows, size_t cols1, size_t cols2)
	                 : data1(EXT::ConstBegin(mobile)), data2(EXT::ConstBegin(immobile)),
	                   result(EXT::CreateZeroMatrix(rows,cols1),EXT::CreateZeroMatrix(rows,cols2)), row(0) {
	 }
	 void apply(const KineticReactionType& r){
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
	 void apply(const KineticReactionType& r){
		r.Add(EXT::BeginRowWise(result,row),r.Dependencies());
	 }
	 StructureOp& operator++() { ++row; return *this; }
	};

 public:
	BasicKineticContainer(const BiochemicalSystem<KineticType>& system) : IKineticContainer<EXT>(Init(system)), problem(system) {}

	const BiochemistryType& Problem() const { return problem; }

	//interface
	VectorType ReactionRatesImpl(const VectorType& all, size_t stage) const {
		VectorType ret =EXT::CreateVector(problem.KineticReactions().size());
		Impl1RateOp tmp(all,ret);
		ReactionLoop(tmp,FullIndex(problem.KineticReactions().size()));
		return ret;
		/*assert(this->CheckSize(all,stage));
		const std::vector<SHARED_RPTR>& r=problem.KineticReactions();
		const size_t s= r.size();
		VectorType ret = EXT::CreateVector(s);
		typename EXT::VectorIteratorType it = EXT::Begin(ret);
		const typename EXT::ConstVectorIteratorType data = EXT::ConstBegin(all);
		for(size_t i=0;i<s;++i, ++it){
			*it = r[i].Rate(data);
		}
		return ret;*/
	}

	VectorType ReactionRatesImpl(const VectorType& mobile, const VectorType& immobile, size_t stage) const {
		VectorType ret = EXT::CreateVector(problem.KineticReactions().size());
		Impl2RateOp tmp(mobile,immobile,ret);
		ReactionLoop(tmp,FullIndex(problem.KineticReactions().size()));
		return ret;
		/*assert(this->CheckSize(all,stage));
		const std::vector<SHARED_RPTR>& r=problem.KineticReactions();
		const size_t s= r.size();
		VectorType ret = EXT::CreateVector(s);
		typename EXT::VectorIteratorType it = EXT::Begin(ret);
		const typename EXT::ConstVectorIteratorType data1 = EXT::ConstBegin(all.Mobile());
		const typename EXT::ConstVectorIteratorType data2 = EXT::ConstBegin(all.Immobile());
		for(size_t i=0;i<s;++i, ++it){
			*it = r[i].Rate(data1,data2);
		}
		return ret;*/
	}


	VectorType SubReactionRatesImpl(const VectorType& all, const std::vector<size_t>& I, size_t stage) const {
		VectorType ret = EXT::CreateVector(I.size());
		Impl1RateOp tmp(all,ret);
		ReactionLoop(tmp,I);
		return ret;
		/*assert(this->CheckSize(all,stage));
		const std::vector<SHARED_RPTR>& r=problem.KineticReactions();
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

	VectorType SubReactionRatesImpl(const VectorType& mobile, const VectorType& immobile, const std::vector<size_t>& I, size_t stage) const {
		VectorType ret = EXT::CreateVector(I.size());
		Impl2RateOp tmp(mobile,immobile,ret);
		ReactionLoop(tmp,I);
		return ret;
		/*assert(this->CheckSize(all,stage));
		const std::vector<SHARED_RPTR>& r=problem.KineticReactions();
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
	
	MatrixType DiffReactionRatesImpl(const VectorType& all, size_t stage) const {
		const size_t s = problem.KineticReactions().size();
		Impl1DiffRateOp tmp(all,s,this->GlobalSpecies());
		return ReactionLoop(tmp,FullIndex(s));
		/*assert(this->CheckSize(all,stage));
		const std::vector<SHARED_RPTR>& r=problem.KineticReactions();
		const size_t s= r.size();
		MatrixType ret = EXT::CreateZeroMatrix(s,this->GlobalSpecies());
		const typename EXT::ConstVectorIteratorType data = EXT::ConstBegin(all);
		for(size_t i=0;i<s;++i){
			r[i].Add(EXT::BeginRowWise(ret,i),r[i].DiffRate(data));
		}
		return ret;*/
	}

	MatrixPair<EXT> DiffReactionRatesImpl(const VectorType& mobile, const VectorType& immobile, size_t stage) const {
		const size_t s = problem.KineticReactions().size();
		Impl2DiffRateOp tmp(mobile,immobile,s,this->MobileSpecies(),this->ImmobileSpecies());
		return ReactionLoop(tmp,FullIndex(s));
		/*assert(this->CheckSize(all,stage));
		const std::vector<SHARED_RPTR>& r=problem.KineticReactions();
		const size_t s= r.size();
		MatrixPair<EXT> ret(EXT::CreateZeroMatrix(s,this->MobileSpecies()),EXT::CreateZeroMatrix(s,this->ImmobileSpecies()));
		const typename EXT::ConstVectorIteratorType data1 = EXT::ConstBegin(all.Mobile());
		const typename EXT::ConstVectorIteratorType data2 = EXT::ConstBegin(all.Immobile());

		for(size_t i=0;i<s;++i){
			r[i].Add(EXT::BeginRowWise(ret.Mobile(),i),EXT::BeginRowWise(ret.Immobile(),i),r[i].DiffRate(data1,data2));
		}
		return ret;*/
	}

	MatrixType DiffSubReactionRatesImpl(const VectorType& all, const std::vector<size_t>& I, size_t stage) const {
		Impl1DiffRateOp tmp(all,I.size(),this->GlobalSpecies());
		return ReactionLoop(tmp,I);
	}

	MatrixPair<EXT> DiffSubReactionRatesImpl(const VectorType& mobile, const VectorType& immobile,
	                                         const std::vector<size_t>& I, size_t stage) const {
		Impl2DiffRateOp tmp(mobile,immobile,I.size(),this->MobileSpecies(),this->ImmobileSpecies());
		return ReactionLoop(tmp,I);
	}

	VectorType ConstSpeciesRatesImpl(const VectorType& all, size_t stage) const {
		return EXT::CreateZeroVector(this->GlobalSpecies(0));
	}

	VectorPair<EXT> ConstSpeciesRatesImpl(const VectorType& mobile, const VectorType& immobile, size_t stage) const {
		return VectorPair<EXT>(EXT::CreateZeroVector(this->MobileSpecies()), EXT::CreateZeroVector(this->ImmobileSpecies()));
	}


	const std::vector<InitializedSpecies>& Participants() const { return problem.Participants(); }

	MatrixType RateStructure(size_t stage, const std::vector<size_t>& I) const {
		assert(stage==0);
		StructureOp tmp(problem.KineticReactions().size(),Participants().size());
		return ReactionLoop(tmp,I);
	}
	MatrixType RateStructure(size_t stage) const {
		assert(stage==0);
		StructureOp tmp(problem.KineticReactions().size(),Participants().size());
		return ReactionLoop(tmp,FullIndex(problem.KineticReactions().size()));
	}		

	MatrixType SpeciesStructure(size_t stage) const {
		StructureOp tmp(problem.KineticReactions().size(),Participants().size());
		return BooleanMatrix<EXT>(this->StoichiometricMatrices()[0])*RateStructure(0);
	}

	MatrixType MobileTransformation() const { return Diag<EXT>(std::vector<sp_scalar>(this->MobileSpecies(),1)); }
	MatrixType ImmobileTransformation() const { return Diag<EXT>(std::vector<sp_scalar>(this->ImmobileSpecies(),1)); }
	MatrixType Transformation() const { return Diag<EXT>(std::vector<sp_scalar>(this->GlobalSpecies(),1)); }

	IKineticContainer<EXT>* copy() const { return new BasicKineticContainer(*this); }
 };

 template<typename EXT>
 using BasicKineticStoichiometry = BasicKineticContainer<EXT>;
} // namespace StoichPack

#endif
