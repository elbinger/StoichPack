#ifndef __H_STOICHPACK_ONESIDED__
#define __H_STOICHPACK_ONESIDED__

#include "StoichPackIHierarchicalLinearStoichiometry.h"
#include <algorithm>

namespace StoichPack{

template<typename EEXT>
class Rearangement{
private:
	std::vector<size_t> I;
	typename EEXT::MatrixType fromBase, toBase;

public:
	Rearangement(const std::vector<size_t>& i, size_t n) : I(i), fromBase(EEXT::CreateZeroMatrix(i.size(),n)),
	                                                       toBase(EEXT::CreateMatrix(n,i.size())) {
		std::sort(I.begin(),I.end());
		for(size_t i=0;i<I.size();++i){
			assert(I[i]<EEXT::cols(fromBase));
			*(EEXT::BeginRowWise(fromBase,i)+I[i])=1;
		}
		toBase=EEXT::Transposed(fromBase);
	}

	const std::vector<size_t>& Index() const { return I; }
	const typename EEXT::MatrixType& FromBase() const { return fromBase; }
	const typename EEXT::MatrixType& ToBase() const { return toBase; }
};
	
template<typename EXT, class BT = IKineticStoichiometry<EXT> >
class OneSidedStoichiometry : public IHierarchicalLinearStoichiometry<EXT,BT > {
 	typedef typename EXT::VectorType VectorType;
	typedef typename EXT::VectorPairType VectorPairType;
	typedef typename EXT::VectorArrayType VectorArrayType;
	typedef typename EXT::VectorArrayPairType VectorArrayPairType;
	typedef typename EXT::MatrixType MatrixType;
	typedef typename EXT::MatrixPairType MatrixPairType;
	typedef typename EXT::MatrixQuadType MatrixQuadType;

	typedef std::vector<size_t> IndexType;

	static IndexType GetSub(const IndexType& x, const IndexType& I){
		const size_t s=I.size();
		IndexType ret(s);
		for(size_t i=0;i<s;++i) ret[i]=x[I[i]];
		return ret;
	}

private:
	using IHierarchicalLinearStoichiometry<EXT,BT >::Base;

	std::vector<IndexType> subreactions, const_subreactions;
	std::vector<MatrixType> stoich_const, stoich_const_mobile, stoich_const_immobile;

	bool Update(typename EXT::RowWiseIteratorType x, typename EXT::ConstRowWiseIteratorType y, size_t n) const {
		bool update=false;
		for(size_t i=0;i<n;++i){
			if(*(y+i)!=0 && *(x+i)==0){
				*(x+i)=1;
				update=true;
			}
		}
		return update;
	}

	MatrixType GetDependencies(size_t substage) const {
		MatrixType structure = Base().SpeciesStructure(substage);
		const size_t n=EXT::rows(structure);
		for(size_t i=0;i<n;++i) *(EXT::BeginColWise(structure,i)+i)=1;
		bool update=true;
		while(update){
			update=false;
			for(size_t i=0;i<n;++i){
				typename EXT::RowWiseIteratorType it = EXT::BeginRowWise(structure,i);
				for(size_t j=0;j<n;++j){
					if(*(it+j)!=0 && i!=j && Update(it,EXT::ConstBeginRowWise(structure,j),n)) update=true;
				}
			}
		}
		return structure;
	}

	IndexType StableSortValues(const MatrixType& M) const {
		class mypair{
			public:
			size_t index, values;
			mypair(size_t i, size_t v) : index(i), values(v) {}
			bool operator<(const mypair& oth) const { return values<oth.values; }
		};
		std::vector<mypair> values;
		const size_t n=EXT::rows(M);
		for(size_t i=0;i<n;++i){
			const size_t tmp = std::count(EXT::ConstBeginRowWise(M,i),EXT::ConstEndRowWise(M,i),sp_scalar(0));
			values.push_back(mypair(i,n-tmp));
		}
		std::sort(values.begin(),values.end());
		std::vector<size_t> ret;
		for(size_t i=0;i<values.size();++i) ret.push_back(values[i].index);
		return ret;
	}

	void ProcessSubstage(size_t substage){
		IndexType known_reactions;
		MatrixType dependencies = GetDependencies(substage);
		IndexType order = StableSortValues(dependencies);

		while(order.size()!=0){
			IndexType I;
			const size_t row = order.front();
			typename EXT::ConstRowWiseIteratorType it=EXT::ConstBeginRowWise(dependencies,row);
			for(size_t i=0;i<EXT::rows(dependencies);++i){
				if(*(it+i)!=0){
					IndexType::iterator pos=std::find(order.begin(),order.end(),i);
					if(pos!=order.end()){
						I.push_back(i);
						order.erase(pos);
					}
				}
			}
			Rearangement<EXT> r(I,EXT::rows(dependencies));
			AddStage(substage,r,known_reactions);
		}
	}

	void AddStage(size_t substage, const Rearangement<EXT>& R, IndexType& known_reactions){
		IndexType mob, immob;
		const size_t basemob = Base().MobileSpecies(substage);

		for(IndexType::const_iterator it=R.Index().begin();it!=R.Index().end();++it){
			if(*it<basemob) mob.push_back(*it);
			else immob.push_back((*it)-basemob);
		}

		Rearangement<EXT> subspecies_mobile(mob,basemob);
		Rearangement<EXT> subspecies_immobile(immob,Base().ImmobileSpecies(substage));

		MatrixType tmp = R.FromBase()*Base().StoichiometricMatrices()[substage];
		IndexType varreactions, constreactions;
		const size_t n_r=EXT::cols(tmp);
		const size_t n_s=EXT::rows(tmp);

		for(size_t j=0;j<n_r;++j){
			if(size_t(std::count(EXT::BeginColWise(tmp,j),EXT::EndColWise(tmp,j),sp_scalar(0)))<n_s){
				if(std::find(known_reactions.begin(),known_reactions.end(),j)==known_reactions.end()){
					known_reactions.push_back(j);
					varreactions.push_back(j);
				} else constreactions.push_back(j);
			}
		}

		subreactions.push_back(varreactions);
		const_subreactions.push_back(constreactions);

		Rearangement<EXT> rvar(varreactions,n_r);
		IHierarchicalLinearStoichiometry<EXT,BT >::AddStage(R.FromBase()*Base().StoichiometricMatrices()[substage]*rvar.ToBase(),substage,true,subspecies_mobile.ToBase(),subspecies_immobile.ToBase(),subspecies_mobile.FromBase(),subspecies_immobile.FromBase());

		Rearangement<EXT> rconst(constreactions,n_r);
		const MatrixType& M=rconst.ToBase();
		stoich_const.push_back(R.FromBase()*Base().StoichiometricMatrices()[substage]*M);
		stoich_const_mobile.push_back(subspecies_mobile.FromBase()*Base().MobileStoichiometricMatrices()[substage]*M);
		stoich_const_immobile.push_back(subspecies_immobile.FromBase()*Base().ImmobileStoichiometricMatrices()[substage]*M);
	}

public:
	OneSidedStoichiometry(const BT* bt) : IHierarchicalLinearStoichiometry<EXT,BT >(bt) {
		for(size_t i=0;i<bt->Stages();++i) ProcessSubstage(i);
		this->Finish();
	}

	VectorType ReactionRatesImpl1(const VectorType& all, size_t stage) const { 
		return Base().SubReactionRatesImpl1(all,subreactions[stage],this->SubStage(stage));
	}
	VectorType ReactionRatesImpl2(const VectorPairType& all, size_t stage) const {
		return Base().SubReactionRatesImpl2(all,subreactions[stage],this->SubStage(stage));
	}
	VectorType SubReactionRatesImpl1(const VectorType& all, const std::vector<size_t>& I, size_t stage) const { 
		return Base().SubReactionRatesImpl1(all,GetSub(subreactions[stage],I),this->SubStage(stage));
	}
	VectorType SubReactionRatesImpl2(const VectorPairType& all, const std::vector<size_t>& I, size_t stage) const {
		return Base().SubReactionRatesImpl2(all,GetSub(subreactions[stage],I),this->SubStage(stage));
	}

	MatrixType DiffReactionRatesImpl1(const VectorType& all, size_t stage) const {
		return Base().DiffSubReactionRatesImpl1(all,subreactions[stage],this->SubStage(stage))*this->toBase(stage);
	}
	MatrixPairType DiffReactionRatesImpl2(const VectorPairType& all, size_t stage) const {
		const MatrixPairType tmp = Base().DiffSubReactionRatesImpl2(all,subreactions[stage],this->SubStage(stage));
		return MatrixPairType(tmp.Mobile()*this->toBaseMobile(stage),tmp.Immobile()*this->toBaseImmobile(stage));
	}

	MatrixType DiffSubReactionRatesImpl1(const VectorType& all, const std::vector<size_t>& I, size_t stage) const {
		return Base().DiffSubReactionRatesImpl1(all,GetSub(subreactions[stage],I),this->SubStage(stage))*this->toBase(stage);
	}
	MatrixPairType DiffSubReactionRatesImpl2(const VectorPairType& all, const std::vector<size_t>& I, size_t stage) const {
		const MatrixPairType tmp = Base().DiffSubReactionRatesImpl2(all,GetSub(subreactions[stage],I),this->SubStage(stage));
		return MatrixPairType(tmp.Mobile()*this->toBaseMobile(stage),tmp.Immobile()*this->toBaseImmobile(stage));
	}

	VectorType ConstSpeciesRatesImpl1(const VectorType& all, size_t stage) const {
		return stoich_const[stage]*Base().SubReactionRatesImpl1(all,const_subreactions[stage],this->SubStage(stage))+
		       this->fromBase(stage)*Base().ConstSpeciesRatesImpl1(all,this->SubStage(stage));
	}

	VectorPairType ConstSpeciesRatesImpl2(const VectorPairType& all, size_t stage) const {
		const size_t substage = this->SubStage(stage);
		const VectorPairType tmp1 = Base().ConstSpeciesRatesImpl2(all,substage);
		const VectorType tmp2 = Base().SubReactionRatesImpl2(all,subreactions[stage],substage);
		return VectorPairType( stoich_const_mobile[stage]*tmp2 + this->fromBaseMobile(stage)*tmp1.Mobile(),
		                       stoich_const_immobile[stage]*tmp2 + this->fromBaseImmobile(stage)*tmp1.Immobile());
	}

	VectorType ConstMobileSpeciesRatesImpl(const VectorPairType& all, size_t stage) const {
		const size_t substage = this->SubStage(stage);
		const VectorType tmp1 = Base().ConstMobileSpeciesRatesImpl(all,substage);
		const VectorType tmp2 = Base().SubReactionRatesImpl2(all,const_subreactions[stage],substage);
		return stoich_const_mobile[stage]*tmp2 + this->fromBaseMobile(stage)*tmp1;
	}
	VectorType ConstImmobileSpeciesRatesImpl(const VectorPairType& all, size_t stage) const {
		const size_t substage = this->SubStage(stage);
		const VectorType tmp1 = Base().ConstImmobileSpeciesRatesImpl(all,substage);
		const VectorType tmp2 = Base().SubReactionRatesImpl2(all,const_subreactions[stage],substage);
		return stoich_const_immobile[stage]*tmp2 + this->fromBaseImmobile(stage)*tmp1;
	}

	bool ApplyMobileCorrectionImpl(VectorArrayType& mobile, const VectorArrayType& allowed) const {
		VectorArrayType tmp = this->ToBaseMobile(mobile);
		const bool ret = Base().ApplyMobileCorrectionImpl(tmp,this->ToBaseMobile(allowed));
		EXT::set(mobile,this->FromBaseMobile(tmp));
		return ret;
	}
	bool ApplyImmobileCorrectionImpl(VectorArrayType& immobile, const VectorArrayType& allowed) const {
		VectorArrayType tmp = this->ToBaseImmobile(immobile);
		const bool ret = Base().ApplyImmobileCorrectionImpl(tmp,this->ToBaseImmobile(allowed));
		EXT::set(immobile, this->FromBaseImmobile(tmp));
		return ret;
	}
		
	bool ApplyCorrectionImpl(VectorArrayType& all, const VectorArrayType& allowed) const {
		VectorArrayType tmp = this->ToBase(all);
		const bool ret = Base().ApplyCorrectionImpl(tmp,this->ToBase(allowed));
		EXT::set(all,this->FromBase(tmp));
		return ret;
	}

	MatrixType SpeciesStructure(size_t stage) const { return EXT::CreateMatrix(1,1); }

};
}

#endif
