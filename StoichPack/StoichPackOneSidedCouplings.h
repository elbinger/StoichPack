#ifndef __H_STOICHPACK_ONESIDED_KINETIC_COUPLINGS__
#define __H_STOICHPACK_ONESIDED_KINETIC_COUPLINGS__

#include "StoichPackIHierarchicalKineticContainer.h"
#include <algorithm>

namespace StoichPack{

template<typename EXT>
class OneSidedCouplingsStage{
public:
	typedef typename EXT::MatrixType MatrixType;

	OneSidedCouplingsStage(const IKineticContainer<EXT>& _Base, const std::vector<size_t>& _SubSpecies,
		const std::vector<size_t>& _ConstReactions, size_t _SubStage) :
			base(_Base), subspecies(_SubSpecies), substage(_SubStage) {
		for(size_t j : RelevantCols()){
			if(std::find(_ConstReactions.begin(),_ConstReactions.end(),j)!=_ConstReactions.end())
				constsubreactions.push_back(j);
		}
	}

	size_t SubStage() const { return substage; }
	const std::vector<size_t>& SubSpecies() const { return subspecies; }
	const std::vector<size_t>& ConstSubReactions() const { return constsubreactions; }
	
	std::vector<size_t> VarSubReactions() const {
		const MatrixType stoich = base.Global(substage).StoichiometricMatrix();
		std::vector<size_t> result;
		for(size_t j : RelevantCols()){
			if(std::find(constsubreactions.begin(),constsubreactions.end(),j)==constsubreactions.end())
				result.push_back(j);
		}
		return result;
	}
		
	MatrixType ToBase() const{
		MatrixType result = EXT::CreateZeroMatrix(base.Global(substage).GlobalSpecies(),subspecies.size());
		for(size_t j=0;j<subspecies.size();++j){
			*(EXT::BeginColWise(result,j)+subspecies[j])=1;
		}
		return result;
	}

	MatrixType FromBaseGlobal() const { return ToBase().transpose(); }
	MatrixType StoichiometricMatrix() const {
		const MatrixType stoich = base.Global(substage).StoichiometricMatrix();
		return SubRows<EXT>(SubCols<EXT>(stoich,VarSubReactions()),subspecies);
	}

	MatrixType ConstStoichiometricMatrix() const{
		const MatrixType stoich = base.Global(substage).StoichiometricMatrix();
		return SubRows<EXT>(SubCols<EXT>(stoich,constsubreactions),subspecies);
	}

	size_t MobileSpecies() const {
		size_t result=0;
		for(size_t s : subspecies) {
			if(s<base.Global(substage).MobileSpecies()) ++result;
		}
		return result;
	}

private:
	const IKineticContainer<EXT>& base;
	const std::vector<size_t> subspecies;
	std::vector<size_t> constsubreactions;
	const size_t substage;

	std::vector<size_t> RelevantCols() const {
		std::vector<size_t> result;
		const MatrixType M = base.Global(substage).StoichiometricMatrix();
		for(size_t j=0;j<EXT::cols(M);++j){
			bool found=false;
			for(size_t i : subspecies){
				if(*(EXT::ConstBeginColWise(M,j)+i)!=0) found=true;
			}
			if(found) result.push_back(j);
		}
		return result;
	}		
 };

 template<typename EXT>
 class SpeciesStructure{
  public:
	class container{
		public:
		container(){}
		template<typename IT>
		container(size_t i, IT begin, IT end){
			for(size_t pos=0; begin!=end; ++begin, ++pos){
				if(*begin!=0) entries.push_back(pos);
			}
			if(std::find(entries.begin(),entries.end(),i)==entries.end()) entries.push_back(i);
			std::sort(entries.begin(),entries.end());
			
		}
		bool operator<(const container& other) const {
			if(entries.size()==other.entries.size()) return entries[0]<other.entries[0];
			else return entries.size()<other.entries.size();
		}
		bool operator==(const container& other) const { return other.entries == entries; }
		const std::vector<size_t>& Entries() const { return entries; }
		void Remove(const container& other){
			std::vector<size_t>::iterator end = entries.end();
			for(size_t x : other.Entries()) end = std::remove(entries.begin(),end,x);
			entries.resize(std::distance(entries.begin(),end));
		}

		private:
		std::vector<size_t> entries;
	};
			
	SpeciesStructure(typename EXT::MatrixType M) {
		while(UpdateEntries(M)){}

		for(size_t i=0;i<EXT::rows(M);++i){
			structure.push_back(container(i,EXT::BeginRowWise(M,i),EXT::EndRowWise(M,i)));
		}

		std::sort(structure.begin(),structure.end());
		typename std::vector<container>::iterator end = std::unique(structure.begin(),structure.end());
		structure.resize(std::distance(structure.begin(),end));

		for(size_t i=0;i<structure.size();++i){
			for(size_t j=i+1;j<structure.size();++j) structure[j].Remove(structure[i]);
			std::sort(structure.begin()+i+1,structure.end());
			end = std::unique(structure.begin()+i,structure.end());
			structure.resize(std::distance(structure.begin(),end));
		}
	}
	
	typename std::vector<container>::const_iterator begin() const { return structure.begin(); }
	typename std::vector<container>::const_iterator end() const { return structure.end(); }

	template<typename IT>
	static bool UpdateEntries(IT begin1, IT end1, IT begin2){
		for(;begin1!=end1;++begin1,++begin2){
			if(*begin2!=0 && *begin1==0){
				*begin1=1;
				return true;
			}
		}
		return false;
	}
	static bool UpdateEntries(typename EXT::MatrixType& M){
		for(size_t i=0;i<EXT::rows(M);++i){
			for(size_t j=0;j<EXT::cols(M);++j){
				if(*(EXT::BeginRowWise(M,i)+j)!=0){
					if(UpdateEntries(EXT::BeginRowWise(M,i),EXT::EndRowWise(M,i),EXT::BeginRowWise(M,j))) return true;
				}
			}
		}
		return false;
	}
 private:
	std::vector<container> structure;
 };
	
 template<typename EXT>
 class OneSidedCouplingsStageArray : public HierarchicalStageInfo<EXT>{
  private:
	const IKineticContainer<EXT>& base;
	std::vector< OneSidedCouplingsStage<EXT> > stages;
	OneSidedCouplingsStageArray(); //FORBID

	void ProcessSubStage(size_t s) {
		SpeciesStructure<EXT> structure(base.SpeciesStructure(s));
		std::vector<size_t> constreactions;
		for(auto x : structure) {
			OneSidedCouplingsStage<EXT> stage(base,x.Entries(),constreactions,s);
			for(size_t r : stage.VarSubReactions()) constreactions.push_back(r);
			stages.push_back(stage);
		}
	}

  public:
	typedef typename EXT::MatrixType MatrixType;
	OneSidedCouplingsStageArray(const IKineticContainer<EXT>& Base) : base(Base){
		for(size_t i=0;i<base.Stages();++i) ProcessSubStage(i);
	}

	size_t Stages() const { return stages.size(); }
	size_t SubStage(size_t stage) const { return stages[stage].SubStage(); }
	std::vector<size_t> SubReactions(size_t stage) const { return stages[stage].VarSubReactions(); }
	MatrixType StoichiometricMatrix(size_t stage) const { return stages[stage].StoichiometricMatrix(); }
	MatrixType ToBase(size_t stage) const { return stages[stage].ToBase(); }
	MatrixType FromBaseGlobal(size_t stage) const { return stages[stage].FromBaseGlobal(); }
	size_t MobileSpecies(size_t stage) const { return stages[stage].MobileSpecies(); }
	bool ForceNoCorrection(size_t stage) const { return false; }
	const IKineticContainer<EXT>& Base() const { return base; }

	typename std::vector<OneSidedCouplingsStage<EXT> >::const_iterator begin() const { return stages.begin(); }
	typename std::vector<OneSidedCouplingsStage<EXT> >::const_iterator end() const { return stages.end(); }
};

 template<typename EXT, class BT = IKineticContainer<EXT> >
 class OneSidedCouplings : public IHierarchicalKineticContainer<EXT,BT > {
 public:
 	typedef typename EXT::VectorType VectorType;
	typedef typename EXT::VectorArrayType VectorArrayType;
	typedef typename EXT::MatrixType MatrixType;

	OneSidedCouplings(const OneSidedCouplingsStageArray<EXT>& info, const BT& bt)
	                  : IHierarchicalKineticContainer<EXT,BT >(info, bt) {
		Load(info);
	}

	OneSidedCouplings(const BT& bt) : IHierarchicalKineticContainer<EXT,BT >(OneSidedCouplingsStageArray<EXT>(bt), bt) {
		Load(OneSidedCouplingsStageArray<EXT>(bt));
	}

	VectorType ConstSpeciesRatesImpl(const VectorType& all, size_t stage) const {
		return stoich_const[stage]*Base().SubReactionRatesImpl(all,const_subreactions[stage],SubStages()[stage])
		       + FromBaseGlobal()[stage]*Base().ConstSpeciesRatesImpl(all,SubStages()[stage]);
	}
	VectorPair<EXT> ConstSpeciesRatesImpl(const VectorType& mobile, const VectorType& immobile, size_t stage) const {
		const VectorType myconst = Base().SubReactionRatesImpl(mobile,immobile,const_subreactions[stage],SubStages()[stage]);
		const VectorPair<EXT> baseconst = Base().ConstSpeciesRatesImpl(mobile,immobile,SubStages()[stage]);
		return VectorPair<EXT>(stoich_const_mobile[stage]*myconst+FromBaseMobile()[stage]*baseconst.Mobile(),
		                       stoich_const_immobile[stage]*myconst+FromBaseImmobile()[stage]*baseconst.Immobile());
	}

	MatrixType RateStructure(size_t stage) const {
		return EXT::CreateMatrix(this->SubReactions()[stage].size(),EXT::rows(stoich_const[stage]),1);
	}
		
	IKineticContainer<EXT>* copy() const { return new OneSidedCouplings(*this); }

 private:
	std::vector< std::vector<size_t> > const_subreactions;
	std::vector<MatrixType> stoich_const, stoich_const_mobile, stoich_const_immobile;

	using IHierarchicalKineticContainer<EXT,BT>::Base;
	using IHierarchicalKineticContainer<EXT,BT>::SubStages;
	using IHierarchicalKineticContainer<EXT,BT>::FromBaseGlobal;
	using IHierarchicalKineticContainer<EXT,BT>::FromBaseMobile;
	using IHierarchicalKineticContainer<EXT,BT>::FromBaseImmobile;

	void Load(const OneSidedCouplingsStageArray<EXT>& info){
		for(auto x : info){
			const_subreactions.push_back(x.ConstSubReactions());
			stoich_const.push_back(x.ConstStoichiometricMatrix());
			MatrixPair<EXT> tmp = DivideRows<EXT>(x.ConstStoichiometricMatrix(),x.MobileSpecies());
			stoich_const_mobile.push_back(tmp.Mobile());
			stoich_const_immobile.push_back(tmp.Immobile());
		}
	}
};

template<typename EXT, class BT = IKineticContainer<EXT> >
using OneSidedStoichiometry = OneSidedCouplings<EXT,BT>;
}

#endif
