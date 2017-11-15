/* File: StoichPackOneSidedCouplings.h
 * Purpose: Provide basic informations about stages.
 * Author: Tobias Elbinger <elbinger@math.fau.de> */

/* One sided couplings:
 * Imagine the following biochemical system:
 * * 3 species, called A, B, C.
 * * a reaction r1 from A to B, the reaction rate of r1 does not depend on B.
 * * a reaction r2 from B to C, the reaction rate of r2 does not depend on C.
 * 
 * We can solve the system in the following way:
 * 1) solve for A (B and C are not needed)
 * 2) solve for B (A is known, C is not needed)
 * 3) solve for C (B is known)
 *
 * In terms of StoichPack, this means creating a reaction container with 3 stages:
 * * stage 0 contains the species A
 * * stage 1 contains the species B
 * * stage 2 contains the species C
 *
 * Remark: In stage 1, the reaction r1 will be treated as constant reaction, i.e. as some kind of source term.
 *         The same happens for reaction r2 in stage 2. */

#ifndef __H_STOICHPACK_ONESIDED_KINETIC_COUPLINGS__
#define __H_STOICHPACK_ONESIDED_KINETIC_COUPLINGS__

#include "StoichPackIHierarchicalKineticContainer.h"

namespace StoichPack{

//Provide basic information about a stage of OneSidedCouplings
template<typename EXT>
class OneSidedCouplingsStage{
private:
	const IKineticContainer<EXT>& base; //the original container
	const std::vector<size_t> subspecies; //all species involved in this stage
	std::vector<size_t> constsubreactions; //all reactions that are considered as constant
	const size_t substage; //the corresponding stage of the original container

	//find all reactions that have in infulence on species involved in this stage
	//Those reactions are called the relevant reactions.
	std::vector<size_t> RelevantReactions() const {
		std::vector<size_t> result;
		const MatrixType M = base.Global(substage).StoichiometricMatrix();
		//iterate over all columns of M
		for(size_t j=0;j<EXT::cols(M);++j){
			bool found=false;
			for(size_t i : subspecies){ //iterate over all involved species
				if(*(EXT::ConstBeginColWise(M,j)+i)!=0) //the species is influenced
					found=true;
			}
			if(found) result.push_back(j); //mark as relevant
		}
		return result;
	}

	OneSidedCouplingsStage(); //FORBID

public:
	typedef typename EXT::MatrixType MatrixType;

	/* Constructor:
	 * Parameters:
	 * * _Base: a reference to the original container.
	 * * _SubSpecies: all involved species.
	 * * _ConstReactions: all reactions that can for sure be treated as constant
	 * * _SubStage: the corresponding stage of the original container. */
	OneSidedCouplingsStage(const IKineticContainer<EXT>& _Base, const std::vector<size_t>& _SubSpecies,
		const std::vector<size_t>& _ConstReactions, size_t _SubStage) :
			base(_Base), subspecies(_SubSpecies), substage(_SubStage) {
		// the only uninitialzed member is constsubreactions
		// we only store relevant reactions that can be treated as constant
		for(size_t j : RelevantReactions()){ // iterate over relevant reactions
			if(std::find(_ConstReactions.begin(),_ConstReactions.end(),j)!=_ConstReactions.end())
				//relevant and constant
				constsubreactions.push_back(j);
		}
	}

	//returns all relevant reactions that can not be treated as constant
	std::vector<size_t> VarSubReactions() const {
		std::vector<size_t> result;
		for(size_t j : RelevantReactions()){ //iterate over relevant reactions
			if(std::find(constsubreactions.begin(),constsubreactions.end(),j)==constsubreactions.end())
				//this reaction cannot be treated as constant
				result.push_back(j);
		}
		return result;
	}

	//Returns a matrix describing the contribution of this stage to the original solution (all species).
	//This is similar to a permuatation matrix.
	MatrixType ToBase() const{
		MatrixType result = EXT::CreateMatrix(base.Global(substage).GlobalSpecies(),subspecies.size(),0);
		for(size_t j=0;j<subspecies.size();++j){
			*(EXT::BeginColWise(result,j)+subspecies[j])=1;
		}
		return result;
	}

	//Returns a matrix describing the contribution of an original solution the values in this stage (all species).
	//This is similar to the inverse permutation of the permutation mentioned in ToBase.
	MatrixType FromBaseGlobal() const { return EXT::Transposed(ToBase()); }

	//The stoichiometric matrix, i.e. the contribution of all relevant reactions that cannot be treated as constant.
	MatrixType StoichiometricMatrix() const {
		const MatrixType stoich = base.Global(substage).StoichiometricMatrix();
		//SubRows<EXT>(...,subspecies): only the rows corresponding to the involved species.
		//SubCols<EXT>(...,VarSubReactions()): only the columns corresponding to relevant reactions
		//                                     that cannot be treated as constant.
		return SubRows<EXT>(SubCols<EXT>(stoich,VarSubReactions()),subspecies);
	}

	//The contibution of relevant reactions that can be treated as constant.
	MatrixType ConstStoichiometricMatrix() const{
		const MatrixType stoich = base.Global(substage).StoichiometricMatrix();
		return SubRows<EXT>(SubCols<EXT>(stoich,constsubreactions),subspecies);
	}

	//Number of Mobile Species in this stage.
	size_t MobileSpecies() const {
		//mobile species are stored in front of immobile species (cf StoichPackSpecies.h)
		// --> count all entries of supspecies that are smaller than the number of mobile species
		//     in the corresponding stage of the original container
		size_t result=0;
		for(size_t s : subspecies) { //iterate over the entries of supspecies
			if(s<base.Global(substage).MobileSpecies()) ++result;
		}
		return result;
	}

	//getters
	size_t SubStage() const { return substage; }
	const std::vector<size_t>& SubSpecies() const { return subspecies; }
	const std::vector<size_t>& ConstSubReactions() const { return constsubreactions; }
 };

 //group species to stages
 template<typename EXT>
 class SpeciesStructure{
  public:
	//a class to store all dependencies of a species
	class container{
		public:
		container(){}

		/* Constructor:
		 * Parameters:
		 * * i: the index of the species for which dependencies are stored
		 * * begin: an iterator to a container that stores dependency information, i.e. if the j-th element
		 * *        of the container is nonzero, species i depends on species j.
		 * * end: iterator to past-end element of the container specified in begin. */
		template<typename IT>
		container(size_t i, IT begin, IT end){
			for(size_t pos=0; begin!=end; ++begin, ++pos){ //iterate over the container
				if(*begin!=0) entries.push_back(pos);
			}
			//species i always depends on species i --> make sure i is contained in entries
			//no element must occur twice!
			if(std::find(entries.begin(),entries.end(),i)==entries.end()) entries.push_back(i);

			//sort entries to assure that the relative oreder is not changed
			std::sort(entries.begin(),entries.end());	
		}

		//helper for sort functions: a container is "less" than another container if it has less entries or
		//if it has the same number of entries but the first species has a lower index
		bool operator<(const container& other) const {
			if(entries.size()==other.entries.size()) return entries[0]<other.entries[0];
			else return entries.size()<other.entries.size();
		}

		bool operator==(const container& other) const { return other.entries == entries; }

		//remove all entries that are equal to an entry of other
		void Remove(const container& other){
			//see stl
			std::vector<size_t>::iterator end = entries.end();
			for(size_t x : other.Entries()) end = std::remove(entries.begin(),end,x);
			entries.resize(std::distance(entries.begin(),end));
		}

		//getter
		const std::vector<size_t>& Entries() const { return entries; }

		private:
		std::vector<size_t> entries; //dependencies
	};

	/* Constructor:
	 * Parameters:
	 * * M: a matrix that specifies whether there is a direct infulence of species j on species i, i.e. the
	 *      entry in the i-th row and j-th column is 0 if there is a reaction that influences species i and is
	 *      influenced by species j. This kind of matrix is returned by IKineticContainer::SpeciesStructure. */
	SpeciesStructure(typename EXT::MatrixType M) {
		while(UpdateEntries(M)){} //create dependency graph

		//M has the following property now:
		//	entry in i-th row and j-th row is 0 <--> species i can be solved independent of species j

		for(size_t i=0;i<EXT::rows(M);++i){
			structure.push_back(container(i,EXT::BeginRowWise(M,i),EXT::EndRowWise(M,i)));
		}

		//Remark: If species i depends on species j, than i depends on all species stored in structure[j],
		//        therefore structure[i]<structure[j] cannot be true.

		//sort and remove multiple entries
		std::sort(structure.begin(),structure.end());
		typename std::vector<container>::iterator end = std::unique(structure.begin(),structure.end());
		structure.resize(std::distance(structure.begin(),end));

		//the species stored in sturcture[0] can be solved independent of all other species.

		for(size_t i=0;i<structure.size();++i){
			//the species in structure[i] only depend on species stored in sturcture[j], j<=i
			//stage i will only contain the species in structure[i]

			//remove all species of stage from remaining entries in structure
			for(size_t j=i+1;j<structure.size();++j) structure[j].Remove(structure[i]);

			//sort and remove multiple entries
			std::sort(structure.begin()+i+1,structure.end());
			end = std::unique(structure.begin()+i,structure.end());
			structure.resize(std::distance(structure.begin(),end));

			//the species stored in structure[i+1] can be solved independent of all species in structure[j], j>i+1
		}
		//each entry of structure represents the species involved in one stage now
	}

	/* UpdateEntries: craeate full dependency graph
	 * Let species i depend on species j and let species j depend on species k --> species i depends on species k. */

	/* Check if updates need to be made in [begin1, end1), where the species corresponding to [begin1, end1) (:=s1)
	 * depends on the species corresponding to [begin2, ...) (:=s2). Update up to one entry if neccessary. */
	template<typename IT>
	static bool UpdateEntries(IT begin1, IT end1, IT begin2){
		for(;begin1!=end1;++begin1,++begin2){
			if(*begin2!=0 && *begin1==0){
				//s2 depends on another species, s1 does not depend on this species
				//--> update and return true
				*begin1=1;
				return true;
			}
		}
		return false; //no updates necessary so far
	}

	/* Check if M needs to be updated in order to represent a full dependency graph, i.e.
	 * 	entry i-th row and j-th column is nonzero <--> species i depends on species j.
	 * Do up to 1 update. If this function returns false, M represents a full dependency graph. */
	static bool UpdateEntries(typename EXT::MatrixType& M){
		for(size_t i=0;i<EXT::rows(M);++i){ // iterate over rows
			for(size_t j=0;j<EXT::cols(M);++j){ // iterate over columns
				if(*(EXT::BeginRowWise(M,i)+j)!=0){
					//species of respective row depends on species of respecive col
					//--> update dependencies
					//return true if update necessary, continue otherwise
					if(UpdateEntries(EXT::BeginRowWise(M,i),EXT::EndRowWise(M,i),EXT::BeginRowWise(M,j))) return true;
				}
			}
		}
		return false; //no update was necessary --> full dependency graph
	}

	//iterators
	typename std::vector<container>::const_iterator begin() const { return structure.begin(); }
	typename std::vector<container>::const_iterator end() const { return structure.end(); }

 private:
	std::vector<container> structure; //involved species for every stage
	SpeciesStructure(); //FORBID
 };

 //Array of OneSidedCouplingsStag in order to provide all information needed for construction of OneSidedCouplings
 template<typename EXT>
 class OneSidedCouplingsStageArray : public IHierarchicalStageInfo<EXT>{
  private:
	const IKineticContainer<EXT>& base; // the original container
	std::vector< OneSidedCouplingsStage<EXT> > stages; // the array
	OneSidedCouplingsStageArray(); //FORBID

	// apply preprocessing to one stage of the original container
	void ProcessSubStage(size_t s) {
		SpeciesStructure<EXT> structure(base.SpeciesStructure(s)); //get substages
		std::vector<size_t> constreactions; //all reactions that can be considered as constant for subsequent stages
		for(auto x : structure) {
			OneSidedCouplingsStage<EXT> stage(base,x.Entries(),constreactions,s); //create stage

			//If a reaction is relevant for the new stage, this means all needed concentrations for this reaction
			//are known after solving this stage and can therefore be treated as constant in subsequent stages.
			//Relevant reactions that can be treated as constant in the new stage are already contained in
			//constreactions, therfore we add all other relevant reactions of the new stage to constreactions. 
			for(size_t r : stage.VarSubReactions()) constreactions.push_back(r); //update constreactions

			stages.push_back(stage); //add stage
		}
	}

  public:
	typedef typename EXT::MatrixType MatrixType;
	OneSidedCouplingsStageArray(const IKineticContainer<EXT>& Base) : base(Base){
		for(size_t i=0;i<base.Stages();++i) ProcessSubStage(i); //preprocess all stages
	}

	/* functions inherited from interface */
	size_t Stages() const { return stages.size(); }
	size_t SubStage(size_t stage) const { return stages[stage].SubStage(); }
	std::vector<size_t> SubReactions(size_t stage) const { return stages[stage].VarSubReactions(); }
	MatrixType StoichiometricMatrix(size_t stage) const { return stages[stage].StoichiometricMatrix(); }
	MatrixType ToBase(size_t stage) const { return stages[stage].ToBase(); }
	MatrixType FromBaseGlobal(size_t stage) const { return stages[stage].FromBaseGlobal(); }
	size_t MobileSpecies(size_t stage) const { return stages[stage].MobileSpecies(); }
	//bool ForceNoCorrection(size_t stage) const { return false; } //not needed in this version
	const IKineticContainer<EXT>& Base() const { return base; }

	//iterators
	typename std::vector<OneSidedCouplingsStage<EXT> >::const_iterator begin() const { return stages.begin(); }
	typename std::vector<OneSidedCouplingsStage<EXT> >::const_iterator end() const { return stages.end(); }
};

 //a reaction container using one sided couplings for preprocessing
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

	/* functions inherited from interface */

	/* Constant reactions:
	 *   reactions that are not constant in the original container but can be treated as constant in this stage
	 * + constant reactions of original container */

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
		//structure of corresponding stage of the original container
		const MatrixType& orig = Base().RateStructure(SubStages()[stage]);

		const MatrixType tmp = orig*EXT::Transposed(FromBaseGlobal()[stage]); // only consider involve species
		return SubRows<EXT>(tmp,this->SubReactions()[stage]); // only consider (non const) relevant reactions
	}
		
	IKineticContainer<EXT>* copy() const { return new OneSidedCouplings(*this); }

 private:
	//reactions that are not not constant in the original container but can be treated as constant in the respective stage.
	std::vector< std::vector<size_t> > const_subreactions;
	//stoichiometric matrices for those reactions
	std::vector<MatrixType> stoich_const, stoich_const_mobile, stoich_const_immobile;

	//access to information that is stored in IHierarchicalKineticContainer
	using IHierarchicalKineticContainer<EXT,BT>::Base;
	using IHierarchicalKineticContainer<EXT,BT>::SubStages;
	using IHierarchicalKineticContainer<EXT,BT>::FromBaseGlobal;
	using IHierarchicalKineticContainer<EXT,BT>::FromBaseMobile;
	using IHierarchicalKineticContainer<EXT,BT>::FromBaseImmobile;

	//initialize members in dependency of a OneSidedCouplingsStageArray
	void Load(const OneSidedCouplingsStageArray<EXT>& info){
		for(auto x : info){ //iterate over stages
			//copy information
			const_subreactions.push_back(x.ConstSubReactions());
			stoich_const.push_back(x.ConstStoichiometricMatrix());
			MatrixPair<EXT> tmp = DivideRows<EXT>(x.ConstStoichiometricMatrix(),x.MobileSpecies());
			stoich_const_mobile.push_back(tmp.Mobile());
			stoich_const_immobile.push_back(tmp.Immobile());
		}
	}
};

//compatibility for older naming convention
template<typename EXT, class BT = IKineticContainer<EXT> >
using OneSidedStoichiometry = OneSidedCouplings<EXT,BT>;
}

#endif
