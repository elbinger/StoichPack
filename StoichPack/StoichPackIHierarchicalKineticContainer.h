/* File: StoichPackIHierarchicalKineticContainer.h
 * Purpose: Define an interface for hierarchical containers.
 * Author: Tobias Elbinger <elbinger@math.fau.de> */


/* Concept of hierarchical containers:
 * As mentioned in REAME_CONCEPTS, hierarchical containers take another container (that may be preprocessed or not) and do (further)
 * preprocessing with this container. This preprocessing is done for each stage of the original container separately in order to
 * benefit from earlier preprocessing steps. Consequently, the following informations a relevant for a general interface for
 * hierarchical containers:
   * the original container itself.
   * mappings from and to concentration values of the original container.
   * a mapping of each stage to the corresponding stage of the original conainer (each stage of the original container may be divided
     into several stages).
   * the set of all relevant reactions.
*/

#ifndef __H_STOICHPACK_HIERARCHICAL_KINETIC_CONTAINER__
#define __H_STOICHPACK_HIERARCHICAL_KINETIC_CONTAINER__

#include "StoichPackIKineticContainer.h"

namespace StoichPack{
 //provide information about stages of hierarchical containers
 template<typename EXT>
 class IHierarchicalStageInfo {
 public:
	typedef typename EXT::MatrixType MatrixType;

	virtual size_t Stages() const = 0; //number of stages

	/* Provide a function that returns the corresponding stage of the original container */
	virtual size_t SubStage(size_t stage) const = 0;

	/* Provide a function that returns all relevant reactions of a given stage */
	virtual std::vector<size_t> SubReactions(size_t stage) const = 0;

	/* Provide a function that returns the stoichiometric matrix of a stage (all species) */
	virtual MatrixType StoichiometricMatrix(size_t stage) const = 0;

	/* Provide a function that returns the contributions of a stage to the corresponding stage of the original container (all species)*/
	virtual MatrixType ToBase(size_t stage) const = 0;

	/* Provide a function that returns the contributions of the corresponding stage of the original container (all species) */
	virtual MatrixType FromBaseGlobal(size_t stage) const = 0;

	/* Mobile species of a stage */
	virtual size_t MobileSpecies(size_t stage) const =0;

	//virtual bool ForceNoCorrection(size_t stage) const =0; //not needed in this version

	/* return a reference to the original container */
	virtual const IKineticContainer<EXT>& Base() const = 0;

	/* return the corresponding StageInfoArray (needed for the constructor of IKineticContainer) */
	StageInfoArray<EXT> GetStageInfoArray() const {
		std::vector<BasicStageInfo<EXT> > result;
		for(size_t i=0;i<Stages();++i){
			//bool correction = (!ForceNoCorrection(i)) && Base().Global(SubStage(i)).Correct();

			//toOriginal: apply ToOriginal on ToBase
			const MatrixType toOriginal = Base().Global(SubStage(i)).ToOriginal()*ToBase(i);

			//fromOriginal: apply FromBase on FromOriginal
			const MatrixType fromOriginal = FromBaseGlobal(i)*Base().Global(SubStage(i)).FromOriginal();

			BasicStageInfo<EXT> tmp(StoichiometricMatrix(i),toOriginal,fromOriginal,MobileSpecies(i),
			                        Base().MobileSpecies()/*,correction*/);
			result.push_back(tmp);
		}
		return StageInfoArray<EXT>(result);
	}

	/* Contributions of the corresponding stage of the original container (mobile species) */
	MatrixType FromBaseMobile(size_t stage) const {
		//mobile species are stored in front of immobile species (cf. StoichPackSpecies.h).
		return SplitBlocks<EXT>(FromBaseGlobal(stage),MobileSpecies(stage),Base().Global(SubStage(stage)).MobileSpecies()).Mobile();
	}

	/* Contributions of the corresponding stage of the original container (immobile species) */
	MatrixType FromBaseImmobile(size_t stage) const {
		//immobile species are stored behind mobile species (cf. StoichPackSpecies.h).
		return SplitBlocks<EXT>(FromBaseGlobal(stage),MobileSpecies(stage),Base().Global(SubStage(stage)).MobileSpecies()).Immobile();
	}

	//dtor : do nothing
	virtual ~IHierarchicalStageInfo() {}
 };

 /* class BaseStorage: stroage for original container. */
 template<typename EXT, typename BT>
 class BaseStorage{
	private:
	const BT storage;
	public:
	BaseStorage(const BT& x) : storage(x) {}
	const BT& get() const { return storage; }
 };

 /* specialization for IKineticContainer<EXT>: storage only possible as pointer */
 template<typename EXT>
 class BaseStorage<EXT,IKineticContainer<EXT> >{
	private:
	const std::shared_ptr<IKineticContainer<EXT> > storage;
	public:
	BaseStorage(const IKineticContainer<EXT>& x) : storage(x.copy()) {}
	const IKineticContainer<EXT>& get() const { return *storage; }
 };

 /* class IHierarchicalKineticContainer: interface for hierarchical containers. */
 template<typename EXT, typename BT = IKineticContainer<EXT> >
 class IHierarchicalKineticContainer : public IKineticContainer<EXT>{
 public:
 	typedef typename EXT::VectorType VectorType;
	typedef typename EXT::VectorArrayType VectorArrayType;
	typedef typename EXT::MatrixType MatrixType;

	/* Functions inherited from interface */

	/* Provide a function to evaluate all reaction rates of a stage, given the values of all concentrations (original presentation). */
	VectorType ReactionRatesImpl(const VectorType& all, size_t stage) const {
		// only use the relevant reactions
		return Base().SubReactionRatesImpl(all,subreactions[stage],substage[stage]);
	}
	/* Provide a funtion that only evaluates the reactions speciefied in I. */
	VectorType SubReactionRatesImpl(const VectorType& all, const std::vector<size_t>& I, size_t stage) const {
		// use a subset of relevant reactions
		return Base().SubReactionRatesImpl(all,GetSub(I,stage),substage[stage]);
	}
	/* Provide a function to evaluate all reaction rates of a stage, given the values of all mobile concentrations
	 * and all immobile concentrations separately (original presentation). */
	VectorType ReactionRatesImpl(const VectorType& mobile, const VectorType& immobile, size_t stage) const {
		// only use relevant reactions
		return Base().SubReactionRatesImpl(mobile,immobile,subreactions[stage],substage[stage]);
	}
	/* Provide a funtion that only evaluates the reactions speciefied in I. */
	VectorType SubReactionRatesImpl(const VectorType& mobile, const VectorType& immobile,
	                                const std::vector<size_t>& I, size_t stage) const {
		// use a subset of relevant reactions
		return Base().SubReactionRatesImpl(mobile,immobile,GetSub(I,stage),substage[stage]);
	}

	/* Provide a function for the Jacobian of ReactionRatesImpl. */
	MatrixType DiffReactionRatesImpl(const VectorType& all, size_t stage) const {
		return Base().DiffSubReactionRatesImpl(all,subreactions[stage],substage[stage]);
	}
	/* Provide a function for the Jacobian of SubReactionRatesImpl. */
	MatrixType DiffSubReactionRatesImpl(const VectorType& all, const std::vector<size_t>& I, size_t stage) const {
		return Base().DiffSubReactionRatesImpl(all,GetSub(I,stage),substage[stage]);
	}
	/* Provide a function for the Jacobian of ReactionRatesImpl. */
	MatrixPair<EXT> DiffReactionRatesImpl(const VectorType& mobile, const VectorType& immobile, size_t stage) const {
		return Base().DiffSubReactionRatesImpl(mobile,immobile,subreactions[stage],substage[stage]);
	}
	/* Provide a function for the Jacobian of SubReactionRatesImpl. */
	MatrixPair<EXT> DiffSubReactionRatesImpl(const VectorType& mobile, const VectorType& immobile,
	                                         const std::vector<size_t>& I, size_t stage) const{
		return Base().DiffSubReactionRatesImpl(mobile,immobile,GetSub(I,stage),substage[stage]);
	}

	const std::vector<InitializedSpecies>& Participants() const {
		//Participants is the same as in the original container
		return Base().Participants();
	}

	/*the following pure virtual functions of IKineticContainer still need to be implemented:
	virtual VectorType ConstSpeciesRatesImpl(const VectorType& all, size_t stage) const =0;
	virtual VectorPair<EXT> ConstSpeciesRatesImpl(const VectorType& mobile, const VectorType& immobile, size_t stage) const =0;
	virtual MatrixType RateStructure(size_t stage) const = 0;
	virtual IKineticContainer<EXT>* copy() const =0*/

 protected:
	const BT& Base() const { return base.get(); } //reference to original container
	const std::vector< std::vector<size_t> >& SubReactions() const { return subreactions; } //relevant reactions
	const std::vector<size_t>& SubStages() const { return substage; } //correponding stage of original container

	/* Contributions of the corresponding stage of the original container (all species) */
	const std::vector<MatrixType>& FromBaseGlobal() const { return fromBaseGlobal; }

	/* Contributions of the corresponding stage of the original container (mobile species) */
	const std::vector<MatrixType>& FromBaseMobile() const { return fromBaseMobile; }

	/* Contributions of the corresponding stage of the original container (immobile species) */
	const std::vector<MatrixType>& FromBaseImmobile() const { return fromBaseImmobile; }

	/* return a subset of relevant reactions for a certain stage */
	std::vector<size_t> GetSub(const std::vector<size_t>& index, size_t stage) const {
		std::vector<size_t> result;
		result.reserve(index.size());
		for(size_t i : index) result.push_back(subreactions[stage][i]);
		return result;
	}

	/* Constructor:
	 * Parameters:
	 * * info: information about the stages.
	 * * Base: the original container. */
	IHierarchicalKineticContainer(const IHierarchicalStageInfo<EXT>& info, const BT& Base)
	                              : IKineticContainer<EXT>(info.GetStageInfoArray()), base(Base) {
		//initialize information for every stage
		for(size_t i=0;i<this->Stages();++i){
			subreactions.push_back(info.SubReactions(i));
			substage.push_back(info.SubStage(i));

			fromBaseGlobal.push_back(info.FromBaseGlobal(i));
			fromBaseMobile.push_back(info.FromBaseMobile(i));
			fromBaseImmobile.push_back(info.FromBaseImmobile(i));
		}
	}
 private:
	const BaseStorage<EXT,BT> base; // storage for original container
	std::vector<std::vector<size_t> > subreactions; // relevant reaction for each stage
	std::vector<size_t> substage; // corresponding stage of the original container for each stage
	std::vector<MatrixType> fromBaseGlobal, fromBaseMobile,fromBaseImmobile; //transformation from original container

	IHierarchicalKineticContainer(); //FORBID
};

 //ensure backward compatibility for older naming convention
 template<typename EXT, typename BT>
 using IHierarchicalStoichiometry = IHierarchicalKineticContainer<EXT,BT>;
}

#endif
