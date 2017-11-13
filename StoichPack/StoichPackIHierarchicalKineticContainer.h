/* File: StoichPackIHierarchicalKineticContainer.h
 * Purpose: Define an interface for hierarchical containers.
 * Author: Tobias Elbinger <elbinger@math.fau.de> */

#ifndef __H_STOICHPACK_HIERARCHICAL_KINETIC_CONTAINER__
#define __H_STOICHPACK_HIERARCHICAL_KINETIC_CONTAINER__

#include "StoichPackIKineticContainer.h"

namespace StoichPack{

 template<typename EXT>
 class HierarchicalStageInfo {
 public:
	typedef typename EXT::MatrixType MatrixType;
	virtual size_t Stages() const = 0;
	virtual size_t SubStage(size_t stage) const = 0;
	virtual std::vector<size_t> SubReactions(size_t stage) const = 0;
	virtual MatrixType StoichiometricMatrix(size_t stage) const = 0;	
	virtual MatrixType ToBase(size_t stage) const = 0;
	virtual MatrixType FromBaseGlobal(size_t stage) const = 0;
	virtual size_t MobileSpecies(size_t stage) const =0;
	virtual bool ForceNoCorrection(size_t stage) const =0;
	virtual const IKineticContainer<EXT>& Base() const = 0;

	StageInfoArray<EXT> GetStageInfoArray() const {
		std::vector<BasicStageInfo<EXT> > result;
		for(size_t i=0;i<Stages();++i){
			bool correction = (!ForceNoCorrection(i)) && Base().Global(SubStage(i)).Correct();
			const MatrixType toOriginal = Base().Global(SubStage(i)).ToOriginal()*ToBase(i);
			const MatrixType fromOriginal = FromBaseGlobal(i)*Base().Global(SubStage(i)).FromOriginal();
			BasicStageInfo<EXT> tmp(StoichiometricMatrix(i),toOriginal,fromOriginal,MobileSpecies(i),Base().MobileSpecies(),correction);
			result.push_back(tmp);
		}
		return StageInfoArray<EXT>(result);
	}

	MatrixType FromBaseMobile(size_t stage) const {
		return SplitBlocks<EXT>(FromBaseGlobal(stage),MobileSpecies(stage),Base().Global(SubStage(stage)).MobileSpecies()).Mobile();
	}
	MatrixType FromBaseImmobile(size_t stage) const {
		return SplitBlocks<EXT>(FromBaseGlobal(stage),MobileSpecies(stage),Base().Global(SubStage(stage)).MobileSpecies()).Immobile();
	}

	virtual ~HierarchicalStageInfo() {}

 private:
	MatrixPair<EXT> BlocksFrom(size_t stage) const {
		return SplitBlocks<EXT>(FromBaseGlobal(stage),MobileSpecies(),Base().Global(SubStage(stage)).MobileSpecies());
	}
 };

 /* class BaseStorage: stroage for underlying container. */
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

	/* Provide a function to evaluate all reaction rates of a stage, given the values of all concentrations (original presentation). */
	VectorType ReactionRatesImpl(const VectorType& all, size_t stage) const {
		return Base().SubReactionRatesImpl(all,subreactions[stage],substage[stage]);
	}
	/* Provide a funtion that only evaluates the reactions speciefied in I. */
	VectorType SubReactionRatesImpl(const VectorType& all, const std::vector<size_t>& I, size_t stage) const {
		return Base().SubReactionRatesImpl(all,GetSub(I,stage),substage[stage]);
	}
	/* Provide a function to evaluate all reaction rates of a stage, given the values of all mobile concentrations
	 * and all immobile concentrations separately (original presentation). */
	VectorType ReactionRatesImpl(const VectorType& mobile, const VectorType& immobile, size_t stage) const {
		return Base().SubReactionRatesImpl(mobile,immobile,subreactions[stage],substage[stage]);
	}
	/* Provide a funtion that only evaluates the reactions speciefied in I. */
	VectorType SubReactionRatesImpl(const VectorType& mobile, const VectorType& immobile,
	                                const std::vector<size_t>& I, size_t stage) const {
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

	const std::vector<InitializedSpecies>& Participants() const { return Base().Participants(); }

	/*the following pure virtual functions need to be implemented:
	virtual VectorType ConstSpeciesRatesImpl(const VectorType& all, size_t stage) const =0;
	virtual VectorPair<EXT> ConstSpeciesRatesImpl(const VectorType& mobile, const VectorType& immobile, size_t stage) const =0;
	virtual MatrixType RateStructure(size_t stage) const = 0;
	virtual IKineticContainer<EXT>* copy() const =0*/

 protected:
	const BT& Base() const { return base.get(); }
	const std::vector< std::vector<size_t> >& SubReactions() const { return subreactions; }
	const std::vector<size_t>& SubStages() const { return substage; }
	const std::vector<MatrixType>& FromBaseGlobal() const { return fromBaseGlobal; }
	const std::vector<MatrixType>& FromBaseMobile() const { return fromBaseMobile; }
	const std::vector<MatrixType>& FromBaseImmobile() const { return fromBaseImmobile; }

	std::vector<size_t> GetSub(const std::vector<size_t>& index, size_t stage) const {
		std::vector<size_t> result;
		result.reserve(index.size());
		for(size_t i : index) result.push_back(subreactions[stage][i]);
		return result;
	}

	IHierarchicalKineticContainer(const HierarchicalStageInfo<EXT>& info, const BT& Base)
	                              : IKineticContainer<EXT>(info.GetStageInfoArray()), base(Base) {
		for(size_t i=0;i<this->Stages();++i){
			subreactions.push_back(info.SubReactions(i));
			substage.push_back(info.SubStage(i));

			fromBaseGlobal.push_back(info.FromBaseGlobal(i));
			fromBaseMobile.push_back(info.FromBaseMobile(i));
			fromBaseImmobile.push_back(info.FromBaseImmobile(i));
		}
	}
 private:
	const BaseStorage<EXT,BT> base;
	std::vector<std::vector<size_t> > subreactions;
	std::vector<size_t> substage;
	std::vector<MatrixType> fromBaseGlobal, fromBaseMobile,fromBaseImmobile;

	IHierarchicalKineticContainer(); //FORBID
};

 //ensure backward compatibility for older naming convention
 template<typename EXT, typename BT>
 using IHierarchicalStoichiometry = IHierarchicalKineticContainer<EXT,BT>;
}

#endif
