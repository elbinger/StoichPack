/* File: StoichPackBasicKineticContainer.h
 * Purpose: Provide a reaction container for the original representation of the problem.
 * Author: Tobias Elbinger <elbinger@math.fau.de> */

#ifndef __H_STOICHPACK_BASIC_KINETIC_CONTAINER__
#define __H_STOICHPACK_BASIC_KINETIC_CONTAINER__

#include "StoichPackIKineticContainer.h"
#include "StoichPackBiochemicalSystem.h"

namespace StoichPack{

 //Implementation of a biochemical system without preprocessing.
 template<typename EXT, typename KineticType = IKineticReaction>
 class BasicKineticContainer : public IKineticContainer<EXT>{
 typedef InitializedBiochemicalSystem<KineticType> BiochemistryType;
 typedef typename BiochemistryType::KineticReactionType KineticReactionType;
 typedef typename EXT::VectorType VectorType;
 typedef typename EXT::MatrixType MatrixType;

 private:
	const BiochemistryType problem; //underlying biochemical problem
	BasicKineticContainer(); //FORBID

	// return StageInfoArray for system (needed for the constructor of IKineticContainer)
	static StageInfoArray<EXT> Init(const BiochemicalSystem<KineticType>& system){
		BiochemistryType problem(system);

		//assemble stoichiometric matrix
		 const std::vector<KineticReactionType>& r=problem.KineticReactions();
		 const size_t n = problem.Count(species_type::mobile)+problem.Count(species_type::immobile); //number of species
		 MatrixType stoich = EXT::CreateMatrix(n,r.size(),0);

		 //iterate over all reactions and set the corresponding values
		 for(size_t i=0;i<r.size();++i){
		 	r[i].Add(EXT::BeginColWise(stoich,i),r[i].Coefficients());
		 }

		//assemble toBase / fromBase (=identity)
		const MatrixType identity = Diag<EXT>(std::vector<sp_scalar>(n,1.));

		const size_t nmob = problem.Count(species_type::mobile); // number of mobile species

		const std::vector<BasicStageInfo<EXT> > tmp(1, BasicStageInfo<EXT>(stoich,identity,identity,nmob,nmob/*,true*/));
		return StageInfoArray<EXT>(tmp);
	}

	// helper class: all reactions
	class FullIndex{
	 private:
		const size_t n;
		FullIndex();
	 public:
		FullIndex(size_t s) : n(s) {}
		size_t size() const { return n; }
		size_t operator[](size_t i) const { return i; }
	};

	//Evaluate the reaction rates of all reactions specified in I. all contains the concentration values of all species.
	template<class IType>
	VectorType RRates(const VectorType& all, const IType& I) const {
		const size_t nr = I.size(); //number of reactions
		VectorType result = EXT::CreateVector(nr);
		auto it = EXT::Begin(result); //iterator over the entries of ret
		for(size_t i=0;i<nr;++i, ++it){
			//evaluate rate
			*it = problem.KineticReactions()[I[i]].Rate(EXT::ConstBegin(all));
		}
		return result;
	}

	//Evaluate the reaction rates of all reactions specified in I. mobile contains the concentration values of the mobile species,
	//immobile contains the concentration values of the immobile species.
	template<class IType>	
	VectorType RRates(const VectorType& mobile, const VectorType& immobile, const IType& I) const {
		const size_t nr = I.size(); //number of reactions
		VectorType ret = EXT::CreateVector(nr);
		auto it = EXT::Begin(ret);
		for(size_t i=0;i<nr;++i, ++it){
			//evaluate rate
			*it = problem.KineticReactions()[I[i]].Rate(EXT::ConstBegin(mobile),EXT::ConstBegin(immobile));
		}
		return ret;
	}

	/* DiffRRates: Jacobians of RRates */

	template<class IType>
	MatrixType DiffRRates(const VectorType& all, const IType& I) const {
		const size_t nr = I.size();
		MatrixType result = EXT::CreateMatrix(nr,EXT::size(all),0); //initialize with 0, since we use additive write
		for(size_t i=0;i<nr;++i) {
			const KineticReactionType& r = problem.KineticReactions()[I[i]];
			//r.Add will take care about writing the entries to the correct column
			//r.DiffRate is just a local derivative evaluation
			r.Add(EXT::BeginRowWise(result,i),r.DiffRate(EXT::ConstBegin(all)));
		}
		return result;
	}

	template<class IType>
	MatrixPair<EXT> DiffRRates(const VectorType& mobile, const VectorType& immobile, const IType& I) const {
		const size_t nr = I.size();
		//initialize with 0, since we use additive write
		MatrixPair<EXT> result(EXT::CreateMatrix(nr,EXT::size(mobile),0),EXT::CreateMatrix(nr,EXT::size(immobile),0));
		for(size_t i=0;i<nr;++i) {
			const KineticReactionType& r = problem.KineticReactions()[I[i]];
			//r.Add will take care about writing the entries to the correct column
			//r.DiffRate is just a local derivative evaluation
			r.Add(EXT::BeginRowWise(result.Mobile(),i),EXT::BeginRowWise(result.Immobile(),i),
			      r.DiffRate(EXT::ConstBegin(mobile),EXT::ConstBegin(immobile)));
		}
		return result;
	}

 public:
	BasicKineticContainer(const BiochemicalSystem<KineticType>& system) : IKineticContainer<EXT>(Init(system)), problem(system) {}

	const BiochemistryType& Problem() const { return problem; }

	/* functions inherited from interface */
	VectorType ReactionRatesImpl(const VectorType& all, size_t stage) const {
		assert(stage==0);
		FullIndex I(problem.KineticReactions().size()); // all reactions
		return RRates(all,I);
	}

	VectorType ReactionRatesImpl(const VectorType& mobile, const VectorType& immobile, size_t stage) const {
		assert(stage==0);
		FullIndex I(problem.KineticReactions().size()); // all reactions
		return RRates(mobile,immobile,I);
	}


	VectorType SubReactionRatesImpl(const VectorType& all, const std::vector<size_t>& I, size_t stage) const {
		assert(stage==0);
		return RRates(all,I);
	}

	VectorType SubReactionRatesImpl(const VectorType& mobile, const VectorType& immobile,
	                                const std::vector<size_t>& I, size_t stage) const {
		assert(stage==0);
		return RRates(mobile,immobile,I);
	}
	
	MatrixType DiffReactionRatesImpl(const VectorType& all, size_t stage) const {
		assert(stage==0);
		FullIndex I(problem.KineticReactions().size()); // all reactions
		return DiffRRates(all,I);
	}

	MatrixPair<EXT> DiffReactionRatesImpl(const VectorType& mobile, const VectorType& immobile, size_t stage) const {
		assert(stage==0);
		FullIndex I(problem.KineticReactions().size()); // all reactions
		return DiffRRates(mobile,immobile,I);
	}

	MatrixType DiffSubReactionRatesImpl(const VectorType& all, const std::vector<size_t>& I, size_t stage) const {
		assert(stage==0);
		return DiffRRates(all,I);
	}

	MatrixPair<EXT> DiffSubReactionRatesImpl(const VectorType& mobile, const VectorType& immobile,
	                                         const std::vector<size_t>& I, size_t stage) const {
		assert(stage==0);
		return DiffRRates(mobile,immobile,I);
	}

	VectorType ConstSpeciesRatesImpl(const VectorType& all, size_t stage) const {
		return EXT::CreateVector(this->GlobalSpecies(0),0);
	}

	VectorPair<EXT> ConstSpeciesRatesImpl(const VectorType& mobile, const VectorType& immobile, size_t stage) const {
		return VectorPair<EXT>(EXT::CreateVector(this->MobileSpecies(),0), EXT::CreateVector(this->ImmobileSpecies(),0));
	}


	const std::vector<InitializedSpecies>& Participants() const { return problem.Participants(); }

	MatrixType RateStructure(size_t stage) const {
		const size_t nr = problem.KineticReactions().size();
		MatrixType result = EXT::CreateMatrix(nr,this->GlobalSpecies(),0);
		for(size_t i=0;i<nr;++i) {
			//r.Add will take care about writing the entries to the correct column
			//r.Dependencies is just a local structure evaluation
			const KineticReactionType& r = problem.KineticReactions()[i];
			r.Add(EXT::BeginRowWise(result,i),r.Dependencies());
		}
		return result;
	}		

	IKineticContainer<EXT>* copy() const { return new BasicKineticContainer(*this); }
 };

 //ensure compatibility for old naming convention
 template<typename EXT>
 using BasicKineticStoichiometry = BasicKineticContainer<EXT>;
} // namespace StoichPack

#endif
