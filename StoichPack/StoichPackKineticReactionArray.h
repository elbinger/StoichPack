#ifndef __H_STOICPACK_KINETICREACTION_ARRAY__
#define __H_STOICPACK_KINETICREACTION_ARRAY__

#include "StoichPackIKineticReaction.h"
#include <memory>

namespace StoichPack{

template<typename ReactionType = IKineticReaction>
class KineticReactionArray{
	typedef ReactionType* RPTR;
	typedef std::shared_ptr<ReactionType> SHARED_RPTR;
	typedef const ReactionType& RREF;
private:
	std::vector<SHARED_RPTR> reactions;
	std::vector<Species> species;
	size_t mobile, immobile;
	bool initialized;

	void set(const KineticReactionArray& other){
		if(!other.Initialized()){
			throw StoichPackException("Cannot copy uninitalized KineticReactionArray!");
		}
		if(Initialized()) throw StoichPackException("Cannot overwrite initalized KineticReactionArray");
		reactions=other.reactions;
		species=other.species;
		mobile=other.mobile;
		immobile=other.immobile;
		initialized=other.initialized;
	}


public:
	KineticReactionArray() : mobile(0), immobile(0), initialized(false) {}
	KineticReactionArray(const KineticReactionArray& other) : initialized(false) { set(other); }
	KineticReactionArray<ReactionType>& operator=(const KineticReactionArray& other) { set(other); }

	void AddSpecies(const std::string& name, species_type type){
		if(Initialized()) throw StoichPackException("Already initialized!");
		if(type==species_type::unknown) throw StoichPackException("Cannot add species with unknown type!");

		Species s(name,type);
		std::vector<Species>::iterator it=std::find(species.begin(),species.end(),s);
		if(it!=species.end()){
			if(it->Type()!=s.Type())
				throw StoichPackException("Species "+s.Name()+" already contained with different type!");
		} else species.push_back(s);
	}

	void AddReaction(RPTR r){
		if(Initialized()) throw StoichPackException("Already initialized!");
		reactions.push_back(SHARED_RPTR(r));
	}

	void Initialize() {	
		if(initialized) throw StoichPackException("Already initialized!");
		if(species.size()==0) throw StoichPackException("No species!");

		initialized=true;

		class predicate{
			species_type type;
			public:
			predicate(species_type t = species_type::unknown) : type(t) {}
			bool operator()(const Species& x, const Species& y) const {
				if(x.Type()==y.Type()) return x.Name().compare(y.Name())<0;
				if(x.Type()==species_type::mobile) return true;
				//if(x.Type()==species_type::immobile) return y.Type()==species_type::constant;
				return false;
			}
			bool operator()(const Species& x) { return x.Type()==type; }
		};

		std::sort(species.begin(),species.end(),predicate());

		assert(std::is_partitioned(species.begin(),species.end(),predicate(species_type::mobile)));
		mobile = std::count_if(species.begin(),species.end(),predicate(species_type::mobile));
		assert(mobile==0 || species[mobile-1].Type()==species_type::mobile);

		assert(std::is_partitioned(species.begin()+mobile,species.end(),predicate(species_type::immobile)));
		immobile = std::count_if(species.begin()+mobile,species.end(),predicate(species_type::immobile));

		if(mobile+immobile==0) throw StoichPackException("No variable species!");

		for(size_t i=0;i<species.size();++i) species[i].Initialize(i,mobile,immobile);

		for(typename std::vector<SHARED_RPTR>::iterator it=reactions.begin();it!=reactions.end();){
			(*it)->Initialize(species);
			/*if(it->AllConst()) it=erase(it);
			else*/ ++it;
		}
	}

	bool Initialized() const { return initialized; }

	size_t MobileSpecies() const { return mobile; }
	size_t ImmobileSpecies() const { return immobile; }

	const std::vector<Species>& Participants() const { return species; }
	const std::vector<SHARED_RPTR>& Reactions() const { return reactions; }
 };
	
} //namespace

#endif

