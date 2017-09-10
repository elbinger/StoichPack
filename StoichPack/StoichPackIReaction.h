#ifndef __H_STOICHPACK_REACTION__
#define __H_STOICHPACK_REACTION__

#include "StoichPackSpecies.h"
#include <vector>
#include <algorithm>

namespace StoichPack{

 class IReaction{
 private:
	std::vector<Species> participants;
	std::vector<sp_scalar> coefficients;

	IReaction& operator=(const IReaction&); // FORBID

 protected:
	void AddSpecies(const std::string& species, sp_scalar weight) {
		if(weight==0) throw StoichPackException("Reaction "+Name()+": Cannot add species "+species+" with weight 0!");
		participants.push_back(Species(species));
		coefficients.push_back(weight);
	}

 public:
	const std::vector<Species>& Participants() const { return participants; }
	const std::vector<sp_scalar>& Coefficients() const { return coefficients; }

	void Initialize(const std::vector<Species>& all) {
		if(participants.size()==0) throw StoichPackException("Reaction "+Name()+": no participants");

		for(std::vector<Species>::iterator part=participants.begin();part!=participants.end();++part){
			std::vector<Species>::const_iterator pos = find(all.begin(),all.end(),*part);
			if(pos==all.end()) throw StoichPackException("Reaction "+Name()+": Unknown species "+part->Name());
			else part->Update(*pos);
		}
	}

	bool Initialized() const {
		if(participants.size()==0) throw StoichPackException("Reaction "+Name()+": no participants");
		return participants[0].Initialized();
	}

	/*bool AllConst() const {
		for(size_t i=0;i<participants.size();++i){
			if(participants[i].Type()!=species_type::constant)
				return false;
			}
		}
		return participants.size()!=0;
	}*/

	template<typename ITV>
	std::vector<sp_scalar> GetValues(ITV v) const {
		const size_t s=participants.size();
		std::vector<sp_scalar> ret(s);
		for(size_t i=0;i<s;++i) ret[i]=participants[i].GetValue(v);
		return ret;
	}

	template<typename ITM, typename ITIMM>
	std::vector<sp_scalar> GetValues(ITM m, ITIMM imm) const {
		const size_t s=participants.size();
		std::vector<sp_scalar> ret(s);
		for(size_t i=0;i<s;++i) ret[i]=participants[i].GetValue(m,imm);
		return ret;
	}

	template<typename ITM, typename ITIMM>
	void Add(ITM itm, ITIMM itimm, const std::vector<sp_scalar>& c) const {
		for(size_t i=0;i<participants.size();++i){
			if(participants[i].Type()==species_type::mobile) *(itm+participants[i].TypePos())+=c[i];
			else if(participants[i].Type()==species_type::immobile) *(itimm+participants[i].TypePos())+=c[i];
		}
	}

	template<typename IT>
	void Add(IT it, const std::vector<sp_scalar>& c) const {
		for(size_t i=0;i<participants.size();++i){
			/*if(participants[i].Type()!=species_type::constant)*/ *(it+participants[i].GlobalPos())+=c[i];
		}
	}

	virtual std::string Name() const { return "<Reaction>"; }

	virtual ~IReaction(){}
 };

}

#endif
