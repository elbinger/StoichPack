#include "MyStoichiometries.h"
#include "MyReactions.h"

void Example1(KineticReactionArray<>& r){
	r.AddReaction(new MichaelisMenten("A","B",1,1));
	r.AddReaction(new MichaelisMenten("B","C",1,1));

	r.AddSpecies("A",species_type::mobile);
	r.AddSpecies("B",species_type::immobile);
	r.AddSpecies("C",species_type::mobile);
}

void Example2(KineticReactionArray<>& r){
	Example1(r);
	r.AddReaction(new MichaelisMenten("C","B",1,1));
}

BasicKineticStoichiometry<EXT>* LoadStoichiometry(const string& name){
	BasicKineticStoichiometry<EXT>* ret = new BasicKineticStoichiometry<EXT>();

	if(name=="Example1") Example1(ret->stoichiometry);
	else if(name=="Example2") Example2(ret->stoichiometry);
	else{
		delete ret;
		throw StoichPackException(string("Cound not find stoichiometry ")+name);
	}

	ret->Initialize();
	return ret;
}

