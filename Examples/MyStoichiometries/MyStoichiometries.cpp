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

void Cycle(KineticReactionArray<>& r, const std::vector<string>& species){
	assert(species.size()>1);
	for(size_t i=0;i<species.size()-1;++i) r.AddReaction(new MichaelisMenten(species[i],species[i+1],1,1));
	r.AddReaction(new MichaelisMenten(species.back(),species.front(),1,1));
}

void Example3(KineticReactionArray<>& r){
	std::vector<string> cyc11, cyc12, cyc13, cyc21, cyc22, cyc3;
	cyc11.push_back("A11"); cyc11.push_back("B11"); cyc11.push_back("C11"); cyc11.push_back("D11");
	cyc12.push_back("A12"); cyc12.push_back("B12"); cyc12.push_back("C12"); cyc12.push_back("D12");
	cyc13.push_back("A13"); cyc13.push_back("B13"); cyc13.push_back("C13"); cyc13.push_back("D13");
	cyc21.push_back("A21"); cyc21.push_back("B21"); cyc21.push_back("C21"); cyc21.push_back("D21");
	cyc22.push_back("A22"); cyc22.push_back("B22"); cyc22.push_back("C22"); cyc22.push_back("D22");
	cyc3.push_back("A3"); cyc3.push_back("B3"); cyc3.push_back("C3"); cyc3.push_back("D3");

	r.AddReaction(new MichaelisMenten("B22","B3",1,1));
	r.AddReaction(new MichaelisMenten("C21","C3",1,1));
	Cycle(r,cyc13);
	r.AddReaction(new MichaelisMenten("A13","A22",1,1));
	Cycle(r,cyc21);
	r.AddReaction(new MichaelisMenten("A11","A21",1,1));
	Cycle(r,cyc3);
	Cycle(r,cyc11);
	r.AddReaction(new MichaelisMenten("B12","B22",1,1));
	r.AddReaction(new MichaelisMenten("C12","C21",1,1));
	Cycle(r,cyc22);
	Cycle(r,cyc12);

	r.AddSpecies("C11",species_type::mobile);
	r.AddSpecies("C12",species_type::immobile);
	r.AddSpecies("C13",species_type::immobile);
	r.AddSpecies("C21",species_type::immobile);
	r.AddSpecies("C22",species_type::immobile);
	r.AddSpecies("C3",species_type::immobile);
	r.AddSpecies("D11",species_type::mobile);
	r.AddSpecies("D12",species_type::mobile);
	r.AddSpecies("D13",species_type::immobile);
	r.AddSpecies("D21",species_type::immobile);
	r.AddSpecies("D22",species_type::immobile);
	r.AddSpecies("D3",species_type::immobile);
	r.AddSpecies("B11",species_type::mobile);
	r.AddSpecies("B12",species_type::mobile);
	r.AddSpecies("B13",species_type::mobile);
	r.AddSpecies("B21",species_type::immobile);
	r.AddSpecies("B22",species_type::immobile);
	r.AddSpecies("B3",species_type::immobile);
	r.AddSpecies("A11",species_type::mobile);
	r.AddSpecies("A12",species_type::mobile);
	r.AddSpecies("A13",species_type::mobile);
	r.AddSpecies("A21",species_type::mobile);
	r.AddSpecies("A22",species_type::immobile);
	r.AddSpecies("A3",species_type::mobile);
}

BasicKineticStoichiometry<EXT> LoadSystem(const string& name){
	BasicKineticStoichiometry<EXT> ret;

	if(name=="Example1") Example1(ret.reactions);
	else if(name=="Example2") Example2(ret.reactions);
	else if(name=="Example3") Example3(ret.reactions);
	else throw StoichPackException(string("Cound not find system ")+name);

	ret.Initialize();
	return ret;
}

BasicKineticStoichiometry<EXT> LoadStoichiometry(const string& name){ return LoadSystem(name); }

