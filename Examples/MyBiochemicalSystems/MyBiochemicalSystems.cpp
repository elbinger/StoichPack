#include "StoichPackBiochemicalSystem.h"
#include "MyReactions.h"
#include "SpatialModifier.h"

BiochemicalSystem<> System1(){
	BiochemicalSystem<> r;
	r.AddKineticReaction(new MichaelisMenten("A","B",1,1));
	r.AddKineticReaction(new MichaelisMenten("B","C",1,1));

	r.AddSpecies({"A", "C"},species_type::mobile);
	r.AddSpecies("B",species_type::immobile);

	return r;
}

std::vector<SpatialModifier> System1Default(){
	return { SpatialModifier("A",1,{false,false,false}) };
}

BiochemicalSystem<> System2(){
	BiochemicalSystem<> r = System1();
	r.AddKineticReaction(new MichaelisMenten("C","B",1,1));
	return r;
}

std::vector<SpatialModifier> System2Default(){
	return System1Default();
}


void Cycle(BiochemicalSystem<>& r, const std::vector<string>& species){
	assert(species.size()>1);
	for(size_t i=0;i<species.size()-1;++i) r.AddKineticReaction(new MichaelisMenten(species[i],species[i+1],1,1));
	r.AddKineticReaction(new MichaelisMenten(species.back(),species.front(),1,1));
}

BiochemicalSystem<> System3(){
	BiochemicalSystem<> r;
	std::vector<string> cyc11, cyc12, cyc13, cyc21, cyc22, cyc3;
	cyc11.push_back("A11"); cyc11.push_back("B11"); cyc11.push_back("C11"); cyc11.push_back("D11");
	cyc12.push_back("A12"); cyc12.push_back("B12"); cyc12.push_back("C12"); cyc12.push_back("D12");
	cyc13.push_back("A13"); cyc13.push_back("B13"); cyc13.push_back("C13"); cyc13.push_back("D13");
	cyc21.push_back("A21"); cyc21.push_back("B21"); cyc21.push_back("C21"); cyc21.push_back("D21");
	cyc22.push_back("A22"); cyc22.push_back("B22"); cyc22.push_back("C22"); cyc22.push_back("D22");
	cyc3.push_back("A3"); cyc3.push_back("B3"); cyc3.push_back("C3"); cyc3.push_back("D3");

	r.AddKineticReaction(new MichaelisMenten("B22","B3",1,1));
	r.AddKineticReaction(new MichaelisMenten("C21","C3",1,1));
	Cycle(r,cyc13);
	r.AddKineticReaction(new MichaelisMenten("A13","A22",1,1));
	Cycle(r,cyc21);
	r.AddKineticReaction(new MichaelisMenten("A11","A21",1,1));
	Cycle(r,cyc3);
	Cycle(r,cyc11);
	r.AddKineticReaction(new MichaelisMenten("B12","B22",1,1));
	r.AddKineticReaction(new MichaelisMenten("C12","C21",1,1));
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

	return r;
}

std::vector<SpatialModifier> System3Default(){
	SpatialModifier A11("A11",1,{false,false,false});
	SpatialModifier A12("A12",1,{true,false,false});
	SpatialModifier A13("A13",1,{false,true,false});
	SpatialModifier A21("A21",1,{true,true,false});
	SpatialModifier A22("A22",1,{false,false,true});
	SpatialModifier A3("A3",1,{true,false,true});
	SpatialModifier C11("C11",0.5,{false,true,true});
	SpatialModifier C12("C12",0.5,{true,true,true});
	SpatialModifier C13("C13",0.5,{false,true,true});
	SpatialModifier C21("C21",0.5,{true,false,true});
	SpatialModifier C22("C22",0.5,{false,false,true});
	SpatialModifier C3("C3",0.5,{true,true,false});
	return {A11,A12,A13,A21,A22,A3,C11,C12,C13,C21,C22,C3};
}

BiochemicalSystem<> LoadSystem(const string& name){
	if(name=="System1") return System1();
	else if(name=="System2") return System2();
	else if(name=="System3") return System3();
	else throw StoichPackException(string("Cound not find system ")+name);
}

vector<SpatialModifier> GetDefaultModifiers(const string& name){
	if(name=="System1") return System1Default();
	else if(name=="System2") return System2Default();
	else if(name=="System3") return System3Default();
	else throw StoichPackException(string("Cound not find system ")+name);
}
	
vector<sp_scalar> GetDefaultValuesImpl(const string& name, const vector<InitializedSpecies>& participants, const vector<sp_scalar>& X){
	vector<sp_scalar> result;

	const vector<SpatialModifier> mod = GetDefaultModifiers(name);

	for(auto s : participants){
		vector<SpatialModifier>::const_iterator it = std::find(mod.begin(),mod.end(),s.Name());
		if(it==mod.end()) result.push_back(0);
		else result.push_back(it->GetValue(X));
	}
	return result;
}

