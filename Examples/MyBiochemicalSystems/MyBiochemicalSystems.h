#ifndef __H_MYSTOICHIOMETRIES__
#define __H_MYSTOICHIOMETRIES__

#include "StoichPackBiochemicalSystem.h"
#include "../Utility/FDMVector.h"

StoichPack::BiochemicalSystem<> LoadSystem(const std::string& name);
std::vector<StoichPack::sp_scalar> GetDefaultValuesImpl(const std::string& name, const StoichPack::InitializedBiochemicalSystem<>& system,
                                                        const std::vector<StoichPack::sp_scalar>& X);

template<typename EXT>
typename EXT::VectorType GetDefaultValues(const std::string& name, const StoichPack::InitializedBiochemicalSystem<>& system){
	const std::vector<StoichPack::sp_scalar>& data = GetDefaultValuesImpl(name,system,std::vector<StoichPack::sp_scalar>(0));
	typename EXT::VectorType result = EXT::CreateVector(data.size());
	std::copy(data.begin(),data.end(),EXT::Begin(result));
	return result;
}

template<typename EXT>
FDMVector GetDefaultValues(const std::string& name,const StoichPack::InitializedBiochemicalSystem<>& system, const FDMMesh& mesh){
	FDMVector u(system.Count(StoichPack::species_type::mobile)+system.Count(StoichPack::species_type::immobile),mesh);
	for(int i=0;i<mesh.Nodes();++i) u.SetSub(i,GetDefaultValues<EXT>,name,system,mesh.Coordinates(i));
	return u;
}

#endif
