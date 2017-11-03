#ifndef __H_MYSTOICHIOMETRIES__
#define __H_MYSTOICHIOMETRIES__

#include "StoichPack.h"
#include "../Utility/FDMVector.h"

StoichPack::BiochemicalSystem<> LoadSystem(const std::string& name);
std::vector<StoichPack::sp_scalar> GetDefaultValuesImpl(const std::string& name,
                                                        const std::vector<StoichPack::InitializedSpecies>& participants,
                                                        const std::vector<StoichPack::sp_scalar>& X);

template<typename EXT>
typename EXT::VectorType GetDefaultValues(const std::string& name, const StoichPack::IKineticContainer<EXT>& S, const std::vector<StoichPack::sp_scalar>& X){
	const std::vector<StoichPack::sp_scalar>& data = GetDefaultValuesImpl(name,S.Participants(),X);
	typename EXT::VectorType result = EXT::CreateVector(data.size());
	std::copy(data.begin(),data.end(),EXT::Begin(result));
	return result;
}

template<typename EXT>
typename EXT::VectorType GetDefaultValues(const std::string& name, const StoichPack::IKineticContainer<EXT>& S){
	return GetDefaultValues<EXT>(name,S,std::vector<StoichPack::sp_scalar>(0));
}


template<typename EXT>
FDMVector GetDefaultValues(const std::string& name, const StoichPack::IKineticContainer<EXT>& S, const FDMMesh& mesh){
	FDMVector u(S.AllSpecies(),mesh);
	for(int i=0;i<mesh.Nodes();++i) u.SetSub(i,GetDefaultValues<EXT>(name,S,mesh.Coordinates(i)));
	return u;
}

#endif
