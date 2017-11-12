/* File:    MyBiochemicalSystems.h
 * Author:  Tobias Elbinger (elbinger@math.fau.de)
 * Purpose: Define functions for loading biochemical systems and default initial values.
 */

#ifndef __H_MYSTOICHIOMETRIES__
#define __H_MYSTOICHIOMETRIES__

#include "StoichPack.h"
#include "../Utility/FDMVector.h"

//Load system with the specified name
StoichPack::BiochemicalSystem<> LoadSystem(const std::string& name);

//Get default initial values for the system with the specified name at coordinates as specified in X.
//The parameter participants defines the storage order of the biochemical species.
std::vector<StoichPack::sp_scalar> GetDefaultValuesImpl(const std::string& name,
                                                        const std::vector<StoichPack::InitializedSpecies>& participants,
                                                        const std::vector<StoichPack::sp_scalar>& X);

//Get default initial values for the system with the specified name that is stored in S at the coordinates specified in X.
template<typename EXT>
typename EXT::VectorType GetDefaultValues(const std::string& name, const StoichPack::IKineticContainer<EXT>& S, const std::vector<StoichPack::sp_scalar>& X){
	const std::vector<StoichPack::sp_scalar>& data = GetDefaultValuesImpl(name,S.Participants(),X);
	//convert std::vector to EXT::VectorType
	typename EXT::VectorType result = EXT::CreateVector(data.size());
	std::copy(data.begin(),data.end(),EXT::Begin(result));
	return result;
}

//Get default initial values for the system with the specified name that is stored in S for ODEProblem.
template<typename EXT>
typename EXT::VectorType GetDefaultValues(const std::string& name, const StoichPack::IKineticContainer<EXT>& S){
	return GetDefaultValues<EXT>(name,S,std::vector<StoichPack::sp_scalar>(0));
}

//Get default initial values for the system with the specified name that is stored in S for PDEProblem on mesh.
template<typename EXT>
FDMVector GetDefaultValues(const std::string& name, const StoichPack::IKineticContainer<EXT>& S, const FDMMesh& mesh){
	FDMVector u(S.GlobalSpecies(),mesh);
	//iterate over all nodes and call GetDefaultValues for the respective position
	for(int i=0;i<mesh.Nodes();++i) u.SetSub(i,GetDefaultValues<EXT>(name,S,mesh.Coordinates(i)));
	return u;
}

#endif
