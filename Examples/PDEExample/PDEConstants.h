#ifndef __H_PDECONSTANTS__
#define __H_PDECONSTANTS__

#include <string>
#include <vector>
#include "StoichPackDefines.h"

using namespace StoichPack;

enum PreprocessingType { NONE , ONESIDED , REDUCED };
class PDEConstants{
private:
	
public:
	static const int max_iterations;
	static const sp_scalar epsilon_mobile;
	static const sp_scalar epsilon_immobile;

	static const sp_scalar linear_epsilon;
	static const int linear_max_iterations;

	static const sp_scalar T;
	static const sp_scalar dt;

	static const int dimension;
	static const int nodes_per_dimension;

	static const std::string system;

	static const std::vector<PreprocessingType> preprocessing;

	static void print() ;
	static std::string FileName();
};

#endif
