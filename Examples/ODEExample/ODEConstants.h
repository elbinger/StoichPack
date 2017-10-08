#ifndef __H_ODECONSTANTS__
#define __H_ODECONSTANTS__

#include <string>
#include "StoichPackDefines.h"

using namespace StoichPack;

enum PreprocessingType { NONE , ONESIDED , REDUCED };
class ODEConstants{
private:
	
public:
	static const int max_iterations;
	static const sp_scalar epsilon;

	static const sp_scalar T;
	static const sp_scalar dt;

	static const std::string system;

	static const PreprocessingType preprocessing;

	static void print() ;
	static std::string FileName();
};

#endif
