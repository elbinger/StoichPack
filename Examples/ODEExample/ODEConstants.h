#ifndef __H_ODECONSTANTS__
#define __H_ODECONSTANTS__

#include <string>
#include "StoichPackDefines.h"

using namespace StoichPack;

class ODEConstants{
private:
	
public:
	static const int max_iterations;
	static const sp_scalar epsilon;

	static const sp_scalar T;
	static const sp_scalar dt;

	static const std::string example;
	static const std::string preprocessing;

	static void print() ;
};

#endif
