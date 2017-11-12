/* File:    ODEConstants.h
 * Author:  Tobias Elbinger (elbinger@math.fau.de)
 * Purpose: "config file" for ODE problems
 */

#ifndef __H_ODECONSTANTS__
#define __H_ODECONSTANTS__

#include <string>
#include "StoichPackDefines.h"

enum PreprocessingType { NONE , ONESIDED , REDUCED };
class ODEConstants{
private:
	
public:
	static const int max_iterations; //max. number of itertation in Newton method
	static const StoichPack::sp_scalar epsilon; //stopping criterion for Newton method

	static const StoichPack::sp_scalar T; //end of time interval (start=0)
	static const StoichPack::sp_scalar dt; //time step width

	static const std::string system; //"SystemX", X = 1 ... 3, the name of the biochemical system that shall be solved.

	/*preprocessing: PreprocessingType::NONE for no preprocessing, PreprocessingType::ONESIDED for exploiting one sided couplings,
	 *PreprocessingType::REDUCED for using linear combinations of species.*/
	static const PreprocessingType preprocessing;

	//print configuration
	static void print() ;
	//file name of the output file
	static std::string FileName();
};

#endif
