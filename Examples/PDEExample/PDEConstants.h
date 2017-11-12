/* File:    PDEConstants.h
 * Author:  Tobias Elbinger (elbinger@math.fau.de)
 * Purpose: "config file" for PDE problems
 */

#ifndef __H_PDECONSTANTS__
#define __H_PDECONSTANTS__

#include <string>
#include <vector>
#include "StoichPackDefines.h"

enum PreprocessingType { NONE , ONESIDED , REDUCED };
class PDEConstants{
private:
	
public:
	static const int max_iterations; //max. number of itertation in Newton method
	static const StoichPack::sp_scalar epsilon_mobile; //stopping criterion for Newton method (mobile species)
	static const StoichPack::sp_scalar epsilon_immobile; //stopping criterion for Newton method (immobile species)

	static const StoichPack::sp_scalar linear_epsilon; //max. iterations for linear solver (0-->Eigen default)
	static const int linear_max_iterations; //stopping criterion for linear solver (0-->Eigen default)

	static const StoichPack::sp_scalar T; //end of time interval (start=0)
	static const StoichPack::sp_scalar dt; //time step width

	static const int dimension; //dimension (1,2 or 3)
	static const int nodes_per_dimension; //discretization points in each direction

	static const std::string system; //"SystemX", X = 1 ... 3, the name of the biochemical system that shall be solved.

	/*preprocessing: PreprocessingType::NONE for no preprocessing, PreprocessingType::ONESIDED for exploiting one sided couplings,
	 *PreprocessingType::REDUCED for using linear combinations of species.*/
	static const std::vector<PreprocessingType> preprocessing;

	//print configuration
	static void print() ;
	//file name of the output file
	static std::string FileName();
};

#endif
