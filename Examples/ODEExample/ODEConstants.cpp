/* File:    ODEConstants.cpp
 * Author:  Tobias Elbinger (elbinger@math.fau.de)
 * Purpose: "config file" for ODE problems
 */

#include "ODEConstants.h"
#include <iostream>
#include <sstream>

using namespace std;
using namespace StoichPack;

const int ODEConstants::max_iterations = 50; //max. number of newton iterations
const sp_scalar ODEConstants::epsilon = 1e-12; //stopping criterion for newton method

const sp_scalar ODEConstants::T = 50; // end of time interval (start=0)
const sp_scalar ODEConstants::dt = 0.1; // time step width

const std::string ODEConstants::system = "System3"; //"SystemX", X = 1 ... 3, the name of the biochemical system that shall be solved.

/*preprocessing: PreprocessingType::NONE for no preprocessing, PreprocessingType::ONESIDED for exploiting one sided couplings,
 *PreprocessingType::REDUCED for using linear combinations of species.*/
const PreprocessingType ODEConstants::preprocessing = PreprocessingType::ONESIDED;

//Print configuration
void ODEConstants::print() {
	cout<<"max_iterations = "<<max_iterations<<endl<<"epsilon = "<<epsilon<<endl<<"T = "<<T<<endl<<"dt = "<<dt<<endl;
	cout<<endl<<"preprocessing = "<<preprocessing<<endl<<"system = "<<system<<endl<<endl;

	//check
	if(max_iterations<=0) throw StoichPackException("ILLEGAL VALUE FOR max_iterations");
	if(epsilon<=0) throw StoichPackException("ILLEGAL VALUE FOR epsilon");
	if(T<=0) throw StoichPackException("ILLEGAL VALUE FOR T");
	if(dt<=0) throw StoichPackException("ILLEGAL VALUE FOR dt");
}

//file name of the output file
std::string ODEConstants::FileName() {
	std::stringstream s;
	s<<"out_"<<system<<"_"<<preprocessing<<"_"<<T;
	return s.str();
}
