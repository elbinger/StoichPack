// "config file" for ODE problems

#include "ODEConstants.h"
#include <iostream>
#include <cassert>
#include <sstream>

using namespace std;

const int ODEConstants::max_iterations = 50; //max. number of newton iterations
const sp_scalar ODEConstants::epsilon = 1e-12; //stopping criterion for newton method

const sp_scalar ODEConstants::T = 50; // end of time interval (start=0)
const sp_scalar ODEConstants::dt = 0.1; // time step width

const std::string ODEConstants::system = "System3"; //"SystemX", X = 1 ... 3

/*preprocessing: PreprocessingType::NONE for no preprocessing, PreprocessingType::ONESIDED for exploiting one sided couplings,
 *PreprocessingType::REDUCED for using linear combinations of species.*/
const PreprocessingType ODEConstants::preprocessing = PreprocessingType::NONE;

void ODEConstants::print() {
	cout<<"max_iterations = "<<max_iterations<<endl<<"epsilon = "<<epsilon<<endl<<"T = "<<T<<endl<<"dt = "<<dt<<endl;
	cout<<endl<<"preprocessing = "<<preprocessing<<endl<<"system = "<<system<<endl<<endl;

	assert(max_iterations>0);
	assert(epsilon>0);
	assert(T>0);
	assert(dt>0);
}

std::string ODEConstants::FileName() {
	std::stringstream s;
	s<<"out_"<<system<<"_"<<preprocessing<<"_"<<T;
	return s.str();
}
