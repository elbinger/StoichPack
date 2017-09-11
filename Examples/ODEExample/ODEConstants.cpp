// "config file"

#include "ODEConstants.h"
#include <iostream>
#include <cassert>

using namespace std;

const int ODEConstants::max_iterations = 50; //max. number of newton iterations
const sp_scalar ODEConstants::epsilon = 1e-12; //stopping criterion for newton method

const sp_scalar ODEConstants::T = 10; // end of time interval (start=0)
const sp_scalar ODEConstants::dt = 0.1; // time step width

const std::string ODEConstants::example = "Example2"; //"Example1" or "Example2"
const std::string ODEConstants::preprocessing = "onesided"; //"onesided" for using one sided couplings, anything else for standard treatment

void ODEConstants::print() {
	cout<<"max_iterations = "<<max_iterations<<endl<<"epsilon = "<<epsilon<<endl<<"T = "<<T<<endl<<"dt = "<<dt<<endl;
	cout<<endl<<"example = "<<example<<endl;

	assert(max_iterations>0);
	assert(epsilon>0);
	assert(T>0);
	assert(dt>0);
}

