// "config file" for ODE problems

#include "PDEConstants.h"
#include <iostream>
#include <cassert>
#include <sstream>

using namespace std;

const int PDEConstants::max_iterations = 50; //max. number of newton iterations
const sp_scalar PDEConstants::epsilon_mobile = 1e-8; //stopping criterion for newton method
const sp_scalar PDEConstants::epsilon_immobile = PDEConstants::epsilon_mobile*1e-4; //stopping criterion for newton method

const int PDEConstants::linear_max_iterations = 0; //max. number of newton iterations
const sp_scalar PDEConstants::linear_epsilon = 0; //stopping criterion for newton method

const sp_scalar PDEConstants::T = 1; // end of time interval (start=0)
const sp_scalar PDEConstants::dt = 0.1; // time step width

const int PDEConstants::dimension = 3;
const int PDEConstants::nodes_per_dimension = 20;

const std::string PDEConstants::system = "System3"; //"SystemX", X = 1 ... 3

/*preprocessing: PreprocessingType::NONE for no preprocessing, PreprocessingType::ONESIDED for exploiting one sided couplings,
 *PreprocessingType::REDUCED for using linear combinations of species.*/
const std::vector<PreprocessingType> PDEConstants::preprocessing = { PreprocessingType::ONESIDED };

void PDEConstants::print() {
	cout<<"max_iterations = "<<max_iterations<<endl;
	cout<<"epsilon_mobile = "<<epsilon_mobile<<endl<<"epsilon_immobile="<<epsilon_immobile<<endl;
	cout<<"T = "<<T<<endl<<"dt = "<<dt<<endl;
	cout<<"dimension = "<<dimension<<endl<<"nodes_per_dimension = "<<nodes_per_dimension<<endl;
	cout<<endl<<"preprocessing = ";
	for(auto p : preprocessing) cout<<p<<" ";
	cout<<endl<<"system = "<<system<<endl<<endl;

	assert(max_iterations>0);
	assert(epsilon_mobile>=epsilon_immobile);
	assert(epsilon_immobile>0);
	if(epsilon_immobile>epsilon_mobile*1e-2) {
		cout<<"Choosing epsilon_immobile > 0.01 * epsilon_mobile may be dangerous!"<<endl;
		cout<<"Small epsilon_immobile will have only small impact on the performance!"<<endl;
		cout<<"There is no reason for choosing huge epsilon_immobile!"<<endl;
	}
	assert(T>0);
	assert(dt>0);
	assert(dimension>=1 && dimension<=3);
	assert(nodes_per_dimension>=2);
}

std::string PDEConstants::FileName() {
	std::stringstream s;
	s<<"out_"<<system<<"_";
	for(auto p : preprocessing) s<<p;
	s<<"_"<<T<<"_"<<dimension<<"_"<<nodes_per_dimension;
	return s.str();
}
