/* File:    PDEConstants.cpp
 * Author:  Tobias Elbinger (elbinger@math.fau.de)
 * Purpose: "config file" for PDE problems
 */

#include "PDEConstants.h"
#include <iostream>
#include <sstream>

using namespace std;
using namespace StoichPack;

const int PDEConstants::max_iterations = 50; //max. number of newton iterations
const sp_scalar PDEConstants::epsilon_mobile = 1e-8; //stopping criterion for newton method (mobile)
const sp_scalar PDEConstants::epsilon_immobile = PDEConstants::epsilon_mobile*1e-4; //stopping criterion for newton method (immobile)

const int PDEConstants::linear_max_iterations = 0; //max. number of iterations for linear solver (0: use Eigen default)
const sp_scalar PDEConstants::linear_epsilon = 0; //stopping criterion for linear solver (0: use Eigen default)

const sp_scalar PDEConstants::T = 1; // end of time interval (start=0)
const sp_scalar PDEConstants::dt = 0.1; // time step width

const int PDEConstants::dimension = 2; //dimension (1,2 or 3)
const int PDEConstants::nodes_per_dimension = 20; //number of discretization points in each direction

const std::string PDEConstants::system = "System3"; //"SystemX", X = 1 ... 3

/*preprocessing: PreprocessingType::NONE for no preprocessing, PreprocessingType::ONESIDED for exploiting one sided couplings,
 *PreprocessingType::REDUCED for using linear combinations of species.*/
const std::vector<PreprocessingType> PDEConstants::preprocessing = { PreprocessingType::ONESIDED };

//print configuration
void PDEConstants::print() {
	cout<<"max_iterations = "<<max_iterations<<endl;
	cout<<"epsilon_mobile = "<<epsilon_mobile<<endl<<"epsilon_immobile="<<epsilon_immobile<<endl;

	cout<<"linear_max_iterations = "<<linear_max_iterations;
	if(linear_max_iterations <=0 ) cout<<" --> using Eigen default";
	cout<<endl<<"linear_epsilon = "<<linear_epsilon;
	if(linear_epsilon<= 0) cout<<" --> using Eigen default";
	cout<<endl;

	cout<<"T = "<<T<<endl<<"dt = "<<dt<<endl;

	cout<<"dimension = "<<dimension<<endl<<"nodes_per_dimension = "<<nodes_per_dimension<<endl;

	cout<<endl<<"preprocessing = ";
	for(auto p : preprocessing) cout<<p<<" ";

	cout<<endl<<"system = "<<system<<endl<<endl;

	//check
	if(max_iterations<=0) throw StoichPackException("ILLEGAL VALUE FOR max_iterations");
	if(epsilon_mobile<epsilon_immobile) throw StoichPackException("ILLEGAL VALUE FOR epsilon_mobile");
	if(epsilon_immobile<=0) throw StoichPackException("ILLEGAL VALUE FOR epsilon_immobile");

	if(epsilon_immobile>epsilon_mobile*1e-2) {
		cout<<"Choosing epsilon_immobile > 0.01 * epsilon_mobile may be dangerous!"<<endl;
		cout<<"Small epsilon_immobile will have only small impact on the performance!"<<endl;
		cout<<"There is no reason for choosing huge epsilon_immobile!"<<endl;
	}

	if(T<=0) throw StoichPackException("ILLEGAL VALUE FOR T");
	if(dt<=0) throw StoichPackException("ILLEGAL VALUE FOR dt");
	if(dimension<1 || dimension>3) throw StoichPackException("ILLEGAL VALUE FOR dimension");
	if(nodes_per_dimension<2) throw StoichPackException("ILLEGAL VALUE FOR nodes_per_dimension");
}

//file name for output file
std::string PDEConstants::FileName() {
	std::stringstream s;
	s<<"out_"<<system<<"_";
	for(auto p : preprocessing) s<<p;
	s<<"_"<<T<<"_"<<dimension<<"_"<<nodes_per_dimension;
	return s.str();
}
