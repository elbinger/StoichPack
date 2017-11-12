/*
File: ODEExemple.cpp
Author: Tobias Elbinger <elbinger@math.fau.de>
Purpose:
 1) Demonstrate the "template"-interface of StoichPack
 2) Demonstrate how to solve ODEs with StoichPack
*/

#include "MyBiochemicalSystems.h"
#include "ODEConstants.h"
#include "../Utility/Newton.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace StoichPack;

//Problem class for newton method
template<typename T>
class ODENewtonProblem{
public:
	const T& S;				//the StoichPack container type for the problem
	EXT::VectorArrayType& allstages;	//values of all stages
	const EXT::VectorType& cold;		//value of this stage in old time step
	size_t stage;				//current stage

	//Residual. Backward euler --> R = (c-c_old)/dt - f(c)
	void Residual(EXT::VectorType& c, EXT::VectorType& R){
		R = (c-cold)/ODEConstants::dt-S.SpeciesRates(allstages,stage);
	}

	//the derivative of R
	void Jacobi(const EXT::VectorType& c, EXT::MatrixType& J){
		J = Diag<EXT>(vector<sp_scalar>(c.rows(),1./ODEConstants::dt)) - S.DiffSpeciesRates(allstages,stage);
	}
	EXT::VectorType Solve(const EXT::MatrixType& J, const EXT::VectorType& R) { return J.fullPivLu().solve(R); }

	//container: the StoichPack container describing the (preorocessed) biochemical problem
	//all: the values of all stages.
	//old: the values of stage[s] in the old time step
	//s: the stage to be solved
	ODENewtonProblem(const T& container, EXT::VectorArrayType& all, const EXT::VectorType& old, size_t s)
	                 : S(container), allstages(all), cold(old), stage(s) {}
};

//solve dt c = f(c), where f is defined by StoichPack container
template<typename T>
EXT::VectorType SolveProblem(const T& S, const EXT::VectorType& c0){
	EXT::VectorArrayType c = S.FromOriginalGlobal(c0); //get transformed variables, i.e. get preprcoessed inital values
	EXT::VectorArrayType cold = c;
	// backward euler
	for(double t = 0;t<=ODEConstants::T;){
		t+=ODEConstants::dt;
		cout<<"| From "<<t-ODEConstants::dt<<" to "<<t<<":"<<endl;
		for(size_t stage=0;stage<S.Stages();++stage){
			//solving each stage
			cout<<"| | Stage "<<stage<<":"<<endl;
			const size_t species = S.GlobalSpecies(stage);

			ODENewtonProblem<T> problem(S,c,cold[stage],stage); //problem class for the newton method
			EXT::MatrixType J(species,species); //storage for the jacobian
			NewtonReturn result=Newton(c[stage],J,problem,
			                           ODEConstants::epsilon,ODEConstants::max_iterations); //newton method

			cout<<"| | r("<<result.iterations<<") = "<<result.residual<<" --> ";
			if(!result.converged){
				cout<<"NOT CONVERGED!"<<endl;
				throw "Could not solve subproblem!";
			}
			cout<<"CONVERGED!"<<endl;
		}
		cold=c;
	}
	return S.ToOriginalGlobal(c);
}
		
int main(){
	try{
		ODEConstants::print(); //print configuration
		const BiochemicalSystem<> system = LoadSystem(ODEConstants::system); //load specified system
		BasicKineticContainer<EXT> S(system);

		EXT::VectorType c =GetDefaultValues<EXT>(ODEConstants::system,S); //load intial values

		cout<<"Loaded and initialized the following system:"<<endl<<S.Problem()<<endl; //print specified problem
		cout<<endl<<"c(0)="<<endl<<c<<endl<<endl; //print intital values

		//preprocess and solve
		if(ODEConstants::preprocessing==PreprocessingType::ONESIDED){
			OneSidedCouplings<EXT,BasicKineticContainer<EXT> > tmp(S); //preprocessing
			c=SolveProblem(tmp,c); //solve preprocessed problem, return values in original form
		} else if(ODEConstants::preprocessing==PreprocessingType::REDUCED) {
			ReducedKineticContainer<EXT,BasicKineticContainer<EXT> > tmp(S); //preprocessing
			c=SolveProblem(tmp,c); //solve preprocessed problem, return values in original form
		} else c=SolveProblem(S,c);

		//print and store results
		cout<<endl<<"c("<<ODEConstants::T<<") = "<<endl<<c<<endl;
		ofstream os(ODEConstants::FileName().c_str());
		assert(os.is_open());
		os<<c;
		os.close();
		return 0;
	} catch(exception& ex){
		cout<<"EXCEPTION!"<<endl;
		cout<<"what() = "<<ex.what()<<endl;
		return -1;
	}
}
