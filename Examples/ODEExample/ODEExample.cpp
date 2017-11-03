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

template<typename T>
class ODENewtonProblem{
public:
	const T& S;
	EXT::VectorArrayType& all;
	const EXT::VectorType& cold;
	size_t stage;
	void Residual(EXT::VectorType& c, EXT::VectorType& R){
		S.ApplyCorrection(all,stage);
		R = (c-cold)/ODEConstants::dt-S.SpeciesRates(all,stage);
	}
	void Jacobi(const EXT::VectorType& c, EXT::MatrixType& J){
		J = - S.DiffSpeciesRates(all,stage);
		for(int j=0;j<c.rows();++j){
			J(j,j)+=1./ODEConstants::dt;
		}
	}
	EXT::VectorType Solve(const EXT::MatrixType& J, const EXT::VectorType& R) { return J.fullPivLu().solve(R); }

	ODENewtonProblem(const T& container, EXT::VectorArrayType& a, const EXT::VectorType& old, size_t s)
	                 : S(container), all(a), cold(old), stage(s) {}
};

template<typename T>
EXT::VectorType SolveProblem(const T& S, const EXT::VectorType& c0){
	EXT::VectorArrayType c = S.FromOriginal(c0);
	EXT::VectorArrayType cold = c;
	for(double t = 0;t<=ODEConstants::T;){
		t+=ODEConstants::dt;
		cout<<"| From "<<t-ODEConstants::dt<<" to "<<t<<":"<<endl;
		for(size_t stage=0;stage<S.Stages();++stage){
			cout<<"| | Stage "<<stage<<":"<<endl;
			const size_t species = S.AllSpecies(stage);
			ODENewtonProblem<T> problem(S,c,cold[stage],stage);
			EXT::MatrixType J(species,species);
			NewtonReturn result=Newton(c[stage],J,problem,
			                           ODEConstants::epsilon,ODEConstants::max_iterations);
			cout<<"| | r("<<result.iterations<<") = "<<result.residual<<" --> ";
			if(!result.converged){
				cout<<"NOT CONVERGED!"<<endl;
				throw "Could not solve subproblem!";
			}
			cout<<"CONVERGED!"<<endl;
		}
		cold=c;
	}
	return S.ToOriginal(c);
}
		
int main(){
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
}
