#include "MyStoichiometries.h"
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

	ODENewtonProblem(const T& stoichiometry, EXT::VectorArrayType& a, const EXT::VectorType& old, size_t s)
	                 : S(stoichiometry), all(a), cold(old), stage(s) {}
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

EXT::VectorType GetInitial(const BasicKineticStoichiometry<EXT>& S){
	if(ODEConstants::example=="Example1" || ODEConstants::example=="Example2"){
		EXT::VectorType c(3);
		c<<1,0,0;
		return c;
	} else if(ODEConstants::example=="Example3"){
		const std::vector<Species>& x = S.Participants();
		EXT::VectorType c(x.size());
		for(size_t i=0;i<x.size();++i) {
			std::string n=x[i].Name().substr(0,1);
			if(n=="A") c[i]=1;
			else if(n=="C") c[i]=0.5;
			else c[i]=0;
		}
		return c;
	} else throw "Unknown example";
}
		
int main(){
	ODEConstants::print();
	BasicKineticStoichiometry<EXT> S = LoadStoichiometry(ODEConstants::example);
	EXT::VectorType c =GetInitial(S);
	cout<<endl<<"c(0)="<<endl<<c<<endl<<endl;
	if(!ODEConstants::virtual_interface){
		if(ODEConstants::preprocessing==PreprocessingType::ONESIDED){
			OneSidedStoichiometry<EXT,BasicKineticStoichiometry<EXT> > tmp(S);
			c=SolveProblem(tmp,c);
		} else c=SolveProblem(S,c);
	} else {
		IKineticStoichiometry<EXT>* tmp = 0;
		if(ODEConstants::preprocessing==PreprocessingType::ONESIDED) tmp=new OneSidedStoichiometry<EXT>(S);
		else tmp=S.copy();
		c=SolveProblem(*tmp,c);
		delete tmp;
	}

	cout<<endl<<"c("<<ODEConstants::T<<") = "<<endl<<c<<endl;
	ofstream os(ODEConstants::FileName().c_str());
	assert(os.is_open());
	os<<c;
	os.close();
	return 0;
}
