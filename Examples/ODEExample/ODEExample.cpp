#include "MyStoichiometries.h"
#include "ODEConstants.h"
#include <iostream>

using namespace std;

bool SolveStage(const IKineticStoichiometry<EXT>& S, const EXT::VectorArrayType& cold_all, EXT::VectorArrayType& c_all, size_t stage){
	assert(stage<S.Stages());

	cout<<"| | Stage "<<stage<<endl;

	EXT::VectorType& c = c_all[stage];
	const EXT::VectorType& cold = cold_all[stage];

	S.SpeciesRates(c_all,stage);

	EXT::VectorType R = (c-cold)/ODEConstants::dt-S.SpeciesRates(c_all,stage);
	EXT::MatrixType J(c.rows(),c.rows());
	double r = R.norm();

	for(int i=0;i<ODEConstants::max_iterations && r==r && r>=ODEConstants::epsilon ;++i){
		cout<<"| | | r("<<i<<")= "<<r<<endl;
		J = - S.DiffSpeciesRates(c_all,stage);
		for(int j=0;j<c.rows();++j){
			J(j,j)+=1./ODEConstants::dt;
		}
		c-=J.fullPivLu().solve(R);
		S.ApplyCorrection(c_all,stage);
		R=(c-cold)/ODEConstants::dt-S.SpeciesRates(c_all,stage);
		r=R.norm();
	}

	cout<<"| | | r(end)= "<<r;
	if(r==r && r<ODEConstants::epsilon){
		cout<<" --> CONVERGED"<<endl;
		return true;
	} else {
		cout<<" --> NOT CONVERGED"<<endl;
		return false;
	}
}

EXT::VectorType SolveProblem(const IKineticStoichiometry<EXT>& S, const EXT::VectorType& c0){
	EXT::VectorArrayType c = S.FromOriginal(c0);
	EXT::VectorArrayType cold = c;
	for(double t = 0;t<=ODEConstants::T;){
		t+=ODEConstants::dt;
		cout<<"| From "<<t-ODEConstants::dt<<" to "<<t<<":"<<endl;
		for(size_t i=0;i<S.Stages();++i){
			if(!SolveStage(S,cold,c,i)) throw "Could not solve subproblem!";
		}
		cold=c;
	}
	return S.ToOriginal(c);
}

IKineticStoichiometry<EXT>* Preprocessing(IKineticStoichiometry<EXT>* S){
	if(ODEConstants::preprocessing == "onesided") return new OneSidedStoichiometry<EXT>(S);
	else return S;
}

EXT::VectorType GetInitial(const IKineticStoichiometry<EXT>& S){
	if(ODEConstants::example=="Example1" || ODEConstants::example=="Example2"){
		EXT::VectorType c(3);
		c<<1,0,0;
		return c;
	} else throw "Unknown example";
}
int main(){
	ODEConstants::print();
	IKineticStoichiometry<EXT>* S = LoadStoichiometry(ODEConstants::example);
	S=Preprocessing(S);
	EXT::VectorType c0 =GetInitial(*S);
	cout<<endl<<"c(0)="<<endl<<c0<<endl<<endl;
	EXT::VectorType c1 = SolveProblem(*S,c0);

	cout<<endl<<"c("<<ODEConstants::T<<") = "<<endl<<c1<<endl;

	delete S;

	return 0;
}
