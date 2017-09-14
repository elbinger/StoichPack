#include "MyStoichiometries.h"
#include "ODEConstants.h"
#include <iostream>
#include <fstream>

using namespace std;

template<typename T>
bool SolveStage(const T& S, const EXT::VectorArrayType& cold_all, EXT::VectorArrayType& c_all, size_t stage){
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

template<typename T>
EXT::VectorType SolveProblem(const T& S, const EXT::VectorType& c0){
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
	BasicKineticStoichiometry<EXT>* S = LoadStoichiometry(ODEConstants::example);
	EXT::VectorType c =GetInitial(*S);
	cout<<endl<<"c(0)="<<endl<<c<<endl<<endl;
	if(!ODEConstants::virtual_interface){
		if(ODEConstants::preprocessing==PreprocessingType::ONESIDED){
			OneSidedStoichiometry<EXT,BasicKineticStoichiometry<EXT> >* tmp=new OneSidedStoichiometry<EXT,BasicKineticStoichiometry<EXT> >(S);
			c=SolveProblem(*tmp,c);
			delete tmp;
		} else {
			c=SolveProblem(*S,c);
			delete S;
		}
	} else {
		IKineticStoichiometry<EXT>* tmp;
		if(ODEConstants::preprocessing==PreprocessingType::ONESIDED) tmp=new OneSidedStoichiometry<EXT>(S);
		else tmp=S;
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
