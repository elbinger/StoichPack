/*
File: PDEExemple.cpp
Author: Tobias Elbinger <elbinger@math.fau.de>
Purpose:
 1) Demonstrate the "virtual"-interface of StoichPack
 2) Demonstrate how to solve PDEs with StoichPack
*/

#include "MyBiochemicalSystems.h"
#include "PDEConstants.h"
#include "../Utility/Newton.h"
#include "../Utility/FDMSparseMatrix.h"
#include "../Utility/FDMVector.h"
#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;

class FDMVectorArrayPair{
private:
	vector<FDMVector> mobile, immobile;
	FDMVectorArrayPair();
public:
	FDMVectorArrayPair(const IKineticContainer<EXT>& S, const FDMVector& c0, const FDMMesh& m){
		mobile.reserve(S.Stages());
		immobile.reserve(S.Stages());
		for(size_t stage=0;stage<S.Stages();++stage){
			mobile.push_back(FDMVector(S.MobileSpecies(stage),m));
			immobile.push_back(FDMVector(S.ImmobileSpecies(stage),m));
		}

		const size_t mobilespec=S.MobileSpecies();
		for(int i=0;i<m.Nodes();++i){
			SetSub(i,S.FromOriginal(Divide<EXT>(c0.GetSub(i),mobilespec)));
		}
	}

	void SetSub(int i, const EXT::VectorArrayPairType& all){
		for(size_t stage=0;stage<mobile.size();++stage) {
			mobile[stage].SetSub(i,all.Mobile()[stage]);
			immobile[stage].SetSub(i,all.Immobile()[stage]);
		}
	}

	EXT::VectorArrayPairType GetSub(int i) const {
		EXT::VectorArrayPairType result(EXT::ReserveVectorArray(mobile.size()),EXT::ReserveVectorArray(mobile.size()));
		for(size_t stage=0;stage<mobile.size();++stage){
			EXT::PushBack(result.Mobile(),mobile[stage].GetSub(i));
			EXT::PushBack(result.Immobile(),immobile[stage].GetSub(i));
		}
		return result;
	}

	vector<FDMVector>& Mobile() { return mobile; }
	vector<FDMVector>& Immobile() { return immobile; }

	const vector<FDMVector>& Mobile() const { return mobile; }
	const vector<FDMVector>& Immobile() const { return immobile; }
};

class PDENewtonSubProblem{
public:
	const IKineticContainer<EXT>& S;
	EXT::VectorArrayPairType& all;
	const EXT::VectorType& cold;
	size_t stage;

	void Residual(EXT::VectorType& c, EXT::VectorType& R){
		//S.ApplyImmobileCorrection(all.Immobile(),stage);
		R=S.ImmobileSpeciesRates(all,stage)-(c-cold)/PDEConstants::dt;
	}
	void Jacobi(const EXT::VectorType& c, EXT::MatrixType& J){
		J=S.DiffImmobileSpeciesRates(all,stage).Immobile()-Diag<EXT>(vector<sp_scalar>(S.ImmobileSpecies(stage),1./PDEConstants::dt));
	}
	EXT::VectorType Solve(const EXT::MatrixType& J, const EXT::VectorType& R) const {
		return J.fullPivLu().solve(R);
	}
	PDENewtonSubProblem(const IKineticContainer<EXT>& container, EXT::VectorArrayPairType& a, const EXT::VectorType& old, size_t s)
	                    : S(container), all(a), cold(old), stage(s) {}
};

class PDENewtonProblem{
public:
	const IKineticContainer<EXT>& S;
	FDMVectorArrayPair& all;
	const FDMVectorArrayPair& allold;
	size_t stage;
	const FDMMesh& mesh;

	EXT::VectorType Stencil(int i, const EXT::VectorArrayPairType& data){
		const sp_scalar tmp = mesh.h*mesh.h;
		EXT::VectorType result = pow(2.,mesh.dimension)/tmp*data.Mobile()[stage];
		const vector<int> neighbours = mesh.Neighbours(i);
		for(int n : neighbours) result-=all.GetSub(n).Mobile()[stage]/tmp;
		return result;
	}

	void Residual(FDMVector& c, FDMVector& R){
		for(int i=0;i<mesh.Nodes();++i){
			EXT::VectorArrayPairType local = all.GetSub(i);
			//S.ApplyMobileCorrection(local.Mobile(),stage);
			EXT::VectorArrayPairType localold = allold.GetSub(i);

			if(S.ImmobileSpecies(stage)!=0){
				PDENewtonSubProblem sub(S,local,localold.Immobile()[stage],stage);
				EXT::MatrixType J(S.ImmobileSpecies(stage),S.ImmobileSpecies(stage));

				NewtonReturn ret = Newton(local.Immobile()[stage],J,sub,PDEConstants::epsilon_immobile,
				                          PDEConstants::max_iterations);
				if(!ret.converged) { cout<<"SUB NOT CONVERGED"<<endl; throw "SUB NOT CONVERGED"; }
			}
			EXT::VectorType r = Stencil(i,local)+(local.Mobile()[stage]-localold.Mobile()[stage])/PDEConstants::dt -
			                    S.MobileSpeciesRates(local,stage);
			R.SetSub(i,r);
			all.SetSub(i,local);
		}
	}

	EXT::MatrixType SolveSmall(const EXT::MatrixType& A, const EXT::MatrixType& b) const {
		if(A.rows()==0){
			assert(A.cols()==0 && b.rows()==0);
			return EXT::MatrixType(0,b.cols());
		} else return A.fullPivLu().solve(b);
	}

	void Jacobi(const FDMVector& c, FDMSparseMatrix& J){
		for(int i=0;i<mesh.Nodes();++i){
			EXT::VectorArrayPairType local = all.GetSub(i);
			const EXT::MatrixQuadType diff = S.DiffSpeciesRates(local,stage);
			const EXT::MatrixType& J11 = diff.Mobile().Mobile();
			const EXT::MatrixType& J12 = diff.Mobile().Immobile();
			const EXT::MatrixType& J21 = diff.Immobile().Mobile();
			const EXT::MatrixType& J22 = diff.Immobile().Immobile();
			const sp_scalar diagentry = pow(2.,mesh.dimension)/(mesh.h*mesh.h)+1./PDEConstants::dt;

			const EXT::MatrixType tmp=Diag<EXT>(vector<sp_scalar>(S.ImmobileSpecies(stage),1./PDEConstants::dt))-J22;
			const EXT::MatrixType afterdiff = J11+J12*SolveSmall(tmp,J21);
			const EXT::MatrixType Jlocal = Diag<EXT>(std::vector<sp_scalar>(S.MobileSpecies(stage),diagentry))-afterdiff;
			J.SetDiagonalEntry(i,Jlocal);
		}
	}
	EXT::VectorType Solve(const FDMSparseMatrix& J, const FDMVector& R) { return J.Solve(R); }

	PDENewtonProblem(const IKineticContainer<EXT>& container, FDMVectorArrayPair& a, const FDMVectorArrayPair& old, size_t s, const FDMMesh& m)
	                 : S(container), all(a), allold(old), stage(s), mesh(m) {}
};

FDMVector FDMToOriginal(const IKineticContainer<EXT>& S, const FDMVectorArrayPair& c, const FDMMesh& mesh){
	FDMVector result(S.AllSpecies(),mesh);
	for(int i=0;i<mesh.Nodes();++i)
		result.SetSub(i,Combine<EXT>(S.ToOriginal(c.GetSub(i))));

	return result;
}

FDMVector SolveProblem(const IKineticContainer<EXT>& S, const FDMVector& c0, const FDMMesh& mesh){
	FDMVectorArrayPair c(S,c0,mesh);
	FDMVectorArrayPair cold = c;

	vector<FDMSparseMatrix> sparse;
	sparse.reserve(S.Stages());
	for(size_t i=0;i<S.Stages();++i) sparse.push_back(FDMSparseMatrix(S.MobileSpecies(i),mesh,PDEConstants::linear_epsilon,PDEConstants::linear_max_iterations));

	time_t t = time(0);
	for(double t = 0;t<PDEConstants::T;){
		t+=PDEConstants::dt;
		cout<<"| From "<<t-PDEConstants::dt<<" to "<<t<<":"<<endl;
		for(size_t stage=0;stage<S.Stages();++stage){
			cout<<"| | Stage "<<stage<<":"<<endl;
			PDENewtonProblem problem(S,c,cold,stage,mesh);
			NewtonReturn result=Newton(c.Mobile()[stage],sparse[stage],problem,
			                           PDEConstants::epsilon_mobile,PDEConstants::max_iterations);
			cout<<"| | r("<<result.iterations<<") = "<<result.residual<<" --> ";
			if(!result.converged){
				cout<<"NOT CONVERGED!"<<endl;
				throw "Could not solve subproblem!";
			}
			cout<<"CONVERGED!"<<endl;
		}
		cold=c;
	}
	cout<<"needed "<<time(0)-t<<" seconds!"<<endl;
	return FDMToOriginal(S,c,mesh);
}

typedef shared_ptr<IKineticContainer<EXT> > SharedContainer;

void LoadPreprocessing(SharedContainer& S){
	for(auto x : PDEConstants::preprocessing) {
		if(x==ONESIDED) SharedContainer(new OneSidedCouplings<EXT>(*S)).swap(S);
		if(x==REDUCED) SharedContainer(new ReducedKineticContainer<EXT>(*S)).swap(S);
	}
}

int main(){
	PDEConstants::print(); //print configuration

	FDMMesh mesh(PDEConstants::nodes_per_dimension,PDEConstants::dimension);

	const BiochemicalSystem<> system = LoadSystem(PDEConstants::system); //load specified system
	SharedContainer S(new BasicKineticContainer<EXT>(system));

	FDMVector c = GetDefaultValues<EXT>(PDEConstants::system,*S,mesh); //load intial values

	LoadPreprocessing(S);
	c=SolveProblem(*S,c,mesh);

	//store results
	ofstream os(PDEConstants::FileName().c_str());
	assert(os.is_open());
	for(auto x : S->Participants()) os<<x.Name()<<" ";
	os<<endl<<c;
	os.close();
	return 0;
}
