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
using namespace StoichPack;

/* Part 1: FDMVectorArrayPair
 * Since StoichPack containers accept arrays of vectors containing concentration values, the vectors representing discrete solutions
 * must also be split into arrays. Since all function evaluations are performed node per node, we have to define functions
 * for reading and writting local arrays, i.e. read/store only the values of one node inside the vector representing the solution
 * on the whole domain (see ../Utility/FDMVector.h for storage orders). Since we will be using implicit elimination for the immobile species,
 * the values for mobile and immobile species have to be stored in different arrays
*/

// A class storing 2 arrays of FDMVectors.
class FDMVectorArrayPair{
private:
	vector<FDMVector> mobile, immobile;
	//mobile: the values of mobile species of each stage
	//immobile: the values of immobile species of each stage

	FDMVectorArrayPair(); //FORBID
public:
	//initialize
	//S: a (preprocessed) StoichPack container
	//c0: initial solution (not preprocessed)
	//m: the mesh that is used
	FDMVectorArrayPair(const IKineticContainer<EXT>& S, const FDMVector& c0, const FDMMesh& m){
		mobile.reserve(S.Stages());
		immobile.reserve(S.Stages());

		//iterate over all stages and append FDMVectors of according size
		for(size_t stage=0;stage<S.Stages();++stage){
			mobile.push_back(FDMVector(S.MobileSpecies(stage),m));
			immobile.push_back(FDMVector(S.ImmobileSpecies(stage),m));
		}

		//preprocess values of c0 and initialize values
		const size_t mobilespec=S.MobileSpecies();
		for(int i=0;i<m.Nodes();++i){
			//Divide: split into mobile and immobile part
			//FromOriginal: preprocessing
			//SetSub: initialize
			SetSub(i,S.FromOriginal(Divide<EXT>(c0.GetSub(i),mobilespec)));
		}
	}

	//store all values belonging to node i
	void SetSub(int i, const VectorArrayPair<EXT>& all){
		//iterate over all stages and store
		for(size_t stage=0;stage<mobile.size();++stage) {
			mobile[stage].SetSub(i,all.Mobile()[stage]);
			immobile[stage].SetSub(i,all.Immobile()[stage]);
		}
	}

	//get all values belonging to node i
	VectorArrayPair<EXT> GetSub(int i) const {
		VectorArrayPair<EXT> result(EXT::ReserveVectorArray(mobile.size()),EXT::ReserveVectorArray(mobile.size()));
		//iterate over all stages and read
		for(size_t stage=0;stage<mobile.size();++stage){
			EXT::PushBack(result.Mobile(),mobile[stage].GetSub(i));
			EXT::PushBack(result.Immobile(),immobile[stage].GetSub(i));
		}
		return result;
	}

	//get mobile parts (all nodes)
	vector<FDMVector>& Mobile() { return mobile; }
	//get immobile parts (all nodes)
	vector<FDMVector>& Immobile() { return immobile; }

	//get mobile parts (const, all nodes)
	const vector<FDMVector>& Mobile() const { return mobile; }
	//get immobile parts (const, all nodes)
	const vector<FDMVector>& Immobile() const { return immobile; }

	//get the unpreprocessed FDMVector for a FDMVectorArrayPair
	FDMVector ToOriginal(const IKineticContainer<EXT>& S, const FDMMesh& mesh){
		FDMVector result(S.GlobalSpecies(),mesh);
		for(int i=0;i<mesh.Nodes();++i)
			result.SetSub(i,Combine<EXT>(S.ToOriginal(GetSub(i))));

		return result;
	}

};

/* Part 2: Newton problem classes
 * We use implicit elimination of immobile species, i.e. we consider the concentrations of the immobile species in time step k c_immob(k)
 * as a function of the mobile species in time step k c_mob(k). This function is defined implicitly via the equation
 * 		(c_immob(k,i)-cimmob(k-1))/dt = f_immob(c_mob(k),c_immob(k,i)) .		(*)
 * For given iteration value c_mob(k,i) of the newton iteration for the mobile species, a newton iteration will have the following form:
 * i) Iterate over all nodes and calculate c_immob(k) using (*) and assuming c_mob(k)=c_mob(k,i) (again with newton method)
 * ii) Calculate residual for c_mob(k,i)
 * iii) Calculate jacobian for c_mob(k,i) using the chain rule.
 * iv) Calculate newton update
 *
 * This is a nested newton iteration, therfore we need to define 2 newton problem classes:
 * i) PDENewtonSubProblem for solving (*)
 * ii) PDENewtonProblem for solving the mobile problem ( dt c_mob - laplace c_mob = f_mob(c_mob, c_immob(c_mob))
 */

//problem class for immobile species
class PDENewtonSubProblem{
public:
	const IKineticContainer<EXT>& S;	//the (preprocessed) StoichPack container
	VectorArrayPair<EXT>& allstages;	//the values of all stages, only for the current node
	const EXT::VectorType& cold;		//values of the old time step, only for the current node/stage, only immobile species
	size_t stage;				//the stage

	void Residual(EXT::VectorType& c, EXT::VectorType& R){
		//R = f_immob(c_mob,c_immob) - (c_immob-c_immob_old)/dt
		R=S.ImmobileSpeciesRates(allstages,stage)-(c-cold)/PDEConstants::dt;
	}
	void Jacobi(const EXT::VectorType& c, EXT::MatrixType& J){
		//derivative of Residual
		J=S.DiffImmobileSpeciesRates(allstages,stage).Immobile()
		  -Diag<EXT>(vector<sp_scalar>(S.ImmobileSpecies(stage),1./PDEConstants::dt));
	}
	EXT::VectorType Solve(const EXT::MatrixType& J, const EXT::VectorType& R) const {
		//solve J*x=R using Gauss with pivotization
		return J.fullPivLu().solve(R);
	}

	//initialize
	//container: the (preprocessed) StoichPack container
	//all: the values of all stages (only current node)
	//old: the values of the last time step (stage s only)
	//s: the stage
	PDENewtonSubProblem(const IKineticContainer<EXT>& container, VectorArrayPair<EXT>& all, const EXT::VectorType& old, size_t s)
	                    : S(container), allstages(all), cold(old), stage(s) {}
};

//problem class for mobile species
class PDENewtonProblem{
public:
	const IKineticContainer<EXT>& S;	//the (preprocessed) StoichPack container
	FDMVectorArrayPair& allstages;		//the values of all stages, all nodes
	const FDMVectorArrayPair& allold;	//values of old time step, all stages/nodes, mobile and immobile species
	size_t stage;				//the stage
	const FDMMesh& mesh;			//the mesh

	//Evaluate -laplace_h c in node i
	EXT::VectorType Stencil(int i, const VectorArrayPair<EXT>& data){
		EXT::VectorType result = pow(2.,mesh.dimension)*data.Mobile()[stage];
		const vector<int> neighbours = mesh.Neighbours(i);
		//iterate over all neighbours
		for(int n : neighbours) result-=allstages.GetSub(n).Mobile()[stage];
		return result/(mesh.h*mesh.h);
	}

	void Residual(FDMVector& c, FDMVector& R){
		//iterate over all nodes
		for(int i=0;i<mesh.Nodes();++i){
			// get values of current iteration value and old time step in node i
			VectorArrayPair<EXT> local = allstages.GetSub(i);
			VectorArrayPair<EXT> localold = allold.GetSub(i);

			if(S.ImmobileSpecies(stage)!=0){
				//we have immobile species in this stage, so we have to solve (*)
				//this done by using Newton for PDENewtonSubProblem
				PDENewtonSubProblem sub(S,local,localold.Immobile()[stage],stage);
				EXT::MatrixType J(S.ImmobileSpecies(stage),S.ImmobileSpecies(stage));

				NewtonReturn ret = Newton(local.Immobile()[stage],J,sub,PDEConstants::epsilon_immobile,
				                          PDEConstants::max_iterations);
				//check for convergence
				if(!ret.converged) { cout<<"SUB NOT CONVERGED"<<endl; throw StoichPackException("SUB NOT CONVERGED"); }
			}
			//values for immobile species are updated now, assemble residual for mobile species
			EXT::VectorType r = Stencil(i,local)+(local.Mobile()[stage]-localold.Mobile()[stage])/PDEConstants::dt
			                    - S.MobileSpeciesRates(local,stage);
			R.SetSub(i,r);
			allstages.SetSub(i,local); //store values for immobile concentrations
		}
	}

	//Solve A*X=b, where X and b are matrices and A is a square matrix. Special treatment if A is empty.
	EXT::MatrixType SolveSmall(const EXT::MatrixType& A, const EXT::MatrixType& b) const {
		assert(A.rows()==A.cols()); //must be square
		if(A.rows()==0){
			//A is empty
			assert(b.rows()==0); //b must be empty
			return EXT::MatrixType(0,b.cols()); //return empty matrix
		} else return A.fullPivLu().solve(b); //Gauss with pivotization
	}

	void Jacobi(const FDMVector& c, FDMSparseMatrix& J){
		//iterate over all nodes
		for(int i=0;i<mesh.Nodes();++i){
			VectorArrayPair<EXT> local = allstages.GetSub(i); //get values in node i, all stages, mobile and immobile species
			const MatrixQuad<EXT> diff = S.DiffSpeciesRates(local,stage); //full jacobian
			const EXT::MatrixType& J11 = diff.Mobile().Mobile();//derivative of mobile species with respect to mobile species
			const EXT::MatrixType& J12 = diff.Mobile().Immobile();//derivative of mobile species with respect to immobile species
			const EXT::MatrixType& J21 = diff.Immobile().Mobile();//...
			const EXT::MatrixType& J22 = diff.Immobile().Immobile();//...
			const sp_scalar diagentry = pow(2.,mesh.dimension)/(mesh.h*mesh.h)+1./PDEConstants::dt;//2^dimension/h^2 + 1/dt

			//application of chain rule
			const EXT::MatrixType tmp=Diag<EXT>(vector<sp_scalar>(S.ImmobileSpecies(stage),1./PDEConstants::dt))-J22;
			//derivative of (*) with respect to c_mob:
			//     dt^(-1) (dc_immob/dc_mob)= (df_immob/dc_mob) + (df_immob/dc_immob) * (dc_immob/dc_cmob)
			//<--> (dt^(-1) - (df_immob/dc_immob) ) * (dc_immob/dc_cmob) = (df_immob/dc_mob)
			//<--> (dt^(-1)- J22) * (dc_immob/dc_mob) = J21
			//<--> tmp * (dc_immob/dc_mob) = J21
	
			// (df_mob(c_mob,c_immob(c_mob))/dc_mob) = d_1 f_mob + d_2 f_mob * (dc_immob/dc_mob) =
			//                                       = J11       + J12 * tmp^(-1)*J21 
			const EXT::MatrixType chain = J11+J12*SolveSmall(tmp,J21);

			const EXT::MatrixType Jlocal = Diag<EXT>(std::vector<sp_scalar>(S.MobileSpecies(stage),diagentry))-chain;
			J.SetDiagonalEntry(i,Jlocal);
		}
	}

	//apply sparse solver
	EXT::VectorType Solve(const FDMSparseMatrix& J, const FDMVector& R) { return J.Solve(R); }

	//initialize
	//container: the (preprocessed) StoichPack container
	//all: the values of all stages (all nodes, all stages, mobile and immobile species)
	//old: the values of the last time step (all nodes, stage s only, mobile and immobile species)
	//s: the stage
	//m: the mesh
	PDENewtonProblem(const IKineticContainer<EXT>& container, FDMVectorArrayPair& all, const FDMVectorArrayPair& old,
	                 size_t s, const FDMMesh& m)
	                 : S(container), allstages(all), allold(old), stage(s), mesh(m) {}
};

//Part 3: the actual solver

FDMVector SolveProblem(const IKineticContainer<EXT>& S, const FDMVector& c0, const FDMMesh& mesh){
	FDMVectorArrayPair c(S,c0,mesh); // get preprocessed values
	FDMVectorArrayPair cold = c;

	//initialize sparse matices
	vector<FDMSparseMatrix> sparse;
	sparse.reserve(S.Stages());
	for(size_t i=0;i<S.Stages();++i)
		sparse.push_back(FDMSparseMatrix(S.MobileSpecies(i),mesh,PDEConstants::linear_epsilon,PDEConstants::linear_max_iterations));

	time_t t = time(0); //start timer

	//backward euler
	for(double t = 0;t<PDEConstants::T;){
		t+=PDEConstants::dt;
		cout<<"| From "<<t-PDEConstants::dt<<" to "<<t<<":"<<endl;
		//iterate over all stages
		for(size_t stage=0;stage<S.Stages();++stage){
			cout<<"| | Stage "<<stage<<":"<<endl;

			//solve problem for mobile species, i.e. apply Newton for PDENewtonProblem
			PDENewtonProblem problem(S,c,cold,stage,mesh);
			NewtonReturn result=Newton(c.Mobile()[stage],sparse[stage],problem,
			                           PDEConstants::epsilon_mobile,PDEConstants::max_iterations);
			cout<<"| | r("<<result.iterations<<") = "<<result.residual<<" --> ";
			if(!result.converged){
				//not good: solver did not converge --> report error and throw
				cout<<"NOT CONVERGED!"<<endl;
				throw StoichPackException("COULD NOT SOLVE TIME STEP!");
			}
			cout<<"CONVERGED!"<<endl;
		}
		cold=c; //update
	}
	cout<<"needed "<<time(0)-t<<" seconds!"<<endl; //stop timer
	return c.ToOriginal(S,mesh); //get original presentation of solution
}

typedef shared_ptr<IKineticContainer<EXT> > SharedContainer;

//apply concatenated preprocessing schemes
void LoadPreprocessing(SharedContainer& S){
	for(auto x : PDEConstants::preprocessing) {
		if(x==ONESIDED) SharedContainer(new OneSidedCouplings<EXT>(*S)).swap(S);
		else if(x==REDUCED) SharedContainer(new ReducedKineticContainer<EXT>(*S)).swap(S);
		else if(x==NONE) { /*do nothing*/ }
		else throw StoichPackException("UNKNOWN PREPROCESSING TYPE!");
	}
}

int main(){
	try{
		PDEConstants::print(); //print configuration

		FDMMesh mesh(PDEConstants::nodes_per_dimension,PDEConstants::dimension); //load mesh

		const BiochemicalSystem<> system = LoadSystem(PDEConstants::system); //load specified system
		SharedContainer S(new BasicKineticContainer<EXT>(system));

		FDMVector c = GetDefaultValues<EXT>(PDEConstants::system,*S,mesh); //load intial values

		LoadPreprocessing(S); //preprocessing
		c=SolveProblem(*S,c,mesh); //solve

		//store results
		ofstream os(PDEConstants::FileName().c_str());
		assert(os.is_open());
		for(auto x : S->Participants()) os<<x.Name()<<" ";
		os<<endl<<c;
		os.close();
		return 0;
	} catch (exception& ex) {
		cout<<"EXCEPTION!"<<endl;
		cout<<"what() = "<<ex.what()<<endl;
		return -1;
	}
}
