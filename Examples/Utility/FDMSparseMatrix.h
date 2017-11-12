/* File: FDMSparseMatrix.h
 * Author: Tobias Elbinger (elbinger@math.fau.de)
 * Purpose: Provide sparse matrix for finite differences (multi-species context).
 */

#ifndef __H_FDM_SPARSE_MATRIX__
#define __H_FDM_SPARSE_MATRIX__

#include "FDMMesh.h"
#include <IterativeLinearSolvers>

/* Storage Order:
 * Vectors will store the value of the i-th species in the point (x*h,y*h,z*h) at the position 
 * 	i+nspec*(x+nodes_per_dimension*(y+nodes_per_dimension*z)),
 * where nspec is the number of species. The FDMSparseMatrix class will use the the same storage order.
 */

class FDMSparseMatrix{
private:
	int species; //number of species
	Eigen::SparseMatrix<double> matrix;
	double epsilon; //stopping criterion for linear solver, epsilon<=0 --> use Eigen default
	int maxiterations; //max. iterations for linear solver, maxiterations<=0 --> use Eigen default

	mutable Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,Eigen::IncompleteLUT<double> > solver; //ILU + BiCGStab

	void InitSolver(double eps, int mit){
		epsilon=eps;
		maxiterations=mit;
		if(epsilon>0) solver.setTolerance(epsilon); // epsilon<=0 --> use Eigen default
		if(maxiterations>0) solver.setMaxIterations(maxiterations); // maxiterations<=0 --> use Eigen default
		solver.analyzePattern(matrix); // initialize memory for preconditioner
	}
public:
	FDMSparseMatrix(int s, const FDMMesh& m, double epsilon=0, int maxiterations=0) : species(s), matrix(m.Nodes()*s,m.Nodes()*s) {
		assert(species>=0);
		if(species!=0) { // species == 0 --> empty matrix
			m.DefaultSparse(matrix,species); // create sparse matrix
			InitSolver(epsilon,maxiterations);
		}
	}

	FDMSparseMatrix(const FDMSparseMatrix& m) : species(m.species), matrix(m.matrix){
		if(species!=0) InitSolver(m.epsilon,m.maxiterations);
	}

	//set a diagonal entry, i.e. a species x species diagonal block
	void SetDiagonalEntry(int i, const Eigen::MatrixXd& m) {
		assert(m.rows()==species && m.cols()==species);
		const int I=species*i;
		for(int j=0;j<species;++j){
			double* entries = &matrix.coeffRef(I,I+j);
			for(int i=0;i<species;++i) *(entries+i)=m(i,j);
		}
	}
	//the offdiagonal blocks are initialized correctly by FDMMesh

	//solve A*x=b using ILU + BiCGStab
	Eigen::VectorXd Solve(const Eigen::VectorXd& b) const {
		if(species==0) // empty matrix --> return empty vector
			return Eigen::VectorXd(0);

		solver.factorize(matrix); //apply preconditioner
		Eigen::VectorXd ret(solver.solve(b));
		if(solver.info()!=Eigen::Success) // check if solver converged
			throw "Solver did not converge!";

		return ret;
	}

	//getters
	size_t Species() const { return species; }
	const Eigen::SparseMatrix<double>& Matrix() const { return matrix; }
};

#endif
