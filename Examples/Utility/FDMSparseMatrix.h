#ifndef __H_FDM_SPARSE_MATRIX__
#define __H_FDM_SPARSE_MATRIX__

#include "FDMMesh.h"
#include <IterativeLinearSolvers>

class FDMSparseMatrix{
private:
	int species;
	Eigen::SparseMatrix<double> matrix;
	double epsilon;
	int maxiterations;
	mutable Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,Eigen::IncompleteLUT<double> > solver;
	void InitSolver(double eps, int mit){
		epsilon=eps;
		maxiterations=mit;
		if(epsilon>0) solver.setTolerance(epsilon);
		if(maxiterations>0) solver.setMaxIterations(maxiterations);
		solver.analyzePattern(matrix);
	}
public:
	FDMSparseMatrix(int s, const FDMMesh& m, double epsilon=0, int maxiterations=0) : species(s), matrix(m.Nodes()*s,m.Nodes()*s) {
		assert(species>=0);
		if(species!=0) {
			m.DefaultSparse(matrix,species);
			InitSolver(epsilon,maxiterations);
		}
	}

	FDMSparseMatrix(const FDMSparseMatrix& m) : species(m.species), matrix(m.matrix){
		if(species!=0) InitSolver(m.epsilon,m.maxiterations);
	}

	void SetDiagonalEntry(int i, const Eigen::MatrixXd& m) {
		assert(m.rows()==species && m.cols()==species);
		const int I=species*i;
		for(int j=0;j<species;++j){
			double* entries = &matrix.coeffRef(I,I+j);
			for(int i=0;i<species;++i) *(entries+i)=m(i,j);
		}
	}

	Eigen::VectorXd Solve(const Eigen::VectorXd& b) const {
		if(species==0) return Eigen::VectorXd(0);
		solver.factorize(matrix);
		Eigen::VectorXd ret(solver.solve(b));
		if(solver.info()!=Eigen::Success) throw "Solver did not converge!";
		return ret;
	}

	size_t Species() const { return species; }
	const Eigen::SparseMatrix<double>& Matrix() const { return matrix; }
};

#endif
