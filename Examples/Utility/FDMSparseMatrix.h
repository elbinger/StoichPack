#ifndef __H_FDM_SPARSE_MATRIX__
#define __H_FDM_SPARSE_MATRIX__

#include "FDMMesh.h"
#include <IterativeLinearSolvers>

class FDMSparseMatrix{
private:
	const int species;
	const FDMMesh& mesh;
	Eigen::SparseMatrix<double> matrix;
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,Eigen::IncompleteLUT<double> > solver;
public:
	FDMSparseMatrix(int s, const FDMMesh& m) : species(s), mesh(m), matrix(m.Nodes()*s,m.Nodes()*s) {
		assert(species>=0);
		if(species!=0) {
			mesh.DefaultSparse(matrix,species);
			solver.analyzePattern(matrix);
		}
	}
	void SetDiagonalEntry(int i, const Eigen::MatrixXd& m) {
		assert(m.rows()==species && m.cols()==species);
		assert(i<mesh.Nodes());
		const int I=species*i;
		for(int j=0;j<species;++j){
			double* entries = &matrix.coeffRef(I,I+j);
			for(int i=0;i<species;++i) *(entries+i)=m(i,j);
		}
	}

	Eigen::VectorXd Solve(const Eigen::VectorXd& b) {
		if(species==0) return Eigen::VectorXd(0);
		solver.factorize(matrix);
		Eigen::VectorXd ret(solver.solve(b));
		if(solver.info()!=Eigen::Success) throw "Solver did not converge!";
		return ret;
	}

	size_t Species() const { return species; }
	const FDMMesh& Mesh() const { return mesh; }
	const Eigen::SparseMatrix<double>& Matrix() const { return matrix; }
};

#endif
