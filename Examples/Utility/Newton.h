#ifndef __H_NEWTON__
#define __H_NEWTON__

#include "../StoichPackEigen.h"

class NewtonReturn{
 public:
 bool converged;
 size_t iterations;
 double residual;
 NewtonReturn(bool conv, size_t it, double r) : converged(conv), iterations(it), residual(r) {}
};

template<typename VectorType, typename JacobiType, typename ProblemType>
NewtonReturn Newton(VectorType& c, JacobiType& J, ProblemType& P, double epsilon, size_t maxiterations){
	VectorType R=c;
	P.Residual(c,R);
	double r=R.norm();
	size_t it=0;
	for(;it<maxiterations && r>=epsilon;++it){
		P.Jacobi(c,J);
		c-=P.Solve(J,R);
		P.Residual(c,R);
		r=R.norm();
	}
	return NewtonReturn(r==r && r<epsilon,it,r);
}
#endif
