/* File: Newton.h
 * Author: Tobias Elbinger (elbinger@math.fau.de)
 * Purpose: Provide a newton method.
 */
#ifndef __H_NEWTON__
#define __H_NEWTON__

#include "../StoichPackEigen.h"

class NewtonReturn{
 public:
 bool converged; //newton method converged?
 size_t iterations; //executed iterations
 double residual; //residual after method
 NewtonReturn(bool conv, size_t it, double r) : converged(conv), iterations(it), residual(r) {}
};

//ProblemType is a class that defines the function Residual and Jacobi,
//i.e. if we want to solve f(c)=0, Residual returns f(c) and Jacobi returns f'(c)
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
	return NewtonReturn(r==r && r<epsilon,it,r); //r==r returns false <--> r=NaN
}
#endif
