/* File: FDMSparseVector.h
 * Author: Tobias Elbinger (elbinger@math.fau.de)
 * Purpose: Provide a vector representing a discrete solution of our problem.
 */

#ifndef __H_FDM_VECTOR__
#define __H_FDM_VECTOR__

#include "FDMMesh.h"

/* Storage Order:
 * Vectors will store the value of the i-th species in the point (x*h,y*h,z*h) at the position 
 * 	i+nspec*(x+nodes_per_dimension*(y+nodes_per_dimension*z)),
 * where nspec is the number of species.
 */

class FDMVector : public Eigen::VectorXd {
private:
	size_t species; //numper of species
public:
	FDMVector(size_t s, const FDMMesh& m) : Eigen::VectorXd(s*m.Nodes()), species(s){}
	size_t Species() const { return species; }

	// Get a vector with the concentrations of all species in node i
	Eigen::VectorXd GetSub(int i) const { 
		return this->block(i*species,0,species,1);
	}

	// Store a vector with the concentrations of all species in node i at the correct position
	void SetSub(int i, const Eigen::VectorXd& x) {
		this->block(i*species,0,species,1)=x;
	}
};

#endif
	
	
