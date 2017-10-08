#ifndef __H_FDM_VECTOR__
#define __H_FDM_VECTOR__

#include "FDMMesh.h"

class FDMVector : public Eigen::VectorXd {
private:
	const size_t species;
	const FDMMesh& mesh;
public:
	FDMVector(size_t s, const FDMMesh& m) : Eigen::VectorXd(s*m.Nodes()), species(s), mesh(m){}
	size_t Species() const { return species; }
	const FDMMesh& Mesh() const { return mesh; }

	Eigen::VectorXd GetSub(int i) const {
		return this->block(i*species,0,species,1);
	}

	void SetSub(int i, const Eigen::VectorXd& x) {
		this->block(i*species,0,species,1)=x;
	}
};

#endif
	
	
