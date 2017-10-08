#ifndef __H_FDMMESH__
#define __H_FDMMESH__

#include <vector>
#include <cassert>
#include <Eigen>
#include <utility>

typedef Eigen::Triplet<double> triplet;

class FDMMesh{
private:
	void AddDiagNode(int i, size_t species, std::vector<triplet>& triplets) const;
	void AddOffDiagNode(int i, int j, size_t specis, std::vector<triplet>& triplets) const;
public:
	const int nodes_per_dimension;
	const int dimension;
	const double h;
	
FDMMesh(size_t npd, size_t d);

int Nodes() const;
int SubNodes(int dim) const;

std::pair<bool,bool> BndNode(int i, int dim) const;

std::vector<int> Neighbours(int i) const ;

void DefaultTriplets(std::vector<triplet>& triplets, size_t species) const;

void DefaultSparse(Eigen::SparseMatrix<double>& A, size_t species) const;

std::vector<double> Coordinates(int i) const;
};

template<typename T>
void MeshIterationImpl(const FDMMesh& mesh, T& data, int dim, int start){
	assert(dim>0);
	assert( (start==0 && dim==mesh.dimension) || dim<mesh.dimension );
	const int subnodes = mesh.SubNodes(dim);
	const int npd = mesh.nodes_per_dimension;
	assert(start%subnodes==0);

	// treat boundary conditions
	for(int i=start;i<start+subnodes;++i) data.boundary(dim,0,mesh,i);
	const int start2 = start+(npd-1)*subnodes;
	for(int i=start2;i<start2+subnodes;++i) data.boundary(dim,1,mesh,i);
	// inner nodes
	if(dim==1){
		for(int i=1;i<npd-1;++i) data.inner(mesh,start+i);
	} else {
		for(int i=1;i<npd-1;++i)
			MeshIterationImpl(mesh,data,dim-1,start+i*subnodes);
	}
}

template<typename T>
void MeshIteration(const FDMMesh& mesh, T& data){
	MeshIterationImpl<T>(mesh,data,mesh.dimension,0);
}


#endif
