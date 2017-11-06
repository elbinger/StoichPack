/* File: FDMMesh.h
 * Author: Tobias Elbinger (elbinger@math.fau.de)
 * Purpose: Provide a mesh class for solving dt c - div(nabla(c)) = f(c) with hom. Neumann boundary in 1,2,3 dimensions with finite differences.
 */

#ifndef __H_FDMMESH__
#define __H_FDMMESH__

#include <vector>
#include <cassert>
#include <Eigen>
#include <utility>

typedef Eigen::Triplet<double> triplet;

class FDMMesh{
private:
	//add an "entry" on the "diagonal", i.e. add a square matrix to the "blockdiagonal" part
	void AddDiagNode(int i, size_t species, std::vector<triplet>& triplets) const;
	//add "offdiagonal" "entry"
	void AddOffDiagNode(int i, int j, size_t specis, std::vector<triplet>& triplets) const;
public:
	//We assume a unit cube/square/interval. We assume ah_x=h_y=h_z.
	const int nodes_per_dimension;
	const int dimension;
	const double h;
	
FDMMesh(size_t nodesperdimension, size_t dim);

int Nodes() const; //number of unknowns for 1 species
int SubNodes(int dim) const; //distance to a neighbour with respect to the dimension

std::pair<bool,bool> BndNode(int i, int dim) const; //node on boundary?
//1st entry: on boundary with entry 0, 2nd entry: on boundary with entry 1

std::vector<int> Neighbours(int i) const ; //Indeces of all neighbours. We are using ghost nodes, nodes might appear twice!

void DefaultTriplets(std::vector<triplet>& triplets, size_t species) const; //Get eigen triplets for creating a sparse matrix

void DefaultSparse(Eigen::SparseMatrix<double>& A, size_t species) const; //Create a sparse matrix

std::vector<double> Coordinates(int i) const; //Coordinates of node i
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
