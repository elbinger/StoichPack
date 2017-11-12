/* File:    FDMMesh.cpp
 * Author:  Tobias Elbinger (elbinger@math.fau.de)
 * Purpose: Provide a mesh class for solving dt c - div(nabla(c)) = f(c) on unit interval/square/cube
 *          with hom. Neumann boundary in 1,2,3 dimensions with finite differences.
 */

#include "FDMMesh.h"
#include <cmath>

/* Storage Order:
 * Vectors will store the value of the i-th species in the point (x*h,y*h,z*h) at the position 
 * 	i+nspec*(x+nodes_per_dimension*(y+nodes_per_dimension*z)),
 * where nspec is the number of species. The FDMMesh class will use the the same storage order.
 */

void FDMMesh::AddDiagNode(int i, size_t species, std::vector<triplet>& triplets) const{
	//Creates triplets for "diagonal" entries. Values are initialized with 0, since diagonal entries change in every Newton iteraton.
	for(size_t s1=0;s1<species;++s1){
		for(size_t s2=0;s2<species;++s2){
			triplets.push_back(triplet(i*species+s1,i*species+s2,0));
		}
	}
}
void FDMMesh::AddOffDiagNode(int i, int j, size_t species, std::vector<triplet>& triplets) const{
	//Creates triplets for "offdiagonal" entries.
	//Values are initialized with the correct value, since they do not need to be updated in every Newton iteration.

	const double coeff = -1./(h*h);
	for(size_t s=0;s<species;++s) triplets.push_back(triplet(i*species+s,j*species+s,coeff));
} 
				
FDMMesh::FDMMesh(size_t nodesperdimension, size_t dim)
		: nodes_per_dimension(nodesperdimension), dimension(dim), h(1./(nodesperdimension-1)){
	assert(nodes_per_dimension>1);
	assert(dimension>=1 && dimension<=3);
}

int FDMMesh::Nodes() const {
	return SubNodes(dimension)*nodes_per_dimension;
}

int FDMMesh::SubNodes(int dim) const {
	assert(dim>0 && dim<=3);
	if(dim==3) return nodes_per_dimension*nodes_per_dimension; //distance to neighbour with respect to the 3rd dimension
	if(dim==2) return nodes_per_dimension; //distance to neighbour with respect to the 2nd dimension
	return 1; //distance to neighbour with respect to the 1st dimension
}

std::pair<bool,bool> FDMMesh::BndNode(int i, int dim) const{
	if(dim>dimension) return std::pair<bool,bool>(false,false);

	//assume node i is (x_1*h,x_2*h,x_3*h) 
	const int tmp = (i-(i%SubNodes(dim)))/SubNodes(dim);
	// dim = 1: tmp = i = x_1+nodes_per_dimension*(...)
	// dim = 2: tmp = x_2 + nodes_per_dimension * x_3
	// dim = 3: tmp = x_3
	// 0<=x_dim<nodes_per_dimension
	const int remainder = tmp%nodes_per_dimension; //this is the value of x_dim
	return std::pair<bool,bool>(remainder == 0 , remainder == nodes_per_dimension-1);
}

std::vector<int> FDMMesh::Neighbours(int i) const {
	std::vector<int> ret;
	ret.reserve(2*dimension);
	for(int d=1;d<=dimension;++d){
		const std::pair<bool,bool> bnd = BndNode(i,d);
		int left = i-SubNodes(d);
		int right = i+SubNodes(d);

		/* ghost nodes:
		 * If i is on the boundary, i.e. i-SubNodes(d) or i+SubNodes(d) is not inside the mesh, we add the other neighbour
		 * twice. We assume that the values of both "neighbour" nodes are equal (approximation of hom. Neumann boundary
		 * condition with central difference quotient in node i), therefore this implementation will lead to the correct
		 * result.
		 */
		if(bnd.first) left=right;
		else if(bnd.second) right=left;

		ret.push_back(left);
		ret.push_back(right);
	}
	return ret;
}

void FDMMesh::DefaultTriplets(std::vector<triplet>& triplets, size_t species) const {
	assert(triplets.size()==0);
	triplets.reserve(Nodes()*(1+2*dimension)*species);

	const int n = Nodes();
	for(int i=0;i<n;++i){ //iterate over all nodes
		AddDiagNode(i,species,triplets); //add a diagonal block
		const std::vector<int> neighbours= Neighbours(i);
		for(size_t j=0;j<neighbours.size();++j) AddOffDiagNode(i,neighbours[j],species,triplets); //add off diagonal block
	}
}

void FDMMesh::DefaultSparse(Eigen::SparseMatrix<double>& A, size_t species) const{
	std::vector<triplet> tmp;
	DefaultTriplets(tmp,species);
	A.setFromTriplets(tmp.begin(), tmp.end());
}

std::vector<double> FDMMesh::Coordinates(int i) const {
	//i=x+nodes_per_dimension*(y+nodes_per_dimension*z)
	assert(i>=0 && i<Nodes());
	std::vector<double> result;
	for(int d=0;d<dimension;++d){
		const int tmp = i%nodes_per_dimension;
		result.push_back(tmp*h);
		i=(i-tmp)/nodes_per_dimension;
	}
	assert(i==0);
	return result;
}

