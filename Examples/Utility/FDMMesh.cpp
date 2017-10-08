#include "FDMMesh.h"
#include <cmath>

void FDMMesh::AddDiagNode(int i, size_t species, std::vector<triplet>& triplets) const{
	const double coeff = pow(2.,dimension)/(h*h);
	for(size_t s1=0;s1<species;++s1){
		for(size_t s2=0;s2<species;++s2){
			const double tmp = (s1==s2) ? coeff : 0.;
			triplets.push_back(triplet(i*species+s1,i*species+s2,tmp));
		}
	}
}
void FDMMesh::AddOffDiagNode(int i, int j, size_t species, std::vector<triplet>& triplets) const{
	const double coeff = -1./(h*h);
	for(size_t s=0;s<species;++s) triplets.push_back(triplet(i*species+s,j*species+s,coeff));
} 
				
FDMMesh::FDMMesh(size_t npd, size_t d)
		: nodes_per_dimension(npd), dimension(d), h(1./(npd-1)){
	assert(nodes_per_dimension>1);
	assert(dimension>=1 && dimension<=3);
}

int FDMMesh::Nodes() const {
	return SubNodes(dimension)*nodes_per_dimension;
}

int FDMMesh::SubNodes(int dim) const {
	assert(dim>0 && dim<=3);
	if(dim==3) return nodes_per_dimension*nodes_per_dimension;
	if(dim==2) return nodes_per_dimension;
	return 1;
}

std::pair<bool,bool> FDMMesh::BndNode(int i, int dim) const{
	if(dim>dimension) return std::pair<bool,bool>(false,false);
	const int tmp = (i-(i%SubNodes(dim)))/SubNodes(dim);
	const int remainder = tmp%nodes_per_dimension;
	return std::pair<bool,bool>(remainder == 0 , remainder == nodes_per_dimension-1);
}

std::vector<int> FDMMesh::Neighbours(int i) const {
	std::vector<int> ret;
	ret.reserve(2*dimension);
	for(int d=1;d<=dimension;++d){
		const std::pair<bool,bool> bnd = BndNode(i,d);
		int left = i-SubNodes(d);
		int right = i+SubNodes(d);
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
	for(int i=0;i<n;++i){
		AddDiagNode(i,species,triplets);
		const std::vector<int> neighbours= Neighbours(i);
		for(size_t j=0;j<neighbours.size();++j) AddOffDiagNode(i,neighbours[j],species,triplets);
	}
}

void FDMMesh::DefaultSparse(Eigen::SparseMatrix<double>& A, size_t species) const{
	std::vector<triplet> tmp;
	DefaultTriplets(tmp,species);
	A.setFromTriplets(tmp.begin(), tmp.end());
}

std::vector<double> FDMMesh::Coordinates(int i) const {
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
	
