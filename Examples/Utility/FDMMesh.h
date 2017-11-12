/* File: FDMMesh.h
 * Author: Tobias Elbinger (elbinger@math.fau.de)
 * Purpose: Provide a mesh class for solving dt c - div(nabla(c)) = f(c) on unit interval/square/cube
 *          with hom. Neumann boundary in 1,2,3 dimensions with finite differences.
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
	const int nodes_per_dimension; //number of discretization points in each dimension
	const int dimension; //dimension of Omega, dimension \in {1, 2, 3}
	const double h; //distance of two discretization points
	
FDMMesh(size_t nodesperdimension, size_t dim);

int Nodes() const; //number of unknowns for 1 species
int SubNodes(int dim) const; //distance to a neighbour with respect to the dimension

std::pair<bool,bool> BndNode(int i, int dim) const; //node on boundary?
//1st entry: on boundary with entry 0, 2nd entry: on boundary with entry 1

std::vector<int> Neighbours(int i) const ; //Indeces of all neighbours. We are using ghost nodes, nodes might appear twice!

void DefaultTriplets(std::vector<triplet>& triplets, size_t species) const; //Get eigen triplets for creating a sparse matrix

//Create a sparse matrix with for the mesh and the given number of species
void DefaultSparse(Eigen::SparseMatrix<double>& A, size_t species) const;

std::vector<double> Coordinates(int i) const; //Coordinates of node i
};

/* Storage Order:
 * Vectors will store the value of the i-th species in the point (x*h,y*h,z*h) at the position 
 * 	i+nspec*(x+nodes_per_dimension*(y+nodes_per_dimension*z)),
 * where nspec is the number of species. The FDMMesh class will use the the same storage order.
 */

#endif
