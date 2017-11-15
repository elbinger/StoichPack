/* File:    StoichPackEigen.h
 * Author:  Tobias Elbinger (elbinger@math.fau.de)
 * Purpose: Make Eigen accesible for StoichPack.
 */

#ifndef __H_STOICHPACK_EIGEN__
#define __H_STOICHPACK_EIGEN__

#include "StoichPack.h"
#include <Eigen>
#include <iterator>

using namespace StoichPack;
using namespace Eigen;

/* Part 1: iterators
 * The data structures in Eigen do not provide iterators, so we write our own iterators.
 * If you use linear algebra packages with iterators, this has not to be done.
 * Nomenclature: rowwise iterator: iterates over all entries in a specified row
 *               columnwise iterator: iterates over all entries in a specified column
 */

//Position of an iterator
class EigenIteratorPos{
public:
	const size_t X; //constant position (row for a rowwise iterator, column for a columnwise iterator)
	size_t Y; //variable position (column for a rowwise iterator, row for a columnwise iterator)

	EigenIteratorPos(size_t x, size_t y) : X(x),Y(y) {}
	
	EigenIteratorPos& operator+=(size_t i) { Y+=i; return *this; }

	EigenIteratorPos operator+(size_t i) const { return EigenIteratorPos(*this)+=i; }

	bool operator==(const EigenIteratorPos& other) const { assert(X==other.X); return Y==other.Y; }
	bool operator!=(const EigenIteratorPos& other) const { return !(*this==other); }

	EigenIteratorPos& operator=(const EigenIteratorPos& other) { assert(X==other.X); Y=other.Y; return *this; }
};

//a const iterator for data structures of type T. rowwsie indicates if it is a rowwise (true) or columnwise (false) iterator
template<typename T, bool rowwise>
class EigenConstIterator : public std::iterator<std::input_iterator_tag,sp_scalar>{
private:
	const T& data;
	EigenIteratorPos pos;
	EigenConstIterator();
public:
	typedef sp_scalar value_type;

	EigenConstIterator(const EigenConstIterator<T,rowwise>& other) : data(other.data), pos(other.pos) {}
	EigenConstIterator<T,rowwise>& operator=(const EigenConstIterator<T,rowwise>& other){
		assert(&data==&other.data);
		pos=other.pos;
		return *this;
	}
	EigenConstIterator(const T& x, const EigenIteratorPos& p) : data(x), pos(p) {}
	sp_scalar operator*() const {
		if(rowwise) return data(pos.X,pos.Y);
		else return data(pos.Y,pos.X);
	}
	const sp_scalar& operator->() const {
		if(rowwise) return data(pos.X,pos.Y);
		else return data(pos.Y,pos.X);
	}

	bool operator==(const EigenConstIterator<T,rowwise>& x) { return pos==x.pos; }
	bool operator!=(const EigenConstIterator<T,rowwise>& x) { return pos!=x.pos; }

	EigenConstIterator<T,rowwise>& operator+=(size_t i) { pos+=i; return *this; }
	EigenConstIterator<T,rowwise>& operator++() { return operator+=(1); }
	EigenConstIterator<T,rowwise> operator++(int) { EigenConstIterator<T,rowwise> ret(*this); operator++(); return ret; }
};

template<typename T, bool rowwise>
EigenConstIterator<T,rowwise> operator+(const EigenConstIterator<T,rowwise>& x, size_t i) { return EigenConstIterator<T,rowwise>(x)+=i; }

//a iterator for data structures of type T. rowwsie indicates if it is a rowwise (true) or columnwise (false) iterator
template<typename T, bool rowwise>
class EigenIterator : public std::iterator<std::output_iterator_tag,sp_scalar>{
private:
	T& data;
	EigenIteratorPos pos;
	EigenIterator();
public:
	typedef sp_scalar value_type;
	EigenIterator(const EigenIterator<T,rowwise>& other) : data(other.data), pos(other.pos) {}
	EigenIterator<T,rowwise>& operator=(const EigenIterator<T,rowwise>& other){
		assert(&data==&other.data);
		pos=other.pos;
		return *this;
	}
	EigenIterator(T& x, const EigenIteratorPos& p) : data(x), pos(p) {}
	sp_scalar& operator*() {
		if(rowwise) return data(pos.X,pos.Y);
		else return data(pos.Y,pos.X);
	}

	sp_scalar operator*() const {
		if(rowwise) return data(pos.X,pos.Y);
		else return data(pos.Y,pos.X);
	}

	sp_scalar& operator->() {
		if(rowwise) return data(pos.X,pos.Y);
		else return data(pos.Y,pos.X);
	}
	const sp_scalar& operator->() const {
		if(rowwise) return data(pos.X,pos.Y);
		else return data(pos.Y,pos.X);
	}

	bool operator==(const EigenIterator<T,rowwise>& x) { return pos==x.pos; }
	bool operator!=(const EigenIterator<T,rowwise>& x) { return pos!=x.pos; }

	EigenIterator<T,rowwise>& operator+=(size_t i) { pos+=i; return *this; }
	EigenIterator<T,rowwise>& operator++() { return operator+=(1); }
	EigenIterator<T,rowwise> operator++(int) { EigenIterator<T,rowwise> ret(*this); operator++(); return ret; }
};

template<typename T, bool rowwise>
EigenIterator<T,rowwise> operator+(const EigenIterator<T,rowwise>& x, size_t i) { return EigenIterator<T,rowwise>(x)+=i; }

/* Part 2: orthogonal decomposition
 * Orthogonal decomposition is used for one preprocessing scheme, so an interface has to be provided.
 * You are free in the implementation of the interface, as long as you assure the following aspects:
 *	- It is constructible from EXT::MatrixType (see below), we call this matrix A.
 *	- It contains a Q that returns an orthogonal matrix.
 *	- It contains a function R that returns a R^(rank(A),cols(A)) matrix.
 *	          [ R() ]
 *	- A = Q()*[     ]
 *	          [  0  ]
 * Remark: R() does not have to return a triangular matrix.
 */

class OrthDecomposition{
private:
	Matrix<sp_scalar,Dynamic,Dynamic> _q,_r;
	OrthDecomposition();
public:
	OrthDecomposition(const Matrix<sp_scalar,Dynamic,Dynamic>& A){
		//use eigen qr-decomposition with pivoting
		ColPivHouseholderQR<Matrix<sp_scalar,Dynamic,Dynamic> > decomp(A);

		_q = decomp.matrixQ();
		Matrix<sp_scalar,Dynamic,Dynamic> r = decomp.matrixR();

		//eliminate rounding errors (this is necessary)
		for(int j=0;j<std::min(r.rows(),r.cols());++j){
			for(int i=j+1;i<r.rows();++i) r(i,j)=0;
		}

		Matrix<sp_scalar,Dynamic,Dynamic> p = decomp.colsPermutation();
		//A*p=_q*r

		const Matrix<sp_scalar,Dynamic,Dynamic> tmp = r*p.transpose(); //A=_q*tmp
		int rank = decomp.rank();
		_r=tmp.block(0,0,rank,tmp.cols()); //truncate
	}
	const Matrix<sp_scalar,Dynamic,Dynamic>& Q() const { return _q; }
	const Matrix<sp_scalar,Dynamic,Dynamic>& R() const { return _r; }
};

/* Part 3: the EXT class
 * This class will be passed as template parameter to StoichPack classes.
 * Make sure, that your implementation has the same functions, typedefs
 * and semantics as this example class.
 */			
class EXT{
public:
	typedef Matrix<sp_scalar,Dynamic,1> VectorType; // (algebraic) vectors are stored like this

	static VectorType CreateVector(size_t s) { return VectorType(s); } //Create a (algebraic) vector of dimension s
	static VectorType CreateVector(size_t s, sp_scalar v) // Create vector of dimension s with all entries initialized as v
		{ return Array<sp_scalar,Dynamic,1>::Constant(s,1,v); }

	static size_t size(const VectorType& x) { return x.rows(); } //dimension of a vector

	static EigenIterator<VectorType,false> Begin(VectorType& x) //iterator to the first entry of a vector
		{ return EigenIterator<VectorType,false>(x,EigenIteratorPos(0,0)); }
	static EigenIterator<VectorType,false> End(VectorType& x) //iterator to the past-end entry of a vector
		{ return EigenIterator<VectorType,false>(x,EigenIteratorPos(0,EXT::size(x))); }
	static EigenConstIterator<VectorType,false> ConstBegin(const VectorType& x) //const iterator to the first entry of a vector
		{ return EigenConstIterator<VectorType,false>(x,EigenIteratorPos(0,0)); }
	static EigenConstIterator<VectorType,false> ConstEnd(const VectorType& x) //const iterator to the past-end entry of a vector
		{ return EigenConstIterator<VectorType,false>(x,EigenIteratorPos(0,size(x))); }

	// arrays of (algebraic) vectors. Matrices are not suitable, since each vector may have different dimension.
	typedef std::vector<VectorType> VectorArrayType;

	static VectorArrayType CreateVectorArray(const VectorType& x) //Create a vector array with 1 entry
		{ return VectorArrayType(1,x); }

	static size_t size(const VectorArrayType& x) { return x.size(); } //length of array

	//iterators and const iterators for VectorArrayType
	static VectorArrayType::iterator Begin(VectorArrayType& x) { return x.begin(); }
	static VectorArrayType::const_iterator ConstBegin(const VectorArrayType& x) { return x.begin(); }
	static VectorArrayType::iterator End(VectorArrayType& x){ return x.end(); }
	static VectorArrayType::const_iterator ConstEnd(const VectorArrayType& x) { return x.end(); }

	static VectorArrayType ReserveVectorArray(size_t s) { //reserve a certain number of entries without creating them
		//this function may also just return an empty array
		//implementation is recommended due to performance issues
		VectorArrayType ret;
		ret.reserve(s);
		return ret;
	}
	static VectorArrayType& PushBack(VectorArrayType& ret, const VectorType& value){ //append value to ret
		ret.push_back(value);
		return ret;
	}

	static void set(VectorArrayType& x, const VectorArrayType& y) { //copy element by element without reallocation
		//this is not x=y !!!!
		const size_t s=y.size();
		assert(x.size()==s);
		for(size_t i=0;i<s;++i) x[i]=y[i];
	}


	typedef Matrix<sp_scalar,Dynamic,Dynamic> MatrixType; //matrices are stored like this

	/* In contrast to vectors, we do not need a matrix array class. */

	static size_t rows(const MatrixType& A) { return A.rows(); } //number of rows
	static size_t cols(const MatrixType& A) { return A.cols(); } //number of columns
	static MatrixType CreateMatrix(size_t r, size_t c) { return MatrixType(r,c); } //Create r x c matrix
	static MatrixType CreateMatrix(size_t n) { return CreateMatrix(n,n); } //Create n x n matrix
	static MatrixType CreateMatrix(size_t r, size_t c, sp_scalar v) //Create r x c matrix with all entries initialized as v
		{return MatrixType::Constant(r,c,v); }
	static MatrixType Transposed(const MatrixType& M) { return M.transpose(); } // transposed matrix

	//iterators and const iterators for matrices
	//ColWise: iterate over all entries in a fixed column
	//RowWise: iterate over all entries in a fixed row
	static EigenIterator<MatrixType,false> BeginColWise(MatrixType& M, size_t col) {
		return EigenIterator<MatrixType,false>(M,EigenIteratorPos(col,0));
	}
	static EigenConstIterator<MatrixType,false> ConstBeginColWise(const MatrixType& M, size_t col) {
		return EigenConstIterator<MatrixType,false>(M,EigenIteratorPos(col,0));
	}
	static EigenIterator<MatrixType,false> EndColWise(MatrixType& M, size_t col) {
		return EigenIterator<MatrixType,false>(M,EigenIteratorPos(col,M.rows()));
	}
	static EigenConstIterator<MatrixType,false> ConstEndColWise(const MatrixType& M, size_t col) {
		return EigenConstIterator<MatrixType,false>(M,EigenIteratorPos(col,M.rows()));
	}
	static EigenIterator<MatrixType,true> BeginRowWise(MatrixType& M, size_t row) {
		return EigenIterator<MatrixType,true>(M,EigenIteratorPos(row,0));
	}
	static EigenConstIterator<MatrixType,true> ConstBeginRowWise(const MatrixType& M, size_t row) {
		return EigenConstIterator<MatrixType,true>(M,EigenIteratorPos(row,0));
	}
	static EigenIterator<MatrixType,true> EndRowWise(MatrixType& M, size_t row) {
		return EigenIterator<MatrixType,true>(M,EigenIteratorPos(row,M.cols()));
	}

	static EigenConstIterator<MatrixType,true> ConstEndRowWise(const MatrixType& M, size_t row) {
		return EigenConstIterator<MatrixType,true>(M,EigenIteratorPos(row,M.cols()));
	}

	typedef OrthDecomposition OrthogonalDecompositionType;
};

#endif
