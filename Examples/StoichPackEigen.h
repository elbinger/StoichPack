#ifndef __H_STOICHPACK_EIGEN__
#define __H_STOICHPACK_EIGEN__

#include "StoichPack.h"
#include <Eigen>
#include <iterator>

using namespace StoichPack;
using namespace Eigen;

class EigenIteratorPos{
public:
	const size_t X;
	size_t Y;

	EigenIteratorPos(size_t x, size_t y) : X(x),Y(y) {}
	
	EigenIteratorPos& operator+=(size_t i) { Y+=i; return *this; }

	EigenIteratorPos operator+(size_t i) const { return EigenIteratorPos(*this)+=i; }

	bool operator==(const EigenIteratorPos& other) const { assert(X==other.X); return Y==other.Y; }
	bool operator!=(const EigenIteratorPos& other) const { return !(*this==other); }

	EigenIteratorPos& operator=(const EigenIteratorPos& other) { assert(X==other.X); Y=other.Y; return *this; }
};

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

#include <iostream>
using namespace std;

class OrthDecomposition{
private:
	Matrix<sp_scalar,Dynamic,Dynamic> _q,_r;
	OrthDecomposition();
public:
	OrthDecomposition(const Matrix<sp_scalar,Dynamic,Dynamic>& A){
		ColPivHouseholderQR<Matrix<sp_scalar,Dynamic,Dynamic> > decomp(A);
		Matrix<sp_scalar,Dynamic,Dynamic> rr = decomp.matrixR();
		for(int j=0;j<min(rr.rows(),rr.cols());++j){
			for(int i=j+1;i<rr.rows();++i) rr(i,j)=0;
		}
		Matrix<sp_scalar,Dynamic,Dynamic> p = decomp.colsPermutation();
		Matrix<sp_scalar,Dynamic,Dynamic> tmp = rr*p.transpose();
		int rank = decomp.rank();
		_r=tmp.block(0,0,rank,tmp.cols());
		_q=decomp.matrixQ();
	}
	const Matrix<sp_scalar,Dynamic,Dynamic>& Q() const { return _q; }
	const Matrix<sp_scalar,Dynamic,Dynamic>& R() const { return _r; }
};
			
class EXT{
public:
	typedef Matrix<sp_scalar,Dynamic,1> VectorType;
	typedef Pair<VectorType> VectorPairType;

	static size_t size(const VectorType& x) { return x.rows(); }
	static VectorType CreateVector(size_t s) { return VectorType(s); }
	static VectorType CreateVector(size_t s, sp_scalar v){ return Array<sp_scalar,Dynamic,1>::Constant(s,1,v); }
	static VectorType CreateZeroVector(size_t s) { return Array<sp_scalar,Dynamic,1>::Zero(s,1); }

	typedef EigenIterator<VectorType,false> VectorIteratorType;
	typedef EigenConstIterator<VectorType,false> ConstVectorIteratorType;

	static VectorIteratorType Begin(VectorType& x) { return VectorIteratorType(x,EigenIteratorPos(0,0)); }
	static VectorIteratorType End(VectorType& x) { return VectorIteratorType(x,EigenIteratorPos(0,EXT::size(x))); }
	static ConstVectorIteratorType ConstBegin(const VectorType& x) { return ConstVectorIteratorType(x,EigenIteratorPos(0,0)); }
	static ConstVectorIteratorType ConstEnd(const VectorType& x) { return ConstVectorIteratorType(x,EigenIteratorPos(0,size(x))); }

	typedef std::vector<VectorType> VectorArrayType;
	typedef Pair<VectorArrayType> VectorArrayPairType;

	static VectorArrayType CreateVectorArray(const VectorType& x) { return VectorArrayType(1,x); }
	static VectorArrayPairType CreateVectorArrayPair(const VectorType& x, const VectorType& y) {
		return VectorArrayPairType(VectorArrayType(1,x),VectorArrayType(1,y));
	}

	static size_t size(const VectorArrayType& x) { return x.size(); }

	typedef std::vector<VectorType>::iterator VectorArrayIteratorType;
	typedef std::vector<VectorType>::const_iterator ConstVectorArrayIteratorType;

	static VectorArrayIteratorType Begin(VectorArrayType& x) { return x.begin(); }
	static ConstVectorArrayIteratorType ConstBegin(const VectorArrayType& x) { return x.begin(); }
	static VectorArrayIteratorType End(VectorArrayType& x){ return x.end(); }
	static ConstVectorArrayIteratorType ConstEnd(const VectorArrayType& x) { return x.end(); }

	static VectorArrayType ReserveVectorArray(size_t s) {
		VectorArrayType ret;
		ret.reserve(s);
		return ret;
	}
	static VectorArrayType& PushBack(VectorArrayType& ret, const VectorType& value){
		ret.push_back(value);
		return ret;
	}

	typedef Matrix<sp_scalar,Dynamic,Dynamic> MatrixType;
	typedef Pair<MatrixType> MatrixPairType;
	typedef Pair<MatrixPairType> MatrixQuadType;

	static size_t rows(const MatrixType& A) { return A.rows(); }
	static size_t cols(const MatrixType& A) { return A.cols(); }
	static MatrixType CreateMatrix(size_t r, size_t c) { return MatrixType(r,c); }
	static MatrixType CreateZeroMatrix(size_t r, size_t c) { return Array<sp_scalar,Dynamic,Dynamic>::Zero(r,c); }
	static MatrixType Transposed(const MatrixType& M) { return M.transpose(); }

	typedef EigenIterator<MatrixType,false> ColWiseIteratorType;
	typedef EigenConstIterator<MatrixType,false> ConstColWiseIteratorType;
	typedef EigenIterator<MatrixType,true> RowWiseIteratorType;
	typedef EigenConstIterator<MatrixType,true> ConstRowWiseIteratorType;

	static ColWiseIteratorType BeginColWise(MatrixType& M, size_t col) { return ColWiseIteratorType(M,EigenIteratorPos(col,0)); }
	static ConstColWiseIteratorType ConstBeginColWise(const MatrixType& M, size_t col) {
		return ConstColWiseIteratorType(M,EigenIteratorPos(col,0));
	}
	static ColWiseIteratorType EndColWise(MatrixType& M, size_t col) {
		return ColWiseIteratorType(M,EigenIteratorPos(col,M.rows()));
	}
	static ConstColWiseIteratorType ConstEndColWise(const MatrixType& M, size_t col) {
		return ConstColWiseIteratorType(M,EigenIteratorPos(col,M.rows()));
	}
	static RowWiseIteratorType BeginRowWise(MatrixType& M, size_t row) { return RowWiseIteratorType(M,EigenIteratorPos(row,0)); }
	static ConstRowWiseIteratorType ConstBeginRowWise(const MatrixType& M, size_t row) {
		return ConstRowWiseIteratorType(M,EigenIteratorPos(row,0));
	}
	static ConstRowWiseIteratorType ConstEndRowWise(const MatrixType& M, size_t row) {
		return ConstRowWiseIteratorType(M,EigenIteratorPos(row,M.cols()));
	}

	static void set(VectorArrayType& x, const VectorArrayType& y) {
		const size_t s=y.size();
		assert(x.size()==s);
		for(size_t i=0;i<s;++i) x[i]=y[i];
	}

	typedef OrthDecomposition OrthogonalDecompositionType;
};

#endif
