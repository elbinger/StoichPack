#ifndef __H_STOICHPACK_UTILITY__
#define __H_STOICHPACK_UTILITY__

#include <vector>
#include <algorithm>
#include <cassert>
#include "StoichPackDefines.h"

namespace StoichPack{
 //some default types and functions. use them if you want to.
 template<typename T>
 class Pair{
 private:
	T mobile, immobile;
 public:
	Pair(const T& x, const T& y) : mobile(x), immobile(y) {}
	T& Mobile() { return mobile; }
	const T& Mobile() const { return mobile; }
	T& Immobile() { return immobile; }
	const T& Immobile() const { return immobile; }
 };
 
 template<typename EXT>
 typename EXT::VectorType Combine(const typename EXT::VectorType& x, const typename EXT::VectorType&  y){
	const size_t sx=EXT::size(x);
	const size_t sy=EXT::size(y);

	typename EXT::VectorType ret = EXT::CreateVector(sx+sy);

	typename EXT::VectorIteratorType it =std::copy(EXT::ConstBegin(x),EXT::ConstEnd(x),EXT::Begin(ret));
	std::copy(EXT::ConstBegin(y),EXT::ConstEnd(y),it);

	return ret;
 }
 
 template<typename EXT>
 typename EXT::VectorPairType Divide(const typename EXT::VectorType& x, size_t n_mobile) {
	assert(EXT::size(x)>=n_mobile);
	typename EXT::VectorPairType ret(EXT::CreateVector(n_mobile),EXT::CreateVector(EXT::size(x)-n_mobile));
	std::copy(EXT::ConstBegin(x),EXT::ConstBegin(x)+n_mobile,EXT::Begin(ret.Mobile()));
	std::copy(EXT::ConstBegin(x)+n_mobile,EXT::ConstEnd(x),EXT::Begin(ret.Immobile()));
	return ret;
 }

 template<typename EXT>
 typename EXT::MatrixType CombineRows(const typename EXT::MatrixType& x, const typename EXT::MatrixType& y){
	assert(EXT::cols(y)==EXT::cols(x));
	const size_t mx=EXT::rows(x);
	const size_t my=EXT::rows(y);

	typename EXT::MatrixType ret= EXT::CreateMatrix(mx+my,EXT::cols(x));

	for(size_t i=0;i<mx;++i) std::copy(EXT::ConstBeginRowWise(x,i),EXT::ConstEndRowWise(x,i),EXT::BeginRowWise(ret,i));
	for(size_t i=0;i<my;++i) std::copy(EXT::ConstBeginRowWise(y,i),EXT::ConstEndRowWise(y,i),EXT::BeginRowWise(ret,i+mx));

	return ret;
 }

 template<typename EXT>
 typename EXT::MatrixType CombineCols(const typename EXT::MatrixType& x, const typename EXT::MatrixType& y){
	assert(EXT::rows(y)==EXT::rows(x));
	const size_t nx=EXT::cols(x);
	const size_t ny=EXT::cols(y);

	typename EXT::MatrixType ret = EXT::CreateMatrix(EXT::rows(x),nx+ny);

	for(size_t i=0;i<nx;++i) std::copy(EXT::ConstBeginColWise(x,i),EXT::ConstEndColWise(x,i),EXT::BeginColWise(ret,i));
	for(size_t i=0;i<ny;++i) std::copy(EXT::ConstBeginColWise(y,i),EXT::ConstEndColWise(y,i),EXT::BeginColWise(ret,i+nx));

	return ret;
 }
 
 template<typename EXT>
 typename EXT::MatrixPairType DivideRows(const typename EXT::MatrixType& x, size_t m_mobile){
	const size_t m=EXT::rows(x);
	assert(m_mobile<=m);
	const size_t n=EXT::cols(x);

	typename EXT::MatrixPairType ret(EXT::CreateMatrix(m_mobile,n),EXT::CreateMatrix(m-m_mobile,n));

	for(size_t i=0;i<m_mobile;++i) std::copy(EXT::ConstBeginRowWise(x,i),EXT::ConstEndRowWise(x,i),EXT::BeginRowWise(ret.Mobile(),i));
	for(size_t i=m_mobile;i<m;++i) std::copy(EXT::ConstBeginRowWise(x,i),EXT::ConstEndRowWise(x,i),
	                                         EXT::BeginRowWise(ret.Immobile(),i-m_mobile));

	return ret;
 }

 template<typename EXT>
 typename EXT::MatrixPairType DivideCols(const typename EXT::MatrixType& x, size_t n_mobile){
	const size_t n=EXT::cols(x);
	assert(n_mobile<=n);
	const size_t m=EXT::rows(x);

	typename EXT::MatrixPairType ret(EXT::CreateMatrix(m,n_mobile),EXT::CreateMatrix(m,n-n_mobile));

	for(size_t i=0;i<n_mobile;++i) std::copy(EXT::ConstBeginColWise(x,i),EXT::ConstEndColWise(x,i),EXT::BeginColWise(ret.Mobile(),i));
	for(size_t i=n_mobile;i<n;++i) std::copy(EXT::ConstBeginColWise(x,i),EXT::ConstEndColWise(x,i),
	                                         EXT::BeginColWise(ret.Immobile(),i-n_mobile));

	return ret;
 }

 template<typename EXT>
 typename EXT::MatrixType SubCols(const typename EXT::MatrixType& M, const std::vector<size_t>& I){
	const size_t s = I.size();
	typename EXT::MatrixType ret = EXT::CreateMatrix(EXT::rows(M),s);
	for(size_t i=0;i<s;++i){
		assert(I[i]<EXT::cols(M));
		std::copy(EXT::ConstBeginColWise(M,I[i]),EXT::ConstEndColWise(M,I[i]),EXT::BeginColWise(ret,i));
	}
	return ret;
 }

 template<typename EXT>
 typename EXT::VectorType SubEntries(const typename EXT::VectorType& x, const std::vector<size_t>& I){
	const size_t s=I.size();
	typename EXT::VectorType ret = EXT::CreateVector(s);
	typename EXT::VectorIteratorType itret = EXT::Begin(ret);
	typename EXT::ConstVectorIteratorType itx =EXT::ConstBegin(x);

	for(size_t i=0;i<s;++i, ++itret) {
		assert(I[i]<EXT::rows(x));
		*itret=*(itx+I[i]);
	}
	return ret;
 }

 template<typename EXT>
 typename EXT::MatrixType BooleanMatrix(const typename EXT::MatrixType& M){
	typename EXT::MatrixType ret(M);
	for(size_t i=0;i<EXT::cols(M);++i){
		for(typename EXT::ColWiseIteratorType it=EXT::BeginColWise(ret,i); it!=EXT::EndColWise(ret,i); ++it){
			if(*it!=0) *it=1;
		}
	}
	return ret;
 }
 template<typename EXT, typename T=sp_scalar>
 typename EXT::MatrixType Diag(const std::vector<T>& values){
	const size_t s=values.size();
	typename EXT::MatrixType ret = EXT::CreateZeroMatrix(s,s);
	for(size_t i=0;i<s;++i) *(EXT::BeginRowWise(ret,i)+i)=values[i];
	return ret;
 }
} // namespace StoichPack

#endif

