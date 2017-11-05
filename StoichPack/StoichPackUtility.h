/* File: StoichPackUtility.h
 * Purpose: Provide some functions needed in several modules of StoichPack.
 * Author: Tobias Elbinger <elbinger@math.fau.de> */

#ifndef __H_STOICHPACK_UTILITY__
#define __H_STOICHPACK_UTILITY__

#include <vector>
#include <algorithm>
#include <cassert>
#include "StoichPackDefines.h"
#include "StoichPackMultiType.h"

namespace StoichPack{

 /* Combine 2 vectors x and y to one vector containing the entries of x and y, starting with the entries of x. */ 
 template<typename EXT>
 typename EXT::VectorType Combine(const typename EXT::VectorType& x, const typename EXT::VectorType&  y){
	const size_t sx=EXT::size(x);
	const size_t sy=EXT::size(y);

	typename EXT::VectorType ret = EXT::CreateVector(sx+sy);

	auto it =std::copy(EXT::ConstBegin(x),EXT::ConstEnd(x),EXT::Begin(ret));
	std::copy(EXT::ConstBegin(y),EXT::ConstEnd(y),it);

	return ret;
 }

 /* Combine a vector pair p=(x,y) to one vector containing the entries of x and y, starting with the entries of x. */ 
 template<typename EXT>
 typename EXT::VectorType Combine(const VectorPair<EXT>& p) { return Combine<EXT>(p.Mobile(),p.Immobile()); }

 /* Divide the vector x in 2 vectors. The first vector (accesible via result.Mobile()) contains the first n_mobile entries.
  * The second vector (accessible via result.Immobile()) contains the rest. */
 template<typename EXT>
 VectorPair<EXT> Divide(const typename EXT::VectorType& x, size_t n_mobile) {
	assert(EXT::size(x)>=n_mobile);
	VectorPair<EXT> ret(EXT::CreateVector(n_mobile),EXT::CreateVector(EXT::size(x)-n_mobile));
	std::copy(EXT::ConstBegin(x),EXT::ConstBegin(x)+n_mobile,EXT::Begin(ret.Mobile()));
	std::copy(EXT::ConstBegin(x)+n_mobile,EXT::ConstEnd(x),EXT::Begin(ret.Immobile()));
	return ret;
 }

 /* Vertical concatination of the matrices x and y: if x has nx rows and y has ny rows, the result has nx+ny rows. The first nx rows
  * of the result are the rows of x, the last ny rows are the rows of y. x and y must have the same number of columns. */
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

 /* Horizontal concatination of the matrices x and y: if x has nx columns and y has ny columns, the result has nx+ny columns.
  * The first nx columns of the result are the columns of x, the last ny columns are the columns of y.
  * x and y must have the same number of rows. */
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

 /* Divide the matrix x in two submatrices. The first one (accesible via result.Mobile()) contains the first n_mobile rows of x,
  * the second one (accessible via result.Immobile()) contains the rest. */ 
 template<typename EXT>
 MatrixPair<EXT> DivideRows(const typename EXT::MatrixType& x, size_t m_mobile){
	const size_t m=EXT::rows(x);
	assert(m_mobile<=m);
	const size_t n=EXT::cols(x);

	MatrixPair<EXT> ret(EXT::CreateMatrix(m_mobile,n),EXT::CreateMatrix(m-m_mobile,n));

	for(size_t i=0;i<m_mobile;++i) std::copy(EXT::ConstBeginRowWise(x,i),EXT::ConstEndRowWise(x,i),EXT::BeginRowWise(ret.Mobile(),i));
	for(size_t i=m_mobile;i<m;++i) std::copy(EXT::ConstBeginRowWise(x,i),EXT::ConstEndRowWise(x,i),
	                                         EXT::BeginRowWise(ret.Immobile(),i-m_mobile));

	return ret;
 }

 /* Divide the matrix x in two submatrices. The first one (accesible via result.Mobile()) contains the first n_mobile columns of x,
  * the second one (accessible via result.Immobile()) contains the rest. */ 
 template<typename EXT>
 MatrixPair<EXT> DivideCols(const typename EXT::MatrixType& x, size_t n_mobile){
	const size_t n=EXT::cols(x);
	assert(n_mobile<=n);
	const size_t m=EXT::rows(x);

	MatrixPair<EXT> ret(EXT::CreateMatrix(m,n_mobile),EXT::CreateMatrix(m,n-n_mobile));

	for(size_t i=0;i<n_mobile;++i) std::copy(EXT::ConstBeginColWise(x,i),EXT::ConstEndColWise(x,i),EXT::BeginColWise(ret.Mobile(),i));
	for(size_t i=n_mobile;i<n;++i) std::copy(EXT::ConstBeginColWise(x,i),EXT::ConstEndColWise(x,i),
	                                         EXT::BeginColWise(ret.Immobile(),i-n_mobile));

	return ret;
 }

 /* Get a matrix that consists of the columns of M specified in I. */
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
 typename EXT::MatrixType SubRows(const typename EXT::MatrixType& M, const std::vector<size_t>& I){
	const size_t s = I.size();
	typename EXT::MatrixType ret = EXT::CreateMatrix(s,EXT::cols(M));
	for(size_t i=0;i<s;++i){
		assert(I[i]<EXT::rows(M));
		std::copy(EXT::ConstBeginRowWise(M,I[i]),EXT::ConstEndRowWise(M,I[i]),EXT::BeginRowWise(ret,i));
	}
	return ret;
 }

 /* Get a vector that consists of the entries of x specified in I. */
 template<typename EXT>
 typename EXT::VectorType SubEntries(const typename EXT::VectorType& x, const std::vector<size_t>& I){
	const size_t s=I.size();
	typename EXT::VectorType ret = EXT::CreateVector(s);
	auto itret = EXT::Begin(ret);
	auto itx =EXT::ConstBegin(x);

	for(size_t i=0;i<s;++i, ++itret) {
		assert(I[i]<EXT::rows(x));
		*itret=*(itx+I[i]);
	}
	return ret;
 }

 /* Get a matrix that speciefies if a entry of M is nonzero, i.e.:
  * result(i,j)=1 <--> M(i,j)!=0
  * result(i,j)=0 <--> M(i,j)=0 */
 template<typename EXT>
 typename EXT::MatrixType BooleanMatrix(const typename EXT::MatrixType& M){
	typename EXT::MatrixType ret(M);
	for(size_t i=0;i<EXT::cols(M);++i){
		for(auto it=EXT::BeginColWise(ret,i); it!=EXT::EndColWise(ret,i); ++it){
			if(*it!=0) *it=1;
		}
	}
	return ret;
 }

 /* Create a diagonal square matrix, containing the entries of values on its diagonal. */
 template<typename EXT, typename T=sp_scalar>
 typename EXT::MatrixType Diag(const std::vector<T>& values){
	const size_t s=values.size();
	typename EXT::MatrixType ret = EXT::CreateZeroMatrix(s,s);
	for(size_t i=0;i<s;++i) *(EXT::BeginRowWise(ret,i)+i)=values[i];
	return ret;
 }

 template<typename EXT, typename IT>
 typename EXT::MatrixType ColCat(IT begin, IT end){
	assert(begin!=end);
	if(begin+1==end) return *begin;
	else return CombineCols<EXT>(*begin,ColCat<EXT,IT>(begin+1,end));
 }

 template<typename EXT, typename IT>
 typename EXT::MatrixType RowCat(IT begin, IT end){
	assert(begin!=end);
	if(begin+1==end) return *begin;
	else return CombineRows<EXT>(*begin,RowCat<EXT,IT>(begin+1,end));
 }

 template<typename EXT>
 MatrixPair<EXT> SplitBlocks(const typename EXT::MatrixType& M, size_t rows, size_t cols){
	const MatrixPair<EXT> tmp = DivideCols<EXT>(M,cols);
	return MatrixPair<EXT>(DivideRows<EXT>(tmp.Mobile(),rows).Mobile(),DivideRows<EXT>(tmp.Immobile(),rows).Immobile());
 }

 template<typename EXT>
 typename EXT::MatrixType CombineBlocks(const typename EXT::MatrixType M11, const typename EXT::MatrixType M22) {
	const typename EXT::MatrixType M12 = EXT::CreateZeroMatrix(EXT::rows(M11),EXT::cols(M22));
	const typename EXT::MatrixType M21 = EXT::CreateZeroMatrix(EXT::rows(M22),EXT::cols(M11));
	return CombineRows<EXT>(CombineCols<EXT>(M11,M12),CombineCols<EXT>(M21,M22));
 }
} // namespace StoichPack

#endif

