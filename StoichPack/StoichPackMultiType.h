#ifndef __H_STOICHPACK_MULTITYPE__
#define __H_STOICHPACK_MULTITYPE__

namespace StoichPack{

 /* A pair class. The beahaviour is the same as for std::pair, except that first becomes Mobile() and second becomes
  * Immobile(). This increases readability of the code. */
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

 template<typename T>
 class Quad : public Pair< Pair<T> >{
 public:
	Quad(const Pair<T>& x, const Pair<T>& y) : Pair< Pair<T> >(x,y) {}
	Quad(const T& x1, const T& x2, const T& y1, const T& y2) : Pair< Pair<T> >( Pair<T>(x1,x2), Pair<T>(y1,y2) ) {}
 };

 template<typename EXT>
 using VectorPair = Pair<typename EXT::VectorType>;

 template<typename EXT>
 using VectorArrayPair = Pair<typename EXT::VectorArrayType >;

 template<typename EXT>
 using MatrixPair = Pair<typename EXT::MatrixType>;

 template<typename EXT>
 using MatrixQuad = Quad<typename EXT::MatrixType>;

}

#endif
