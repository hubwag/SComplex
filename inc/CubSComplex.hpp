#ifndef CUB_SCOMPLEX_HPP
#define CUB_SCOMPLEX_HPP

#include <boost/assert.hpp>
#include <boost/optional.hpp>
#include <utility>


#include <capd/auxil/CRef.h>
#include <capd/homologicalAlgebra/embeddingDim.h>
#include <capd/vectalg/MatrixSlice.h>
#include <capd/matrixAlgorithms/intMatrixAlgorithms.hpp>
#include <capd/homologicalAlgebra/homologicalAlgebra.hpp>
#include <capd/homologicalAlgebra/homAlgFunctors.hpp>
#include <capd/homologicalAlgebra/cubSetFunctors.hpp>
#include <capd/homologicalAlgebra/ReducibleFreeChainComplex.hpp>
#include <capd/homologicalAlgebra/readCubCellSet.hpp>

#include <capd/bitSet/CubCellSetT.hpp>
#include <capd/bitSet/CubSetT.hpp>
#include <capd/bitSet/EuclBitSetT.hpp>





template<int DIM>
class CubSComplex {
  class CellImpl;
  typedef unsigned long int cluster;
  typedef BitSetT<BitmapT<cluster> > BitSet;
  typedef EuclBitSetT<BitSet, DIM> EuclBitSet;
  
public:
  typedef CubSetT<EuclBitSet> BCubSet;
  typedef CubCellSetT<EuclBitSet> BCubCellSet;

  class BitCoordCellImpl;
  class DynamicBitCoordCellImpl;
  class BitCoordPtrCellImpl;

  template<typename, typename Enable=void>
  class CellProxy;

  typedef CellProxy<BitCoordCellImpl> Cell;
  typedef CellProxy<BitCoordPtrCellImpl> CellRef;
  typedef CellProxy<DynamicBitCoordCellImpl> DynamicCell;
  
  typedef int Dim;
  typedef int Color;
  //typedef std::pair<const BCubCellSet::BitCoordIterator::WordIterator, int> Id;
  typedef long Id;
  
private:
  
  template<typename NumeratorT, bool isConst>
  class IteratorProvider;

  
  class BasicNumerator;
  class CellNumerator;
  class CellDimNumerator;
  class CbdNumerator;
  class BdNumerator;

  template<bool isConst>
  class IteratorsImpl;

  template<bool isConst>
  class ColoredIteratorsImpl {

  public:
	 typedef IteratorsImpl<isConst> Iterators;

	 template<typename CubSComplex::Color color>
	 class Color {
	 public:
		typedef IteratorsImpl<isConst> Iterators;
	 };
  };

public:
  
  typedef IteratorsImpl<false> Iterators;
  typedef IteratorsImpl<true> ConstIterators;


  typedef ColoredIteratorsImpl<false> ColoredIterators;
  typedef ColoredIteratorsImpl<true> ColoredConstIterators;
  
  CubSComplex();
  CubSComplex(const int* A_w, bool clear=false);
  explicit CubSComplex(CRef<BCubCellSet> _bCubCellSet);

  size_t cardinality() {  return bCubCellSet.cardinality(); }
  size_t size() const { return bCubCellSet.getBmpSizeInBits(); }

  Dim getDim() { return bCubCellSet.embDim(); }
  
   ConstIterators iterators() const;
   Iterators iterators();

  typename ColoredConstIterators::Iterators iterators(const Color& color) const;
  typename ColoredIterators::Iterators iterators(const Color& color);

  template<Color color>
  typename ColoredIterators::template Color<color>::Iterators iterators();

  template<Color color>
  typename ColoredConstIterators::template Color<color>::Iterators iterators() const;

  template<typename ImplT>
  bool getUniqueFace(const CellProxy<ImplT>& cell, DynamicCell& coface) const;

  template<typename ImplT>
  bool getUniqueCoFace(const CellProxy<ImplT>& cell, DynamicCell& coface) const;

  template<typename ImplT1, typename ImplT2>
  int coincidenceIndex(const CellProxy<ImplT1> &a, const CellProxy<ImplT2> &b) const;
  
protected:
  CRef<BCubCellSet> _bCubCellSetCR;
  BCubCellSet& bCubCellSet;
  
  template<bool isConst>
  friend class IteratorsImpl;
  
};

#include "CubSComplex_Cell.hpp"
#include "CubSComplex_IteratorProvider.hpp"
#include "CubSComplex_Iterators.hpp"
#include "CubSComplex_ColoredIterators.hpp"
#include "CubSComplex_Numerators.hpp"

// inline CubSComplex::CubSComplex():
//   bCubCellSetCR(new BCubCellSet()), baseDimension(0)
// {}

// inline CubSComplex::CubSComplex(const int* A_w, bool clear):
//   bCubCellSetCR(new BCubCelSet(A_w,clear)),baseDimension(0)
// {}

template<int DIM>
inline CubSComplex<DIM>::CubSComplex(CRef<BCubCellSet> _bCubCellSetCR): _bCubCellSetCR(_bCubCellSetCR), bCubCellSet(_bCubCellSetCR()) {
  bCubCellSet.addEmptyCollar();
}

template<int DIM>
inline typename CubSComplex<DIM>::Iterators CubSComplex<DIM>::iterators() {
  throw std::logic_error("Not implemented yet."); // How cn I iterate over removed elements ?

  return Iterators(*this);
}

template<int DIM>
inline typename CubSComplex<DIM>::ConstIterators CubSComplex<DIM>::iterators() const {
  throw std::logic_error("Not implemented yet."); // How cn I iterate over removed elements ?

  return ConstIterators(*this);
}

template<int DIM>
inline typename CubSComplex<DIM>::ColoredConstIterators::Iterators CubSComplex<DIM>::iterators(const Color& color) const {
  BOOST_ASSERT(color == 1);
  return typename ColoredConstIterators::Iterators(*this);
}

template<int DIM>
inline typename CubSComplex<DIM>::ColoredIterators::Iterators CubSComplex<DIM>::iterators(const Color& color) {
  BOOST_ASSERT(color == 1);
  return typename ColoredIterators::Iterators(*this);
}

template<int DIM>
template<typename CubSComplex<DIM>::Color color>
inline typename CubSComplex<DIM>::ColoredIterators::template Color<color>::Iterators CubSComplex<DIM>::iterators() {
  return typename ColoredIterators::template Color<color>::Iterators(*this);
}

template<int DIM>
template<typename CubSComplex<DIM>::Color color>
inline typename CubSComplex<DIM>::ColoredConstIterators::template Color<color>::Iterators CubSComplex<DIM>::iterators() const {
  return typename ColoredConstIterators::template Color<color>::Iterators(*this);
}

template<int DIM>
template<typename ImplT>
inline bool CubSComplex<DIM>::getUniqueCoFace(const CellProxy<ImplT>& cell, DynamicCell& coface) const {
  if (bCubCellSet.isFreeFace(const_cast<typename CubSComplex::BCubCellSet::BitCoordIterator&>(cell.getBitCoordIt()), coface.getBitCoordIt())) {
  	 return true;
  } else {
  	 return false;
  }
}

template<int DIM>
template<typename ImplT>
inline bool CubSComplex<DIM>::getUniqueFace(const CellProxy<ImplT>& cell, DynamicCell& coface) const {
  if (bCubCellSet.isFreeCoFace(const_cast<typename CubSComplex::BCubCellSet::BitCoordIterator&>(cell.getBitCoordIt()), coface.getBitCoordIt())) {
  	 return true;
  } else {
  	 return false;
  }
}

template<int DIM>
template<typename ImplT1, typename ImplT2>
inline int CubSComplex<DIM>::coincidenceIndex(const CellProxy<ImplT1> &_a, const CellProxy<ImplT2> &_b) const {
  int res = 0;
  int sgn = 1;

  const typename CubSComplex::BCubCellSet::BitCoordIterator& a = _a.getBitCoordIt();
  const typename CubSComplex::BCubCellSet::BitCoordIterator& b = _b.getBitCoordIt();
  
  for (size_t i = 0, end = bCubCellSet.embDim(); i < end; ++i) {
  	 if (! (a[i]/2 == b[i]/2 || a[i]/2 + (a[i]%2) == b[i]/2)) {
  		return 0; // b[i] left side doesn't intersect a[i] interval
  	 }
	 
  	 if (a[i] % 2 == 1 && b[i] % 2 == 1) {
  		sgn = -sgn; //both nondegenerated, go to next
  	 } else if (a[i] % 2 == 0 && b[i] % 2 == 0) {
  		// both degenerated
  	 } else if (b[i] % 2 == 1) {		
  		return 0; // a[i] is inside b[i]
  	 } else { // b[i] is inside a[i]
  		if (res != 0) {
  		  return 0; // second time, so not proper face
  		}
  		res = sgn;
  		if (a[i] / 2 != b[i]/2) { // b[i] is the right face of a[i]
  		  res = -res;
  		}
  	 }
  }
  return res;
}


template<int DIM>
boost::shared_ptr<CubSComplex<DIM> > readCubSComplex(std::string fileName) {
  return boost::shared_ptr<CubSComplex<DIM> >(new CubSComplex<DIM>(readCubCellSet<typename CubSComplex<DIM>::BCubSet, typename CubSComplex<DIM>::BCubCellSet>(fileName.c_str())));
}

#endif

