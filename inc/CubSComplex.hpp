#ifndef CUB_SCOMPLEX_HPP
#define CUB_SCOMPLEX_HPP


#include <capd/repSet/ElementaryCell.h>
#include <capd/auxil/CRef.h>
#include <capd/bitSet/CubCellSetT.hpp>
#include <capd/bitSet/CubSetT.hpp>
#include <capd/homologicalAlgebra/readCubCellSet.hpp>


#include <boost/assert.hpp>
#include <boost/optional.hpp>


class CubSComplex {
public:
  
  class Cell;
  typedef size_t Dim;
  typedef int Color;

private:

  template<typename NumeratorT, bool isConst>
  class IteratorProvider;

  class BasicNumerator;
  class CellNumerator;
  class CellDimNumerator;
  class CbdNumerator;
  class BdNumerator;

  typedef unsigned long int cluster;
  typedef BitSetT<BitmapT<cluster> > BitSet;
  typedef EuclBitSetT<BitSet,embeddingDim> EuclBitSet;
  typedef CubSetT<EuclBitSet> BCubSet;
  typedef CubCellSetT<EuclBitSet> BCubCellSet;


  template<bool isConst>
  class IteratorsImpl;

  template<bool isConst>
  class ColoredIteratorsImpl {

  public:
	 typedef IteratorsImpl<isConst> Iterators;

	 template<CubSComplex::Color color>
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

  size_t cardinality();

   ConstIterators iterators() const;
   Iterators iterators();

  ColoredConstIterators::Iterators iterators(const Color& color) const;
  ColoredIterators::Iterators iterators(const Color& color);

  template<Color color>
  typename ColoredIterators::Color<color>::Iterators iterators();

  template<Color color>
  typename ColoredConstIterators::Color<color>::Iterators iterators() const;

  bool getUniqueCoFace(const Cell& cell, Cell& coface) const;

  bool getUniqueFace(const Cell& cell, Cell& coface) const;
  
  Dim getBaseDimension() const;

  int coincidenceIndex(const Cell &a, const Cell &b) const;
  
protected:
  CRef<BCubCellSet> bCubCellSetCR;
  Dim baseDimension;
  
  template<bool isConst>
  friend class IteratorsImpl;
  
};

inline CubSComplex::Dim CubSComplex::getBaseDimension() const {
  return baseDimension;
}

inline size_t CubSComplex::cardinality() {
  return bCubCellSetCR().cardinality();
}

#include "CubSComplex_Cell.hpp"
#include "CubSComplex_IteratorProvider.hpp"
#include "CubSComplex_Iterators.hpp"
#include "CubSComplex_ColoredIterators.hpp"
#include "CubSComplex_Numerators.hpp"

inline CubSComplex::CubSComplex():
  bCubCellSetCR(new BCubCellSet()), baseDimension(0)
{}

inline CubSComplex::CubSComplex(const int* A_w, bool clear):
  bCubCellSetCR(new BCubCelSet(A_w,clear)),baseDimension(0)
{}

inline CubSComplex::CubSComplex(CRef<BCubCellSet> _bCubCellSet):baseDimension(0){
  bCubCellSetCR = _bCubCellSet;
  bCubCellSetCR().addEmptyCollar();
}

inline CubSComplex::Iterators CubSComplex::iterators() {
  throw std::logic_error("Not implemented yet."); // How cn I iterate over removed elements ?

  return Iterators(*this);
}

inline CubSComplex::ConstIterators CubSComplex::iterators() const {
  throw std::logic_error("Not implemented yet."); // How cn I iterate over removed elements ?

  return ConstIterators(*this);
}

inline CubSComplex::ColoredConstIterators::Iterators CubSComplex::iterators(const Color& color) const {
  BOOST_ASSERT(color == 1);
  return ColoredConstIterators::Iterators(*this);
}

inline CubSComplex::ColoredIterators::Iterators CubSComplex::iterators(const Color& color) {
  BOOST_ASSERT(color == 1);
  return ColoredIterators::Iterators(*this);
}

template<>
inline CubSComplex::ColoredIterators::Color<1>::Iterators CubSComplex::iterators<1>() {
  return ColoredIterators::Color<1>::Iterators(*this);
}

template<>
inline CubSComplex::ColoredConstIterators::Color<1>::Iterators CubSComplex::iterators<1>() const {
  return ColoredConstIterators::Color<1>::Iterators(*this);
}

inline bool CubSComplex::getUniqueCoFace(const Cell& cell, Cell& coface) const {
  if (bCubCellSetCR().isFreeFace(const_cast<Cell&>(cell), coface)) {
	 return true;
  } else {
	 return false;
  }
}

inline bool CubSComplex::getUniqueFace(const Cell& cell, Cell& coface) const {
  if (bCubCellSetCR().isFreeCoFace(const_cast<Cell&>(cell), coface)) {
	 return true;
  } else {
	 return false;
  }
}

inline int CubSComplex::coincidenceIndex(const Cell &a, const Cell &b) const {
  if (bCubCellSetCR().embDim() != bCubCellSetCR().embDim()) {
	 return 0;
  }

  int res = 0;
  int sgn = 1;
  for (size_t i = 0, end = bCubCellSetCR().embDim(); i < end; ++i) {
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

#endif

