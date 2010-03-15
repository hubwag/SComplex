#ifndef CUB_SCOMPLEX_ITERATORS_HPP
#define CUB_SCOMPLEX_ITERATORS_HPP

#include "CubSComplex.hpp"
#include "util/Iterators.hpp"

#include <boost/mpl/if.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/range.hpp>
#include <boost/shared_ptr.hpp>
#include <functional>


template<bool isConst>
class CubSComplex::IteratorsImpl {
  typedef typename boost::mpl::if_c<isConst, const CubSComplex&, CubSComplex&>::type SComplexRef;
public:
  typedef CubSComplex::IteratorProvider<CellNumerator, isConst> AllCells;
  //typedef CubSComplex::IteratorProvider<CellDimNumerator, isConst> DimCells;
  typedef CubSComplex::IteratorProvider<BdNumerator, isConst> BdCells;  
  typedef CubSComplex::IteratorProvider<CbdNumerator, isConst> CbdCells;

  typedef CellProxy<BitCoordPtrCellImpl> CellType;


  template<typename Derived>
  class CoordIterator: public boost::iterator_facade<Derived, CellType,
																	  boost::forward_traversal_tag, CellType>
  {
  
  public:
	 typedef BCubCellSet::BitCoordIterator BitCoordIterator;
	 typedef BCubCellSet::BitIterator BitIterator;
  
	 CoordIterator(): coordIt(*((EuclBitSet*)NULL)) {}

	 CoordIterator(CubSComplex& s): coordIt(new BitCoordIterator(s.bCubCellSet)) {}
  
  protected:
	 friend class boost::iterator_core_access;
	 BitCoordIterator* coordIt;
  
	 bool equal(const Derived& other) const {
		return (BitIterator&)(*coordIt) == ((BitIterator&)(*other.coordIt));
	 }

	 CellType dereference() const {
		return CellType(coordIt);
	 }
  };

  class DimCoordIterator: public  CoordIterator<DimCoordIterator> {

  public:
	 DimCoordIterator(): CoordIterator<DimCoordIterator>(), dim(std::numeric_limits<Dim>::max()) {}
  
	 DimCoordIterator(CubSComplex& s, const Dim _dim): CoordIterator<DimCoordIterator>(s), dim(_dim) { //begin iterator
		findDim();
	 }

	 DimCoordIterator(CubSComplex& s): CoordIterator<DimCoordIterator>(s), dim(std::numeric_limits<Dim>::max()) { //end constructor
		coordIt->wIt = const_cast<BCubCellSet::BitCoordIterator::WordIterator>(coordIt->getBitmap().getBitmapEnd());
	 }
  
  private:
	 using CoordIterator<DimCoordIterator>::coordIt;
	 
	 friend class boost::iterator_core_access;

	 void findDim() {
		while (coordIt->wIt < coordIt->getBitmap().getBitmapEnd()) {
		  coordIt->moveToFirstPixel();
		  if(coordIt->ownDim()==dim) break;
		  ++(*coordIt);
		}	 
	 }
  
	 void increment() {
		++(*coordIt);
		findDim();
	 }

	 Dim dim;  
  };

  typedef typename boost::iterator_range<DimCoordIterator > DimCells;


  IteratorsImpl(SComplexRef _scomplex): scomplex(_scomplex) {}
  
  AllCells allCells() const;
  DimCells dimCells(const Dim& dim);
  BdCells bdCells(const Cell& cell) const;
  CbdCells cbdCells(const Cell& cell) const;
  
private:
  SComplexRef scomplex;
};

template<bool isConst>
inline typename CubSComplex::IteratorsImpl<isConst>::AllCells CubSComplex::IteratorsImpl<isConst>::allCells() const {
  return AllCells(CellNumerator(scomplex));
}

template<bool isConst>
inline typename CubSComplex::IteratorsImpl<isConst>::DimCells CubSComplex::IteratorsImpl<isConst>::dimCells(const Dim& dim) {
  //return DimCells(CellDimNumerator(scomplex, dim));
  return DimCells(DimCoordIterator(scomplex, dim), DimCoordIterator(scomplex));
}

template<bool isConst>
inline typename CubSComplex::IteratorsImpl<isConst>::CbdCells CubSComplex::IteratorsImpl<isConst>::cbdCells(const Cell& cell) const {
  return CbdCells(CbdNumerator(scomplex, cell));
}

template<bool isConst>
inline typename CubSComplex::IteratorsImpl<isConst>::BdCells CubSComplex::IteratorsImpl<isConst>::bdCells(const Cell& cell) const {
  return BdCells(BdNumerator(scomplex, cell));
}

#endif
