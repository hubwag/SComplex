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
  //typedef CubSComplex::IteratorProvider<CbdNumerator, isConst> CbdCells;

  typedef CubCellProxy<BitCoordPtrCellImpl> CellType;


  template<typename Derived>
  class CoordIterator: public boost::iterator_facade<Derived, CellType,
																	  boost::forward_traversal_tag, CellType>
  {
  
  public:
	 typedef BCubCellSet::BitCoordIterator BitCoordIterator;
	 typedef BCubCellSet::BitIterator BitIterator;
  
	 CoordIterator(): coordIt(NULL) {}
	 CoordIterator(CubSComplex& s): coordIt(new (coordItMem) BitCoordIterator(s.bCubCellSet)) {}

	 CoordIterator(const CoordIterator& other) {
		if (other.coordIt) {
		  this->coordIt = new (coordItMem) BitCoordIterator(*other.coordIt);
		} else {
		  this->coordIt = NULL;
		}
	 }
	 
  protected:
	 friend class boost::iterator_core_access;
	 mutable BitCoordIterator* coordIt;
	 char coordItMem[sizeof(BitCoordIterator)];
	 
	 bool equal(const Derived& other) const {
		return (coordIt == NULL && other.coordIt == NULL) ||
		  ( (coordIt != NULL && other.coordIt != NULL)
			 && (BitIterator&)(*coordIt) == ((BitIterator&)(*other.coordIt)));
	 }

	 CellType dereference() const {
		return CellType(&(*coordIt));
	 }
  };

  class DimCoordIterator: public  CoordIterator<DimCoordIterator> {

  public:
	 DimCoordIterator(): CoordIterator<DimCoordIterator>(), dim(std::numeric_limits<Dim>::max()) {}
  
	 DimCoordIterator(CubSComplex& s, const Dim _dim): CoordIterator<DimCoordIterator>(s), dim(_dim) { //begin iterator
		findDim();
	 }

	 DimCoordIterator(CubSComplex& s): CoordIterator<DimCoordIterator>(), dim(std::numeric_limits<Dim>::max()) { //end constructor
		//coordIt->wIt = const_cast<BCubCellSet::BitCoordIterator::WordIterator>(coordIt->getBitmap().getBitmapEnd());
	 }
  
  private:
	 using CoordIterator<DimCoordIterator>::coordIt;
	 
	 friend class boost::iterator_core_access;

	 void findDim() {
		while (coordIt->wIt < coordIt->getBitmap().getBitmapEnd()) {
		  coordIt->moveToFirstPixel();
		  if(coordIt->getBit() && coordIt->ownDim()==dim)
			 return;
		  ++(*coordIt);
		}
		coordIt = NULL;
	 }
  
	 void increment() {
		++(*coordIt);
		findDim();
	 }

	 Dim dim;  
  };

  class CbdCoordIterator: public  CoordIterator<CbdCoordIterator> {
  public:
	 CbdCoordIterator(): CoordIterator<CbdCoordIterator>() {}

	 template<typename ImplT>
	 CbdCoordIterator(CubSComplex& s, const CubCellProxy<ImplT>& c): CoordIterator<CbdCoordIterator>(s),
																							i(-1), dim(s.bCubCellSet.embDim()) { //begin iterator
		*coordIt = c.getBitCoordIt();
		findCbd();
	 }

	 CbdCoordIterator(CubSComplex& s): CoordIterator<CbdCoordIterator>(s),
												  i(0), dim(s.bCubCellSet.embDim()) { //end constructor
		toEnd();
	 }
  
  private:
	 using CoordIterator<CbdCoordIterator>::coordIt;
	 
	 friend class boost::iterator_core_access;

	 void toEnd() {
		//coordIt->wIt = const_cast<BCubCellSet::BitCoordIterator::WordIterator>(coordIt->getBitmap().getBitmapEnd());
		coordIt = NULL;
	 }
	 
	 void findCbd() {
		if (i > -1 && i/2 < dim) {
		  if (i%2 == 0) {
			 coordIt->incInDir(i/2);
		  } else {
			 coordIt->decInDir(i/2);
		  }
		}

		++i;
		while(i/2 < dim) {
		  while(i/2 < dim && coordIt->odd(i/2)) {
			 i += 2;
		  }
		  if (i/2 < dim) {
			 if (i%2 == 0) {
				coordIt->decInDir(i/2);
			 } else {
				coordIt->incInDir(i/2);
			 }
			 if (coordIt->getBit()) {
				return;
			 } else {
				if (i%2 == 0) {
				  coordIt->incInDir(i/2);
				} else {
				  coordIt->decInDir(i/2);
				}
				++i;
			 }
		  }
		}
		toEnd();
	 }
  
	 void increment() {
		findCbd();
	 }
	 
	 int i;
	 Dim dim;
  };

  typedef typename boost::iterator_range<DimCoordIterator > DimCells;
  typedef typename boost::iterator_range<CbdCoordIterator > CbdCells;


  IteratorsImpl(SComplexRef _scomplex): scomplex(_scomplex) {}
  
  AllCells allCells() const;
  DimCells dimCells(const Dim& dim);
  BdCells bdCells(const Cell& cell) const;

  template<typename ImplT>
  CbdCells cbdCells(const CubCellProxy<ImplT>& cell) const {
	 return CbdCells(CbdCoordIterator(scomplex, cell), CbdCoordIterator(scomplex));
  }
  
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
inline typename CubSComplex::IteratorsImpl<isConst>::BdCells CubSComplex::IteratorsImpl<isConst>::bdCells(const Cell& cell) const {
  return BdCells(BdNumerator(scomplex, cell));
}

#endif
