#ifndef CUB_SCOMPLEX_ITERATORS_HPP
#define CUB_SCOMPLEX_ITERATORS_HPP

#include "CubSComplex.hpp"
#include "../../util/Iterators.hpp"

#include <boost/mpl/if.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/range.hpp>
#include <boost/shared_ptr.hpp>
#include <functional>


template<int DIM>
template<bool isConst>
class CubSComplex<DIM>::IteratorsImpl {
  typedef typename boost::mpl::if_c<isConst, const CubSComplex&, CubSComplex&>::type SComplexRef;
public:

  typedef CellProxy<BitCoordPtrCellImpl> CellType;


  template<typename Derived>
  class CoordIterator: public boost::iterator_facade<Derived, CellType,
																	  boost::forward_traversal_tag, CellType>
  {
  
  public:
	 typedef typename CubSComplex::BCubCellSet::BitCoordIterator BitCoordIterator;
	 typedef typename CubSComplex::BCubCellSet::BitIterator BitIterator;
  
	 CoordIterator():  initializedCoordIt(false), coordIt(NULL) {}
	 
	 CoordIterator(const CoordIterator& other) {
		if (other.coordIt) {
		  this->coordIt = new (coordItMem) BitCoordIterator(*other.coordIt);
		} else {
		  this->coordIt = NULL;
		}
		this->initializedCoordIt = other.initializedCoordIt;
	 }

  protected:
	 BitCoordIterator* getCoordIt() {
		if (! initializedCoordIt) {
		  initializedCoordIt = true;
		  ((Derived*)this)->initCoordIt((BitCoordIterator*)coordItMem);
		}
		return coordIt;
	 }

	 BitCoordIterator* getCoordIt() const {
	  	return const_cast<CoordIterator*>(this)->getCoordIt();
	 }
	 
  protected:
	 friend class boost::iterator_core_access;
	 bool initializedCoordIt;
	 mutable BitCoordIterator* coordIt;
	 char coordItMem[sizeof(BitCoordIterator)];
	 
	 
	 bool equal(const Derived& other) const {
		BitCoordIterator* it1 = this->getCoordIt();
		BitCoordIterator* it2 = other.getCoordIt();
		return (it1 == NULL && it2 == NULL) ||
		  ( (it1 != NULL && it2 != NULL)
			 && (BitIterator&)(*it1) == ((BitIterator&)(*it2)));
	 }

	 CellType dereference() const {
		return CellType(getCoordIt());
	 }
  };

  class DimCoordIterator: public  CoordIterator<DimCoordIterator> {
	 typedef typename CubSComplex::BCubCellSet::BitCoordIterator BitCoordIterator;
  public:
  
	 DimCoordIterator(const CubSComplex& s, const Dim& _dim): CoordIterator<DimCoordIterator>(), complex(&s), dim(_dim) { //begin iterator
	 }

	 DimCoordIterator(const CubSComplex& s): CoordIterator<DimCoordIterator>(), complex(NULL), dim(std::numeric_limits<Dim>::max()) { //end constructor
	 }
  
  private:	 
	 friend class boost::iterator_core_access;
	 friend class CoordIterator<DimCoordIterator>;
	 
	 void initCoordIt(BitCoordIterator* allocatedMemory) {
		if (complex == NULL) {
		  coordIt = NULL;
		} else {
		  coordIt = new (allocatedMemory) BitCoordIterator(complex->bCubCellSet);
		  findDim();
		}
	 }
	 
	 void findDim() {
		while (coordIt->wIt < coordIt->getBitmap().end().wIt) {
		  if (!coordIt->findPoint())
		    break;
		  if(coordIt->getBit() && coordIt->ownDim()==dim)
			 return;
		  ++(*coordIt);
		}
		coordIt = NULL;
	 }
  
	 void increment() {
		++(*CoordIterator<DimCoordIterator>::getCoordIt());
		findDim();
	 }

	 using CoordIterator<DimCoordIterator>::coordIt;	 
	 const CubSComplex* complex;
	 Dim dim;  
  };

  class AllCellsCoordIterator: public  CoordIterator<AllCellsCoordIterator> {
	 typedef typename CubSComplex::BCubCellSet::BitCoordIterator BitCoordIterator;
  public:
  
    AllCellsCoordIterator(const CubSComplex& s, int): CoordIterator<AllCellsCoordIterator>(), complex(&s) { //begin iterator
	 }

    AllCellsCoordIterator(const CubSComplex& s): CoordIterator<AllCellsCoordIterator>(), complex(NULL) { //end constructor
    }
  
  private:	 
	 friend class boost::iterator_core_access;
	 friend class CoordIterator<AllCellsCoordIterator>;
	 
	 void initCoordIt(BitCoordIterator* allocatedMemory) {
		if (complex == NULL) {
		  coordIt = NULL;
		} else {
		  coordIt = new (allocatedMemory) BitCoordIterator(complex->bCubCellSet);
		  findFirst();
		}
	 }
	 
	 void findFirst() {
		while (coordIt->wIt < coordIt->getBitmap().end().wIt) {
		  if (!coordIt->findPoint())
		    break;
		  if(coordIt->getBit())
		    return;
		  ++(*coordIt);
		}
		coordIt = NULL;
	 }
  
	 void increment() {
		++(*CoordIterator<AllCellsCoordIterator>::getCoordIt());
		findFirst();
	 }

	 using CoordIterator<AllCellsCoordIterator>::coordIt;	 
	 const CubSComplex* complex;
  };

  class CbdCoordIterator: public  CoordIterator<CbdCoordIterator> {
	 typedef typename CubSComplex::BCubCellSet::BitCoordIterator BitCoordIterator;
  public:
	 CbdCoordIterator(): CoordIterator<CbdCoordIterator>() {}

	 template<typename ImplT>
	 CbdCoordIterator(CubSComplex& s, const CellProxy<ImplT>& c): CoordIterator<CbdCoordIterator>(), 
								      sourceCoordIt(&c.getBitCoordIt()),
								      i(-1), dim(s.bCubCellSet.embDim()) { //begin iterator
	 }

	 CbdCoordIterator(CubSComplex& s): CoordIterator<CbdCoordIterator>(), sourceCoordIt(NULL),
					   i(0), dim(0) { //end constructor
	   toEnd();
	 }
  
  private:	 
	 friend class boost::iterator_core_access;
	 friend class  CoordIterator<CbdCoordIterator>;
	 
	 void initCoordIt(BitCoordIterator* allocatedMemory) {
		if (sourceCoordIt == NULL) {
		  coordIt = NULL;
		} else {
		  coordIt = new (allocatedMemory) BitCoordIterator(*sourceCoordIt);
		  findCbd();
		}
	 }

	 void toEnd() {
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
		CoordIterator<CbdCoordIterator>::getCoordIt();
		findCbd();
	 }

	 using CoordIterator<CbdCoordIterator>::coordIt;
	 const BitCoordIterator* sourceCoordIt;	 
	 int i;
	 Dim dim;
  };

  class BdCoordIterator: public  CoordIterator<BdCoordIterator> {
	 typedef typename CubSComplex::BCubCellSet::BitCoordIterator BitCoordIterator;
  public:
	 BdCoordIterator(): CoordIterator<BdCoordIterator>() {}

	 template<typename ImplT>
	 BdCoordIterator(CubSComplex& s, const CellProxy<ImplT>& c): CoordIterator<BdCoordIterator>(), 
								      sourceCoordIt(&c.getBitCoordIt()),
								      i(-1), dim(s.bCubCellSet.embDim()) { //begin iterator
	 }

	 BdCoordIterator(CubSComplex& s): CoordIterator<BdCoordIterator>(), sourceCoordIt(NULL),
					   i(0), dim(0) { //end constructor
	   toEnd();
	 }
  
  private:	 
	 friend class boost::iterator_core_access;
	 friend class  CoordIterator<BdCoordIterator>;
	 
	 void initCoordIt(BitCoordIterator* allocatedMemory) {
		if (sourceCoordIt == NULL) {
		  coordIt = NULL;
		} else {
		  coordIt = new (allocatedMemory) BitCoordIterator(*sourceCoordIt);
		  findBd();
		}
	 }

	 void toEnd() {
		coordIt = NULL;
	 }
	 
	 void findBd() {
		if (i > -1 && i/2 < dim) {
		  if (i%2 == 0) {
			 coordIt->incInDir(i/2);
		  } else {
			 coordIt->decInDir(i/2);
		  }
		}

		++i;
		while(i/2 < dim) {
		  while(i/2 < dim && !coordIt->odd(i/2)) {
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
		CoordIterator<BdCoordIterator>::getCoordIt();
		findBd();
	 }

	 using CoordIterator<BdCoordIterator>::coordIt;
	 const BitCoordIterator* sourceCoordIt;	 
	 int i;
	 Dim dim;
  };

  typedef typename boost::iterator_range<AllCellsCoordIterator > AllCells;
  typedef typename boost::iterator_range<DimCoordIterator > DimCells;
  typedef typename boost::iterator_range<CbdCoordIterator > CbdCells;
  typedef typename boost::iterator_range<BdCoordIterator > BdCells;

  IteratorsImpl(SComplexRef _scomplex): scomplex(_scomplex) {}
  

  AllCells allCells() {return AllCells(AllCellsCoordIterator(scomplex, 0), AllCellsCoordIterator(scomplex));}
  DimCells dimCells(const Dim& dim) {return DimCells(DimCoordIterator(scomplex, dim), DimCoordIterator(scomplex));}


  template<typename ImplT>
  CbdCells cbdCells(const CellProxy<ImplT>& cell) const {
	 return CbdCells(CbdCoordIterator(scomplex, cell), CbdCoordIterator(scomplex));
  }

  template<typename ImplT>
  BdCells bdCells(const CellProxy<ImplT>& cell) const {
	 return BdCells(BdCoordIterator(scomplex, cell), BdCoordIterator(scomplex));
  }
  
private:
  SComplexRef scomplex;
};




  


#endif
