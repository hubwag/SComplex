#ifndef CUB_SCOMPLEX_CELL_HPP
#define CUB_SCOMPLEX_CELL_HPP

#include "CubSComplex.hpp"

class CubSComplex::Cell {
  friend class CellNumerator;
  friend class CbdNumerator;
  friend class CellDimNumerator;
  friend class CubSComplex;
  typedef BCubCellSet::BitCoordIterator::WordIterator WordIterator;

  BCubCellSet::BitIterator bitIt;
  BCubCellSet::BitCoordIterator* externalBitCoordIt;
  BCubCellSet::BitCoordIterator* internalBitCoordIt;
  Dim dim;
  
public:
  explicit Cell(const CubSComplex& s, const Dim& _dim = std::numeric_limits<Dim>::max()): bitIt(s.bCubCellSet),
																														externalBitCoordIt(NULL), internalBitCoordIt(NULL), dim(_dim) {
  }

  explicit Cell(BCubCellSet::BitCoordIterator* b): bitIt(*b), externalBitCoordIt(b), internalBitCoordIt(NULL), dim(b->ownDim())  {
  }

  explicit Cell(const BCubCellSet::BitIterator& b, const Dim& _dim): bitIt(b), externalBitCoordIt(NULL), internalBitCoordIt(NULL), dim(_dim)  {}
  
  BCubCellSet::BitCoordIterator& getBitCoordIt() {	
	 if ( !internalBitCoordIt && (externalBitCoordIt && bitIt == *((BCubCellSet::BitIterator*)externalBitCoordIt))) {
		return *externalBitCoordIt;
	 }
	 if (!internalBitCoordIt || (bitIt == *((BCubCellSet::BitIterator*)internalBitCoordIt))) {
		if (internalBitCoordIt) {
		  delete internalBitCoordIt;
		}
		internalBitCoordIt = new BCubCellSet::BitCoordIterator(bitIt);
		std::cout << "Created new internalBitCoordIt" << std::endl;
	 }
	 return *internalBitCoordIt;
  }

  ~Cell() {
	 if (internalBitCoordIt) {
		delete internalBitCoordIt;
		internalBitCoordIt = NULL;
	 }
  }
  
  Color getColor() const{
	 return bitIt.getBit() ? 1 : 2;
  }

  template<Color color>
  void setColor();

  void setColor(const Color& color) {
	 BOOST_ASSERT(color == 2);
	 bitIt.clearBit();
  }

	 
  Dim getDim() const {
	 return dim;
  }

  bool operator<(const Cell& b) const {
	 return bitIt < b.bitIt;
  }
  
};

template<>
inline void CubSComplex::Cell::setColor<2>() {
  setColor(2);
}

#endif
