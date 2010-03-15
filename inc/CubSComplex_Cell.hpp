#ifndef CUB_SCOMPLEX_CELL_HPP
#define CUB_SCOMPLEX_CELL_HPP

#include "CubSComplex.hpp"
#include <boost/shared_ptr.hpp>

template<typename CellImplT>
class CubSComplex::CellProxy {
  CellImplT impl;
public:
  CellProxy(const CellImplT& _impl): impl(_impl) {}
  
  Color getColor() const{
	 return impl.getColor();
  }

  template<Color color>
  void setColor() {
	 impl.template setColor<color>();
  }

  void setColor(const Color& color) {
	 impl.setColor(color);
  }
	 
  Dim getDim() const {
	 return impl.getDim();
  }

  const BCubCellSet::BitCoordIterator& getBitCoordIt() const {
	 return impl.getBitCoordIt();
  }

  BCubCellSet::BitCoordIterator& getBitCoordIt() {
	 return impl.getBitCoordIt();
  }

  bool operator<(const CellProxy& b) const {
	 return impl < b.impl;
  }
};

class CubSComplex::BitCoordCellImpl {
  Dim dim;

protected:
  BCubCellSet::BitCoordIterator bitIt;
  
public:
  BitCoordCellImpl(const BCubCellSet::BitCoordIterator& b): bitIt(b), dim(b.ownDim()) {}

  BitCoordCellImpl(const CubSComplex& s): bitIt(s.bCubCellSet), dim(std::numeric_limits<Dim>::max()) {}
  
  Color getColor() const{
	 return bitIt.getBit() ? 1 : 2;
  }

  template<Color color>
  void setColor() {
	 setColor(color);
  }

  void setColor(const Color& color) {
	 BOOST_ASSERT(color == 2);
	 BOOST_ASSERT(getColor() == 1);
	 bitIt.clearBit();
  }
	 
  Dim getDim() const {
	 return dim;
  }

  const BCubCellSet::BitCoordIterator& getBitCoordIt() const {
	 return bitIt;
  }

  bool operator<(const BitCoordCellImpl& b) const {
	 return bitIt < b.bitIt;
  }
  
};

class CubSComplex::DynamicBitCoordCellImpl: public BitCoordCellImpl  {

public:
 DynamicBitCoordCellImpl(const BCubCellSet::BitCoordIterator& b): BitCoordCellImpl(b) {}

 DynamicBitCoordCellImpl(const CubSComplex& s): BitCoordCellImpl(s.bCubCellSet) {}
  
  Dim getDim() const {
	 return bitIt.ownDim();
  }

  BCubCellSet::BitCoordIterator& getBitCoordIt() {
	 return bitIt;
  }

};

class CubSComplex::BitCoordPtrCellImpl {
protected:
  BCubCellSet::BitCoordIterator* bitIt;
  
public:
  BitCoordPtrCellImpl(BCubCellSet::BitCoordIterator* b): bitIt(b) {}

  BitCoordPtrCellImpl(const CubSComplex& s): bitIt(NULL) {}

  template<typename ImplT>
  BitCoordPtrCellImpl(CellProxy<ImplT>& b): bitIt(& b.getBitCoordIt()) {}

  Color getColor() const{
	 return bitIt->getBit() ? 1 : 2;
  }

  template<Color color>
  void setColor() {
	 setColor(color);
  }

  void setColor(const Color& color) {
	 BOOST_ASSERT(color == 2);
	 BOOST_ASSERT(getColor() == 1);
	 bitIt->clearBit();
  }
	 
  Dim getDim() const {
	 return bitIt->ownDim();
  }

  const BCubCellSet::BitCoordIterator& getBitCoordIt() const {
	 return *bitIt;
  }

  BCubCellSet::BitCoordIterator& getBitCoordIt() {
	 return *bitIt;
  }

  bool operator<(const BitCoordPtrCellImpl& b) const {
	 return *bitIt < *b.bitIt;
  }
  
};

#endif
