#ifndef CUB_SCOMPLEX_CELL_HPP
#define CUB_SCOMPLEX_CELL_HPP

#include "CubSComplex.hpp"
#include "CellProxy.hpp"
#include <boost/shared_ptr.hpp>

template<typename CellImplT>
class CubSComplex::CubCellProxy: public CellProxy<CellImplT> {
protected:
  using CellProxy<CellImplT>::impl;
public:
  CubCellProxy(const CellImplT& _impl): CellProxy<CellImplT>(_impl) {}

  template<typename CubCellImplT2>
  CubCellProxy(const CubCellProxy<CubCellImplT2>& other): CellProxy<CellImplT>(other.impl) {}
  
  const BCubCellSet::BitCoordIterator& getBitCoordIt() const {
	 return impl.getBitCoordIt();
  }

  BCubCellSet::BitCoordIterator& getBitCoordIt() {
	 return impl.getBitCoordIt();
  }

};

class CubSComplex::CellImpl {
public:
  typedef CubSComplex::Dim Dim;
  typedef CubSComplex::Color Color;

};

class CubSComplex::BitCoordCellImpl: public CellImpl {
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

class CubSComplex::BitCoordPtrCellImpl: public CellImpl {
protected:
  BCubCellSet::BitCoordIterator* bitIt;
  
public:
  BitCoordPtrCellImpl(BCubCellSet::BitCoordIterator* b): bitIt(b) {}

  BitCoordPtrCellImpl(const CubSComplex& s): bitIt(NULL) {}

  template<typename ImplT>
  BitCoordPtrCellImpl(CubCellProxy<ImplT>& b): bitIt(& b.getBitCoordIt()) {}

  operator BitCoordCellImpl () const {
	 return BitCoordCellImpl(*bitIt);
  }
  
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
