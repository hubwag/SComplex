#ifndef CUB_SCOMPLEX_CELL_HPP
#define CUB_SCOMPLEX_CELL_HPP

#include "CubSComplex.hpp"
#include "CellProxy.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>
#include <boost/type_traits.hpp>

template<typename CellImplT, typename Enable>
class CubSComplex::CubCellProxy: public Enable {
};

template<typename CellImplT>
class CubSComplex::CubCellProxy<CellImplT,
										  typename boost::enable_if<typename boost::is_base_of<CubSComplex::CellImpl,
																												 typename boost::remove_pointer<CellImplT>::type >::type
																			 >::type > : public CellProxy<CellImplT> {
protected:
  using CellProxy<CellImplT>::impl;
public:
  template<typename CubCellImplT2>
  CubCellProxy(const CubCellImplT2& _impl): CellProxy<CellImplT>(_impl) {}

  template<typename CubCellImplT2>
  CubCellProxy(const CubCellProxy<CubCellImplT2>& other): CellProxy<CellImplT>((CellImplT)other.impl) {}
  
  const BCubCellSet::BitCoordIterator& getBitCoordIt() const {
	 return CellProxy<CellImplT>::getImpl()->getBitCoordIt();
  }

  BCubCellSet::BitCoordIterator& getBitCoordIt() {
  	 return CellProxy<CellImplT>::getImpl()->getBitCoordIt();
  }

};


class CubSComplex::CellImpl {
public:
  typedef CubSComplex::Dim Dim;
  typedef CubSComplex::Color Color;

};

class CubSComplex::BitCoordCellImpl: public CellImpl {
protected:
  BCubCellSet::BitCoordIterator bitIt;
private:
  //Dim dim;
public:
  
  BitCoordCellImpl(const BCubCellSet::BitCoordIterator& b): bitIt(b)/*, dim(b.ownDim())*/ {}

  BitCoordCellImpl(const CubSComplex& s): bitIt(s.bCubCellSet)/*, dim(std::numeric_limits<Dim>::max())*/ {}
  
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
	 //return dim;
	 return bitIt.ownDim();
  }

  const BCubCellSet::BitCoordIterator& getBitCoordIt() const {
	 return bitIt;
  }

  BCubCellSet::BitCoordIterator& getBitCoordIt() {
	 return bitIt;
  }

  bool operator<(const BitCoordCellImpl& b) const {
	 return bitIt < b.bitIt;
  }
  
};

class CubSComplex::BitCoordPtrCellImpl: public CellImpl {
protected:
  mutable BCubCellSet::BitCoordIterator* bitIt;
  
public:
  BitCoordPtrCellImpl(BCubCellSet::BitCoordIterator* b): bitIt(b) {}

  explicit BitCoordPtrCellImpl(const CubSComplex& s): bitIt(NULL) {}

  template<typename ImplT>
  BitCoordPtrCellImpl(CubCellProxy<ImplT>& b): bitIt(& b.getBitCoordIt()) {}

  operator const BCubCellSet::BitCoordIterator& () const {
	 return *bitIt;
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

#endif
