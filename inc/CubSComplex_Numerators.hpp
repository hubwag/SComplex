#ifndef CUB_SCOMPLEX_NUMERATORS_HPP
#define CUB_SCOMPLEX_NUMERATORS_HPP

#include "CubSComplex.hpp"
#include <limits>
#include <boost/shared_ptr.hpp>

class CubSComplex::BasicNumerator {
protected:
  explicit BasicNumerator(const CubSComplex& s): bitCoordIt(s.bCubCellSet) {}
  
public:
  void toEnd() {
	 this->bitCoordIt.wIt=const_cast<BCubCellSet::BitCoordIterator::WordIterator>(this->bitCoordIt.getBitmap().getBitmapEnd());
  }
		
  CubCellProxy<BitCoordPtrCellImpl> Current(){
	 return CubCellProxy<BitCoordPtrCellImpl>(&bitCoordIt);
  }

  bool operator==(const BasicNumerator& o) const {
	 return bitCoordIt.wIt == o.bitCoordIt.wIt;
  }

  bool operator!=(const BasicNumerator& o) const {
	 return this->bitCoordIt.wIt != o.bitCoordIt.wIt;
  }

protected:
  BCubCellSet::BitCoordIterator bitCoordIt;
  
};

class CubSComplex::CellNumerator: public BasicNumerator{
  
public:
	 typedef Cell value_type;
  
  CellNumerator(const CubSComplex& s): BasicNumerator(s) {
  }
  
	 bool MoveNext(){
		++(bitCoordIt);
		bitCoordIt.moveToFirstPixel();
		return bitCoordIt.wIt < bitCoordIt.getBitmap().getBitmapEnd();
	 }
};

class CubSComplex::CellDimNumerator: public BasicNumerator {
public:
	 typedef Cell value_type;
		
  CellDimNumerator(const CubSComplex& s,int d): BasicNumerator(s),dim(d) {
	 --(bitCoordIt);
  }

  bool MoveNext(){
	 for(;;){
		++(bitCoordIt);
		bitCoordIt.moveToFirstPixel();
		if(bitCoordIt.wIt == bitCoordIt.getBitmap().getBitmapEnd()) return false;
		if(bitCoordIt.ownDim()==dim) return true;
	 }
  }
	 
protected:
  int dim;
};

class CubSComplex::CbdNumerator: public BasicNumerator {
public:
  typedef Cell value_type;

  template<typename ImplT>
  CbdNumerator(const CubSComplex& s, const CubCellProxy<ImplT>& c):BasicNumerator(s), center(c.getBitCoordIt()),i(0),downDir(true), dim(s.bCubCellSet.embDim()){
  }
		
  bool MoveNext(){
	 while(i < dim){
		// process only directions in which cell is degenerate
		if(!downDir || !center.odd(i)){
		  ((BCubCellSet::BitCoordIterator&)(bitCoordIt)) = center;
		  // First check the bottom face
		  if(downDir){
			 bitCoordIt.decInDir(i);
			 downDir=false;
			 // and now go to the top face
		  }else{
			 bitCoordIt.incInDir(i);;
			 downDir=true;
			 ++i;
		  }
		  return true;
		}else{
		  ++i;
		}
	 }
	 //cCell=false;
	 toEnd();
	 return false;
  }
protected:
  const BCubCellSet::BitCoordIterator center;
  int i;
  bool downDir;
  const int dim;
};

class CubSComplex::BdNumerator: public BasicNumerator{
public:
	 typedef Cell value_type;

  template<typename ImplT>
  BdNumerator(const CubSComplex& s, const CubCellProxy<ImplT>& c): BasicNumerator(s), center(c.getBitCoordIt()),i(0),downDir(true), dim(s.bCubCellSet.embDim()) {
	 }
		
  bool MoveNext(){
	 while(i < dim){
		// process only directions in which cell is degenerate
		if(!downDir || center.odd(i)){
		  ((BCubCellSet::BitCoordIterator&)(bitCoordIt)) = center;
		  // First check the bottom face
		  if(downDir){
			 bitCoordIt.decInDir(i);
			 downDir=false;
			 // and now go to the top face
		  }else{
			 bitCoordIt.incInDir(i);;
			 downDir=true;
			 ++i;
		  }
		  return true;
		}else{
		  ++i;
		}
	 }
	 //cCell=false;
	 toEnd();
	 return false;
  }

protected:
  const BCubCellSet::BitCoordIterator center;
  int i;
  bool downDir;
  const int dim;
};

#endif // CUB_SCOMPLEX_NUMERATORS_HPP
