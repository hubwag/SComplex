#ifndef CUB_SCOMPLEX_NUMERATORS_HPP
#define CUB_SCOMPLEX_NUMERATORS_HPP

#include "CubSComplex.hpp"

class CubSComplex::BasicNumerator: public BCubCellSet::BitCoordIterator {
protected:
  explicit BasicNumerator(const CubSComplex& s): BCubCellSet::BitCoordIterator(s.bCubCellSetCR()), current(s) {}
  
public:
  void toEnd() {
	 this->wIt=const_cast<WordIterator>(this->getBitmap().getBitmapEnd());
  }
		
  Cell& Current(){
	 return (current = Cell(*this));
  }

  bool operator==(const BasicNumerator& o) const {
	 return this->wIt == o.wIt;
  }

  bool operator!=(const BasicNumerator& o) const {
	 return this->wIt != o.wIt;
  }

private:
  Cell current; // to remove, for walid reference only
};

class CubSComplex::CellNumerator: public BasicNumerator{
  
public:
	 typedef Cell value_type;
  
  CellNumerator(const CubSComplex& s): BasicNumerator(s) {
  }
  
	 bool MoveNext(){
		++(*this);
		moveToFirstPixel();
		return wIt < getBitmap().getBitmapEnd();
	 }
};

class CubSComplex::CellDimNumerator: public BasicNumerator {
public:
	 typedef Cell value_type;
		
  CellDimNumerator(const CubSComplex& s,int d): BasicNumerator(s),dim(d) {
	 --(*this);
  }

  bool MoveNext(){
	 for(;;){
		++(*this);
		moveToFirstPixel();
		if(wIt == getBitmap().getBitmapEnd()) return false;
		if(ownDim()==dim) return true;
	 }
  }
	 
protected:
  int dim;
};

class CubSComplex::CbdNumerator: public BasicNumerator {
public:
  typedef Cell value_type;
		
  CbdNumerator(const CubSComplex& s, const Cell& c):BasicNumerator(s), center(c),i(0),downDir(true), dim(c.embDim()){
  }
		
  bool MoveNext(){
	 while(i < dim){
		// process only directions in which cell is degenerate
		if(!downDir || !center.odd(i)){
		  ((BCubCellSet::BitCoordIterator&)(*this)) = center;
		  // First check the bottom face
		  if(downDir){
			 decInDir(i);
			 downDir=false;
			 // and now go to the top face
		  }else{
			 incInDir(i);;
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
  const Cell& center;
  int i;
  bool downDir;
  const int dim;
};

class CubSComplex::BdNumerator: public BasicNumerator{
public:
	 typedef Cell value_type;
		
  BdNumerator(const CubSComplex& s, const Cell& c): BasicNumerator(s), center(c),i(0),downDir(true), dim(c.embDim()) {
	 }
		
  bool MoveNext(){
	 while(i < dim){
		// process only directions in which cell is degenerate
		if(!downDir || center.odd(i)){
		  ((BCubCellSet::BitCoordIterator&)(*this)) = center;
		  // First check the bottom face
		  if(downDir){
			 decInDir(i);
			 downDir=false;
			 // and now go to the top face
		  }else{
			 incInDir(i);;
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
  const Cell& center;
  int i;
  bool downDir;
  const int dim;
};

#endif // CUB_SCOMPLEX_NUMERATORS_HPP
