#ifndef _SCOMPLEX_HPP
#define _SCOMPLEX_HPP

#include <algorithm>
#include <vector>
#include <list>
#include <boost/ref.hpp>
#include <boost/mpl/if.hpp>

#include "util/Iterators.hpp"


template<typename CellT>
class ColorListNeighbourModel {
  typedef CellT Cell;
  typedef typename std::list<typename Cell::CellRef> CellsList;
  typedef typename std::vector<CellsList> CellsByColor;  
  typedef typename std::vector<typename Cell::CellRef> CellsVect;
  
public:

  typedef typename CellsList::iterator NeighbourId;
  typedef typename CellsVect::iterator NeighboursIterator;

  
  explicit ColorListNeighbourModel(size_t colors): neighboursByColor(colors) {}

  NeighboursIterator begin() {
	 return neighbours.begin();
  }

  NeighboursIterator end() {
	 return neighbours.end();
  }

  NeighbourId add(Cell& cell, const typename Cell::Color& color) {
	 return neighboursByColor[color].insert(neighboursByColor[color].end(), cell);	 
  }

  void changeColor(NeighbourId it, const typename Cell::Color& oldColor, const typename Cell::Color& newColor) {
	 neighboursByColor[newColor].splice(neighboursByColor[newColor].end(), neighboursByColor[oldColor], it);
  }

private:
  CellsVect neighbours;
  CellsByColor neighboursByColor;
};

template<typename NeighboursModelT>
class SComplex {
public:
  typedef int Color;
  typedef int Dim;
  typedef int Id;
  typedef NeighboursModelT NeighboursModel;

  class Cell;
private:
  
  template<bool isConst>
  class IteratorsImpl {	 
  public:

	 typedef Util::Iterators::CollectionBeginEnd<NeighboursModel> AllCells;
	 typedef Util::Iterators::CollectionBeginEnd<NeighboursModel> DimCells;
	 typedef Util::Iterators::CollectionBeginEnd<NeighboursModel> BdCells;  
	 typedef Util::Iterators::CollectionBeginEnd<NeighboursModel> CbdCells;

	 explicit IteratorsImpl(SComplex& _complex): complex(_complex) {}
	 
	 AllCells allCells() const;
	 DimCells dimCells(const Dim& dim) const;
	 
	 BdCells bdCells(const Cell& cell) const {
		return Util::Iterators::CollectionBeginEnd<NeighboursModel>(complex.boundaries[cell.getId()]);
	 }
	 
	 CbdCells cbdCells(const Cell& cell) const {
		return Util::Iterators::CollectionBeginEnd<NeighboursModel>(complex.coboundaries[cell.getId()]);
	 }

  private:
	 SComplex& complex;
	 
  };
  
public:
  
  typedef IteratorsImpl<false> Iterators;
  typedef IteratorsImpl<true> ConstIterators;

  typedef IteratorsImpl<false> ColoredIterators;
  typedef IteratorsImpl<true> ColoredConstIterators;
  
  class Cell {

  public:
	 typedef typename SComplex::Color Color;
	 typedef boost::reference_wrapper<Cell> CellRef;
	 
	 explicit Cell(const Id& _id, size_t colors): id(_id) {}

	 void setColor(const Color& newColor);
	 
  private:
	 Color color;
	 const Id id;
  };


  ConstIterators iterators() const { return Iterators(*this); }
  Iterators iterators() { return ConstIterators(*this); }

private:
  typedef std::vector<Cell> Cells;

  Cells cells; // sorted by dim
  std::vector<std::pair<typename Cells::iterator, typename Cells::iterator> > cellsByDim; // pairs of begin/end in cells

  std::vector<NeighboursModel> boundary, coboundary; // (co)boundaries by cell id
};

template<typename NeighboursModelT>
void SComplex<NeighboursModelT>::Cell::setColor(const Color& newColor) {
  Color oldColor = this->color;

  std::for_each(iterators().bdCells(*this).begin(), iterators().bdCells(*this).end(), 1);
			  
  this->color = newColor;
}


#endif // _SCOMPLEX_HPP

