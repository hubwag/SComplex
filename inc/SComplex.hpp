#ifndef _SCOMPLEX_HPP
#define _SCOMPLEX_HPP

#include <algorithm>
#include <vector>
#include <list>
#include <functional>

#include <boost/ref.hpp>
#include <boost/mpl/if.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/construct.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/function.hpp>

#include "util/Iterators.hpp"



template<template<typename Object> class NeighboursModelT>
class SComplex {
public:
  typedef int Color;
  typedef int Dim;
  typedef int Id;
  class Cell;
  
  typedef std::vector<Cell*> Cells;
  
  class Cell {

  public:
	 typedef typename SComplex::Color Color;
	 
	 explicit Cell(const Id& _id, size_t colors): id(_id) {}

	 void setColor(const Color& newColor);
	 
  private:
	 Color color;
	 const Id id;
  };

  typedef NeighboursModelT<Cell> NeighboursModel;

private:
  
  template<bool isConst>
  class IteratorsImpl {	 
  public:

	 typedef Util::Iterators::CollectionBeginEnd<Cells,
																boost::function<boost::indirect_iterator<typename Cells::iterator> (typename Cells::iterator) > > AllCells;
	 
	 typedef Util::Iterators::CollectionBeginEnd<NeighboursModel> DimCells;
	 typedef Util::Iterators::CollectionBeginEnd<NeighboursModel> BdCells;  
	 typedef Util::Iterators::CollectionBeginEnd<NeighboursModel> CbdCells;

	 explicit IteratorsImpl(SComplex& _complex): complex(_complex) {}
	 
	 AllCells allCells() const {
		return AllCells(complex.cells,
							 boost::lambda::bind(boost::lambda::constructor<boost::indirect_iterator<typename Cells::iterator> >(), boost::lambda::_1));
	 }

	 
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

  SComplex(size_t colors, size_t _size, std::vector<std::pair<Id, Id> > _cells): cells(_size) {
	 using namespace boost::lambda;
	 using boost::ref;
	 using boost::lambda::_1;
	 using boost::lambda::bind;
		
	 Id id = 0;
	 std::for_each(cells.begin(), cells.end(), _1 = bind(new_ptr<Cell>(), ref(id)++, colors));
  }

  size_t cardinality() const { return cells.size(); }
  
  Iterators iterators() const { return Iterators(*this); }
  ConstIterators iterators() { return ConstIterators(*this); }

private:
  Cells cells; // sorted by dim
  std::vector<std::pair<typename Cells::iterator, typename Cells::iterator> > cellsByDim; // pairs of begin/end in cells

  std::vector<NeighboursModel> boundary, coboundary; // (co)boundaries by cell id
};

template<template<typename Object> class NeighboursModelT>
void SComplex<NeighboursModelT>::Cell::setColor(const Color& newColor) {
  Color oldColor = this->color;

  std::for_each(iterators().bdCells(*this).begin(), iterators().bdCells(*this).end(), 1);
			  
  this->color = newColor;
}


#endif // _SCOMPLEX_HPP

