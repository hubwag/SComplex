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
#include <boost/tuple/tuple.hpp>
#include <boost/foreach.hpp>

#include "util/Iterators.hpp"



template<template<typename Object, typename Color> class NeighboursModelT>
class SComplex {
public:
  typedef size_t Color;
  typedef size_t Dim;
  typedef size_t Id;
  class Cell;
  
  typedef std::vector<Cell> Cells;
  
  class Cell {

  public:
	 typedef typename SComplex::Color Color;
	 
	 explicit Cell(const Id& _id, const Color& _color): id(_id), color(_color) {}

	 void setColor(const Color& newColor);
	 const Color& getColor() const;

	 const Id& getId() const;
	 
  private:
	 Id id;
	 Color color;
  };

  typedef NeighboursModelT<boost::reference_wrapper<Cell>, Color> NeighboursModel;

private:
  
  template<bool isConst>
  class IteratorsImpl {	 
  public:

	 // typedef Util::Iterators::CollectionBeginEnd<Cells,
	 // 															boost::function<boost::indirect_iterator<typename Cells::iterator> (typename Cells::iterator) > > AllCells;
	 typedef Util::Iterators::CollectionBeginEnd<Cells> AllCells;

	 
	 typedef Util::Iterators::CollectionBeginEnd<NeighboursModel> DimCells;
	 typedef Util::Iterators::CollectionBeginEnd<NeighboursModel> BdCells;  
	 typedef Util::Iterators::CollectionBeginEnd<NeighboursModel> CbdCells;

	 explicit IteratorsImpl(SComplex& _complex): complex(_complex) {}
	 
	 AllCells allCells() const {
		// return AllCells(complex.cells,
		// 					 boost::lambda::bind(boost::lambda::constructor<boost::indirect_iterator<typename Cells::iterator> >(), boost::lambda::_1));
		return AllCells(complex.cells);
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

  typedef std::vector<boost::tuple<Id, Id, int> > KappaMap;
  
  SComplex(size_t colors, size_t _size, const KappaMap& kappaMap): boundaries(_size), coboundaries(_size) {
	 using namespace boost::lambda;
	 using boost::ref;
	 using boost::lambda::_1;
	 using boost::lambda::bind;
	 using boost::get;

	 cells.reserve(_size);
	 // Id id = 0;
	 // std::for_each(cells.begin(), cells.end(), _1 = bind(new_ptr<Cell>(), ref(id)++, 0));
	 for (Id id = 0; id < _size; ++id) {
		cells.push_back(Cell(id, 0));
	 }
		
	 //TODO call init with correct size
	 std::for_each(boundaries.begin(), boundaries.end(), bind(&NeighboursModel::init, _1, colors, 0));
	 std::for_each(coboundaries.begin(), coboundaries.end(), bind(&NeighboursModel::init, _1, colors, 0));

	 BOOST_FOREACH(KappaMap::value_type kappa, kappaMap) {
		Id coface = get<0>(kappa);
		Id face = get<1>(kappa);

		boundaries[coface].add(boost::ref(cells[face]), cells[face].getColor());
		coboundaries[face].add(boost::ref(cells[coface]), cells[coface].getColor());
	 }
	 
  }

  size_t cardinality() const { return cells.size(); }

  Cell& operator[](const Id id) {
	 return cells[id];
  }
  
  Iterators iterators() const { return Iterators(*this); }
  ConstIterators iterators() { return ConstIterators(*this); }

private:
  Cells cells;
  std::vector<std::pair<typename Cells::iterator, typename Cells::iterator> > cellsByDim; // pairs of begin/end in cells

  std::vector<NeighboursModel> boundaries, coboundaries; // (co)boundaries by cell id
};

template<template<typename Object, typename Color> class NeighboursModelT>
void SComplex<NeighboursModelT>::Cell::setColor(const Color& newColor) {
  Color oldColor = this->color;

  std::for_each(iterators().bdCells(*this).begin(), iterators().bdCells(*this).end(), 1);
			  
  this->color = newColor;
}

template<template<typename Object, typename Color> class NeighboursModelT>
inline const typename SComplex<NeighboursModelT>::Color& SComplex<NeighboursModelT>::Cell::getColor() const {
  return color;
}

template<template<typename Object, typename Color> class NeighboursModelT>
inline const typename SComplex<NeighboursModelT>::Id& SComplex<NeighboursModelT>::Cell::getId() const {
  return id;
}

#endif // _SCOMPLEX_HPP

