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
	 
	 explicit Cell(SComplex* _complex, const Id& _id, const Color& _color): id(_id), color(_color), complex(_complex) {}

	 void setColor(const Color& newColor);
	 const Color& getColor() const;

	 const Id& getId() const;
	 
  private:
	 Id id;
	 Color color;
	 SComplex* complex;
  };

  typedef NeighboursModelT<boost::reference_wrapper<Cell>, Color> NeighboursModel;

private:
  
  template<bool isConst>
  class IteratorsImpl {
	 
	 struct NeighboursModelCellExtractor: public std::unary_function<typename boost::mpl::if_c<isConst,
																															 const typename NeighboursModel::NeighbourLink&, 
																															 typename NeighboursModel::NeighbourLink&
																															 >::type,
																						  typename boost::mpl::if_c<isConst,
																															 const Cell&,
																															 Cell&
																															 >::type>
	 {
		typedef std::unary_function<typename boost::mpl::if_c<isConst,
																				const typename NeighboursModel::NeighbourLink&, 
																				typename NeighboursModel::NeighbourLink&
																				>::type,
											 typename boost::mpl::if_c<isConst,
																				const Cell&,
																				Cell&
																				>::type> BaseClass;
		typedef typename BaseClass::result_type result_type;
		typedef typename BaseClass::argument_type argument_type;
		
		result_type operator()(const typename NeighboursModel::NeighbourLink&  link) const {
		  return link.objectRef.get();
		}


	 };

  public:

	 typedef typename boost::mpl::if_c<isConst, const SComplex&, SComplex&>::type Complex;
	 
	 typedef Util::Iterators::CollectionBeginEnd<isConst, Cells, std::_Identity<Cell> > AllCells;	 
	 typedef Util::Iterators::CollectionBeginEnd<isConst, NeighboursModel, NeighboursModelCellExtractor> DimCells;
	 typedef Util::Iterators::CollectionBeginEnd<isConst, NeighboursModel, NeighboursModelCellExtractor> BdCells;  
	 typedef Util::Iterators::CollectionBeginEnd<isConst, NeighboursModel, NeighboursModelCellExtractor> CbdCells;

	 explicit IteratorsImpl(Complex _complex): complex(_complex) {}
	 
	 AllCells allCells() const {
		return AllCells(complex.cells);
	 }
	 
	 DimCells dimCells(const Dim& dim) const;
	 
	 BdCells bdCells(const Cell& cell) const {
		return BdCells(complex.boundaries[cell.getId()]);
	 }
	 
	 CbdCells cbdCells(const Cell& cell) const {
		return CbdCells(complex.coboundaries[cell.getId()]);
	 }

  private:
	 Complex complex;
	 
  };

  template<bool isConst>
  class ColoredIteratorsImpl {

  public:
	 typedef IteratorsImpl<isConst> Iterators;

	 template<Color color>
	 class Color {
	 public:
		typedef IteratorsImpl<isConst> Iterators;
	 };
  };

public:
  
  typedef IteratorsImpl<false> Iterators;
  typedef IteratorsImpl<true> ConstIterators;

  
  typedef ColoredIteratorsImpl<false> ColoredIterators;
  typedef ColoredIteratorsImpl<true> ColoredConstIterators;

  typedef std::vector<boost::tuple<Id, Id, int> > KappaMap;
  
  SComplex(size_t colors, size_t _size, const KappaMap& kappaMap = KappaMap()): boundaries(_size), coboundaries(_size) {
	 using namespace boost::lambda;
	 using boost::ref;
	 using boost::lambda::_1;
	 using boost::lambda::bind;
	 using boost::get;

	 cells.reserve(_size);
	 // Id id = 0;
	 // std::for_each(cells.begin(), cells.end(), _1 = bind(new_ptr<Cell>(), ref(id)++, 0));
	 for (Id id = 0; id < _size; ++id) {
		cells.push_back(Cell(this, id, 0));
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
  
  Iterators iterators() { return Iterators(*this); }
  ConstIterators iterators() const { return ConstIterators(*this); }

  typename ColoredConstIterators::Iterators iterators(const Color& color) const;
  typename ColoredIterators::Iterators iterators(const Color& color);

  template<Color color>
  typename ColoredIterators::template Color<color>::Iterators iterators();

  template<Color color>
  typename ColoredConstIterators::template Color<color>::Iterators iterators() const;

private:
  Cells cells;
  
  std::vector<std::pair<typename Cells::iterator, typename Cells::iterator> > cellsByDim; // pairs of begin/end in cells

  std::vector<NeighboursModel> boundaries, coboundaries; // (co)boundaries by cell id
};

template<template<typename Object, typename Color> class NeighboursModelT>
inline void SComplex<NeighboursModelT>::Cell::setColor(const Color& newColor) {
  Color oldColor = this->color;

  BOOST_FOREACH(typename SComplex<NeighboursModelT>::Iterators::BdCells::iterator::value_type bd,
					 complex->iterators().bdCells(*this)) {
	 
  }
			  
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

template<template<typename Object, typename Color> class NeighboursModelT>
inline typename SComplex<NeighboursModelT>::ColoredConstIterators::Iterators SComplex<NeighboursModelT>::iterators(const Color& color) const {

}

template<template<typename Object, typename Color> class NeighboursModelT>
inline typename SComplex<NeighboursModelT>::ColoredIterators::Iterators SComplex<NeighboursModelT>::iterators(const Color& color) {

}

template<template<typename Object, typename Color> class NeighboursModelT>
template<typename SComplex<NeighboursModelT>::Color color>
inline typename SComplex<NeighboursModelT>::ColoredIterators::template Color<color>::Iterators SComplex<NeighboursModelT>::iterators() {

}

template<template<typename Object, typename Color> class NeighboursModelT>
template<typename SComplex<NeighboursModelT>::Color color>
inline typename SComplex<NeighboursModelT>::ColoredConstIterators::template Color<color>::Iterators SComplex<NeighboursModelT>::iterators() const {

}


#endif // _SCOMPLEX_HPP

