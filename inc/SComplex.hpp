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

#include "util/ColorListModel.hpp"

template<template<typename Object, typename Color> class NeighboursModelT>
class SComplex {
public:
  typedef size_t Color;
  typedef size_t Dim;
  typedef size_t Id;
  class Cell_impl;
  
  
  class Cell_impl {

  public:
	 typedef typename SComplex::Color Color;
	 
	 explicit Cell_impl(const Id& _id, const Color& _color): id(_id), color(_color) {}

	 const Color& getColor() const { return color; }

	 void setColor(const Color& _color) { this->color = _color; }
	 
	 const Id& getId() const { return id; }
	 
  private:
	 Id id;
	 Color color;
  };


  class ConstCell {	 
  public:
	 explicit ConstCell(const SComplex* _complex, Cell_impl* _cell): complex(_complex), cell(_cell) {}
	 
	 const Color& getColor() const { return cell->getColor(); }

	 const Id& getId() const { return cell->getId(); }

  protected:
	 const SComplex* complex;
	 Cell_impl* cell;
  };

  class Cell: public ConstCell {
  public:
	 explicit Cell(SComplex* _complex, Cell_impl* _cell): ConstCell(_complex, _cell), complex(_complex) {}

	 void setColor(const Color& newColor) {
		Color oldColor = cell->getColor();

		// BOOST_FOREACH(typename SComplex<NeighboursModelT>::Iterators::BdCells::iterator::value_type bd,
		// 				  complex->iterators().bdCells(*this)) {
	 
		// }
			  
		cell->setColor(newColor);
	 }
	 
  private:
	 using ConstCell::cell;
	 SComplex* complex;
  };
	 
  struct NeighbourLink;
  typedef Util::Neighbours::ColorListModel<NeighbourLink, Color> NeighboursModel;

  struct NeighbourLink {
	 Cell_impl* cell;
	 typename NeighboursModel::ObjectPtrsIterator neighbourLinkPtrsIterator; //an iterator to NeighbourLinkPtrs not in this, but in an instance of the class for a neighbour.

	 explicit NeighbourLink(Cell_impl* _cell) : cell(_cell) {}
  };

  //typedef NeighboursModelT<boost::reference_wrapper<Cell_impl>, Color> NeighboursModel;
  //typedef std::vector<Cell_impl*> Cells;
  typedef Util::Neighbours::ColorListModel<Cell_impl*, Color> Cells;

private:

  template<typename CellType>
  struct CellFromNeighbourLinkExtractor: public std::unary_function<const NeighbourLink&, CellType>
  {
	 SComplex* complex;
	 explicit CellFromNeighbourLinkExtractor(SComplex* _complex): complex(_complex) {}
		
	 CellType operator()(const NeighbourLink&  link) const {
		return CellType(complex, link.cell);
	 }
  };

  template<typename CellType>
  struct CellFromCellsExtractor: public std::unary_function<Cell_impl*, CellType>
  {
	 SComplex* complex;
	 explicit CellFromCellsExtractor(SComplex* _complex): complex(_complex) {}
		
	 CellType operator()(Cell_impl* cell) const {
		return CellType(complex, cell);
	 }
  };

  
  template<bool isConst>
  class IteratorsImpl {
	 typedef typename boost::mpl::if_c<isConst, ConstCell, Cell>::type CellType;	 
  public:
	 
	 typedef Util::Iterators::RangeTransform<const typename Cells::AllObjects, CellFromCellsExtractor<CellType> > AllCells;	 
	 typedef Util::Iterators::RangeTransform<const typename Cells::AllObjects, CellFromCellsExtractor<CellType> > DimCells;
	 typedef Util::Iterators::RangeTransform<const typename NeighboursModel::AllObjects, CellFromNeighbourLinkExtractor<CellType> > BdCells, CbdCells;
	 
	 explicit IteratorsImpl(SComplex* _complex): complex(_complex) {}
	 
	 AllCells allCells() {
		return AllCells(complex->cells.allObjects(), CellFromCellsExtractor<CellType>(complex)); 
	 }
	 
	 DimCells dimCells(const Dim& dim);
	 
	 BdCells bdCells(const ConstCell& cell) {
		return BdCells(complex->boundaries[cell.getId()].allObjects(), CellFromNeighbourLinkExtractor<CellType>(complex));
	 }
	 
	 CbdCells cbdCells(const ConstCell& cell) {
		return CbdCells(complex->coboundaries[cell.getId()].allObjects(), CellFromNeighbourLinkExtractor<CellType>(complex));
	 }

  private:
	 SComplex* complex;	 
  };

  template<bool isConst>
  class ColoredIteratorsImpl {

	 class IteratorsImpl {
		typedef typename boost::mpl::if_c<isConst, ConstCell, Cell>::type CellType;	 
	 public:
	 
		typedef Util::Iterators::RangeTransform<Cells, CellFromCellsExtractor<CellType> > AllCells;	 
		typedef Util::Iterators::RangeTransform<Cells, CellFromCellsExtractor<CellType> > DimCells;
		typedef Util::Iterators::RangeTransform<typename NeighboursModel::ObjectsInColor, CellFromNeighbourLinkExtractor<CellType> > BdCells, CbdCells;
	 
		explicit IteratorsImpl(SComplex* _complex, const Color& _color): complex(_complex), color(_color)  {}
	 
		AllCells allCells() {
		  return AllCells(complex->cells, CellFromCellsExtractor<CellType>(complex)); 
		}
	 
		DimCells dimCells(const Dim& dim);
	 
		BdCells bdCells(const ConstCell& cell) {
		  return BdCells(complex->boundaries[cell.getId()].allNeighbours(), CellFromNeighbourLinkExtractor<CellType>(complex));
		}
	 
		CbdCells cbdCells(const ConstCell& cell) {
		  return CbdCells(complex->coboundaries[cell.getId()].allNeighbours(), CellFromNeighbourLinkExtractor<CellType>(complex));
		}

	 private:
		SComplex* complex;
		Color color;
	 };
	 
  public:
	 typedef IteratorsImpl Iterators;

	 template<Color color>
	 class Color {
	 public:
		typedef IteratorsImpl Iterators;
	 };
  };

public:
  
  typedef IteratorsImpl<false> Iterators;
  typedef IteratorsImpl<true> ConstIterators;

  
  typedef ColoredIteratorsImpl<false> ColoredIterators;
  typedef ColoredIteratorsImpl<true> ColoredConstIterators;

  typedef std::vector<boost::tuple<Id, Id, int> > KappaMap;
  
  SComplex(size_t colors, size_t _size, const KappaMap& kappaMap = KappaMap()): boundaries(_size), coboundaries(_size), nonConstThis(*this) {
	 using namespace boost::lambda;
	 using boost::ref;
	 using boost::lambda::_1;
	 using boost::lambda::bind;
	 using boost::get;

	 cells.init(_size, colors);
	 // Id id = 0;
	 // std::for_each(cells.begin(), cells.end(), _1 = bind(new_ptr<Cell>(), ref(id)++, 0));
	 for (Id id = 0; id < _size; ++id) {
		//cells.push_back(new Cell_impl(id, 0));
		Color color = 0;
		cells.add(new Cell_impl(id, color), color);
	 }
		
	 //TODO call init with correct size
	 std::for_each(boundaries.begin(), boundaries.end(), bind(&NeighboursModel::init, _1, colors, 0));
	 std::for_each(coboundaries.begin(), coboundaries.end(), bind(&NeighboursModel::init, _1, colors, 0));

	 BOOST_FOREACH(KappaMap::value_type kappa, kappaMap) {
		Id coface = get<0>(kappa);
		Id face = get<1>(kappa);

		boundaries[coface].add(NeighbourLink(cells.allObjects()[face]), cells.allObjects()[face]->getColor());
		coboundaries[face].add(NeighbourLink(cells.allObjects()[coface]), cells.allObjects()[coface]->getColor());
	 }
	 
  }

  size_t cardinality() const { return cells.allObjects().size(); }

  Cell operator[](const Id id) {
	 return Cell(this, cells.allObjects()[id]);
  }

  ConstCell operator[](const Id id) const {
	 return ConstCell(this, cells.allObjects()[id]);
  }
  
  Iterators iterators() { return Iterators(nonConstThis.get_pointer()); }
  ConstIterators iterators() const { return ConstIterators(nonConstThis.get_pointer()); }

  typename ColoredConstIterators::Iterators iterators(const Color& color) const;
  typename ColoredIterators::Iterators iterators(const Color& color);

  template<Color color>
  typename ColoredIterators::template Color<color>::Iterators iterators();

  template<Color color>
  typename ColoredConstIterators::template Color<color>::Iterators iterators() const;

private:
  mutable Cells cells;
  //  std::vector<std::pair<typename Cells::iterator, typename Cells::iterator> > cellsByDim; // pairs of begin/end in cells
  std::vector<NeighboursModel> boundaries, coboundaries; // (co)boundaries by cell id
  boost::reference_wrapper<SComplex> nonConstThis;  
};

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

