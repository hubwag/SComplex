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
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

#include "util/Iterators.hpp"

#include "util/ColorListModel.hpp"

template<template<typename Object, typename Color> class NeighboursModelT>
class SComplex {
public:
  typedef size_t Color;
  typedef size_t Dim;
  typedef size_t Id;

  class Cell_impl;  
  typedef Util::Neighbours::ColorListModel<boost::shared_ptr<Cell_impl>, Color> Cells;
  
  class Cell_impl {

  public:
	 typedef typename SComplex::Color Color;

	 typedef typename Cells::ObjectPtrsIterator IteratorInCells;

	 
	 explicit Cell_impl(const Id& _id, const Color& _color): id(_id), color(_color) {}

	 const Color& getColor() const { return color; }

	 void setColor(const Color& _color) { this->color = _color; }
	 
	 const Id& getId() const { return id; }
	 
	 const IteratorInCells& getIteratorInCells() const { return iteratorInCells; }

	 void setIteratorInCells(const IteratorInCells& _iteratorInCells) { this->iteratorInCells = _iteratorInCells; }
	 
  private:
	 Id id;
	 Color color;
	 IteratorInCells iteratorInCells;
  };

  typedef boost::shared_ptr<Cell_impl> CellImplPtr;

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

		// BOOST_FOREACH(typename SComplex<NeighboursModelT>::Iterators::BdCells::iterator::value_type bd,
		// 				  complex->iterators().bdCells(*this)) {
	 
		// }
		complex->cells.changeColor(cell->getIteratorInCells(), cell->getColor(), newColor);
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

	 explicit NeighbourLink(const CellImplPtr& _cell) : cell(_cell.get()) {}
  };

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
  struct CellFromCellsExtractor: public std::unary_function<const CellImplPtr&, CellType>
  {
	 SComplex* complex;
	 explicit CellFromCellsExtractor(SComplex* _complex): complex(_complex) {}
		
	 CellType operator()(const CellImplPtr& cell) const {
		return CellType(complex, cell.get());
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
	 
		typedef Util::Iterators::RangeTransform<const typename Cells::ObjectsInColor, CellFromCellsExtractor<CellType> > AllCells;	 
		typedef Util::Iterators::RangeTransform<const typename Cells::ObjectsInColor, CellFromCellsExtractor<CellType> > DimCells;
		typedef Util::Iterators::RangeTransform<const typename NeighboursModel::ObjectsInColor, CellFromNeighbourLinkExtractor<CellType> > BdCells, CbdCells;
	 
		explicit IteratorsImpl(SComplex* _complex, const Color& _color): complex(_complex), color(_color)  {}
	 
		AllCells allCells() {
		  return AllCells(complex->cells.objectsInColor(color), CellFromCellsExtractor<CellType>(complex)); 
		}
	 
		DimCells dimCells(const Dim& dim);
	 
		BdCells bdCells(const ConstCell& cell) {
		  return BdCells(complex->boundaries[cell.getId()].objectsInColor(color), CellFromNeighbourLinkExtractor<CellType>(complex));
		}
	 
		CbdCells cbdCells(const ConstCell& cell) {
		  return CbdCells(complex->coboundaries[cell.getId()].objectsInCOlor(color), CellFromNeighbourLinkExtractor<CellType>(complex));
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
		CellImplPtr cell(new Cell_impl(id, color));
		cell->setIteratorInCells(cells.add(cell, color));
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
	 return Cell(this, cells.allObjects()[id].get());
  }

  ConstCell operator[](const Id id) const {
	 return ConstCell(this, cells.allObjects()[id].get());
  }
  
  Iterators iterators() { return Iterators(nonConstThis.get_pointer()); }
  ConstIterators iterators() const { return ConstIterators(nonConstThis.get_pointer()); }

  typename ColoredConstIterators::Iterators iterators(const Color& color) const { return typename ColoredConstIterators::Iterators(nonConstThis.get_pointer(), color); }
  typename ColoredIterators::Iterators iterators(const Color& color) { return typename ColoredIterators::Iterators(nonConstThis.get_pointer(), color); }

  template<Color color>
  typename ColoredIterators::template Color<color>::Iterators iterators() { return iterators(color); }

  template<Color color>
  typename ColoredConstIterators::template Color<color>::Iterators iterators() const { return iterators(color); }

private:
  mutable Cells cells;
  //  std::vector<std::pair<typename Cells::iterator, typename Cells::iterator> > cellsByDim; // pairs of begin/end in cells
  std::vector<NeighboursModel> boundaries, coboundaries; // (co)boundaries by cell id
  boost::reference_wrapper<SComplex> nonConstThis;  
};


#endif // _SCOMPLEX_HPP

