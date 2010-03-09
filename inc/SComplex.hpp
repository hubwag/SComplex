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
#include <boost/utility.hpp>

#include "util/Iterators.hpp"

#include "util/ColorListModel.hpp"

template<template<typename Object, typename Color> class NeighboursModelT>
class SComplex: boost::noncopyable {
public:
  typedef size_t Color;
  typedef size_t Dim;
  typedef size_t Id;

private:
  class CellImpl;  
  typedef Util::Neighbours::ColorListModel<boost::shared_ptr<CellImpl>, Color> Cells;
  typedef std::vector<Util::Neighbours::ColorListModel<boost::reference_wrapper<CellImpl>, Color> > CellsByDim;
  
  class CellImpl: boost::noncopyable {
  public:
	 typedef typename SComplex::Color Color;

	 typename Cells::ObjectPtrsIterator iteratorInCells;
	 typename CellsByDim::value_type::ObjectPtrsIterator iteratorInCellsByDim;
	 
	 explicit CellImpl(const Id& _id, const Dim& _dim, const Color& _color): id(_id), dim(_dim), color(_color) {}

	 const Color& getColor() const { return color; }	 
	 const Id& getId() const { return id; }
	 const Dim& getDim() const { return dim; }
	 void setColor(const Color& _color) { this->color = _color; }
	 
  private:
	 const Id id;
	 const Dim dim;
	 Color color;
  };

  typedef boost::shared_ptr<CellImpl> CellImplPtr;

public:
  
  class ConstCell {	 
  public:
	 explicit ConstCell(const SComplex* _complex, CellImpl* _cell): complex(_complex), cell(_cell) {}
	 
	 const Color& getColor() const { return cell->getColor(); }
	 const Id& getId() const { return cell->getId(); }
	 const Dim& getDim() const { return cell->getDim(); }
  protected:
	 const SComplex* complex;
	 CellImpl* cell;
  };

  class Cell: public ConstCell {
  public:
	 explicit Cell(SComplex* _complex, CellImpl* _cell): ConstCell(_complex, _cell), complex(_complex) {}

	 void setColor(const Color& newColor) {
		for (typename NeighboursModel::AllObjects::iterator it = complex->boundaries[this->getId()].allObjects().begin(),
				 end = complex->boundaries[this->getId()].allObjects().end(); it != end; ++it) {
		  complex->coboundaries[it->cell->getId()].changeColor(it->iteratorInNeighbour, this->getColor(), newColor);
		}

		for (typename NeighboursModel::AllObjects::iterator it = complex->coboundaries[this->getId()].allObjects().begin(),
				 end = complex->coboundaries[this->getId()].allObjects().end(); it != end; ++it) {
		  complex->boundaries[it->cell->getId()].changeColor(it->iteratorInNeighbour, this->getColor(), newColor);		  
		}

		complex->cells.changeColor(cell->iteratorInCells, cell->getColor(), newColor);
		complex->cellsByDim[cell->getDim()].changeColor(cell->iteratorInCellsByDim, cell->getColor(), newColor);
		cell->setColor(newColor);
	 }
	 
  private:
	 using ConstCell::cell;
	 SComplex* complex;
  };
	 
  struct NeighbourLink;
  typedef Util::Neighbours::ColorListModel<NeighbourLink, Color> NeighboursModel;

  struct NeighbourLink {
	 typedef typename NeighboursModel::ObjectPtrsIterator IteratorInNeighbour;
	 CellImpl* cell;	 
	 IteratorInNeighbour iteratorInNeighbour; //an iterator to NeighbourLinkPtrs not in this, but in an instance of the class for a neighbour.

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

  template<typename CellType>
  struct CellFromCellsByDimExtractor: public std::unary_function<const boost::reference_wrapper<CellImpl>&, CellType>
  {
	 SComplex* complex;
	 explicit CellFromCellsByDimExtractor(SComplex* _complex): complex(_complex) {}
		
	 CellType operator()(const boost::reference_wrapper<CellImpl>& cell) const {
		return CellType(complex, cell.get_pointer());
	 }
  };
  
  template<bool isConst>
  class IteratorsImpl {
	 typedef typename boost::mpl::if_c<isConst, ConstCell, Cell>::type CellType;	 
  public:
	 
	 typedef Util::Iterators::RangeTransform<const typename Cells::AllObjects, CellFromCellsExtractor<CellType> > AllCells;	 
	 typedef Util::Iterators::RangeTransform<const typename CellsByDim::value_type::AllObjects, CellFromCellsByDimExtractor<CellType> > DimCells;
	 typedef Util::Iterators::RangeTransform<const typename NeighboursModel::AllObjects, CellFromNeighbourLinkExtractor<CellType> > BdCells, CbdCells;
	 
	 explicit IteratorsImpl(SComplex* _complex): complex(_complex) {}
	 
	 AllCells allCells() {
		return AllCells(complex->cells.allObjects(), CellFromCellsExtractor<CellType>(complex)); 
	 }
	 
	 DimCells dimCells(const Dim& dim) {
		return DimCells(complex->cellsByDim[dim].allObjects(), CellFromCellsByDimExtractor<CellType>(complex)); 
	 }
	 
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
		typedef Util::Iterators::RangeTransform<const typename CellsByDim::value_type::ObjectsInColor, CellFromCellsByDimExtractor<CellType> > DimCells;
		typedef Util::Iterators::RangeTransform<const typename NeighboursModel::ObjectsInColor, CellFromNeighbourLinkExtractor<CellType> > BdCells, CbdCells;
	 
		explicit IteratorsImpl(SComplex* _complex, const Color& _color): complex(_complex), color(_color)  {}
	 
		AllCells allCells() {
		  return AllCells(complex->cells.objectsInColor(color), CellFromCellsExtractor<CellType>(complex)); 
		}
	 
		DimCells dimCells(const Dim& dim) {
		  return DimCells(complex->cellsByDim[dim].objectsInColor(color), CellFromCellsByDimExtractor<CellType>(complex)); 
		}
	 
		BdCells bdCells(const ConstCell& cell) {
		  return BdCells(complex->boundaries[cell.getId()].objectsInColor(color), CellFromNeighbourLinkExtractor<CellType>(complex));
		}
	 
		CbdCells cbdCells(const ConstCell& cell) {
		  return CbdCells(complex->coboundaries[cell.getId()].objectsInColor(color), CellFromNeighbourLinkExtractor<CellType>(complex));
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
  typedef std::vector<Dim> Dims;
  
  SComplex(size_t colors, size_t _size, const Dims& dims, const KappaMap& kappaMap = KappaMap()): boundaries(_size), coboundaries(_size), nonConstThis(*this) {
	 using boost::get;
	 BOOST_ASSERT(_size == dims.size());
	 

	 // init dimension related structures
	 Dim maxDim = _size > 0 ? *std::max_element(dims.begin(), dims.end()) : 0;	 
	 std::vector<size_t> cellsInDim(maxDim + 1);

	 cellsByDim.resize(maxDim + 1);
	 
	 for (Id id = 0; id < _size; ++id) {
		++cellsInDim[dims[id]];
	 }
	 
	 for (Dim dim = 0; dim <= maxDim; ++dim) {
		cellsByDim[dim].init(colors, cellsInDim[dim]);
	 }

	 //init cells
	 cells.init(_size, colors);
	 
	 for (Id id = 0; id < _size; ++id) {
		Color color = 0;
		CellImplPtr cell(new CellImpl(id, dims[id], color));
		cell->iteratorInCells = cells.add(cell, color);
		cell->iteratorInCellsByDim = cellsByDim[cell->getDim()].add(boost::ref(*(cell.get())), color);
	 }

	 // init neighbours
	 std::vector<size_t> boundariesSize(_size);
	 std::vector<size_t> coboundariesSize(_size);
	 BOOST_FOREACH(KappaMap::value_type kappa, kappaMap) {
		++boundariesSize[get<0>(kappa)];
		++coboundariesSize[get<1>(kappa)];
	 }

	 for (size_t i = 0; i < _size; ++i) {
		boundaries[i].init(colors, boundariesSize[i]);
		coboundaries[i].init(colors, coboundariesSize[i]);
	 }

	 BOOST_FOREACH(KappaMap::value_type kappa, kappaMap) {
		Id coface = get<0>(kappa);
		Id face = get<1>(kappa);

		BOOST_ASSERT(coface != face);
		typename NeighbourLink::IteratorInNeighbour bdIteratorInNeighbour = boundaries[coface].add(NeighbourLink(cells.allObjects()[face]), cells.allObjects()[face]->getColor());
		typename NeighbourLink::IteratorInNeighbour cbdIteratorInNeighbour = coboundaries[face].add(NeighbourLink(cells.allObjects()[coface]), cells.allObjects()[coface]->getColor());
		boundaries[coface].allObjects().back().iteratorInNeighbour = cbdIteratorInNeighbour;
		coboundaries[face].allObjects().back().iteratorInNeighbour = bdIteratorInNeighbour;
	 }
	 
  }

  size_t cardinality() const { return nonConstThis.get().cells.allObjects().size(); }

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
  Cells cells;
  CellsByDim cellsByDim;
  std::vector<NeighboursModel> boundaries, coboundaries; // (co)boundaries by cell id
  boost::reference_wrapper<SComplex> nonConstThis;  
};


#endif // _SCOMPLEX_HPP

