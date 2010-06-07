#ifndef SCOMPLEX_ALGS_HPP_
#define SCOMPLEX_ALGS_HPP_


#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/optional.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>


#include <utility>
#include <algorithm>
#include <map>
#include <set>

#include "Coreduction.hpp"
#include "Shave.hpp"
#include "strategy/DefaultReduceStrategy.hpp"
#include "../RedHomCAPD.h"


template<typename SComplexT, typename ReducibleFreeChainComplexT>
class ReducibleFreeChainComplexOverZFromSComplexAlgorithm {

public:
  typedef SComplexT SComplex;
  typedef typename DefaultReduceStrategy<SComplexT>::Cell Cell;

  typedef ReducibleFreeChainComplexT ReducibleFreeChainComplex;

  class SComplexChainCell {
	 SComplex& complex;
	 const Cell cell;
	 const int embededDim;

  public:
	 SComplexChainCell(SComplex& _complex, const Cell& _cell, int _embededDim): complex(_complex), cell(_cell), embededDim(_embededDim) {

	 }

	 int embDim() const {
		return this->embededDim;
	 }

	 int ownDim() const {
	   return cell.getDim();
	 }

	 void boundary(std::map<SComplexChainCell,int>& A_boundary) const {
	   typename SComplex::ColoredIterators::Iterators::BdCells bdCells = complex.iterators(1).bdCells(this->cell);

	   for (typename SComplex::ColoredIterators::Iterators::BdCells::iterator it = bdCells.begin(),
		  end = bdCells.end(); it != end; ++it) {
	     A_boundary.insert(std::make_pair(SComplexChainCell(complex, *it, embededDim),
					      complex.coincidenceIndex(this->cell, *it) ));
	   }
	 }

	 bool operator<(const SComplexChainCell& b) const {
		return this->cell < b.cell;
	 }
  };

public:


  ReducibleFreeChainComplexOverZFromSComplexAlgorithm(SComplex& _s): s(_s) {}
  CRef<ReducibleFreeChainComplex> operator()();

private:
  SComplex& s;
};

template<typename SComplexT, typename ReducibleFreeChainComplexT>
inline CRef<ReducibleFreeChainComplexT> ReducibleFreeChainComplexOverZFromSComplexAlgorithm<SComplexT, ReducibleFreeChainComplexT>::operator()(){

  std::set<SComplexChainCell> cells;

  size_t maxDim = (DefaultReduceStrategy<SComplexT>(s)).getMaxDim(); // TODO add strategy as a member

  for (typename SComplex::ColoredIterators::Iterators::AllCells::iterator it = s.iterators(1).allCells().begin(),
			end = s.iterators(1).allCells().end();
		 it != end; ++it) {
	 cells.insert(SComplexChainCell(s, *it, maxDim));
  }

  CRef<ReducibleFreeChainComplex> rfccCR( new ReducibleFreeChainComplex(cells));

  return rfccCR;
}


#endif // SCOMPLEX_ALGS_HPP_
