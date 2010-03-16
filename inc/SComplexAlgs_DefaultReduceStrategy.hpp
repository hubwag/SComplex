#ifndef SCOMPLEX_ALGS_DEFAULT_REDUCE_STRATEGY_HPP_
#define SCOMPLEX_ALGS_DEFAULT_REDUCE_STRATEGY_HPP_

#include "CellProxy.hpp"
#include <boost/optional.hpp>
#include <boost/foreach.hpp>
#include <utility>

template<typename SComplexT>
class ReduceStrategyTraits {
public:

  template<typename ImplT>
  class ReductionPair {
  public:
	 typedef CellProxy<typename SComplexT::Cell> first_type;
	 typedef CellProxy<typename SComplexT::Cell> second_type;
	 typedef std::pair<first_type, second_type > type;
  };
};

template<typename SComplexT>
class DefaultReduceStrategyBase {

public:
  typedef SComplexT SComplex;
  typedef ReduceStrategyTraits<SComplex> Traits;
  typedef typename SComplex::Cell Cell;


  DefaultReduceStrategyBase(SComplex& _complex): complex(_complex), dummyCell1(_complex), dummyCell2(_complex),  dummyCell3(_complex) {}
  
  SComplex& getComplex() const {
	 return complex;
  }

  // Cell& getFace(CoreductionPair& coRedPair) {
  // 	 return coRedPair.first;
  // }
  
  // bool reduced(const Cell& cell) const {
  // 	 return cell.getColor() == 2;
  // }

  // void coreduce(CoreductionPair& coRedPair) const {
  // 	 coRedPair.first.template setColor<2>();
  // 	 coRedPair.second.template setColor<2>();
  // }

  // template<typename ImplT>
  // void reduce(typename Traits::template ReductionPair<ImplT>::type& redPair) const {
  // 	 redPair.first.template setColor<2>();
  // 	 redPair.second.template setColor<2>();
  // }

  template<typename ImplT>
  void reduce(CellProxy<ImplT>& cell) {
	 cell.template setColor<2>();
  }
  
  // boost::optional<Cell&> extract() {
  // 	 typename SComplex::ColoredIterators::Iterators::DimCells::iterator end = complex.iterators(1).dimCells(0).end(),
  // 		it = complex.iterators(1).dimCells(0).begin();

  // 	 if (it != end) { 
  // 		dummyCell1 = *it;
  // 		return dummyCell1;
  // 	 }
  // 	 return boost::optional<Cell&>();	 
  // }
  
  // boost::optional<CoreductionPair> forceCoreductionPair() {
  // 	 return boost::optional<CoreductionPair>();
  // }

  // boost::optional<CoreductionPair> getCoreductionPair(Cell& cell) {
  // 	 int times = 0;
  // 	 BOOST_FOREACH(typename SComplex::ColoredIterators::Iterators::BdCells::iterator::value_type v,
  // 						complex.iterators(1).bdCells(cell)) {
  // 		if (times == 0) {
  // 		  dummyCell3 = v;
  // 		}
  // 		++times;
  // 		if (times == 2) {
  // 		  return boost::optional<ReductionPair>();
  // 		}
  // 	 }

  // 	 if (times == 1) {
  // 		return std::make_pair(dummyCell3, cell);
  // 	 }
  // 	 return boost::optional<ReductionPair>();
  // }

  // boost::optional<ReductionPair> getReductionPair(Cell& cell) {
  // 	 int times = 0;
  // 	 BOOST_FOREACH(typename SComplex::ColoredIterators::Iterators::CbdCells::iterator::value_type v,
  // 						complex.iterators(1).cbdCells(cell)) {
  // 		if (times == 0) {
  // 		  dummyCell2 = v;
  // 		}
  // 		++times;
  // 		if (times == 2) {
  // 		  return boost::optional<ReductionPair>();
  // 		}
  // 	 }

  // 	 if (times == 1) {
  // 		return std::make_pair(cell, dummyCell2);
  // 	 }
  // 	 return boost::optional<ReductionPair>();
  // }


  size_t getMaxDim() {
	 size_t maxDim = 0;
	 for (typename SComplex::ColoredIterators::Iterators::AllCells::iterator it = complex.template iterators<1>().allCells().begin(),
			  end = complex.template iterators<1>().allCells().end();
			it != end; ++it) {

		maxDim = std::max(maxDim, it->getDim());
	 }
	 
	 return maxDim;
  }

protected:
  SComplex& complex;
  Cell dummyCell1, dummyCell2, dummyCell3;
};

template<typename SComplexT>
class DefaultReduceStrategy: public DefaultReduceStrategyBase<SComplexT> {

public:
  DefaultReduceStrategy(SComplexT& _complex): DefaultReduceStrategyBase<SComplexT>(_complex) {}
};

#endif

