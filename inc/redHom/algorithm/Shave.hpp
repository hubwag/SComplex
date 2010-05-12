#ifndef _SCOMPLEX_ALGS_SHAVE_HPP
#define _SCOMPLEX_ALGS_SHAVE_HPP

#include "strategy/DefaultReduceStrategy.hpp"
#include "strategy/DefaultReduceStrategy_CubSComplex.hpp"

template<typename StrategyT>
class ShaveAlgorithm {

public:
  typedef StrategyT Strategy;
  typedef typename Strategy::SComplex SComplex;
  typedef typename SComplex::Cell Cell;

  ShaveAlgorithm(Strategy* _strategy): strategy(_strategy) {}

  ~ShaveAlgorithm() {
	 delete strategy;
  }

  void operator()();

private:
  Strategy* strategy;
};

class ShaveAlgorithmFactory {

public:
  template<typename SComplex>
  static ShaveAlgorithm<DefaultReduceStrategy<SComplex> > createDefault(SComplex& s) {
	 return ShaveAlgorithm<DefaultReduceStrategy<SComplex> >(new DefaultReduceStrategy<SComplex>(s));
  }
};

template<typename StrategyT>
inline void ShaveAlgorithm<StrategyT>::operator()(){

  for(typename SComplex::Dim d = strategy->getMaxDim() - 1; d >= 0; --d){
	 typedef typename SComplex::ColoredIterators::Iterators::DimCells::iterator DimIt;

	 typename SComplex::ColoredIterators::Iterators::DimCells dimCells = strategy->getComplex().template iterators<1>().dimCells(d);
	 for (DimIt it = dimCells.begin(),
	 		  end = dimCells.end();
	 		it != end; ++it) {
	 // BOOST_FOREACH(typename DimIt::value_type v,
	 // 					dimCells) {
	 //strategy->reduceIfPossible(v);
		typename DimIt::value_type v = *it;
		typename StrategyT::Traits::template GetReductionPair<typename DimIt::value_type>::result_type reductionPair =
		  strategy->getReductionPair(v);
		if (reductionPair) {
		  strategy->reduce(*reductionPair);
		  strategy->reduce(v);
		}
	 }
  }
}

#endif
