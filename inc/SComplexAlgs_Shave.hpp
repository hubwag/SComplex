#ifndef _SCOMPLEX_ALGS_SHAVE_HPP
#define _SCOMPLEX_ALGS_SHAVE_HPP

#include "SComplexAlgs_DefaultReduceStrategy.hpp"

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

  for(int d=2/* strategy->getMaxDim()*/;d>=0;--d){
	 typedef typename SComplex::ColoredIterators::Iterators::DimCells::iterator DimIt;

	 for (DimIt it = strategy->getComplex().template iterators<1>().dimCells(d).begin(),
			  end = strategy->getComplex().template iterators<1>().dimCells(d).end();
			it != end; ++it) {
		typename DimIt::value_type v = *it;
		//strategy->reduceIfPossible(v);
		typename StrategyT::Traits::template GetReductionPair<typename DimIt::value_type>::result_type reductionPair =
		  strategy->template getReductionPair<CubSComplex::BitCoordPtrCellImpl>(v);
		if (reductionPair) {
		  strategy->reduce(*reductionPair);
		  strategy->reduce(Strategy::Traits::makeProxy(v));
		}
	 }
  }
}

#endif
