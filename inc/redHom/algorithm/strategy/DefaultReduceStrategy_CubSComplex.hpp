#ifndef SCOMPLEX_ALGS_CUBSCOMPLEX_REDUCE_STRATEGY_HPP_
#define SCOMPLEX_ALGS_CUBSCOMPLEX_REDUCE_STRATEGY_HPP_


#include "../../RedHomCAPD.h"

#include "DefaultReduceStrategy.hpp"
#include "DefaultReduceStrategyTraits_CubSComplex.hpp"
#include "redHom/complex/cubical/CubSComplex.hpp"

#if 0

template<int DIM>
class DefaultReduceStrategy<CubSComplex<DIM> >: public DefaultReduceStrategyBase<CubSComplex<DIM> > {

  typename CubSComplex<DIM>::DynamicCell dynamicCell;
  //using typename DefaultReduceStrategyBase<CubSComplex<DIM> >::Traits;

  using DefaultReduceStrategyBase<CubSComplex<DIM> >::complex;

public:
  typedef typename DefaultReduceStrategyBase<CubSComplex<DIM> >::Traits Traits;

  DefaultReduceStrategy(CubSComplex<DIM>& complex): DefaultReduceStrategyBase<CubSComplex<DIM> >(complex), dynamicCell(complex) {}


  /*
  template<typename ArgT>
  typename Traits::template GetCoreductionPair<ArgT>::result_type
  getCoreductionPair(const ArgT& cell)
  {
  	 if (complex.getUniqueFace(cell, dynamicCell)) {
		return typename Traits::template GetCoreductionPair<ArgT>::result_type(dynamicCell.getImpl());
  	 } else {
		return typename Traits::template GetCoreductionPair<ArgT>::result_type();
  	 }
  }

  template<typename ArgT>
  typename Traits::template GetReductionPair<ArgT>::result_type
  getReductionPair(const ArgT& cell)
  {
	 if (complex.getUniqueCoFace(cell, dynamicCell)) {
		return typename Traits::template GetReductionPair<ArgT>::result_type(dynamicCell.getImpl());
  	 } else {
		return typename Traits::template GetReductionPair<ArgT>::result_type();
  	 }
  } */

  typedef CubSComplex<DIM> SComplex;

  template<typename ArgT>
    typename Traits::template GetCoreductionPair<ArgT>::result_type
    getCoreductionPair(const ArgT& cell)
    {
        int times = 0;
        BOOST_FOREACH(typename SComplex::ColoredIterators::Iterators::BdCells::iterator::value_type v,
                      this->complex.iterators(1).bdCells(cell))
        {
            if (times == 0)
            {
                this->dummyCell3 = v;
            }
            ++times;
            if (times == 2)
            {
                break;
            }
        }

        if (times == 1)
        {
            return typename Traits::template GetCoreductionPair<ArgT>::result_type(this->dummyCell3);
        }
        return typename Traits::template GetCoreductionPair<ArgT>::result_type();
    }

    template<typename ArgT>
    typename Traits::template GetReductionPair<ArgT>::result_type
    getReductionPair(const ArgT &cell)
    {
        int times = 0;
        BOOST_FOREACH(typename SComplex::ColoredIterators::Iterators::CbdCells::iterator::value_type v,
                      this->complex.iterators(1).cbdCells(cell))
        {
            if (times == 0)
            {
                this->dummyCell2 = v;
            }
            ++times;
            if (times == 2)
            {
                break;
            }
        }

        if (times == 1)
        {
            return typename Traits::template GetReductionPair<ArgT>::result_type(this->dummyCell2);
        }
        return typename Traits::template GetReductionPair<ArgT>::result_type();
    }


size_t getMaxDim__() {
	 typename SComplex::Dim maxDim = 0;

	 for (typename SComplex::ColoredIterators::Iterators::AllCells::iterator it = this->complex.iterators(1).allCells().begin(),
			  end = this->complex.iterators(1).allCells().end();
			it != end; ++it) {
		maxDim = std::max(maxDim, it->getDim());
	 }

	 return maxDim;
  }


  size_t getMaxDim() {
  	cout << "GETTING SIZE OF CUB COMPLEX!" << endl;
	return complex.getDim();
  }

};

#endif

#endif
