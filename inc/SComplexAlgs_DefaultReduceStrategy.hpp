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
  class Proxy {
  public:
	 typedef CellProxy<ImplT> type;
  };
  
  template<typename ImplT>
  class ReductionPair {
  public:
	 typedef Proxy<typename SComplexT::Cell> first_type;
	 typedef Proxy<typename SComplexT::Cell> second_type;
	 typedef std::pair<first_type, second_type > type;
  };
};

template<typename SComplexT>
class DefaultReduceStrategyBase {

public:
  typedef SComplexT SComplex;
  typedef ReduceStrategyTraits<SComplex> Traits;
  typedef typename SComplex::Cell Cell;


  DefaultReduceStrategyBase(SComplex& _complex): complex(_complex), dummyCell2(_complex),  dummyCell3(_complex) {}
  
  SComplex& getComplex() const {
	 return complex;
  }
  
  template<typename ImplT>
  static bool reduced(const typename Traits::template Proxy<ImplT>& cell) {
  	 return cell.getColor() == 2;
  }

  template<typename ImplT1, typename ImplT2>	 
  static void coreduce(const typename Traits::template Proxy<ImplT1>& a, const typename Traits::template Proxy<ImplT2>& b)  {
  	 a.template setColor<2>();
  	 b.template setColor<2>();
  }

  template<typename ImplT1, typename ImplT2>	 
  static void reduce(const typename Traits::template Proxy<ImplT1>& a, const typename Traits::template Proxy<ImplT2>& b)  {
  	 a.template setColor<2>();
  	 b.template setColor<2>();
  }

  template<typename ImplT>
  static void reduce(const typename Traits::template Proxy<ImplT>& cell) {
	 cell.template setColor<2>();
  }
  
  typename Traits::Extract::result_type extract() {
	 typename SComplex::ColoredIterators::Iterators::DimCells dimCells = complex.iterators(1).dimCells(0);
  	 typename SComplex::ColoredIterators::Iterators::DimCells::iterator end = dimCells.end(),
  		it = dimCells.begin();

  	 if (it != end) { 
  		return typename Traits::Extract::result_type::value_type(*it);
  	 }
  	 return typename Traits::Extract::result_type();	 
  }
  
  static typename Traits::ForceCoreduction::result_type forceCoreductionPair() {
  	 return typename Traits::ForceCoreduction::result_type();
  }

  template<typename ImplT>
  typename Traits::template GetCoreductionPair<typename Traits::template Proxy<ImplT> >::result_type
  getCoreductionPair(typename Traits::template GetCoreductionPair<typename Traits::template Proxy<ImplT> >::argument_type cell)
  {
  	 int times = 0;
  	 BOOST_FOREACH(typename SComplex::ColoredIterators::Iterators::BdCells::iterator::value_type v,
  						complex.iterators(1).bdCells(cell)) {
  		if (times == 0) {
  		  dummyCell3 = v;
  		}
  		++times;
  		if (times == 2) {
		  break;
  		}
  	 }

  	 if (times == 1) {
  		return typename Traits::template GetCoreductionPair<typename Traits::template Proxy<ImplT> >::result_type(dummyCell3);
  	 }
  	 return typename Traits::template GetCoreductionPair<typename Traits::template Proxy<ImplT> >::result_type();
  }

  template<typename ImplT>
  typename Traits::template GetReductionPair<typename Traits::template Proxy<ImplT> >::result_type
  getReductionPair(typename Traits::template GetReductionPair<typename Traits::template Proxy<ImplT> >::argument_type cell)
  {
  	 int times = 0;
  	 BOOST_FOREACH(typename SComplex::ColoredIterators::Iterators::CbdCells::iterator::value_type v,
  						complex.iterators(1).cbdCells(cell)) {
  		if (times == 0) {
  		  dummyCell2 = v;
  		}
  		++times;
  		if (times == 2) {
  		  break;
  		}
  	 }

  	 if (times == 1) {
  		return typename Traits::template GetReductionPair<typename Traits::template Proxy<ImplT> >::result_type(dummyCell2);
  	 }
  	 return typename Traits::template GetReductionPair<typename Traits::template Proxy<ImplT> >::result_type();
  }


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
  Cell dummyCell2, dummyCell3;
};

template<typename SComplexT>
class DefaultReduceStrategy: public DefaultReduceStrategyBase<SComplexT> {

public:
  DefaultReduceStrategy(SComplexT& _complex): DefaultReduceStrategyBase<SComplexT>(_complex) {}
};

#endif

