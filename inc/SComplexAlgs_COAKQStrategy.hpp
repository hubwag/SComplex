#ifndef SCOMPLEX_ALGS_COAKQ_STRATEGY_HPP_
#define SCOMPLEX_ALGS_COAKQ_STRATEGY_HPP_

#include <boost/optional.hpp>
#include <boost/foreach.hpp>
#include <utility>

template<typename SComplexT>
class COAKQStrategyTraits {
public:
  
  template<typename ArgT>
  struct GetReductionPair: public std::unary_function<const ArgT&,
						      boost::optional<typename SComplexT::Cell> > {};

  struct ForceReduction {
    typedef boost::optional<std::pair<typename SComplexT::Cell,
				      typename SComplexT::Cell> > result_type;
  };

  struct Extract {
    typedef boost::optional<typename SComplexT::Cell >  result_type;
  };
};

template<typename SComplexT>
class COAKQStrategyBase {

public:
  typedef SComplexT SComplex;
  typedef COAKQStrategyTraits<SComplex> Traits;
  typedef typename SComplex::Cell Cell;


  COAKQStrategyBase(SComplex& _complex): complex(_complex), dummyCell2(_complex),  dummyCell3(_complex) {}
  
  SComplex& getComplex() const {
	 return complex;
  }
  
  template<typename ImplT>
  static bool reduced(const typename SComplex::template CellProxy<ImplT>& cell) {
  	 return cell.getColor() == 2;
  }

  template<typename ImplT1, typename ImplT2>	 
  static void coreduce(const typename SComplex::template CellProxy<ImplT1>& a, const typename SComplex::template CellProxy<ImplT2>& b)  {
  	 a.template setColor<2>();
  	 b.template setColor<2>();
  }

  template<typename ImplT1, typename ImplT2>	 
  static void reduce(const typename SComplex::template CellProxy<ImplT1>& a, const typename SComplex::template CellProxy<ImplT2>& b)  {
  	 a.template setColor<2>();
  	 b.template setColor<2>();
  }

  template<typename ImplT>
  static void reduce(const typename SComplex::template CellProxy<ImplT>& cell) {
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
  
  static typename Traits::ForceReduction::result_type forceReductionPair() {
  	 return typename Traits::ForceReduction::result_type();
  }

  template<typename ArgT>
  typename Traits::template GetReductionPair<ArgT>::result_type
  getReductionPair(const ArgT& cell)
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
  		return typename Traits::template GetReductionPair<ArgT>::result_type(dummyCell3);
  	 }
  	 return typename Traits::template GetReductionPair<ArgT>::result_type();
  }

  template<typename ImplT>
  typename Traits::template GetReductionPair<typename SComplex::template CellProxy<ImplT> >::result_type
  getReductionPair(typename Traits::template GetReductionPair<typename SComplex::template CellProxy<ImplT> >::argument_type cell)
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
  		return typename Traits::template GetReductionPair<typename SComplex::template CellProxy<ImplT> >::result_type(dummyCell2);
  	 }
  	 return typename Traits::template GetReductionPair<typename SComplex::template CellProxy<ImplT> >::result_type();
  }


  size_t getMaxDim() {
	 typename SComplex::Dim maxDim = 0;
	 for (typename SComplex::ColoredIterators::Iterators::AllCells::iterator it = complex.template iterators<1>().allCells().begin(),
			  end = complex.template iterators<1>().allCells().end();
			it != end; ++it) {

		maxDim = std::max(maxDim, (*it).getDim());
	 }
	 
	 return maxDim;
  }

protected:
  SComplex& complex;
  Cell dummyCell2, dummyCell3;
};

template<typename SComplexT>
class COAKQStrategy: public COAKQStrategyBase<SComplexT> {

public:
  COAKQStrategy(SComplexT& _complex): COAKQStrategyBase<SComplexT>(_complex) {}
};

#endif

