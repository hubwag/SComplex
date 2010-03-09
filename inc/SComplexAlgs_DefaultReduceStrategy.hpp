#ifndef SCOMPLEX_ALGS_DEFAULT_REDUCE_STRATEGY_HPP_
#define SCOMPLEX_ALGS_DEFAULT_REDUCE_STRATEGY_HPP_

#include <boost/optional.hpp>
#include <boost/foreach.hpp>


template<typename SComplexT>
class DefaultReduceStrategyBase {

public:
  typedef SComplexT SComplex;
  typedef typename SComplex::Cell Cell;
  typedef std::pair<boost::reference_wrapper<Cell>, boost::reference_wrapper<Cell> > CoreductionPair;
  typedef std::pair<boost::reference_wrapper<Cell>, boost::reference_wrapper<Cell> > ReductionPair;  


  DefaultReduceStrategyBase(SComplex& _complex): complex(_complex), dummyCell1(_complex), dummyCell2(_complex),  dummyCell3(_complex) {}
  
  SComplex& getComplex() const {
	 return complex;
  }

  Cell& getFace(const CoreductionPair& coRedPair) {
	 return boost::unwrap_ref(coRedPair.first);
  }
  
  bool reduced(const Cell& cell) const {
	 return cell.getColor() == 2;
  }

  void coreduce(const CoreductionPair& coRedPair) const {
	 boost::unwrap_ref(coRedPair.first).template setColor<2>();
	 boost::unwrap_ref(coRedPair.second).template setColor<2>();
  }

  void reduce(const ReductionPair& redPair) const {
	 boost::unwrap_ref(redPair.first).template setColor<2>();
	 boost::unwrap_ref(redPair.second).template setColor<2>();
  }

  void reduce(Cell& cell) {
	 cell.template setColor<2>();
  }
  
  boost::optional<Cell&> extract() {
	 typename SComplex::ColoredIterators::Iterators::DimCells::iterator end = complex.iterators(1).dimCells(0).end(),
		it = complex.iterators(1).dimCells(0).begin();

	 if (it != end) { 
		dummyCell1 = *it;
		return dummyCell1;
	 }
	 return boost::optional<Cell&>();	 
  }
  
  boost::optional<CoreductionPair> forceCoreductionPair() {
	 return boost::optional<CoreductionPair>();
  }

  boost::optional<CoreductionPair> getCoreductionPair(Cell& cell) {
	 int times = 0;
	 BOOST_FOREACH(typename SComplex::ColoredIterators::Iterators::BdCells::iterator::value_type v,
						complex.iterators(1).bdCells(cell)) {
		if (times == 0) {
		  dummyCell3 = v;
		}
		++times;
		if (times == 2) {
		  return boost::optional<ReductionPair>();
		}
	 }

	 if (times == 1) {
		return std::make_pair(boost::ref(dummyCell3), boost::ref(cell));
	 }
	 return boost::optional<ReductionPair>();
  }
  
  boost::optional<ReductionPair> getReductionPair(Cell& cell) {
	 int times = 0;
	 BOOST_FOREACH(typename SComplex::ColoredIterators::Iterators::CbdCells::iterator::value_type v,
						complex.iterators(1).cbdCells(cell)) {
		if (times == 0) {
		  dummyCell2 = v;
		}
		++times;
		if (times == 2) {
		  return boost::optional<ReductionPair>();
		}
	 }

	 if (times == 1) {
		return std::make_pair(boost::ref(cell), boost::ref(dummyCell2));
	 }
	 return boost::optional<ReductionPair>();
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

