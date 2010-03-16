#ifndef SCOMPLEX_ALGS_CUBSCOMPLEX_REDUCE_STRATEGY_HPP_
#define SCOMPLEX_ALGS_CUBSCOMPLEX_REDUCE_STRATEGY_HPP_


#include <capd/homologicalAlgebra/embeddingDim.h>

#include <capd/vectalg/MatrixSlice.h>
#include <capd/matrixAlgorithms/intMatrixAlgorithms.hpp>

#include <capd/homologicalAlgebra/homologicalAlgebra.hpp>
#include <capd/homologicalAlgebra/homAlgFunctors.hpp>
#include <capd/homologicalAlgebra/cubSetFunctors.hpp>
#include <capd/homologicalAlgebra/ReducibleFreeChainComplex.hpp>

#include <capd/matrixAlgorithms/matrixAlgorithmsLib.h>
#include <capd/homologicalAlgebra/homologicalAlgebra.hpp>
#include <capd/homologicalAlgebra/homAlgFunctors.hpp>
#include <capd/homologicalAlgebra/cubSetFunctors.hpp>
#include <capd/homologicalAlgebra/ReducibleFreeChainComplex.hpp>

#include "CellProxy.hpp"
#include "SComplexAlgs_DefaultReduceStrategy.hpp"
#include "CubSComplex.hpp"
#include <capd/vectalg/MatrixSlice.h>

template<>
class ReduceStrategyTraits<CubSComplex> {
public:
   template<typename>
   class ReductionPair;

  template<typename ImplT>
  class ReductionPair<CubSComplex::CubCellProxy<ImplT> > {
  public:
	 typedef CellProxy<ImplT*> first_type;
	 typedef CellProxy<CubSComplex::DynamicCell::Impl*> second_type;
	 typedef std::pair<first_type, second_type> type;
  };

  template<typename ImplT>
  class CoreductionPair: public ReductionPair<ImplT> {
  };  
};

template<>
class DefaultReduceStrategy<CubSComplex>: public DefaultReduceStrategyBase<CubSComplex> {

  CubSComplex::DynamicCell dynamicCell;  
public:
  
  DefaultReduceStrategy(CubSComplex& complex): DefaultReduceStrategyBase<CubSComplex>(complex), dynamicCell(complex) {}
  
  // template<typename ImplT>
  // boost::optional<typename Traits::ReductionPair<ImplT>::type> getCoreductionPair(const CubSComplex::CubCellProxy<ImplT>& cell) {
  // 	 if (complex.getUniqueFace(cell, dynamicCell)) {
  // 		//	  	return std::make_pair(boost::ref(dummyCell3), boost::ref(cell));
  // 		return std::make_pair(dummyCell3, cell);
  // 	 } else {
  // 		return boost::optional<ReductionPair>();
  // 	 }
  // }

  template<typename ImplT>
  boost::optional<typename Traits::ReductionPair<CubSComplex::CubCellProxy<ImplT> >::second_type> getReductionPair(const CubSComplex::CubCellProxy<ImplT>& cell) {
  	 if (complex.getUniqueCoFace(cell, dynamicCell)) {
  		// return std::make_pair(typename Traits::ReductionPair<CubSComplex::CubCellProxy<ImplT> >::first_type(cell.getImpl()),
		// 							 typename Traits::ReductionPair<CubSComplex::CubCellProxy<ImplT> >::second_type(dynamicCell.getImpl()));
		return typename Traits::ReductionPair<CubSComplex::CubCellProxy<ImplT> >::second_type(dynamicCell.getImpl());
  	 } else {
		return boost::optional<typename Traits::ReductionPair<CubSComplex::CubCellProxy<ImplT> >::second_type>();
  		//return typename Traits::ReductionPair<CubSComplex::CubCellProxy<ImplT> >::type();
  	 }
  }

  template<typename ImplT>
  void reduceIfPossible(CubSComplex::CubCellProxy<ImplT>& cell) {
  	 if (complex.getUniqueCoFace(cell, dynamicCell)) {
  		reduce(cell);
  		reduce(dynamicCell);
  		//return dummyCell2;
  	 } else {
  		//return boost::optional<Cell>();
  	 }
  }
  
};

#endif
