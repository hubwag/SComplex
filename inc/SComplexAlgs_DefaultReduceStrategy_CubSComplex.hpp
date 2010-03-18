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

#include "SComplexAlgs_DefaultReduceStrategyTraits_CubSComplex.hpp"

template<>
class DefaultReduceStrategy<CubSComplex>: public DefaultReduceStrategyBase<CubSComplex> {

  CubSComplex::DynamicCell dynamicCell;  
public:
  
  DefaultReduceStrategy(CubSComplex& complex): DefaultReduceStrategyBase<CubSComplex>(complex), dynamicCell(complex) {}
  
  template<typename ArgT>
  typename Traits::GetCoreductionPair<ArgT>::result_type
  getCoreductionPair(const ArgT& cell)
  {  
  	 if (complex.getUniqueFace(cell, dynamicCell)) {
		return typename Traits::GetCoreductionPair<ArgT>::result_type(dynamicCell.getImpl());
  	 } else {
		return typename Traits::GetCoreductionPair<ArgT>::result_type();
  	 }
  }

  template<typename ArgT>
  typename Traits::GetReductionPair<ArgT>::result_type
  getReductionPair(const ArgT& cell)
  {
	 if (complex.getUniqueCoFace(cell, dynamicCell)) {
		return typename Traits::GetReductionPair<ArgT>::result_type(dynamicCell.getImpl());
  	 } else {
		return typename Traits::GetReductionPair<ArgT>::result_type();
  	 }
  }

    
};

#endif
