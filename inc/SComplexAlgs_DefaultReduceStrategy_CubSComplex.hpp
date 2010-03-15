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

#include "SComplexAlgs_DefaultReduceStrategy.hpp"
#include "CubSComplex.hpp"
#include <capd/vectalg/MatrixSlice.h>


template<>
class DefaultReduceStrategy<CubSComplex>: public DefaultReduceStrategyBase<CubSComplex> {

  CubSComplex::DynamicCell dynamicCell;
  
public:
  DefaultReduceStrategy(CubSComplex& complex): DefaultReduceStrategyBase<CubSComplex>(complex), dynamicCell(complex) {}
  
  // boost::optional<CoreductionPair> getCoreductionPair(Cell& cell) {
  // 	 if (complex.getUniqueFace(cell, dummyCell3)) {
  // 		//	  	return std::make_pair(boost::ref(dummyCell3), boost::ref(cell));
  // 		return std::make_pair(dummyCell3, cell);
  // 	 } else {
  // 		return boost::optional<ReductionPair>();
  // 	 }
  // }
  
  // boost::optional<ReductionPair> getReductionPair( Cell& cell) {
  // 	 if (complex.getUniqueCoFace(cell, dummyCell2)) {
  // 		//	  	return std::make_pair(boost::ref(cell), boost::ref(dummyCell2));
  // 		return std::make_pair(cell, dummyCell2);
  // 	 } else {
  // 		return boost::optional<ReductionPair>();
  // 	 }
  // }


  void reduce(CubSComplex::TempCell cell) {
	 cell.setColor<2>();
  }

  void reduceIfPossible(CubSComplex::TempCell cell) {
	 if (complex.getUniqueCoFace(cell, dynamicCell)) {
		reduce(cell);
		reduce(CubSComplex::TempCell(dynamicCell));
		//return dummyCell2;
	 } else {
		//return boost::optional<Cell>();
	 }
  }
  
};

#endif
