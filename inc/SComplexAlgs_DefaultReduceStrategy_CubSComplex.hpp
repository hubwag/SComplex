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

#include <functional>

template<>
class ReduceStrategyTraits<CubSComplex> {
public:

  template<typename ImplT>
  class Proxy: public CubSComplex::CubCellProxy<ImplT> {
  public:	 
	 Proxy(const CubSComplex::CubCellProxy<ImplT>& o): CubSComplex::CubCellProxy<ImplT>(o) {}
	 Proxy(const ImplT& impl): CubSComplex::CubCellProxy<ImplT>(impl) {}
  };

  template<typename ImplT>
  static Proxy<ImplT*> makeProxy(const CubSComplex::CubCellProxy<ImplT>& impl) {
	 return Proxy<ImplT*>(impl.getImpl());
  }
  
  template<typename>
  class GetReductionPair;  
  template<typename>
  class GetCoreductionPair;

  template<typename ImplT>
  struct GetReductionPair<typename CubSComplex::CubCellProxy<ImplT> >:  std::unary_function<const typename CubSComplex::CubCellProxy<ImplT>&,
																																 boost::optional<Proxy<CubSComplex::DynamicCell::Impl*> > > {};
  template<typename ImplT>
  struct GetCoreductionPair<typename CubSComplex::CubCellProxy<ImplT> >:  std::unary_function<const typename CubSComplex::CubCellProxy<ImplT>&,
																																	boost::optional<Proxy<CubSComplex::DynamicCell::Impl*> > > {};  
  struct ForceCoreduction {
	 typedef boost::optional<std::pair<Proxy<CubSComplex::DynamicCell::Impl*>,
												  Proxy<CubSComplex::DynamicCell::Impl*> > > result_type;
  };

  struct Extract {
	 typedef boost::optional<Proxy<CubSComplex::Cell> >  result_type;
  };

};

template<>
class DefaultReduceStrategy<CubSComplex>: public DefaultReduceStrategyBase<CubSComplex> {

  CubSComplex::DynamicCell dynamicCell;  
public:
  
  DefaultReduceStrategy(CubSComplex& complex): DefaultReduceStrategyBase<CubSComplex>(complex), dynamicCell(complex) {}
  
  template<typename ImplT>
  typename Traits::GetCoreductionPair<typename CubSComplex::CubCellProxy<ImplT> >::result_type
  //  getCoreductionPair(typename Traits::GetCoreductionPair<typename CubSComplex::CubCellProxy<ImplT> >::argument_type cell)
  getCoreductionPair(const typename CubSComplex::CubCellProxy<ImplT>& cell)
  {  
  	 if (complex.getUniqueFace(cell, dynamicCell)) {
		return typename Traits::GetCoreductionPair<typename CubSComplex::CubCellProxy<ImplT> >::result_type(dynamicCell.getImpl());
  	 } else {
		return typename Traits::GetCoreductionPair<typename CubSComplex::CubCellProxy<ImplT> >::result_type();
  	 }
  }

  template<typename ImplT>
  typename Traits::GetReductionPair<typename CubSComplex::CubCellProxy<ImplT> >::result_type
  getReductionPair(typename Traits::GetReductionPair<typename CubSComplex::CubCellProxy<ImplT> >::argument_type cell)
  {
	 if (complex.getUniqueCoFace(cell, dynamicCell)) {
		return typename Traits::GetReductionPair<typename CubSComplex::CubCellProxy<ImplT> >::result_type(dynamicCell.getImpl());
  	 } else {
		return typename Traits::GetReductionPair<typename CubSComplex::CubCellProxy<ImplT> >::result_type();
  	 }
  }

    
};

#endif
