#ifndef SCOMPLEX_ALGS_HPP_
#define SCOMPLEX_ALGS_HPP_


#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/optional.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>


#include <utility>
#include <algorithm>

#include "SComplexAlgs_Coreduction.hpp"
#include "SComplexAlgs_Shave.hpp"


template<typename SComplexT, typename ReducibleFreeChainComplexT>
class ReducibleFreeChainComplexOverZFromSComplexAlgorithm {

public:
  typedef SComplexT SComplex;
  typedef typename SComplex::Cell Cell;

  typedef ReducibleFreeChainComplexT ReducibleFreeChainComplex;
  
  class SComplexChainCell {
	 SComplex& complex;
	 const Cell cell;
	 const int embededDim;
	 
  public:
	 SComplexChainCell(SComplex& _complex, const Cell& _cell, int _embededDim): complex(_complex), cell(_cell), embededDim(_embededDim) {
		
	 }

	 int embDim() const {
		return this->embededDim;
	 }

	 int ownDim() const {
		return cell.getDim();
	 }

	 void boundary(std::map<SComplexChainCell,int>& A_boundary) const {
		for (typename SComplex::ColoredIterators::Iterators::BdCells::iterator it = complex.template iterators<1>().bdCells(this->cell).begin(),
				 end = complex.template iterators<1>().bdCells(this->cell).end();
			  it != end; ++it) {
		  A_boundary.insert(std::make_pair(SComplexChainCell(complex, *it, embededDim),
													  //complex.coincidenceIndex(*it, this->cell)
													  complex.coincidenceIndex(this->cell, *it)
													  ));
		}
	 }

	 bool operator<(const SComplexChainCell& b) const {
		return this->cell < b.cell;
	 }
  };
  
public:
  

  ReducibleFreeChainComplexOverZFromSComplexAlgorithm(SComplex& _s): s(_s) {}
  CRef<ReducibleFreeChainComplex> operator()();
  
private:
  SComplex& s;
};


template<typename SComplexT, typename ReducibleFreeChainComplexT>
inline CRef<ReducibleFreeChainComplexT> ReducibleFreeChainComplexOverZFromSComplexAlgorithm<SComplexT, ReducibleFreeChainComplexT>::operator()(){

  Stopwatch sw;
  std::set<SComplexChainCell> cells;

  size_t maxDim = 0;
  for (typename SComplex::ColoredIterators::Iterators::AllCells::iterator it = s.template iterators<1>().allCells().begin(),
			end = s.template iterators<1>().allCells().end();
		 it != end; ++it) {

	 maxDim = std::max(maxDim, it->getDim());
  }
  
  for (typename SComplex::ColoredIterators::Iterators::AllCells::iterator it = s.template iterators<1>().allCells().begin(),
			end = s.template iterators<1>().allCells().end();
		 it != end; ++it) {
	 cells.insert(SComplexChainCell(s, *it, maxDim));
  }

  CRef<ReducibleFreeChainComplex> rfccCR( new ReducibleFreeChainComplex(cells));

  fcout << "Reducible chain complex (over Z) construction of CubCelSet completed in " << sw  << std::endl;
  return rfccCR;
}


#endif // SCOMPLEX_ALGS_HPP_
