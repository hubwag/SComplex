#ifndef _SCOMPLEX_ALGS_COREDUCTION_HPP
#define _SCOMPLEX_ALGS_COREDUCTION_HPP

#include "SComplexAlgs_DefaultReduceStrategy.hpp"
#include <deque>

template<typename StrategyT>
class CoreductionAlgorithm {

public:
  typedef StrategyT Strategy;

  typedef typename Strategy::SComplex SComplex;
  typedef typename Strategy::Cell Cell;  
  
  CoreductionAlgorithm(Strategy* _strategy): strategy(_strategy) {}

  ~CoreductionAlgorithm() {
	 if (strategy)
		delete strategy;
  }
  
  int operator()();

private:
  template<typename ImplT>
  void storeGenerator(const typename Strategy::Traits::template Proxy<ImplT>& a) {
	 collectedHomGenerators.push_back(a);
  }

  template<typename ImplT1, typename ImplT2>
  void storeCoreductionPair(const typename Strategy::Traits::template Proxy<ImplT1>& a, const typename Strategy::Traits::template Proxy<ImplT2>& b) {}
  
  bool coreduceNextPair();  

  template<typename ImplT>
  void addCellsToProcess(const typename Strategy::Traits::template Proxy<ImplT>& sourceFace) {
	 // Finally, put all present cofaces of the source face
	 // into the queue
	 typename SComplex::ColoredIterators::Iterators::CbdCells cbdCells = strategy->getComplex().template iterators<1>().cbdCells(sourceFace);
	 for (typename SComplex::ColoredIterators::Iterators::CbdCells::iterator cbdn = cbdCells.begin(),
			  end = cbdCells.end(); cbdn != end; ++cbdn) {
		cellsToProcess.push_back(*cbdn);
	 }
  }


  template<typename ImplT1, typename ImplT2>
  void doCoreduction(const typename Strategy::Traits::template Proxy<ImplT1>& a, const typename Strategy::Traits::template Proxy<ImplT2>& b) {
	 storeCoreductionPair(a, b);
	 addCellsToProcess(a);		
	 strategy->coreduce(a, b);
  }
  
  Strategy* strategy;
  
  std::vector<Cell> collectedHomGenerators;
  std::deque<Cell> cellsToProcess;
};

class CoreductionAlgorithmFactory {

public:
  template<typename SComplex>
  static CoreductionAlgorithm<DefaultReduceStrategy<SComplex> > createDefault(SComplex& s) {
	 return CoreductionAlgorithm<DefaultReduceStrategy<SComplex> >(new DefaultReduceStrategy<SComplex>(s));
  }

};


template<typename StrategyT>
inline bool CoreductionAlgorithm<StrategyT>::coreduceNextPair() {
  
  while (! cellsToProcess.empty() ) {
	 Cell& cell = cellsToProcess.front();
	 
	 if (! strategy->reduced(cell)) {
		typename StrategyT::Traits::template GetCoreductionPair<Cell>::result_type coreductionPair = strategy->getCoreductionPair(cell);
		
		if (coreductionPair) {		
		  doCoreduction(*coreductionPair, cell);
		  cellsToProcess.pop_front();
		  return true;
		} else {
		  addCellsToProcess(cell);
		}
	 }
	 cellsToProcess.pop_front();	 
  }

  // Originally there are no candidates in the queue
  // They may also disappear when a connected component
  // is exhausted
  // If we know that a coreduction may be there,
  // For instance when treating a non-compact set
  typename StrategyT::Traits::ForceCoreduction::result_type force = strategy->forceCoreductionPair();
  if (force) {
	 doCoreduction(force->first, force->second);
	 return true;
  }
  return false;
}


template<typename StrategyT>
inline int CoreductionAlgorithm<StrategyT>::operator()(){
  int cnt=0;

  for(;;){
	 if (coreduceNextPair()) {
		++cnt;++cnt;
	 } else {
		// If the search failed or when we even did not try to search
		// and we know that a cell of lowest dimension is always
		// a homology generator like in the case of a vertex in
		// a compact set, we just pick up such a cell and
		// remove it from the complex
		typename StrategyT::Traits::Extract::result_type sourceFace = strategy->extract();

		if(sourceFace){
		  storeGenerator(*sourceFace);		  
		  addCellsToProcess(*sourceFace);
		  
		  strategy->reduce(*sourceFace);
		  ++cnt;
		}else{
		  break; // no base face left: quit any further processing
		}
	 } 
  }
  
  return cnt; // the number of cells removed
}

#endif
