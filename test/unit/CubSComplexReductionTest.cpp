#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <iostream>

#include <redHom/complex/cubical/CubSComplex.hpp>
#include <redHom/complex/cubical/CubSComplexReader.hpp>
#include <redHom/algorithm/Algorithms.hpp>
#include <redHom/algorithm/strategy/DefaultReduceStrategy_CubSComplex.hpp>


typedef ElementaryCell ElementaryCellType;
typedef int ScalarType;
typedef FreeModule<int,capd::vectalg::Matrix<int,0,0> > FreeModuleType;
typedef FreeChainComplex<FreeModuleType> FreeChainComplexType;
typedef ReducibleFreeChainComplex<FreeModuleType,int> ReducibleFreeChainComplexType;

BOOST_AUTO_TEST_SUITE(CubSComplex_reductions)


template<typename SComplex>
boost::tuple<int, int, int, std::string>  CrHomS() {
	 boost::tuple<int, int, int, std::string> result;
	 using boost::tuples::get;
	 
  Stopwatch swTot;
  BOOST_TEST_MESSAGE(" --- Reading cubical cellular set ");

  CubSComplexReader<8> reader;
  boost::shared_ptr<SComplex> complex = reader(PROJECT_SOURCE_DIR"data/cubical/qtorus.cub"); 

  get<0>(result) = complex->cardinality();

  Stopwatch swComp,swRed;

  //(ShaveAlgorithmFactory::createDefault(*complex))();  
  BOOST_TEST_MESSAGE(" --- Shave reduced the size to " << complex->cardinality() << " in " << swRed);
  get<1>(result) = complex->cardinality();
  
  Stopwatch swCoRed;

  (*CoreductionAlgorithmFactory::createDefault(*complex))();

  BOOST_TEST_MESSAGE(" --- Coreduction reduced the size to " << complex->cardinality() << " in " << swCoRed);
  get<2>(result) = complex->cardinality();
  
  CRef<ReducibleFreeChainComplexType> RFCComplexCR=
	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<SComplex, ReducibleFreeChainComplexType>(*complex))();
  BOOST_TEST_MESSAGE(" --- RFCC constructed  ");

  CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);
  BOOST_TEST_MESSAGE(" --- Computation completed in " << swComp);
  BOOST_TEST_MESSAGE(" --- Computed homology is: \n\n" << homSignCR());
  get<3>(result) = homSignCR().bettiVector();
  
  BOOST_TEST_MESSAGE(" --- Total computation time is: " << swTot);

  return result;
}


BOOST_AUTO_TEST_CASE(reduction_test) {
  BOOST_CHECK_EQUAL(CrHomS<CubSComplex<8> >(), boost::make_tuple(168, 168, 33, "0,2,1"));
}

BOOST_AUTO_TEST_SUITE_END()
