#include <SComplexAlgs_COAKQ.hpp>

#include <SComplex.hpp>
#include <SComplexDefaultTraits.hpp>

#include <SComplexAlgs.hpp>
#include <SComplexBuilderFromSimplices.hpp>

#include <SimplexSubdivision.hpp>


#include <boost/test/unit_test.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/assign/list_of.hpp>

#include <algorithm>

using namespace boost;
using namespace boost::assign;

typedef ElementaryCell ElementaryCellType;
typedef int ScalarType;
typedef FreeModule<int,capd::vectalg::Matrix<int,0,0> > FreeModuleType;
typedef FreeChainComplex<FreeModuleType> FreeChainComplexType;
typedef ReducibleFreeChainComplex<FreeModuleType,int> ReducibleFreeChainComplexType;


BOOST_AUTO_TEST_SUITE(COAKQSuite)


BOOST_AUTO_TEST_CASE(init) {
  typedef SComplex<SComplexDefaultTraits> Complex;
  
  std::vector<std::set<int> > simplices = makeSpaceFromWelds(makeTorusWelds());

  for (int i = 0; i < 3; i++) {
  	 simplices = subdivide3(simplices);
  }

  SComplexBuilderFromSimplices<long, SComplexDefaultTraits> builder(300);
  boost::shared_ptr<Complex> complex = builder(simplices, 3, 1);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR=
    (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(*complex))();
  CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);
  std::string betti = homSignCR().bettiVector();

  BOOST_CHECK_EQUAL(betti, "0,2,1");


}

BOOST_AUTO_TEST_SUITE_END()

