#include <iostream>
using namespace std;

#include <capd/auxil/Stopwatch.h>
#include <SComplex.hpp>
#include <util/ColorListNeighbourModel.hpp>
#include <SComplexAlgs_Coreduction.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/assign/list_of.hpp>

#include <algorithm>



using namespace boost;
using namespace boost::assign;

BOOST_AUTO_TEST_SUITE(SComplexSuite)


BOOST_AUTO_TEST_CASE(coreduction) {
  typedef SComplex<Util::Neighbours::ColorListNeighbourModel> Complex;
  const int size = 300;
  const int colors = 2;

  Complex::Dims dims(size);

  Complex complex(colors, size, dims);

  Stopwatch swCoRed;

  (CoreductionAlgorithmFactory::createDefault(complex))();
  BOOST_TEST_MESSAGE(" --- Coreduction reduced the size to " << complex.cardinality() << " in " << swCoRed);

}


BOOST_AUTO_TEST_SUITE_END()

