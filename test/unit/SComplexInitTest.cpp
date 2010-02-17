#include <iostream>
using namespace std;

#include <SComplex.hpp>
#include <util/ColorListNeighbourModel.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/assign/list_of.hpp>

#include <algorithm>

using namespace boost;
using namespace boost::assign;

BOOST_AUTO_TEST_SUITE(SComplexSuite)

BOOST_AUTO_TEST_CASE(init) {
  typedef SComplex<Util::Neighbours::ColorListNeighbourModel> Complex;
  const int colors = 3;
  
  std::vector<std::pair<Complex::Id, Complex::Id> > pairs = pair_list_of(0, 1);  
  Complex complex(colors, pairs.size() * 10, pairs);

  BOOST_CHECK_EQUAL(complex.cardinality(), pairs.size() * 10);

  size_t tmpSize = 0;
  for (Complex::Iterators::AllCells::iterator it= complex.iterators().allCells().begin(),
			end = complex.iterators().allCells().end(); it != end; ++it) {
	 ++tmpSize;
  }
  
  BOOST_CHECK_EQUAL(complex.cardinality(), tmpSize);  
}

BOOST_AUTO_TEST_SUITE_END()
