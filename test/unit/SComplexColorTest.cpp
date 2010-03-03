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


BOOST_AUTO_TEST_CASE(setColor) {
  typedef SComplex<Util::Neighbours::ColorListNeighbourModel> Complex;
  const int colors = 3;
    
  Complex complex(colors, 3);

  complex[0].setColor(0);
  complex[1].setColor(1);
  complex[2].setColor(2);

  BOOST_CHECK_EQUAL(complex[0].getColor(), (Complex::Color)0);
  BOOST_CHECK_EQUAL(complex[1].getColor(), (Complex::Color)1);
  BOOST_CHECK_EQUAL(complex[2].getColor(), (Complex::Color)2);
}

BOOST_AUTO_TEST_CASE(setColorCheckNeighbours) {
  typedef SComplex<Util::Neighbours::ColorListNeighbourModel> Complex;
  const int colors = 3;
  
  Complex::KappaMap kappaMap = tuple_list_of(0, 1, 1)(0, 2, 1)(1, 2, 1)(0, 3, 1);
  std::vector<int> bdSizes = list_of(3)(1)(0)(0);
  
  Complex complex(colors, kappaMap.size(), kappaMap);

  boost::reference_wrapper<Complex::Cell> cells[] = {boost::ref(complex[0]), boost::ref(complex[1]), boost::ref(complex[2]), boost::ref(complex[3])};

}

BOOST_AUTO_TEST_SUITE_END()
