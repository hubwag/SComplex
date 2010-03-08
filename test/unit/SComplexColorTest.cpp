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
  const int size = 300;
  const int colors = size;
    
  Complex complex(colors, size);

  for (int i = 0; i < size; ++i) {
	 complex[i].setColor(i);
  }

  for (int i = 0; i < size; ++i) {
	 BOOST_CHECK_EQUAL(complex[i].getColor(), (Complex::Color)i);
  }

  int i = 0;
  BOOST_FOREACH(Complex::Iterators::AllCells::iterator::value_type v,
					 complex.iterators().allCells()) {
	 BOOST_CHECK_EQUAL(v.getColor(), (Complex::Color) i);
	 ++i;
  }
  
  for (int i = 0; i < size; ++i) {
	 BOOST_CHECK_EQUAL( complex.iterators((Complex::Color)i).allCells().begin()->getId(), i);
  }
}

BOOST_AUTO_TEST_CASE(setColorCheckNeighbours) {
  typedef SComplex<Util::Neighbours::ColorListNeighbourModel> Complex;
  const int colors = 3;
  
  Complex::KappaMap kappaMap = tuple_list_of(0, 1, 1)(0, 2, 1)(1, 2, 1)(0, 3, 1);
  std::vector<int> bdSizes = list_of(3)(1)(0)(0);
  
  Complex complex(colors, kappaMap.size(), kappaMap);

  const int size = 4;
  Complex::Cell cells[size] = {complex[0], complex[1], complex[2], complex[3]};

  for (int i = 0; i < size; ++i) {
	 cells[i].setColor(i);
  }

  BOOST_FOREACH(Complex::Iterators::BdCells::iterator::value_type t,
					 complex.iterators().bdCells(cells[0])) {
	 //	 BOOST_CHECK_EQUAL(
  }
    
}

BOOST_AUTO_TEST_SUITE_END()
