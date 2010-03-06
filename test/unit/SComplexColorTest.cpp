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

  Complex::Cell cells[] = {complex[0], complex[1], complex[2], complex[3]};


  BOOST_FOREACH( Complex::Cell cell, cells) {
	 BOOST_FOREACH(Complex::Iterators::BdCells::iterator::value_type t,
						complex.iterators().bdCells(cell)) {
		t.setColor(1);
	 }
  }

  const Complex& refComplex = complex;
 
  for (Complex::ConstIterators::AllCells::const_iterator it= refComplex.iterators().allCells().begin(),
			end = refComplex.iterators().allCells().end(); it != end; ++it) {
	 //it->setColor(0);
  }

  
}

BOOST_AUTO_TEST_SUITE_END()
