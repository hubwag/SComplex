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
  
  Complex::KappaMap kappaMap = tuple_list_of(0, 1, 1);  
  Complex complex(colors, kappaMap.size() * 10, kappaMap);

  BOOST_CHECK_EQUAL(complex.cardinality(), kappaMap.size() * 10);

  size_t tmpSize = 0;
  for (Complex::Iterators::AllCells::iterator it= complex.iterators().allCells().begin(),
			end = complex.iterators().allCells().end(); it != end; ++it) {
	 ++tmpSize;
  }
  BOOST_CHECK_EQUAL(complex.cardinality(), tmpSize);
  
  const Complex& refComplex = complex;
  
  tmpSize = 0;
  for (Complex::Iterators::AllCells::const_iterator it= refComplex.iterators().allCells().begin(),
			end = refComplex.iterators().allCells().end(); it != end; ++it) {
	 ++tmpSize;
  }
  BOOST_CHECK_EQUAL(refComplex.cardinality(), tmpSize);

}

BOOST_AUTO_TEST_CASE(coboundarySize) {
  typedef SComplex<Util::Neighbours::ColorListNeighbourModel> Complex;
  const int colors = 3;
  
  Complex::KappaMap kappaMap = tuple_list_of(0, 1, 1)(0, 2, 1)(1, 2, 1)(0, 3, 1);
  std::vector<int> cbdSizes = list_of(0)(1)(2)(1);
  
  Complex complex(colors, kappaMap.size(), kappaMap);

  Complex::Cell cells[] = {complex[0], complex[1], complex[2], complex[3]} ;

  std::vector<int> tmpCbdSizes;
  BOOST_FOREACH( Complex::Cell cell, cells) {
	 size_t tmpSize = 0;
	 BOOST_FOREACH(Complex::Iterators::BdCells::iterator::value_type t,
						complex.iterators().cbdCells(cell)) {
		tmpSize++;
	 }
	 tmpCbdSizes.push_back(tmpSize);
  }

  BOOST_CHECK_EQUAL_COLLECTIONS(tmpCbdSizes.begin(), tmpCbdSizes.end(), cbdSizes.begin(), cbdSizes.end());

}


BOOST_AUTO_TEST_CASE(boundarySize) {
  typedef SComplex<Util::Neighbours::ColorListNeighbourModel> Complex;
  const int colors = 3;
  
  Complex::KappaMap kappaMap = tuple_list_of(0, 1, 1)(0, 2, 1)(1, 2, 1)(0, 3, 1);
  std::vector<int> bdSizes = list_of(3)(1)(0)(0);
  
  Complex complex(colors, kappaMap.size(), kappaMap);

  Complex::Cell cells[] = {complex[0], complex[1], complex[2], complex[3]};

  std::vector<int> tmpBdSizes;
  BOOST_FOREACH( Complex::Cell cell, cells) {
	 size_t tmpSize = 0;
	 BOOST_FOREACH(Complex::Iterators::BdCells::iterator::value_type t,
						complex.iterators().bdCells(cell)) {
		tmpSize++;
	 }
	 tmpBdSizes.push_back(tmpSize);
  }

  BOOST_CHECK_EQUAL_COLLECTIONS(tmpBdSizes.begin(), tmpBdSizes.end(), bdSizes.begin(), bdSizes.end());

  tmpBdSizes = std::vector<int>();
  BOOST_FOREACH( Complex::Cell cell, cells) {
	 size_t tmpSize = 0;
	 BOOST_FOREACH(Complex::ConstIterators::BdCells::const_iterator::value_type t,
						((const Complex&)complex).iterators().bdCells(cell)) {
		tmpSize++;
	 }
	 tmpBdSizes.push_back(tmpSize);
  }

  BOOST_CHECK_EQUAL_COLLECTIONS(tmpBdSizes.begin(), tmpBdSizes.end(), bdSizes.begin(), bdSizes.end());

}

BOOST_AUTO_TEST_SUITE_END()
