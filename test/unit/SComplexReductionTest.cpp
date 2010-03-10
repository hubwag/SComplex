#include <iostream>
using namespace std;

#include <capd/auxil/Stopwatch.h>
#include <SComplex.hpp>
#include <SComplexDefaultTraits.hpp>
#include <SComplexAlgs_Coreduction.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/assign/list_of.hpp>

#include <algorithm>



using namespace boost;
using namespace boost::assign;

BOOST_AUTO_TEST_SUITE(SComplexSuite)


BOOST_AUTO_TEST_CASE(coreduction_line) {
  typedef SComplex<SComplexDefaultTraits> Complex;
  
  Complex::Dims dims = list_of(0)(1)(0)(1)(0);
  Complex::KappaMap kappaMap = tuple_list_of(1, 0, 1)(1, 2, 1)(3, 2, 1)(3, 4, 1); // .-.-.
  
  Complex complex(3, dims.size(), dims, kappaMap, 1);

  (CoreductionAlgorithmFactory::createDefault(complex))();

  BOOST_FOREACH(Complex::Iterators::AllCells::iterator::value_type v,
					 complex.iterators().allCells()) {
	 BOOST_CHECK_EQUAL(v.getColor(), (Complex::Color)2);
  }

  BOOST_CHECK(complex.iterators(1).allCells().begin() == complex.iterators(1).allCells().end());  
}

BOOST_AUTO_TEST_CASE(coreduction_emptyTriangle) {
  typedef SComplex<SComplexDefaultTraits> Complex;
  
  Complex::Dims dims = list_of(0)(1)(0)(1)(0)(1);
  Complex::KappaMap kappaMap = tuple_list_of(1, 0, 1)(1, 2, 1)(3, 2, 1)(3, 4, 1)(5, 4, 1)(5, 0, 1);
  
  Complex complex(3, dims.size(), dims, kappaMap, 1);

  (CoreductionAlgorithmFactory::createDefault(complex))();

  BOOST_CHECK( ++(complex.iterators(1).dimCells((Complex::Dim)1).begin()) == complex.iterators(1).dimCells((Complex::Dim)1).end());
  BOOST_CHECK( complex.iterators(1).dimCells((Complex::Dim)0).begin() == complex.iterators(1).dimCells((Complex::Dim)0).end());
}

BOOST_AUTO_TEST_CASE(coreduction_fullTriangle) {
  typedef SComplex<SComplexDefaultTraits> Complex;
  
  Complex::Dims dims = list_of(0)(1)(0)(1)(0)(1)(2);
  Complex::KappaMap kappaMap = tuple_list_of(1, 0, 1)(1, 2, 1)(3, 2, 1)(3, 4, 1)(5, 4, 1)(5, 0, 1)
	 (6, 1, 1)(6, 3, 1)(6, 5, 1); 
  
  Complex complex(3, dims.size(), dims, kappaMap, 1);

  (CoreductionAlgorithmFactory::createDefault(complex))();

  BOOST_FOREACH(Complex::Iterators::AllCells::iterator::value_type v,
					 complex.iterators().allCells()) {
	 BOOST_CHECK_EQUAL(v.getColor(), (Complex::Color)2);
  }

  BOOST_CHECK(complex.iterators(1).allCells().begin() == complex.iterators(1).allCells().end());  
}

BOOST_AUTO_TEST_SUITE_END()

