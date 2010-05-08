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


BOOST_AUTO_TEST_CASE(line) {
  typedef SComplex<SComplexDefaultTraits> Complex;

  Complex::Dims dims = list_of(0)(1)(0)(1)(0);
  Complex::KappaMap kappaMap = tuple_list_of(1, 0, 1)(1, 2, 1)(3, 2, 1)(3, 4, 1); // .-.-.
  
  Complex complex(3, dims, kappaMap, 1);

  COAKQAlgorithm<COAKQStrategy<Complex, Complex> > algorithm(new COAKQStrategy<Complex, Complex>(complex));

  algorithm();

  BOOST_FOREACH(Complex::Iterators::AllCells::iterator::value_type v,
					 complex.iterators().allCells()) {
	 BOOST_CHECK_EQUAL(v.getColor(), (Complex::Color)2);
  }

  BOOST_CHECK(complex.iterators(1).allCells().begin() == complex.iterators(1).allCells().end());  

  BOOST_CHECK_EQUAL(algorithm.getStrategy().getOutputComplex().size(), 1);
}


BOOST_AUTO_TEST_CASE(coreduction_emptyTriangle) {
  typedef SComplex<SComplexDefaultTraits> Complex;
  
  Complex::Dims dims = list_of(0)(1)(0)(1)(0)(1);
  Complex::KappaMap kappaMap = tuple_list_of(1, 0, 1)(1, 2, -1)(3, 2, 1)(3, 4, -1)(5, 4, 1)(5, 0, -1);
  
  Complex complex(3, dims, kappaMap, 1);

  COAKQAlgorithm<COAKQStrategy<Complex, Complex> > algorithm(new COAKQStrategy<Complex, Complex>(complex));

  algorithm();

  BOOST_FOREACH(Complex::Iterators::AllCells::iterator::value_type v,
					 complex.iterators().allCells()) {
	 BOOST_CHECK_EQUAL(v.getColor(), (Complex::Color)2);
  }

  BOOST_CHECK(complex.iterators(1).allCells().begin() == complex.iterators(1).allCells().end());  

  Complex& coAKQ = algorithm.getStrategy().getOutputComplex();
  BOOST_CHECK_EQUAL(coAKQ.size(), 2);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(coAKQ))();
  CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);

  std::string betti = homSignCR().bettiVector();

  BOOST_CHECK_EQUAL(betti, "1,1");

}

BOOST_AUTO_TEST_CASE(coreduction_fullTriangle) {
  typedef SComplex<SComplexDefaultTraits> Complex;


  Complex::Dims dims = list_of(0)(1)(0)(1)(0)(1)(2);
  Complex::KappaMap kappaMap = tuple_list_of(1, 0, -1)(1, 2, 1)(3, 2, -1)(3, 4, 1)(5, 4, -1)(5, 0, 1)
	 (6, 1, -1)(6, 3, 1)(6, 5, -1); 
  
  Complex complex(3, dims, kappaMap, 1);

  COAKQAlgorithm<COAKQStrategy<Complex, Complex> > algorithm(new COAKQStrategy<Complex, Complex>(complex));

  algorithm();

  BOOST_FOREACH(Complex::Iterators::AllCells::iterator::value_type v,
					 complex.iterators().allCells()) {
	 BOOST_CHECK_EQUAL(v.getColor(), (Complex::Color)2);
  }

  BOOST_CHECK(complex.iterators(1).allCells().begin() == complex.iterators(1).allCells().end());  

  Complex& coAKQ = algorithm.getStrategy().getOutputComplex();
  BOOST_CHECK_EQUAL(coAKQ.size(), 2);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(coAKQ))();
  CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);

  std::string betti = homSignCR().bettiVector();

  BOOST_CHECK_EQUAL(betti, "");

}


BOOST_AUTO_TEST_SUITE_END()

