#include <iostream>
using namespace std;

#include <capd/auxil/Stopwatch.h>
#include <SComplex.hpp>
#include <SComplexDefaultTraits.hpp>
#include <SComplexAlgs.hpp>
#include <SComplexBuilderFromSimplices.hpp>

#include <SimplexSubdivision.hpp>


#include <boost/test/unit_test.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/assign/list_inserter.hpp>


#include <algorithm>


#include <capd/auxil/Stopwatch.h>
#include <capd/auxil/CRef.h>
#include <capd/homologicalAlgebra/embeddingDim.h>

#include <capd/vectalg/MatrixSlice.h>
#include <capd/matrixAlgorithms/intMatrixAlgorithms.hpp>

#include <capd/homologicalAlgebra/homologicalAlgebra.hpp>
#include <capd/homologicalAlgebra/homAlgFunctors.hpp>
#include <capd/homologicalAlgebra/cubSetFunctors.hpp>
#include <capd/homologicalAlgebra/ReducibleFreeChainComplex.hpp>


using namespace boost;
using namespace boost::assign;


typedef ElementaryCell ElementaryCellType;
typedef int ScalarType;
typedef FreeModule<int,capd::vectalg::Matrix<int,0,0> > FreeModuleType;
typedef FreeChainComplex<FreeModuleType> FreeChainComplexType;
typedef ReducibleFreeChainComplex<FreeModuleType,int> ReducibleFreeChainComplexType;

BOOST_AUTO_TEST_SUITE(SComplexSuite)


BOOST_AUTO_TEST_CASE(coreduction_line) {
  typedef SComplex<SComplexDefaultTraits> Complex;
  
  Complex::Dims dims = list_of(0)(1)(0)(1)(0);
  Complex::KappaMap kappaMap = tuple_list_of(1, 0, 1)(1, 2, 1)(3, 2, 1)(3, 4, 1); // .-.-.
  
  Complex complex(3, dims, kappaMap, 1);

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
  
  Complex complex(3, dims, kappaMap, 1);

  (CoreductionAlgorithmFactory::createDefault(complex))();

  BOOST_CHECK( ++(complex.iterators(1).dimCells((Complex::Dim)1).begin()) == complex.iterators(1).dimCells((Complex::Dim)1).end());
  BOOST_CHECK( complex.iterators(1).dimCells((Complex::Dim)0).begin() == complex.iterators(1).dimCells((Complex::Dim)0).end());
}

BOOST_AUTO_TEST_CASE(coreduction_fullTriangle) {
  typedef SComplex<SComplexDefaultTraits> Complex;
  
  Complex::Dims dims = list_of(0)(1)(0)(1)(0)(1)(2);
  Complex::KappaMap kappaMap = tuple_list_of(1, 0, 1)(1, 2, 1)(3, 2, 1)(3, 4, 1)(5, 4, 1)(5, 0, 1)
	 (6, 1, 1)(6, 3, 1)(6, 5, 1); 
  
  Complex complex(3, dims, kappaMap, 1);

  (CoreductionAlgorithmFactory::createDefault(complex))();

  BOOST_FOREACH(Complex::Iterators::AllCells::iterator::value_type v,
					 complex.iterators().allCells()) {
	 BOOST_CHECK_EQUAL(v.getColor(), (Complex::Color)2);
  }

  BOOST_CHECK(complex.iterators(1).allCells().begin() == complex.iterators(1).allCells().end());  
}

BOOST_AUTO_TEST_CASE(coreduction_simplicialFullTriangle) {
  typedef SComplex<SComplexDefaultTraits> Complex;
  
  SComplexBuilderFromSimplices<long, SComplexDefaultTraits> builder(3);
  std::set<std::vector<int> > simplices;

  insert(simplices)( list_of(0)(1)(2) );

	 
  boost::shared_ptr<Complex> complex = builder(simplices, 3, 1);

  BOOST_CHECK_EQUAL(complex->cardinality(), (size_t)7);

  (CoreductionAlgorithmFactory::createDefault(*complex))();

  BOOST_FOREACH(Complex::Iterators::AllCells::iterator::value_type v,
  					 complex->iterators().allCells()) {
  	 BOOST_CHECK_EQUAL(v.getColor(), (Complex::Color)2);
  }

  BOOST_CHECK(complex->iterators(1).allCells().begin() == complex->iterators(1).allCells().end());  
}

template<typename TraitsT>
std::string reduction(SComplex<TraitsT>& complex) {
  (CoreductionAlgorithmFactory::createDefault(complex))();
  
  CRef<ReducibleFreeChainComplexType> RFCComplexCR=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<SComplex<TraitsT>, ReducibleFreeChainComplexType>(complex))();
  CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);

  return homSignCR().bettiVector();
}

template<typename T>
std::string reduction(const T& simplices) {
  typedef SComplex<SComplexDefaultTraits> Complex;  
  SComplexBuilderFromSimplices<long, SComplexDefaultTraits> builder(3);

  boost::shared_ptr<Complex> complex = builder(simplices, 3, 1);
  return reduction(*complex);
}

BOOST_AUTO_TEST_CASE(coreduction_simplicialEmptyTetrahedron) {
  typedef SComplex<SComplexDefaultTraits> Complex;  
  SComplexBuilderFromSimplices<long, SComplexDefaultTraits> builder(3);
  std::set<std::vector<int> > simplices;

  insert(simplices)( list_of(0)(1)(2)(3) );
	 
  boost::shared_ptr<Complex> complex = builder(simplices, 3, 1);
  complex->iterators().dimCells(3).begin()->setColor(0);
  
  BOOST_CHECK_EQUAL(complex->cardinality(), (size_t) 15);
  BOOST_CHECK(complex->iterators(1).dimCells(3).begin() == complex->iterators(1).dimCells(3).end());
  BOOST_CHECK_EQUAL(reduction(*complex), "0,0,1");
}

BOOST_AUTO_TEST_CASE(coreduction_simplicialTorus) {
  typedef SComplex<SComplexDefaultTraits> Complex;
  
  SComplexBuilderFromSimplices<long, SComplexDefaultTraits> builder(3);
  std::vector<std::set<int> > simplices = makeSpaceFromWelds(makeTorusWelds());

  for (int i = 0; i < 3; i++) {
  	 simplices = subdivide3(simplices);
  }

  BOOST_CHECK_EQUAL(reduction(simplices), "0,2,1");
}

BOOST_AUTO_TEST_SUITE_END()

