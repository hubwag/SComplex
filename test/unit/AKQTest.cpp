#define BOOST_TEST_MAIN

#include <capd/auxil/ofstreamcout.h>
ofstreamcout fcout;

#include <redHom/algorithm/Coreduction.hpp>
#include <redHom/algorithm/AKQStrategy.hpp>
#include <redHom/algorithm/Algorithms.hpp>

#include <redHom/complex/scomplex/SComplex.hpp>
#include <redHom/complex/scomplex/SComplexDefaultTraits.hpp>
#include <redHom/complex/scomplex/SComplexBuilderFromSimplices.hpp>

#include <redHom/complex/cubical/CubSComplex.hpp>
#include <redHom/complex/cubical/CubSComplexReader.hpp>


#include <redHom/complex/simplicial/SimplexIO.hpp>
#include <redHom/SimplexSubdivision.hpp>

#include <boost/test/included/unit_test.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/assign/list_of.hpp>

#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>

using namespace boost;
using namespace boost::assign;
using namespace std;

typedef ElementaryCell ElementaryCellType;
typedef int ScalarType;
typedef FreeModule<int,capd::vectalg::Matrix<int,0,0> > FreeModuleType;
typedef FreeChainComplex<FreeModuleType> FreeChainComplexType;
typedef ReducibleFreeChainComplex<FreeModuleType,int> ReducibleFreeChainComplexType;

BOOST_AUTO_TEST_SUITE(COAKQSuite)

template<typename Complex>
void performTest(Complex &complex, const string &correctBettiVector, int correctSize)
{
    CoreductionAlgorithm<AKQReduceStrategy<Complex> > algorithm(new AKQReduceStrategy<Complex>(complex));

    algorithm();

    BOOST_FOREACH(typename Complex::Iterators::AllCells::iterator::value_type v,
                  complex.iterators().allCells())
    {
        BOOST_CHECK_EQUAL(v.getColor(), (typename Complex::Color)2);
    }

    BOOST_CHECK(complex.iterators(1).allCells().begin() == complex.iterators(1).allCells().end());

    Complex* coAKQ = algorithm.getStrategy()->getOutputComplex();

    BOOST_CHECK_EQUAL(coAKQ->size(), correctSize);

    CRef<ReducibleFreeChainComplexType> RFCComplexCR =
        (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(*coAKQ))();
    CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);

    std::string betti = homSignCR().bettiVector();

    BOOST_CHECK_EQUAL(betti, correctBettiVector);
}

BOOST_AUTO_TEST_CASE(line)
{
    typedef SComplex<SComplexDefaultTraits> Complex;

    Complex::Dims dims = list_of(0)(1)(0)(1)(0);
    Complex::KappaMap kappaMap = tuple_list_of(1, 0, 1)(1, 2, 1)(3, 2, 1)(3, 4, 1); // .-.-.

    Complex complex(3, dims, kappaMap, 1);

    performTest(complex, "1", 1);
}

BOOST_AUTO_TEST_CASE(emptyTriangle)
{
    typedef SComplex<SComplexDefaultTraits> Complex;

    Complex::Dims dims = list_of(0)(1)(0)(1)(0)(1);
    Complex::KappaMap kappaMap = tuple_list_of(1, 0, 1)(1, 2, -1)(3, 2, 1)(3, 4, -1)(5, 4, 1)(5, 0, -1);

    Complex complex(3, dims, kappaMap, 1);
    performTest(complex, "1,1", 2);
}


BOOST_AUTO_TEST_CASE(fullTriangle)
{
    typedef SComplex<SComplexDefaultTraits> Complex;

    Complex::Dims dims = list_of(0)(1)(0)(1)(0)(1)(2);
    Complex::KappaMap kappaMap = tuple_list_of(1, 0, -1)(1, 2, 1)(3, 2, -1)(3, 4, 1)(5, 4, -1)(5, 0, 1)
                                 (6, 1, -1)(6, 3, 1)(6, 5, -1);

    Complex complex(3, dims, kappaMap, 1);

    performTest(complex, "1", 1);
}


BOOST_AUTO_TEST_CASE(torus)
{
    typedef SComplex<SComplexDefaultTraits> Complex;

    SComplexBuilderFromSimplices<long, SComplexDefaultTraits> builder(3);
    std::vector<std::set<int> > simplices = makeSpaceFromWelds(makeTorusWelds());

    for (int i = 0; i < 3; i++)
    {
        simplices = subdivide3(simplices);
    }

    boost::shared_ptr<Complex> complex = builder(simplices, 3, 1);

    CRef<ReducibleFreeChainComplexType> RFCComplexCR_orginal=
        (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(*complex))();
    CRef<HomologySignature> homSignCR_orginal=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR_orginal);

    CoreductionAlgorithm<AKQReduceStrategy<Complex> > algorithm(new AKQReduceStrategy<Complex>(*complex));

    algorithm();


  Complex* AKQ = algorithm.getStrategy()->getOutputComplex();
  BOOST_CHECK_EQUAL(AKQ->size(), 4);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(*AKQ))();
  CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);

  std::ostringstream signature;
  signature << "AKQ: " << homSignCR()<< " | org: " << homSignCR_orginal();
  std::string sig = signature.str();
  std::replace(sig.begin(), sig.end(), '\n', '#');
  BOOST_CHECK_EQUAL(sig, "AKQ:   H^0 = Z^1#  H^1 = Z^2#  H^2 = Z^1# | org:   H^0 = Z^1#  H^1 = Z^2#  H^2 = Z^1#");
}


BOOST_AUTO_TEST_CASE(cubical_projpln) {
  typedef CubSComplex<5> Complex;

  CubSComplexReader<5> reader;
  boost::shared_ptr<Complex> complex = reader(PROJECT_SOURCE_DIR"data/cubical/qprojpln.cub"); 
  
  CRef<ReducibleFreeChainComplexType> RFCComplexCR_orginal=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(*complex))();
  CRef<HomologySignature> homSignCR_orginal=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR_orginal);

  CoreductionAlgorithm<AKQReduceStrategy<Complex> > algorithm(new AKQReduceStrategy<Complex>(*complex));

  algorithm();


  SComplex<SComplexDefaultTraits>* AKQ = algorithm.getStrategy()->getOutputComplex();
  BOOST_CHECK(AKQ != NULL);
  BOOST_CHECK_EQUAL(AKQ->size(), 3);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<SComplex<SComplexDefaultTraits>,
	  ReducibleFreeChainComplexType>(*AKQ))();
  CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);

  std::ostringstream signature;
  signature << "AKQ: " << homSignCR()<< " | org: " << homSignCR_orginal();
  std::string sig = signature.str();
  std::replace(sig.begin(), sig.end(), '\n', '#');
  BOOST_CHECK_EQUAL(sig, "AKQ:   H^0 = Z^1#  H^1 = Z/2# | org:   H^0 = Z^1#  H^1 = Z/2#");
}


/*
BOOST_AUTO_TEST_CASE(klein) {
  typedef SComplex<SComplexDefaultTraits> Complex;

  SComplexBuilderFromSimplices<long, SComplexDefaultTraits> builder(3);
  std::vector<std::set<int> > simplices = makeSpaceFromWelds(makeKleinWelds());

  for (int i = 0; i < 3; i++) {
    //simplices = subdivide3(simplices);
  }

  boost::shared_ptr<Complex> complex = builder(simplices, 3, 1);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR_orginal=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(*complex))();
  CRef<HomologySignature> homSignCR_orginal=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR_orginal);

  CoreductionAlgorithm<AKQReduceStrategy<Complex> > algorithm(new AKQReduceStrategy<Complex, Complex>(*complex));

  algorithm();

  Complex& coAKQ = algorithm.getStrategy().getOutputComplex();
  BOOST_CHECK_EQUAL(coAKQ.size(), 4);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(coAKQ))();
  CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);

  std::ostringstream signature;
  signature << "coAKQ: " << homSignCR()<< " | org: " << homSignCR_orginal();
  std::string sig = signature.str();
  std::replace(sig.begin(), sig.end(), '\n', '#');
  BOOST_CHECK_EQUAL(sig, "coAKQ:   H^0 = Z^1#  H^1 = Z^1 + Z/2# | org:   H^0 = Z^1#  H^1 = Z^1 + Z/2#");
}


BOOST_AUTO_TEST_CASE(projective) {
  typedef SComplex<SComplexDefaultTraits> Complex;

  fstream file(PROJECT_SOURCE_DIR"/data/spaces/projective.dat");

  std::vector<std::set<int> > simplices = parseDat(file);

  SComplexBuilderFromSimplices<long, SComplexDefaultTraits> builder(3);
  boost::shared_ptr<Complex> complex = builder(simplices, 3, 1);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR_orginal=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(*complex))();
  CRef<HomologySignature> homSignCR_orginal=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR_orginal);

  CoreductionAlgorithm<AKQReduceStrategy<Complex> > algorithm(new AKQReduceStrategy<Complex, Complex>(*complex));

  algorithm();

  BOOST_FOREACH(Complex::Iterators::AllCells::iterator::value_type v,
					 complex->iterators().allCells()) {
	 BOOST_CHECK_EQUAL(v.getColor(), (Complex::Color)2);
  }

  BOOST_CHECK(complex->iterators(1).allCells().begin() == complex->iterators(1).allCells().end());

  Complex& coAKQ = algorithm.getStrategy().getOutputComplex();
  BOOST_CHECK_EQUAL(coAKQ.size(), 3);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(coAKQ))();
  CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);

  std::ostringstream signature;
  signature << "coAKQ: " << homSignCR()<< " | org: " << homSignCR_orginal();
  std::string sig = signature.str();
  std::replace(sig.begin(), sig.end(), '\n', '#');
  BOOST_CHECK_EQUAL(sig, "coAKQ:   H^0 = Z^1#  H^1 = Z/2# | org:   H^0 = Z^1#  H^1 = Z/2#");
}


BOOST_AUTO_TEST_CASE(bing) {
  typedef SComplex<SComplexDefaultTraits> Complex;

  fstream file(PROJECT_SOURCE_DIR"/data/spaces/bing.dat");
  std::vector<std::set<int> > simplices = parseDat(file);

  SComplexBuilderFromSimplices<long, SComplexDefaultTraits> builder(3);
  boost::shared_ptr<Complex> complex = builder(simplices, 3, 1);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR_orginal=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(*complex))();
  CRef<HomologySignature> homSignCR_orginal=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR_orginal);

  CoreductionAlgorithm<AKQReduceStrategy<Complex> > algorithm(new AKQReduceStrategy<Complex, Complex>(*complex));

  algorithm();

  BOOST_FOREACH(Complex::Iterators::AllCells::iterator::value_type v,
					 complex->iterators().allCells()) {
	 BOOST_CHECK_EQUAL(v.getColor(), (Complex::Color)2);
  }

  BOOST_CHECK(complex->iterators(1).allCells().begin() == complex->iterators(1).allCells().end());

  Complex& coAKQ = algorithm.getStrategy().getOutputComplex();
  BOOST_CHECK_EQUAL(coAKQ.size(), 1);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(coAKQ))();
  CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);

  std::ostringstream signature;
  signature << "coAKQ: " << homSignCR()<< " | org: " << homSignCR_orginal();
  std::string sig = signature.str();
  std::replace(sig.begin(), sig.end(), '\n', '#');
  BOOST_CHECK_EQUAL(sig, "coAKQ:   H^0 = Z^1# | org:   H^0 = Z^1#");
}


BOOST_AUTO_TEST_CASE(dunce_hat) {
  typedef SComplex<SComplexDefaultTraits> Complex;

  fstream file(PROJECT_SOURCE_DIR"/data/spaces/dunce_hat.dat");
  std::vector<std::set<int> > simplices = parseDat(file);

  SComplexBuilderFromSimplices<long, SComplexDefaultTraits> builder(3);
  boost::shared_ptr<Complex> complex = builder(simplices, 3, 1);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR_orginal=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(*complex))();
  CRef<HomologySignature> homSignCR_orginal=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR_orginal);

  CoreductionAlgorithm<AKQReduceStrategy<Complex> > algorithm(new AKQReduceStrategy<Complex, Complex>(*complex));

  algorithm();

  BOOST_FOREACH(Complex::Iterators::AllCells::iterator::value_type v,
					 complex->iterators().allCells()) {
	 BOOST_CHECK_EQUAL(v.getColor(), (Complex::Color)2);
  }

  BOOST_CHECK(complex->iterators(1).allCells().begin() == complex->iterators(1).allCells().end());

  Complex& coAKQ = algorithm.getStrategy().getOutputComplex();
  BOOST_CHECK_EQUAL(coAKQ.size(), 3);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(coAKQ))();
  CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);

  std::ostringstream signature;
  signature << "coAKQ: " << homSignCR()<< " | org: " << homSignCR_orginal();
  std::string sig = signature.str();
  std::replace(sig.begin(), sig.end(), '\n', '#');
  BOOST_CHECK_EQUAL(sig, "coAKQ:   H^0 = Z^1# | org:   H^0 = Z^1#");
}


BOOST_AUTO_TEST_CASE(nc_sphere) {
  typedef SComplex<SComplexDefaultTraits> Complex;

  fstream file(PROJECT_SOURCE_DIR"/data/spaces/nc_sphere.dat");
  std::vector<std::set<int> > simplices = parseDat(file);

  SComplexBuilderFromSimplices<long, SComplexDefaultTraits> builder(3);
  boost::shared_ptr<Complex> complex = builder(simplices, 3, 1);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR_orginal=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(*complex))();
  CRef<HomologySignature> homSignCR_orginal=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR_orginal);

  CoreductionAlgorithm<AKQReduceStrategy<Complex> > algorithm(new AKQReduceStrategy<Complex, Complex>(*complex));

  algorithm();

  BOOST_FOREACH(Complex::Iterators::AllCells::iterator::value_type v,
					 complex->iterators().allCells()) {
	 BOOST_CHECK_EQUAL(v.getColor(), (Complex::Color)2);
  }

  BOOST_CHECK(complex->iterators(1).allCells().begin() == complex->iterators(1).allCells().end());

  Complex& coAKQ = algorithm.getStrategy().getOutputComplex();
  BOOST_CHECK_EQUAL(coAKQ.size(), 4);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(coAKQ))();
  CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);

  std::ostringstream signature;
  signature << "coAKQ: " << homSignCR()<< " | org: " << homSignCR_orginal();
  std::string sig = signature.str();
  std::replace(sig.begin(), sig.end(), '\n', '#');
  BOOST_CHECK_EQUAL(sig, "coAKQ:   H^0 = Z^1#  H^1 = 0#  H^2 = 0#  H^3 = Z^1# | org:   H^0 = Z^1#  H^1 = 0#  H^2 = 0#  H^3 = Z^1#");
}



BOOST_AUTO_TEST_CASE(poincare) {
  typedef SComplex<SComplexDefaultTraits> Complex;

  fstream file(PROJECT_SOURCE_DIR"/data/spaces/poincare.dat");
  std::vector<std::set<int> > simplices = parseDat(file);

  SComplexBuilderFromSimplices<long, SComplexDefaultTraits> builder(3);
  boost::shared_ptr<Complex> complex = builder(simplices, 3, 1);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR_orginal=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(*complex))();
  CRef<HomologySignature> homSignCR_orginal=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR_orginal);

  CoreductionAlgorithm<AKQReduceStrategy<Complex> > algorithm(new AKQReduceStrategy<Complex, Complex>(*complex));

  algorithm();

  BOOST_FOREACH(Complex::Iterators::AllCells::iterator::value_type v,
					 complex->iterators().allCells()) {
	 BOOST_CHECK_EQUAL(v.getColor(), (Complex::Color)2);
  }

  BOOST_CHECK(complex->iterators(1).allCells().begin() == complex->iterators(1).allCells().end());

  Complex& coAKQ = algorithm.getStrategy().getOutputComplex();
  BOOST_CHECK_EQUAL(coAKQ.size(), 6);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(coAKQ))();
  CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);

  std::ostringstream signature;
  signature << "coAKQ: " << homSignCR()<< " | org: " << homSignCR_orginal();
  std::string sig = signature.str();
  std::replace(sig.begin(), sig.end(), '\n', '#');
  BOOST_CHECK_EQUAL(sig, "coAKQ:   H^0 = Z^1#  H^1 = 0#  H^2 = 0#  H^3 = Z^1# | org:   H^0 = Z^1#  H^1 = 0#  H^2 = 0#  H^3 = Z^1#");
}



BOOST_AUTO_TEST_CASE(nonpl_sphere) {
  typedef SComplex<SComplexDefaultTraits> Complex;

  fstream file(PROJECT_SOURCE_DIR"/data/spaces/nonpl_sphere.dat");
  std::vector<std::set<int> > simplices = parseDat(file);

  SComplexBuilderFromSimplices<long, SComplexDefaultTraits> builder(3);
  boost::shared_ptr<Complex> complex = builder(simplices, 3, 1);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR_orginal=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(*complex))();
  CRef<HomologySignature> homSignCR_orginal=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR_orginal);

  CoreductionAlgorithm<AKQReduceStrategy<Complex> > algorithm(new AKQReduceStrategy<Complex, Complex>(*complex));

  algorithm();

  BOOST_FOREACH(Complex::Iterators::AllCells::iterator::value_type v,
					 complex->iterators().allCells()) {
	 BOOST_CHECK_EQUAL(v.getColor(), (Complex::Color)2);
  }

  BOOST_CHECK(complex->iterators(1).allCells().begin() == complex->iterators(1).allCells().end());

  Complex& coAKQ = algorithm.getStrategy().getOutputComplex();
  BOOST_CHECK_EQUAL(coAKQ.size(), 6);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(coAKQ))();
  CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);

  std::ostringstream signature;
  signature << "coAKQ: " << homSignCR()<< " | org: " << homSignCR_orginal();
  std::string sig = signature.str();
  std::replace(sig.begin(), sig.end(), '\n', '#');
  BOOST_CHECK_EQUAL(sig, "coAKQ:   H^0 = Z^1#  H^1 = 0#  H^2 = 0#  H^3 = 0#  H^4 = 0#  H^5 = Z^1# | org:   H^0 = Z^1#  H^1 = 0#  H^2 = 0#  H^3 = 0#  H^4 = 0#  H^5 = Z^1#");
}


BOOST_AUTO_TEST_CASE(knot) {
  typedef SComplex<SComplexDefaultTraits> Complex;

  fstream file(PROJECT_SOURCE_DIR"/data/spaces/knot.dat");
  std::vector<std::set<int> > simplices = parseDat(file);

  SComplexBuilderFromSimplices<long, SComplexDefaultTraits> builder(3);
  boost::shared_ptr<Complex> complex = builder(simplices, 3, 1);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR_orginal=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(*complex))();
  CRef<HomologySignature> homSignCR_orginal=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR_orginal);

  CoreductionAlgorithm<AKQReduceStrategy<Complex> > algorithm(new AKQReduceStrategy<Complex, Complex>(*complex));

  algorithm();

  BOOST_FOREACH(Complex::Iterators::AllCells::iterator::value_type v,
					 complex->iterators().allCells()) {
	 BOOST_CHECK_EQUAL(v.getColor(), (Complex::Color)2);
  }

  BOOST_CHECK(complex->iterators(1).allCells().begin() == complex->iterators(1).allCells().end());

  Complex& coAKQ = algorithm.getStrategy().getOutputComplex();
  BOOST_CHECK_EQUAL(coAKQ.size(), 3);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(coAKQ))();
  CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);

  std::ostringstream signature;
  signature << "coAKQ: " << homSignCR()<< " | org: " << homSignCR_orginal();
  std::string sig = signature.str();
  std::replace(sig.begin(), sig.end(), '\n', '#');
  BOOST_CHECK_EQUAL(sig, "coAKQ:   H^0 = Z^1# | org:   H^0 = Z^1#");
}


void checkFile(const std::string& fileName) {

}


BOOST_AUTO_TEST_CASE(bjorner) {
  typedef SComplex<SComplexDefaultTraits> Complex;

  fstream file(PROJECT_SOURCE_DIR"/data/spaces/bjorner.dat");
  std::vector<std::set<int> > simplices = parseDat(file);

  SComplexBuilderFromSimplices<long, SComplexDefaultTraits> builder(3);
  boost::shared_ptr<Complex> complex = builder(simplices, 3, 1);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR_orginal=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(*complex))();
  CRef<HomologySignature> homSignCR_orginal=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR_orginal);

  CoreductionAlgorithm<AKQReduceStrategy<Complex> > algorithm(new AKQReduceStrategy<Complex, Complex>(*complex));

  algorithm();

  BOOST_FOREACH(Complex::Iterators::AllCells::iterator::value_type v,
					 complex->iterators().allCells()) {
	 BOOST_CHECK_EQUAL(v.getColor(), (Complex::Color)2);
  }

  BOOST_CHECK(complex->iterators(1).allCells().begin() == complex->iterators(1).allCells().end());

  Complex& coAKQ = algorithm.getStrategy().getOutputComplex();
  BOOST_CHECK_EQUAL(coAKQ.size(), 4);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(coAKQ))();
  CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);

  std::ostringstream signature;
  signature << "coAKQ: " << homSignCR()<< " | org: " << homSignCR_orginal();
  std::string sig = signature.str();
  std::replace(sig.begin(), sig.end(), '\n', '#');
  BOOST_CHECK_EQUAL(sig, "coAKQ:   H^0 = Z^1#  H^1 = 0#  H^2 = Z^1# | org:   H^0 = Z^1#  H^1 = 0#  H^2 = Z^1#");
}


// BOOST_AUTO_TEST_CASE(randomSimplices) {
//   typedef SComplex<SComplexDefaultTraits> Complex;

//   fstream file(PROJECT_SOURCE_DIR"/data/randomSimplices.dat");
//   std::vector<std::set<int> > simplices = parseDat(file);

//   SComplexBuilderFromSimplices<long, SComplexDefaultTraits> builder(3);
//   boost::shared_ptr<Complex> complex = builder(simplices, 3, 1);

//   CRef<ReducibleFreeChainComplexType> RFCComplexCR_orginal=
//   	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(*complex))();
//   CRef<HomologySignature> homSignCR_orginal=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR_orginal);

//   CoreductionAlgorithm<AKQReduceStrategy<Complex> > algorithm(new AKQReduceStrategy<Complex, Complex>(*complex));

//   algorithm();

//   BOOST_FOREACH(Complex::Iterators::AllCells::iterator::value_type v,
// 					 complex->iterators().allCells()) {
// 	 BOOST_CHECK_EQUAL(v.getColor(), (Complex::Color)2);
//   }

//   BOOST_CHECK(complex->iterators(1).allCells().begin() == complex->iterators(1).allCells().end());

//   Complex& coAKQ = algorithm.getStrategy().getOutputComplex();
//   BOOST_CHECK_EQUAL(coAKQ.size(), 174);

//   CRef<ReducibleFreeChainComplexType> RFCComplexCR=
//   	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(coAKQ))();
//   CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);

//   std::ostringstream signature;
//   signature << "coAKQ: " << homSignCR()<< " | org: " << homSignCR_orginal();
//   std::string sig = signature.str();
//   std::replace(sig.begin(), sig.end(), '\n', '#');
//   BOOST_CHECK_EQUAL(sig, "coAKQ:   H^0 = Z^1#  H^1 = Z^39#  H^2 = Z^84# | org:   H^0 = Z^1#  H^1 = Z^39#  H^2 = Z^84#");
// }

*/


BOOST_AUTO_TEST_SUITE_END()

