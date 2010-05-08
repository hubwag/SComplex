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
#include <string>
#include <sstream>
#include<fstream>

using namespace boost;
using namespace boost::assign;
using namespace std;

typedef ElementaryCell ElementaryCellType;
typedef int ScalarType;
typedef FreeModule<int,capd::vectalg::Matrix<int,0,0> > FreeModuleType;
typedef FreeChainComplex<FreeModuleType> FreeChainComplexType;
typedef ReducibleFreeChainComplex<FreeModuleType,int> ReducibleFreeChainComplexType;



vector<set<int> > parseDat(istream &stream)
{
  
  vector<set<int> > simplices;

  for (string s; std::getline(stream, s); )
    {
      if (s.size() == 0 || s[0] == '#')
	continue;

      stringstream ss(s);

      set<int> simp;

      for (int v; ss >> v;)
	simp.insert(v);

      simplices.push_back(simp);
    }

  return simplices;
}

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


BOOST_AUTO_TEST_CASE(emptyTriangle) {
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

BOOST_AUTO_TEST_CASE(fullTriangle) {
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
  BOOST_CHECK_EQUAL(coAKQ.size(), 1);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(coAKQ))();
  CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);

  std::string betti = homSignCR().bettiVector();

  BOOST_CHECK_EQUAL(betti, "1");

}


BOOST_AUTO_TEST_CASE(torus) {
  typedef SComplex<SComplexDefaultTraits> Complex;

  SComplexBuilderFromSimplices<long, SComplexDefaultTraits> builder(3);
  std::vector<std::set<int> > simplices = makeSpaceFromWelds(makeTorusWelds());

  for (int i = 0; i < 3; i++) {
  	 simplices = subdivide3(simplices);
  }

  boost::shared_ptr<Complex> complex = builder(simplices, 3, 1);

  CRef<ReducibleFreeChainComplexType> RFCComplexCR_orginal=
  	 (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(*complex))();
  CRef<HomologySignature> homSignCR_orginal=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR_orginal);

  COAKQAlgorithm<COAKQStrategy<Complex, Complex> > algorithm(new COAKQStrategy<Complex, Complex>(*complex));

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
  BOOST_CHECK_EQUAL(sig, "coAKQ:   H^0 = Z^1#  H^1 = Z^2#  H^2 = Z^1# | org:   H^0 = Z^1#  H^1 = Z^2#  H^2 = Z^1#");
}

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

  COAKQAlgorithm<COAKQStrategy<Complex, Complex> > algorithm(new COAKQStrategy<Complex, Complex>(*complex));

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

  COAKQAlgorithm<COAKQStrategy<Complex, Complex> > algorithm(new COAKQStrategy<Complex, Complex>(*complex));

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

BOOST_AUTO_TEST_SUITE_END()

