#include <iostream>
using namespace std;

#include <boost/test/unit_test.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <algorithm>
#include <set>
#include <queue>

#include <capd/auxil/Stopwatch.h>
#include <capd/auxil/CRef.h>
#include <capd/homologicalAlgebra/embeddingDim.h>

#include <capd/vectalg/MatrixSlice.h>
#include <capd/matrixAlgorithms/intMatrixAlgorithms.hpp>

#include <capd/homologicalAlgebra/homologicalAlgebra.hpp>
#include <capd/homologicalAlgebra/homAlgFunctors.hpp>
#include <capd/homologicalAlgebra/cubSetFunctors.hpp>
#include <capd/homologicalAlgebra/ReducibleFreeChainComplex.hpp>

//#include "CubSComplex.hpp"
#include "SComplexAlgs.hpp"

typedef ElementaryCell ElementaryCellType;
typedef int ScalarType;
typedef FreeModule<int,capd::vectalg::Matrix<int,0,0> > FreeModuleType;
typedef FreeChainComplex<FreeModuleType> FreeChainComplexType;
typedef ReducibleFreeChainComplex<FreeModuleType,int> ReducibleFreeChainComplexType;

typedef ElementaryCell ElementaryCellType;
typedef int ScalarType;
typedef FreeModule<int,capd::vectalg::Matrix<int,0,0> > FreeModuleType;
typedef FreeChainComplex<FreeModuleType> FreeChainComplexType;
typedef ReducibleFreeChainComplex<FreeModuleType,int> ReducibleFreeChainComplexType;

#include "Simplex.hpp"
#include "SimplexSComplex.hpp"

#include "SComplexAlgs.hpp"

BOOST_AUTO_TEST_SUITE(SimplexSComplex_basic)

std::set<int> make_int_set(int a, int b = -1, int c = -1, int d = -1)
{
	std::set<int> s;

	if (a >= 0)
		s.insert(a);
	if (b >= 0)
		s.insert(b);
	if (c >= 0)
		s.insert(c);
	if (d >= 0)
		s.insert(d);

	return s;
}

BOOST_AUTO_TEST_CASE(add_get_simple_test_2d) {
	SimplexSComplex comp;

	comp.add_simplex(make_int_set(0,1));
	BOOST_CHECK(comp.get_simplex(make_int_set(0,1)) != 0);
	BOOST_CHECK(comp.get_simplex(make_int_set(0)) != 0);
	BOOST_CHECK(comp.get_simplex(make_int_set(1)) != 0);
}

BOOST_AUTO_TEST_CASE(add_get_simple_test_3d) {
	SimplexSComplex comp;

	comp.add_simplex(make_int_set(0,1,2));
	BOOST_CHECK(comp.get_simplex(make_int_set(0,1,2)) != 0);
	BOOST_CHECK(comp.get_simplex(make_int_set(0,1)) != 0);
	BOOST_CHECK(comp.get_simplex(make_int_set(1,2)) != 0);
	BOOST_CHECK(comp.get_simplex(make_int_set(0,2)) != 0);
	BOOST_CHECK(comp.get_simplex(make_int_set(0)) != 0);
	BOOST_CHECK(comp.get_simplex(make_int_set(1)) != 0);
	BOOST_CHECK(comp.get_simplex(make_int_set(2)) != 0);
}

BOOST_AUTO_TEST_CASE(shave_and_coreduce_test) {
	SimplexSComplex comp;

	// An example from: COREDUCTION HOMOLOGY ALGORITHM
	// MARIAN MROZEK AND BOGDAN BATKO

	comp.add_simplex(make_int_set(1,2,3));
	comp.add_simplex(make_int_set(0,1));
	comp.add_simplex(make_int_set(0,3));
	comp.add_simplex(make_int_set(0,4));
	comp.add_simplex(make_int_set(3,4));
	comp.add_simplex(make_int_set(0,5));

	(ShaveAlgorithmFactory::createDefault(comp))();
	BOOST_CHECK_EQUAL(9, comp.cardinality());

	(CoreductionAlgorithmFactory::createDefault(comp))();
	BOOST_CHECK_EQUAL(2, comp.cardinality());
}

BOOST_AUTO_TEST_SUITE_END()
