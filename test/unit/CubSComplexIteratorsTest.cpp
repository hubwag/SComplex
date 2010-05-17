#include <boost/test/unit_test.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/shared_ptr.hpp>
#include <algorithm>

#include <iostream>

#include <redHom/complex/cubical/CubSComplex.hpp>
#include <redHom/complex/cubical/CubSComplexReader.hpp>

using namespace std;

BOOST_AUTO_TEST_SUITE(CubSComplex_iterators)


boost::shared_ptr<CubSComplex<3> > createComplex() {
  CubSComplexReader<3> reader;
  
  return reader(PROJECT_SOURCE_DIR"data/cubical/torus.cub");
}

  
BOOST_AUTO_TEST_CASE(iterators_AllCells_init_test) {
	 boost::shared_ptr<CubSComplex<3> > s = createComplex();
	 
	 BOOST_CHECK(s->iterators(1).allCells().begin() != s->iterators(1).allCells().end());
	 BOOST_CHECK(s->iterators(1).allCells().begin() == s->iterators(1).allCells().begin());
	 BOOST_CHECK(s->iterators(1).allCells().end() == s->iterators(1).allCells().end());
	 BOOST_CHECK( s->iterators(1).allCells().begin() != s->iterators(1).allCells().end());
	 BOOST_CHECK( s->iterators(1).allCells().begin() == s->iterators(1).allCells().begin());
	 BOOST_CHECK( s->iterators(1).allCells().end() == s->iterators(1).allCells().end());
}

BOOST_AUTO_TEST_CASE(iterators_DimCells_init_test) {
  boost::shared_ptr<CubSComplex<3> > s = createComplex();

	 BOOST_CHECK(s->iterators(1).dimCells(0).begin() != s->iterators(1).dimCells(0).end());
	 BOOST_CHECK(s->iterators(1).dimCells(0).begin() == s->iterators(1).dimCells(0).begin());
	 BOOST_CHECK(s->iterators(1).dimCells(0).end() == s->iterators(1).dimCells(0).end());
	 BOOST_CHECK( s->iterators(1).dimCells(0).begin() != s->iterators(1).dimCells(0).end());
	 BOOST_CHECK( s->iterators(1).dimCells(0).begin() == s->iterators(1).dimCells(0).begin());
	 BOOST_CHECK( s->iterators(1).dimCells(0).end() == s->iterators(1).dimCells(0).end());
}

BOOST_AUTO_TEST_CASE(iterators_AllCells_forEach_test) {
  typedef CubSComplex<3> SComplex;
  boost::shared_ptr<SComplex> s = createComplex();
	 const int initSize = s->cardinality();
	 
	 BOOST_CHECK_EQUAL(s->cardinality(), 3456);

	 int cnt = 0;
	 for (SComplex::Iterators::AllCells::iterator it = s->iterators(1).allCells().begin();
			it != s->iterators(1).allCells().end(); ++it) {
		if ((*it).getColor() == 1) {
		  ++cnt;
		}
	 }
	 BOOST_CHECK_EQUAL(cnt, initSize);
	 
	 std::for_each(s->iterators(1).allCells().begin(), s->iterators(1).allCells().end(), boost::bind(&SComplex::Cell::setColor, _1, 2));
	 BOOST_CHECK_EQUAL(s->cardinality(), 0);

	 cnt = 0;
	 std::for_each(s->iterators(1).allCells().begin(), s->iterators(1).allCells().end(), boost::lambda::var(cnt)++);
	 //	 BOOST_CHECK_EQUAL(cnt, initSize);	 
}

BOOST_AUTO_TEST_SUITE_END()
