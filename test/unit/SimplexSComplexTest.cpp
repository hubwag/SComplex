#include <iostream>
using namespace std;

#include <boost/test/unit_test.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <algorithm>
#include <set>

#include "SimplexSComplex.hpp"
#include "Simplex.hpp"

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

BOOST_AUTO_TEST_CASE(iterators_AllCells_init_test) {
	SimplexSComplex comp;

	comp.add_simplex(make_int_set(0,1));
	BOOST_CHECK(comp.get_simplex(make_int_set(0,1)) != 0);
	BOOST_CHECK(comp.get_simplex(make_int_set(0)) != 0);
	BOOST_CHECK(comp.get_simplex(make_int_set(1)) != 0);
}


BOOST_AUTO_TEST_SUITE_END()
