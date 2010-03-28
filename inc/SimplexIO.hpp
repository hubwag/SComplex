#ifndef SIMPLEX_IO_HPP
#define SIMPLEX_IO_HPP

#include <string>
#include <sstream>
#include <set>

class SimplexSComplex;

void parseObj(istream &stream, SimplexSComplex &comp)
{
	using std::string;
	using std::stringstream;

	for (string s; getline(stream, s); )
	{
		stringstream ss(s);

		string what;
		ss >> what;

		if (what != "f")
			continue;

		set<int> simp;

		for (int v; ss >> v;)
			simp.insert(v);

		comp.addSimplex(simp);
	}
}

#endif
