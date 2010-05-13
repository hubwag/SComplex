#ifndef SIMPLEX_IO_HPP
#define SIMPLEX_IO_HPP

#include <string>
#include <sstream>
#include <set>

class SimplexSComplex;

void parseDat(istream &stream, SimplexSComplex &comp, int subdivs = 0)
{
	using std::string;
	using std::stringstream;

	vector<set<int> > simplices;

	for (string s; getline(stream, s); )
	{
		if (s.size() == 0 || s[0] == '#')
			continue;

		stringstream ss(s);

		set<int> simp;

		for (int v; ss >> v;)
			simp.insert(v);

		simplices.push_back(simp);

		// comp.addSimplex(simp);
	}

	for (int i = 0; i < subdivs; i++)
		simplices = subdivide3(simplices);

	for (int i = 0; i < simplices.size(); i++)
		comp.addSimplex(simplices[i]);
}

void parseObj(istream &stream, SimplexSComplex &comp, int subdivs = 0)
{
	using std::string;
	using std::stringstream;

	vector<set<int> > simplices;

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

		simplices.push_back(simp);

		// comp.addSimplex(simp);
	}

	for (int i = 0; i < subdivs; i++)
		simplices = subdivide3(simplices);

	for (int i = 0; i < simplices.size(); i++)
		comp.addSimplex(simplices[i]);
}

#endif
