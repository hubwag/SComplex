#include <iterator>
#include <iostream>
#include <vector>
#include <set>
#include <cmath>

using namespace std;

/*
0 1 2 3
4 5 6 7
8 9 A B
C D E F

torus:
0 1 2 0
4 5 6 4
8 9 A 8
0 1 2 0

klein bottle:
0 1 2 0
4 5 6 4
8 9 A 8
0 2 1 0

projective space:
0 1 2 3
4 5 6 7
7 9 A 4
3 2 1 0
*/

set<int> make_int_set(int a = -1, int b = -1, int c = -1, int d = -1)
{
	set<int> s;
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

void print_set(const set<int> &s)
{
	// copy(s.begin(), s.end(), ostream_iterator<int>(cout, ", "));
	for (set<int>::const_iterator it = s.begin(), end = s.end(); it != end; ++it)
	{
		cout << hex << *it << ",";
	}
	cout << dec << endl;
}

int find_set(const vector<int> &v, int id)
{
	while(v[id] != id)
		id = v[id];
	return id;
}

vector<int> makeTorusWelds()
{
	const int n = 16;
	vector<int> x(n);

	for (int i = 0; i < n; i++)
		x[i] = i;

	for (int i = 0; i < 4; i++)
		x[0xC + i] = i;

	for (int i = 3; i < 16; i+=4)
		x[max(i, i-3)] = min(i, i-3);

	return x;
}

vector<int> makeKleinWelds()
{
	const int n = 16;
	vector<int> x(n);

	for (int i = 0; i < n; i++)
		x[i] = i;

	for (int i = 0; i < 4; i++)
		x[max(i, 0xF - i)] = min(i, 0xF-i);

	for (int i = 3; i < 16; i+=4)
		x[max(i, i - 3)] = min(i, i - 3);

	return x;
}

vector<int> makeProjectiveSpaceWelds()
{
	const int n = 16;
	vector<int> x(n);

	for (int i = 0; i < n; i++)
		x[i] = i;

	for (int i = 0; i < 4; i++)
		x[0xF - i] = i;

	for (int i = 3; i < 16; i+=4)
		x[max(i, 0xF-i)] = min(i, 0xF - i);

	return x;
}

vector<set<int> > makeSpaceFromWelds(const vector<int> &welds)
{
	const vector<int> &x = welds;
	const int n = welds.size();
	const int side = static_cast<int>(round(sqrt(n)));

	vector<set<int> > tris;

	for (int i = 0; i+1 < side; i++)
	for (int j = 0; j+1 < side; j++)
	{
		int a = side * i + j; // origin
		int b = side * i + j + 1; // horizontal
		int c = side * i + j + side; // vertical
		int d = side * i + j + side + 1; // diagonal

		tris.push_back(make_int_set(find_set(x,a), find_set(x,b), find_set(x,d)));
		tris.push_back(make_int_set(find_set(x,a), find_set(x,c), find_set(x,d)));
	}

	return tris;
}
// barycentric
vector<set<int> > subdivide6(const vector<set<int> > &v)
{
	assert(v[0].size() == 3); // Dim == 2 only, for now!
	typedef vector<set<int> > simp_vec;

	int mx = 0;
	for (simp_vec::const_iterator it = v.begin(); it != v.end(); ++it)
	{
		mx = max(mx, *max_element(it->begin(), it->end()));
	}

	int next_free_triangle_middle = mx + 1;
	int next_free_edge_middle = 3 * next_free_triangle_middle + 2; // ?

	map<set<int>, int> middle;

	vector<set<int> > ret;

	for (simp_vec::const_iterator it = v.begin(); it != v.end(); ++it)
	{
		vector<int> cycle(it->begin(), it->end());
		cycle.push_back(cycle.front());
		for (size_t i = 0; i + 1 < cycle.size(); ++i)
		{
			int a = cycle[i];
			int b = cycle[i+1];

			int m = 0;
			if (middle.count(make_int_set(a,b)))
				m = middle[make_int_set(a,b)];
			else m = middle[make_int_set(a,b)] = next_free_edge_middle++;

			ret.push_back(make_int_set(a,m,next_free_triangle_middle));
			ret.push_back(make_int_set(b,m,next_free_triangle_middle));
		}
		++next_free_triangle_middle;
	}

	assert(ret.size() == 6 * v.size());

	return ret;
}

vector<set<int> > subdivide3(const vector<set<int> > &v)
{
	typedef vector<set<int> > simp_vec;
	int next = 0;
	for (simp_vec::const_iterator it = v.begin(); it != v.end(); ++it)
	{
		next = max(next, *max_element(it->begin(), it->end()));
	}

	++next;

	vector<set<int> > ret;

	for (simp_vec::const_iterator it = v.begin(); it != v.end(); ++it)
	{
		for (set<int>::const_iterator curr = it->begin(); curr != it->end(); ++curr)
		{
			ret.push_back(*it);
			ret.back().erase(*curr);
			ret.back().insert(next);
		}

		++next;
	}

	return ret;
}

int ___main()
{
	vector<set<int> > comp = makeSpaceFromWelds(makeKleinWelds());

	cout << comp.size() << endl;

	for (int i = 0; i < 3; i++)
	{
		comp = subdivide3(comp);
		cout << (dec) << comp.size() << endl;
	}

	return 0;
}
