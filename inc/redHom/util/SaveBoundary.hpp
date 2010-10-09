#ifndef SAVE_BOUNDARY_HPP
#define SAVE_BOUNDARY_HPP

#include <sstream>
#include <map>
#include <algorithm>
#include <fstream>

template<typename SComplex>
void saveBoundary(SComplex &complex, int d, const std::string &name_prefix)
{
    std::stringstream ss(name_prefix);
    ss << name_prefix << ".b" << d;
    std::ofstream out(ss.str().c_str());

	// some space for the header, to avoid rewriting the whole file!
    out << "                                     \n";

    int mxa = 0, mxb = 0;

    std::map<int, int> firstMap;
    std::map<int, int> secondMap;

    for (typename SComplex::ColoredIterators::Iterators::DimCells::iterator cell = complex.iterators(1).dimCells(d).begin(),
            endc = complex.iterators(1).dimCells(d).end();
            cell != endc; ++cell)
    {
        std::vector<std::pair<int,int> > v;

        for (typename SComplex::ColoredIterators::Iterators::BdCells::iterator it = complex.iterators(1).bdCells(*cell).begin(),
                end = complex.iterators(1).bdCells(*cell).end();
                it != end; ++it)
        {
        	int c = it->getId();
        	if (firstMap.count(c) == 0)
        	{
        		int next = firstMap.size();
				firstMap[c] = next;
        	}

			int mapped = firstMap[c] + 1;

            v.push_back(make_pair(mapped, complex.coincidenceIndex(*cell, *it)));
        }

        std::sort(v.begin(), v.end());
        int id = cell->getId();

		if (secondMap.count(id) == 0)
		{
			int next = secondMap.size();
			secondMap[id] = next;
		}

		id = secondMap[id] + 1;

        mxa = max(mxa, id);

        for (size_t i = 0; i < v.size(); i++)
        {
            out << id << ' ' << v[i].first << ' ' << v[i].second << ' ' << '\n';
            mxb = max(mxb, v[i].first);
        }
    }

    out << "0 0 0\n"; // for linbox and gap/Homology

    out.seekp(0);

    out << mxa << " " << mxb << " M";

    std::cout << "saved bd of dim " <<  d << std::endl;
    std::cout << "matrix dimension: " << mxa << ' ' << mxb << std::endl;
}

#endif
