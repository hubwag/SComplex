#ifndef SIMPLESET_H
#define SIMPLESET_H

#include <vector>
#include <boost/iterator/filter_iterator.hpp>

using namespace std;

template<typename t>
class simple_set
{
    vector<t> elems;
    vector<int> present;

public:

    template<typename iter_t>
    simple_set(iter_t _begin, iter_t _end) : elems(_begin, _end)
    {
        present.resize(elems.size(), true);
    }

    long getInt() const
	{
		long r = 0;
		for (vector<int>::const_iterator it = present.begin(), end = present.end(); it != end; ++it)
		{
			r<<=1;
			r|=*it;
		}

		return r;
	}

    simple_set()
    {
    }

    size_t size() const
    {
        return distance(begin(), end()); // :)
    }

    void insert(const t &val)
    {
        elems.push_back(val);
        present.push_back(true);
    }

    struct is_present
    {
		const simple_set &ss;

		is_present(const simple_set &s) : ss(s) {}

    	bool operator()(const t &x) const
    	{
    		if (ss.elems.size() == 0)
				return false;
    		int off = &x - &ss.elems[0];
    		return ss.present[off] == true;
    	}
    };

    typedef boost::filter_iterator<is_present, typename vector<t>::iterator> iterator;
    typedef boost::filter_iterator<is_present, typename vector<t>::const_iterator> const_iterator;

    void insert_index(int p, const t&)
    {
    	present[p] = true;
    }

	void erase_index(int p, const t&)
    {
    	present[p] = false;
    }

    iterator begin()
    {
        return iterator(is_present(*this), elems.begin(), elems.end());
    }

    iterator end()
    {
        return iterator(is_present(*this), elems.end(), elems.end());
    }

    const_iterator begin() const
    {
        return const_iterator(is_present(*this), elems.begin(), elems.end());
    }

    const_iterator end() const
    {
        return const_iterator(is_present(*this), elems.end(), elems.end());
    }
};



#endif
