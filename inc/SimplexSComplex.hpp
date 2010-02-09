#ifndef SIMPLEXSCOMPLEX_HPP
#define SIMPLEXSCOMPLEX_HPP

#include <vector>
#include <map>
#include <boost/iterator/transform_iterator.hpp>

#include "SimplexCell.hpp"
#include "simple_set.h"

using namespace std;

class Simplex;

class SimplexSComplex
{
public:

    int cardinality()
    {
        return distance(all_begin(), all_end()); // slow
    }

    typedef SimplexCell Cell;

// to be changed later!!
#include "SimplexInternalColoredIterators.hpp"
// to be changed later!!

    ColoredIterators::Iterators iterators(int color)
    {
        return ColoredIterators::Iterators(*this, color);
    }

    template<int color>
    ColoredIterators::Iterators iterators()
    {
        return ColoredIterators::Iterators(*this, color);
    }

private:

    typedef map<vector<int>, Simplex*> map_type;
    typedef vector<Simplex*> per_dimension_storage_type;

    // typedef set<int> vertex_set;
    typedef simple_set<int> vertex_set;

    map<int, per_dimension_storage_type> per_dimension;
    per_dimension_storage_type all_simplices;

    static const int MAX_VERTICES = 100000;
    vector<Simplex*> base;

    void simplex_added_event(Simplex &s)
    {
        per_dimension[s.getDim()].push_back(&s);
        all_simplices.push_back(&s);
    }

    void simplex_removed_event(Simplex &s)
    {
    }

    Simplex* make_base_simplex(vertex_set &s)
    {
        assert(s.size() == 1);

        const int v = *s.begin();

        if (base[v] == 0)
        {
            Simplex *new_simplex = new Simplex(s.begin(), s.end());
            simplex_added_event(*new_simplex);
            base[v] = new_simplex;

            return new_simplex;
        }

        return base[v];
    }

    Simplex* create_simplex_hierarchy(vertex_set &s);

public:

    SimplexSComplex() : base(MAX_VERTICES, static_cast<Simplex*>(0))
    {
    }

    ~SimplexSComplex()
    {
        for (per_dimension_storage_type::const_iterator it = all_simplices.begin(), end = all_simplices.end(); it != end; ++it)
            delete *it;
    }

    typedef vector<Simplex*> memo_type;
    memo_type memo;

    long get_int(const simple_set<int> &s)
    {
        return s.getInt();
    }

    long get_int(const set<int> &s)
    {
        return -1;
    }

    template<typename iterable_t>
    Simplex* get_simplex(const iterable_t &s);

    inline bool getUniqueCoFace(const Cell& cell, Cell& coface) const
    {
    	// could be optimized...
        int cnt = std::distance(cell.getImpl().coborder_begin(1), cell.getImpl().coborder_end(1));

        if (cnt == 1)
        {
            coface = *(cell.getImpl().coborder_begin(1));
            return true;
        }

        return false;
    }

    inline bool getUniqueFace(const Cell& cell, Cell& face) const
    {
        int cnt = std::distance(cell.getImpl().border_begin(1), cell.getImpl().border_end(1));

        if (cnt == 1)
        {
            face = *(cell.getImpl().border_begin(1));
            return true;
        }

        return false;
    }

    typedef Simplex::colored_iterator iterator;

    iterator all_begin(int color = 1)
    {
        return Simplex::colored_iterator(Simplex::colored_iterator_inner(Simplex::has_color(color), all_simplices.begin(), all_simplices.end()) );
    }

    iterator all_end(int color = 1)
    {
        return Simplex::colored_iterator(Simplex::colored_iterator_inner(Simplex::has_color(color), all_simplices.end(), all_simplices.end()) );
    }

    iterator dim_begin(int dim, int color = 1)
    {
        return Simplex::colored_iterator(Simplex::colored_iterator_inner(Simplex::has_color(color), per_dimension[dim].begin(), per_dimension[dim].end()) );
    }

    iterator dim_end(int dim, int color = 1)
    {
        return Simplex::colored_iterator(Simplex::colored_iterator_inner(Simplex::has_color(color), per_dimension[dim].end(), per_dimension[dim].end()) );
    }

private:

    Simplex* add_simplex(vertex_set &s)
    {
        Simplex *found = get_simplex(s);

        if (found)
            return found;

        return create_simplex_hierarchy(s);
    }

public:

    Simplex* add_simplex(const set<int> &s)
    {
        memo.resize(1u<<s.size(), 0);

        vertex_set faster(s.begin(), s.end());
        Simplex *ret = add_simplex(faster);

        memo.clear();

        return ret;
    }
};

Simplex* SimplexSComplex::create_simplex_hierarchy(vertex_set &s)
{
    assert(s.size() > 0);

    if (s.size() == 1)
    {
        return make_base_simplex(s);
    }

    Simplex* root_simplex = new Simplex(s.begin(), s.end());
    simplex_added_event(*root_simplex);

    // set<int> probably can't be used as a drop-in replacement...
    for (vertex_set::iterator it = s.begin(), end = s.end(); it != end; )
    {
        const int &val = *it;
        vertex_set::iterator new_it = it;
        ++new_it;

        s.erase(it);
        Simplex *sub_simplex = add_simplex(s);
        s.insert(it, val);

        root_simplex->add_to_border(*sub_simplex);
        sub_simplex->add_to_coborder(*root_simplex, val);
        // {val} = root \ sub

        it = new_it;
    }

    return root_simplex;
}

template<typename iterable_t>
Simplex* SimplexSComplex::get_simplex(const iterable_t &s)
{
    // assert(s.size() > 0);

    long cv = -1;

    if (memo.size())
    {
        cv = get_int(s);

        if (cv >= 0 && memo[cv])
        {
            return memo[cv];
        }
    }
    // vector<int> v(s.begin(), s.end()); // sorted

    for (typename iterable_t::const_iterator it = s.begin(), end = s.end(); it != end; ++it)
    {
        int vert = *it;

        if (base[vert] == 0)
            return 0;
    }

    typename iterable_t::const_iterator it = s.begin();
    int vert = *it++;

    Simplex *simplex = base[vert];

    const size_t n = s.size();

    for (size_t i = 1u; i < n; i++)
    {
        // assert(simplex->added_cells.size() == simplex->coborder.size());
        const int new_vert = *it++;

        Simplex *coface = 0;

        const size_t added_cells_size = simplex->added_cells.size();

        for (size_t j = 0u; j < added_cells_size; j++)
        {
            if (simplex->added_cells[j] == new_vert)
            {
                coface = simplex->coborder[j];
                break;
            }
        }

        if (coface == 0)
        {
            return 0;
        }

        simplex = coface;
    }

    if (cv >= 0)
        memo[cv] = simplex;

    return simplex;
}

#endif
