#ifndef SIMPLEXSCOMPLEX_HPP
#define SIMPLEXSCOMPLEX_HPP

#include <vector>
#include <map>
#include <boost/iterator/transform_iterator.hpp>

#include "SimplexCell.hpp"
#include "simple_set.h"
#include "BasicCellProxy.hpp"

using namespace std;

class SimplexSComplex
{
public:

    template<typename ImplT>
	class CellProxy: public BasicCellProxy<ImplT>
    {
    public:
        CellProxy(const ImplT& impl): BasicCellProxy<ImplT>(impl) {}
        CellProxy(): BasicCellProxy<ImplT>(0) {} // needed by map

        CellProxy(SimplexSComplex&): BasicCellProxy<ImplT>(0) {}
    };

    typedef int Dim;
    typedef int Id;

    struct Simplex;
    typedef CellProxy<Simplex*> Cell;

    const Dim& getMaxDim() const
    {
    	return per_dimension.rbegin()->first;
    }

    struct Simplex
    {
        typedef char Color;
        Color color;

        typedef int Dim;
        typedef int Id;

        Id id;
		static Id nextId;

        Id getId() const
        {
            return id;
        }

        vector<int> nrs; // not necessary?? can be retreived from border (recursively)?
        vector<Simplex*> border;

        vector<Simplex*> coborder;
        vector<int> added_cells; // could use radix_map<int, Simplex*>  or radix_map<int, int>
        // which would be useful if codorders are big

        template<typename iter_t>
        Simplex(iter_t b, iter_t e) : color(1), nrs(b, e) { id = nextId++; }

        Simplex() : color(1) { id = nextId++; }

        Simplex(SimplexSComplex &) : color(1) { id = nextId++; }

        // border, coborder operations

        void addToBorder(Simplex &to_be_added)
        {
            border.push_back(&to_be_added);
        }

        void addToCoborder(Simplex &to_be_added, int added_cell = -1)
        {
            coborder.push_back(&to_be_added);

            if (added_cell != -1)
                added_cells.push_back(added_cell);
        }

        // colors

        int getColor() const
        {
            return color;
        }

        template<Color col>
        void setColor()
        {
            color = col;
        }

        void setColor(Color col)
        {
            color = col;
        }

        struct has_color
        {
            Color required_color;
            has_color(Color color) : required_color(color) {}

            bool operator()(const Simplex *s) const
            {
                return required_color == -1 || s->color == required_color;
            }
        };

        int getDim() const
        {
            return -1 + static_cast<int>(nrs.size());
        }

        bool operator<(const Simplex& b) const
        {
            return this->nrs < b.nrs;
        }

        // iterators...

        typedef boost::filter_iterator<has_color, vector<Simplex*>::iterator> colored_iterator_inner;
        typedef boost::filter_iterator<has_color, vector<Simplex*>::iterator> border_iterator_inner;
        typedef boost::filter_iterator<has_color, vector<Simplex*>::iterator> coborder_iterator_inner;

        class to_cell : public unary_function<Simplex*, Cell>
        {

        public:
            Cell operator()(Simplex *s) const
            {
                return Cell(s);
            }
        };

        typedef boost::transform_iterator<to_cell, colored_iterator_inner> colored_iterator;
        typedef boost::transform_iterator<to_cell, colored_iterator_inner> border_iterator;
        typedef boost::transform_iterator<to_cell, colored_iterator_inner> coborder_iterator;

        typedef boost::filter_iterator<has_color, vector<Simplex*>::const_iterator> const_colored_iterator_inner;
        typedef boost::filter_iterator<has_color, vector<Simplex*>::const_iterator> const_border_iterator_inner;
        typedef boost::filter_iterator<has_color, vector<Simplex*>::const_iterator> const_coborder_iterator_inner;

        typedef boost::transform_iterator<to_cell, const_colored_iterator_inner> const_colored_iterator;
        typedef boost::transform_iterator<to_cell, const_colored_iterator_inner> const_border_iterator;
        typedef boost::transform_iterator<to_cell, const_colored_iterator_inner> const_coborder_iterator;

        border_iterator border_begin(int color = 1)
        {
            return border_iterator( border_iterator_inner(has_color(color), border.begin(), border.end()), to_cell() );
        }

        border_iterator border_end(int color = 1)
        {
            return border_iterator( border_iterator_inner(has_color(color), border.end(), border.end()), to_cell() );
        }

        const_border_iterator border_begin(int color = 1) const
        {
            return const_border_iterator( const_border_iterator_inner(has_color(color), border.begin(), border.end()), to_cell() );
        }

        const_border_iterator border_end(int color = 1) const
        {
            return const_border_iterator( const_border_iterator_inner(has_color(color), border.end(), border.end()), to_cell() );
        }

        coborder_iterator coborder_begin(Color color = 1)
        {
            return coborder_iterator( coborder_iterator_inner(has_color(color), coborder.begin(), coborder.end()), to_cell() );
        }

        coborder_iterator coborder_end(Color color = 1)
        {
            return coborder_iterator( coborder_iterator_inner(has_color(color), coborder.end(), coborder.end()), to_cell() );
        }

        const_coborder_iterator coborder_begin(Color color = 1) const
        {
            return const_coborder_iterator( const_coborder_iterator_inner(has_color(color), coborder.begin(), coborder.end()), to_cell() );
        }

        const_coborder_iterator coborder_end(Color color = 1) const
        {
            return const_coborder_iterator( const_coborder_iterator_inner(has_color(color), coborder.end(), coborder.end()), to_cell() );
        }
    };

    int size()
    {
        return distance(all_begin(), all_end()); // slow
    }

    int cardinality()
    {
        return distance(all_begin(), all_end()); // slow
    }

    // typedef SimplexCell Cell;


// to be changed later!!
#include "SimplexInternalColoredIterators.hpp"
// to be changed later!!

	Iterators iterators()
    {
        return Iterators(*this);
    }

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

    static const int MAX_VERTICES = 10000;
    vector<Simplex*> base;

    void simplexAddedEvent(Simplex &s)
    {
        per_dimension[s.getDim()].push_back(&s);
        all_simplices.push_back(&s);
    }

    void resizeBase(int v)
    {
    	while(v >= base.size())
			base.push_back(0);
    }

    Simplex* makeBaseSimplex(vertex_set &s)
    {
        assert(s.size() == 1);

        const int v = *s.begin();

        if (v >= base.size())
        {
        	resizeBase(v);
        }

        if (base[v] == 0)
        {
            Simplex *new_simplex = new Simplex(s.begin(), s.end());
            simplexAddedEvent(*new_simplex);
            base[v] = new_simplex;

            return new_simplex;
        }

        return base[v];
    }

    Simplex* createSimplexHierarchy(vertex_set &s);

public:

    SimplexSComplex() : base(MAX_VERTICES, static_cast<Simplex*>(0))
    {
    	Simplex::nextId = 0;
    }

    ~SimplexSComplex()
    {
        for (per_dimension_storage_type::const_iterator it = all_simplices.begin(), end = all_simplices.end(); it != end; ++it)
            delete *it;
    }

    long getInt(const simple_set<int> &s)
    {
        return s.getInt();
    }

    long getInt(const set<int> &s)
    {
        return -1;
    }

    template<typename iterable_t>
    Simplex* getSimplex(const iterable_t &s);

    inline bool getUniqueCoFace(const Cell& cell, Cell& coface) const
    {
        // could be optimized...
        int cnt = std::distance(cell.getImpl()->coborder_begin(1), cell.getImpl()->coborder_end(1));

        if (cnt == 1)
        {
            coface = *(cell.getImpl()->coborder_begin(1));
            return true;
        }

        return false;
    }

    inline bool getUniqueFace(const Cell& cell, Cell& face) const
    {
        int cnt = std::distance(cell.getImpl()->border_begin(1), cell.getImpl()->border_end(1));

        if (cnt == 1)
        {
            face = *(cell.getImpl()->border_begin(1));
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

    Simplex* addSimplex(vertex_set &s)
    {
        Simplex *found = getSimplex(s);

        if (found)
            return found;

        return createSimplexHierarchy(s);
    }

public:

    int coincidenceIndex(const Cell &a, const Cell &b) const
    {
        const Simplex * const sigma = a.getImpl();
        const Simplex * const tau = b.getImpl();

        if (sigma->getDim() != 1 + tau->getDim())
        {
            return 0;
        }

        int inter = -1;
        set_difference(sigma->nrs.begin(), sigma->nrs.end(), tau->nrs.begin(), tau->nrs.end(), &inter);
        assert(inter != -1);

        int i = distance(sigma->nrs.begin(),
                         find(sigma->nrs.begin(), sigma->nrs.end(), inter));

        return (i % 2 == 0) ? 1 : -1; // (-1)**i
    }

    Simplex* addSimplex(const set<int> &s)
    {
        vertex_set faster(s.begin(), s.end());
        Simplex *ret = addSimplex(faster);

        return ret;
    }
};

SimplexSComplex::Simplex* SimplexSComplex::createSimplexHierarchy(vertex_set &s)
{
    assert(s.size() > 0);

    if (s.size() == 1)
    {
        return makeBaseSimplex(s);
    }

    Simplex* root_simplex = new Simplex(s.begin(), s.end());
    simplexAddedEvent(*root_simplex);

    // set<int> probably can't be used as a drop-in replacement...
    for (vertex_set::iterator it = s.begin(), end = s.end(); it != end; )
    {
        const int &val = *it;
        vertex_set::iterator new_it = it;
        ++new_it;

        s.erase(it);
        Simplex *sub_simplex = addSimplex(s);
        s.insert(it, val);

        root_simplex->addToBorder(*sub_simplex);
        sub_simplex->addToCoborder(*root_simplex, val);
        // {val} = root \ sub

        it = new_it;
    }

    return root_simplex;
}

template<typename iterable_t>
SimplexSComplex::Simplex* SimplexSComplex::getSimplex(const iterable_t &s)
{
    // assert(s.size() > 0);

    for (typename iterable_t::const_iterator it = s.begin(), end = s.end(); it != end; ++it)
    {
        int vert = *it;

        if (vert >= base.size() || base[vert] == 0)
            return 0;
    }

    typename iterable_t::const_iterator it = s.begin();
    int vert = *it++;

    if (vert >= base.size())
		resizeBase(vert);

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

    return simplex;
}
SimplexSComplex::Simplex::Id SimplexSComplex::Simplex::nextId = 0;

#endif
