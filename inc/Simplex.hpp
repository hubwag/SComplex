#ifndef SIMPLEXHPP
#define SIMPLEXHPP

#include <vector>
#include <numeric>
#include <algorithm>

#include <boost/iterator/filter_iterator.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/optional.hpp>

#include <boost/ref.hpp>

#include "SimplexCell.hpp"

using namespace std;

class SimplexSComplex;
class SimplexCell;


struct Simplex
{
    typedef char Color;
    Color color;

    vector<int> nrs; // not necessary?? can be retreived from border (recursively)?
    vector<Simplex*> border;

    vector<Simplex*> coborder;
    vector<int> added_cells;

    template<typename iter_t>
    Simplex(iter_t b, iter_t e) : color(1), nrs(b, e) {}

    Simplex() : color(1) {}

    Simplex(SimplexSComplex &) : color(1) {}

    void debug_output__()
    {
    	// int border_size = distance(border_begin(), border_end());
    	// int coborder_size = distance(coborder_begin(), coborder_end());

    	// cout << "border: " << border_size << " coborder: " << coborder_size << endl;
    }

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
            // cout << "c: " << (int)s->color << " ";
            return s->color == required_color;
        }
    };

    int getDim() const
    {
        return -1 + static_cast<int>(nrs.size());
    }

/*
    bool getFaceCompanion(Cell& companion)   // should be const, requires changes in isFreeCoFace to be const
    {
    	int cnt = std::distance(coborder_begin(1), coborder_end(1));
    	if (cnt == 1)
    	{
    		companion = **(coborder.begin());
    		return true;
    	}
        return false;
    }
    */


    /*Cell& getElementaryCell()
    {
    	return *this;
    }*/
/*
	template<typename ComplexT>
    inline boost::optional<Cell> getUniqueCoFace(const ComplexT& complex)
    {
        int cnt = std::distance(coborder_begin(1), coborder_end(1));

    	if (cnt == 1)
    	{
    		return boost::optional<Cell>(**(coborder.begin()));
    	}

        return boost::optional<Cell>();
    }

    template<typename ComplexT>
    inline boost::optional<Cell> getUniqueFace(const ComplexT& complex)
    {
        int cnt = std::distance(border_begin(1), border_end(1));

    	if (cnt == 1)
    	{
    		return boost::optional<Cell>(**(border.begin()));
    	}

        return boost::optional<Cell>();
    }
    */

    typedef boost::filter_iterator<has_color, vector<Simplex*>::iterator> colored_iterator_inner;
    typedef boost::filter_iterator<has_color, vector<Simplex*>::iterator> border_iterator_inner;
    typedef boost::filter_iterator<has_color, vector<Simplex*>::iterator> coborder_iterator_inner;

	class to_cell : public unary_function<Simplex*, SimplexCell&>
	{
		mutable SimplexCell ret_cell;
		public:
		SimplexCell& operator()(Simplex *s) const
		{
			return ret_cell = SimplexCell(s);
		}
	};

	typedef boost::transform_iterator<to_cell, colored_iterator_inner> colored_iterator;
	typedef boost::transform_iterator<to_cell, colored_iterator_inner> border_iterator;
	typedef boost::transform_iterator<to_cell, colored_iterator_inner> coborder_iterator;

    /*typedef boost::indirect_iterator<colored_iterator_inner> colored_iterator;
    typedef boost::indirect_iterator<border_iterator_inner> border_iterator;
    typedef boost::indirect_iterator<coborder_iterator_inner> coborder_iterator;*/

    typedef boost::filter_iterator<has_color, vector<Simplex*>::const_iterator> const_colored_iterator_inner;
    typedef boost::filter_iterator<has_color, vector<Simplex*>::const_iterator> const_border_iterator_inner;
    typedef boost::filter_iterator<has_color, vector<Simplex*>::const_iterator> const_coborder_iterator_inner;

    typedef boost::transform_iterator<to_cell, const_colored_iterator_inner> const_colored_iterator;
	typedef boost::transform_iterator<to_cell, const_colored_iterator_inner> const_border_iterator;
	typedef boost::transform_iterator<to_cell, const_colored_iterator_inner> const_coborder_iterator;

    /*typedef boost::indirect_iterator<const_colored_iterator_inner> const_colored_iterator;
    typedef boost::indirect_iterator<const_border_iterator_inner> const_border_iterator;
    typedef boost::indirect_iterator<const_coborder_iterator_inner> const_coborder_iterator;*/

    border_iterator border_begin(int color = 1)
    {
        return border_iterator(border_iterator_inner(has_color(color), border.begin(), border.end()), to_cell() );
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

template<typename simplex_t>
void add_to_border(simplex_t &s, simplex_t &to_be_added)
{
    s.border.push_back(&to_be_added);
}

template<typename simplex_t>
void add_to_coborder(simplex_t &s, simplex_t &to_be_added, int added_cell = -1)
{
    s.coborder.push_back(&to_be_added);

    if (added_cell != -1)
        s.added_cells.push_back(added_cell);
}

/*
template<typename simplex_t>
void remove_from_border(simplex_t &s, simplex_t &to_be_removed)
{
    // s.border.erase(find(s.border.begin(), s.border.end(), &to_be_removed));
    swap(*s.border.back(), to_be_removed);
    s.border.pop_back();
}

template<typename simplex_t>
void remove_from_coborder(simplex_t &s, simplex_t &to_be_removed)
{
    swap(*s.coborder.back(), to_be_removed);
    s.coborder.pop_back();
}
*/

/*

template<typename simplex_t>
void remove_simplex(simplex_t &s)
{
    for (int i = 0; i < s.border.size(); i++)
        remove_from_coborder(*s.border[i], s);

    for (int i = 0; i < s.coborder.size(); i++)
        remove_from_border(*s.coborder[i], s);
}
*/

// typedef Simplex Cell; // ?

#endif
