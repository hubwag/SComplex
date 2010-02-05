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

    // border, coborder operations

    void add_to_border(Simplex &to_be_added)
    {
        border.push_back(&to_be_added);
    }

    void add_to_coborder(Simplex &to_be_added, int added_cell = -1)
    {
        coborder.push_back(&to_be_added);

        if (added_cell != -1)
            added_cells.push_back(added_cell);
    }

    void debug_output__()
    {
        // int border_size = distance(border_begin(), border_end());
        // int coborder_size = distance(coborder_begin(), coborder_end());

        // cout << "border: " << border_size << " coborder: " << coborder_size << endl;
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
            return s->color == required_color;
        }
    };

    int getDim() const
    {
        return -1 + static_cast<int>(nrs.size());
    }

    /*Cell& getElementaryCell()
    {
    	return *this;
    }*/

	// iterators...

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

    typedef boost::filter_iterator<has_color, vector<Simplex*>::const_iterator> const_colored_iterator_inner;
    typedef boost::filter_iterator<has_color, vector<Simplex*>::const_iterator> const_border_iterator_inner;
    typedef boost::filter_iterator<has_color, vector<Simplex*>::const_iterator> const_coborder_iterator_inner;

    typedef boost::transform_iterator<to_cell, const_colored_iterator_inner> const_colored_iterator;
    typedef boost::transform_iterator<to_cell, const_colored_iterator_inner> const_border_iterator;
    typedef boost::transform_iterator<to_cell, const_colored_iterator_inner> const_coborder_iterator;

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

#endif
