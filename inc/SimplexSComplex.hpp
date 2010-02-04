#ifndef SIMPLEXSCOMPLEX_HPP
#define SIMPLEXSCOMPLEX_HPP

#include "simplex_iterators.hpp"

#include <boost/iterator/transform_iterator.hpp>

class SimplexSComplex
{
public:

    int cardinality()
    {
        return distance(all_begin(), all_end());
    }

    typedef Simplex::Cell Cell;

    struct ColoredIterators
    {
        struct Iterators
        {
            const int color;
            SimplexSComplex &comp;

            Iterators(SimplexSComplex &cmplx, int c) : color(c), comp(cmplx)
            {
            }

            struct DimCells
            {
                typedef Simplex::colored_iterator iterator;

                SimplexSComplex &comp;
                int color;
                int dim;

                explicit DimCells(SimplexSComplex &cmplx, int col, int d) : comp(cmplx), color(col), dim(d) {}

                iterator begin()
                {
                    return comp.dim_begin(dim, color);
                }

                iterator end()
                {
                    return comp.dim_end(dim, color);
                }
            };

            DimCells dimCells(int dim)
            {
                return DimCells(comp, color, dim);
            }

            struct BdCells
            {
                // typedef vector<> iterator;
            };

            BdCells bdCells();

            struct AllCells
            {
                typedef Simplex::colored_iterator iterator;
                SimplexSComplex &comp;
                int color;

                explicit AllCells(SimplexSComplex &cmplx, int col) : comp(cmplx), color(col) {}

                iterator begin()
                {
                    return comp.all_begin(color);
                }

                iterator end()
                {
                    return comp.all_end(color);
                }
            };

            AllCells allCells()
            {
                return AllCells(comp, color);
            }

            struct CbdCells
            {
                typedef Simplex::coborder_iterator iterator;
                typedef Simplex::const_coborder_iterator const_iterator;

                // typedef boost::filter_iterator<Simplex::has_color, vector<Simplex*>::iterator> iterator;
                // typedef boost::filter_iterator<Simplex::has_color, vector<Simplex*>::const_iterator> const_iterator;
                SimplexSComplex &comp;
                const Cell &cell;
                int color;

                explicit CbdCells(SimplexSComplex &cmplx, const Cell &c, int col) : comp(cmplx), cell(c), color(col) {}

                iterator begin()
                {
                    return const_cast<Cell&>(cell).coborder_begin(color);
                }

                iterator end()
                {
                    return const_cast<Cell&>(cell).coborder_end(color);
                }

                const_iterator begin() const
                {
                    return cell.coborder_begin(color);
                }

                const_iterator end() const
                {
                    return cell.coborder_end(color);
                }
            };

            CbdCells cbdCells(const Cell &c)
            {
                return CbdCells(comp, c, color);
            }
        };

    };


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

    typedef set<int> vertex_set;

    map<int, per_dimension_storage_type> per_dimension;
    per_dimension_storage_type all_simplices;

    static const int MAX_VERTICES = 100000;
    vector<Simplex*> base;

    void simplex_added_event(Simplex &s)
    {
        per_dimension[s.get_dim()].push_back(&s);
        all_simplices.push_back(&s);
    }

    void simplex_removed_event(Simplex &s)
    {
        // per_dimension[s.get_dim()].push_back(&s);
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

    Simplex* create_simplex_hierarchy(vertex_set &s)
    {
        assert(s.size() > 0);

        if (s.size() == 1)
        {
            return make_base_simplex(s);
        }

        Simplex* root_simplex = new Simplex(s.begin(), s.end());
        simplex_added_event(*root_simplex);

        set<int> tmp(s.begin(), s.end());

        for (auto it = s.begin(); it != s.end(); ++it)
        {
            int val = *it;

            tmp.erase(*it);
            Simplex *sub_simplex = add_simplex(tmp);
            tmp.insert(*it);

            add_to_border(*root_simplex, *sub_simplex);
            add_to_coborder(*sub_simplex, *root_simplex, val);
        }

        return root_simplex;
    }

public:

    SimplexSComplex() : base(MAX_VERTICES, static_cast<Simplex*>(0))
    {
    }

    ~SimplexSComplex()
    {
        for (auto it = per_dimension.begin(); it != per_dimension.end(); ++it)
            for (auto inner = it->second.begin(); inner != it->second.end(); ++inner)
                delete *inner;
    }

    Simplex* get_simplex(vertex_set &s)
    {
        assert(s.size() > 0);
        // vector<int> v(s.begin(), s.end()); // sorted

        for (auto it = s.begin(); it != s.end(); ++it)
        {
            int vert = *it;

            if (base[vert] == 0)
                return 0;
        }

        auto it = s.begin();
        int vert = *it++;

        Simplex *simplex = base[vert];

        for (size_t i = 1; i < s.size(); i++)
        {
            assert(simplex->added_cells.size() == simplex->coborder.size());
            const int new_vert = *it++;

            Simplex *coface = 0;

            for (size_t j = 0; j < simplex->added_cells.size(); j++)
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

        // cerr << "found!!";

        return simplex;
    }

    inline bool getUniqueCoFace(const Cell& cell, Cell& coface) const
    {
        int cnt = std::distance(cell.coborder_begin(1), cell.coborder_end(1));

        if (cnt == 1)
        {
            coface = **(cell.coborder.begin());
            return true;
        }

        return false;
    }

    inline bool getUniqueFace(const Cell& cell, Cell& face) const
    {
        int cnt = std::distance(cell.border_begin(1), cell.border_end(1));

        if (cnt == 1)
        {
            face = **(cell.border.begin());
            return true;
        }

        return false;
    }

    typedef Simplex::colored_iterator iterator;

    iterator all_begin(int color = 1)
    {
        return Simplex::colored_iterator( Simplex::colored_iterator_inner(Simplex::has_color(color), all_simplices.begin(), all_simplices.end()) );
    }

    iterator all_end(int color = 1)
    {
        return Simplex::colored_iterator( Simplex::colored_iterator_inner(Simplex::has_color(color), all_simplices.end(), all_simplices.end()) );
    }

    iterator dim_begin(int dim, int color = 1)
    {
        return Simplex::colored_iterator( Simplex::colored_iterator_inner(Simplex::has_color(color), per_dimension[dim].begin(), per_dimension[dim].end()) );
    }

    iterator dim_end(int dim, int color = 1)
    {
        return Simplex::colored_iterator( Simplex::colored_iterator_inner(Simplex::has_color(color), per_dimension[dim].end(), per_dimension[dim].end()) );
        // return Simplex::colored_iterator(Simplex::has_color(color), per_dimension[dim].end(), per_dimension[dim].end());
    }

    Simplex* add_simplex(vertex_set &s)
    {
        Simplex *found = get_simplex(s);

        if (found)
            return found;

        return create_simplex_hierarchy(s);
    }
};

#endif
