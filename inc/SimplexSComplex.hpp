#ifndef SIMPLEXSCOMPLEX_HPP
#define SIMPLEXSCOMPLEX_HPP

#include "simplex_iterators.hpp"
#include "SimplexCell.hpp"

#include <boost/iterator/transform_iterator.hpp>

class SimplexSComplex
{
public:

    int cardinality()
    {
        return distance(all_begin(), all_end());
    }

	typedef SimplexCell Cell;

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
                const Simplex &cell;
                int color;

                explicit CbdCells(SimplexSComplex &cmplx, const Simplex &c, int col) : comp(cmplx), cell(c), color(col) {}

                iterator begin()
                {
                    return const_cast<Simplex&>(cell).coborder_begin(color);
                }

                iterator end()
                {
                    return const_cast<Simplex&>(cell).coborder_end(color);
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

            CbdCells cbdCells(const Simplex &c)
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

        vertex_set tmp(s.begin(), s.end());

        int i = 0;

        for (vertex_set::const_iterator it = s.begin(); it != s.end(); ++it, ++i)
        {
            const int &val = *it;

            tmp.erase_index(i, val);
            Simplex *sub_simplex = add_simplex(tmp);
            tmp.insert_index(i, val);

            root_simplex->add_to_border(*sub_simplex);
            sub_simplex->add_to_coborder(*root_simplex, val);
            // {val} = root \ sub
        }

        return root_simplex;
    }

public:

    SimplexSComplex() : base(MAX_VERTICES, static_cast<Simplex*>(0))
    {
    }

    ~SimplexSComplex()
    {
        for (per_dimension_storage_type::const_iterator it = all_simplices.begin(); it != all_simplices.end(); ++it)
			delete *it;
    }

    typedef vector<Simplex*> memo_type;
    memo_type memo;

    long get_int(simple_set<int> &s)
    {
    	return s.getInt();
    }

    long get_int(set<int> &s)
    {
    	return -1;
    }

	template<typename iterable_t>
    Simplex* get_simplex(iterable_t &s)
    {
        // assert(s.size() > 0);

        long cv = -1;

        if (memo.size())
        {
        	cv = get_int(s);

			// cout << "got: " << cv << " ";
			if (cv >= 0 && memo[cv])
			{
				// cout << "found memo: " << cv << " ";
				// return memo[cv];
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

    inline bool getUniqueCoFace(const Cell& cell, Cell& coface) const
    {
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

    Simplex* add_simplex(set<int> &s)
    {
    	memo.resize(1u<<s.size(), 0);

    	vertex_set faster(s.begin(), s.end());
        Simplex *ret = add_simplex(faster);

        memo.clear();

        return ret;
    }
};

#endif
