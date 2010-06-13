#include "redHom/complex/simplicial/SimplexSComplex.hpp"

SimplexSComplex::Simplex* SimplexSComplex::createSimplexHierarchy(vertex_set &s)
{
    assert(s.size() > 0);

    if (s.size() == 1)
    {
        return makeBaseSimplex(s);
    }

    Simplex* rootSimplex = new Simplex(s.begin(), s.end(), nextId++);
    simplexAddedEvent(*rootSimplex);

    // set<int> probably can't be used as a drop-in replacement...
    for (vertex_set::iterator it = s.begin(), end = s.end(); it != end; )
    {
        const int &val = *it;
        vertex_set::iterator new_it = it;
        ++new_it;

        s.erase(it);
        Simplex *subSimplex = addSimplex(s);
        s.insert(it, val);

        rootSimplex->addToBorder(*subSimplex);
        subSimplex->addToCoborder(*rootSimplex, val);
        // {val} = root \ sub

        it = new_it;
    }

    return rootSimplex;
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
        // assert(simplex->addedCells.size() == simplex->coborder.size());
        const int new_vert = *it++;

        Simplex *coface = 0;

        const size_t addedCellsSize = simplex->addedCells.size();

        for (size_t j = 0u; j < addedCellsSize; j++)
        {
            if (simplex->addedCells[j] == new_vert)
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
