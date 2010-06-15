#ifndef ALGS_AKQ_REDUCE_STRATEGY_HPP_
#define ALGS_AKQ_REDUCE_STRATEGY_HPP_

#include <boost/optional.hpp>
#include <boost/foreach.hpp>
#include <utility>
#include <stack>
#include <set>
using std::set;

#include "redHom/complex/scomplex/SComplex.hpp"
#include "redHom/complex/scomplex/SComplexDefaultTraits.hpp"

#include "redHom/algorithm/strategy/DefaultReduceStrategy.hpp"
#include "redHom/algorithm/strategy/DefaultReduceStrategy_CubSComplex.hpp"
#include "redHom/algorithm/strategy/DefaultReduceStrategyTraits_CubSComplex.hpp"

// #define INLINE_CRITICAL_CODE __attribute__ ((noinline))
#define INLINE_CRITICAL_CODE


template<typename SComplexT>
class AKQReduceStrategy: public DefaultReduceStrategy<SComplexT>
{
protected:
    using DefaultReduceStrategyBase<SComplexT>::complex;
public:
    enum AKQType {UNSET, KING, QUEEN, ACE};
    typedef SComplexT SComplex;
    typedef DefaultReduceStrategyTraits<SComplexT> Traits;
    typedef typename SComplex::Cell Cell;

    typedef ::SComplex<SComplexDefaultTraits> OutputComplexType;

    AKQReduceStrategy(SComplex& _complex) :
		DefaultReduceStrategy<SComplex>(_complex),
		extractIt(_complex.iterators(1).dimCells(0).begin()),
		extractEnd(_complex.iterators(1).dimCells(0).end())
    {
      max_d = complex.getDim();

      size_t complexSize = _complex.size();
        morse.resize(complexSize);
        akq.resize(complexSize);
        kerKing.resize(complexSize, Cell(complex));

        extractDim = 0; // bottom up
    }

    bool markNullPaths(Cell curr, Cell from)
    {
    	if (akq[curr.getId()] == ACE && curr.getId() != from.getId())
			return true;

    	if (akq[curr.getId()] == UNSET)
    		return false;

        int ourMorseValue = morse[curr.getId()];

		bool ok = false;

        // case 1: to ace
		BOOST_FOREACH(Cell to, complex.iterators().bdCells(curr))
		{
			if (akq[to.getId()] == ACE && morse[to.getId()] < ourMorseValue)
			{
				ok = true;
			}
		}
		// case 2: to king
		BOOST_FOREACH(Cell bd, complex.iterators().bdCells(curr))
		{
			if (akq[bd.getId()] != QUEEN)
				continue;

			Cell to = kerKing[bd.getId()];

			if (akq[to.getId()] != KING)
			{
				continue;
			}

			BOOST_ASSERT(akq[to.getId()] == KING);

			if (morse[to.getId()] < ourMorseValue)
			{
				// S.push(std::make_pair(to, accumulatedWeight * toKingCoeff(curr, to, bd)));
				if (markNullPaths(to, from))
				{
					ok = true;
				}
			}
		}

		if (!ok && akq[curr.getId()] != ACE)
		{
			akq[curr.getId()] = UNSET;
		}

		return ok;
    }


    template<typename ImplT>
    static bool reduced(const typename SComplex::template CellProxy<ImplT>& cell)
    {
        return cell.getColor() == 2;
    }

    template<typename cellT1>
	int calcMorseValue(const cellT1 &c)
    {
        int val = 0;
        BOOST_FOREACH(Cell el, complex.iterators().bdCells(c))
        {
            val = max(val, morse[el.getId()]);
        }

        return val + 1;
    }

    template<typename T1, typename T2>
	void kingGetsMarried(const T1 &king, const T2 &queen)
    {
    	int kId = king.getId();
    	int qId = queen.getId();

        akq[kId] = KING;
        akq[qId] = QUEEN;
        int v = calcMorseValue(king);
        morse[kId] = v;
        morse[qId] = v;
        kerKing[kId] = king;
    }

    template<typename ImplT1, typename ImplT2>
    void coreduce(const typename SComplex::template CellProxy<ImplT1>& a, const typename SComplex::template CellProxy<ImplT2>& b)
    {
        BOOST_ASSERT(a.getDim() != b.getDim());

        if (a.getDim() > b.getDim())
        {
            kingGetsMarried(a, b);
        }
        else
        {
            kingGetsMarried(b, a);
        }

        a.template setColor<2>();
        b.template setColor<2>();
    }

	INLINE_CRITICAL_CODE int toAceCoeff(Cell x, Cell y)
    {
        BOOST_ASSERT(akq[y.getId()] == ACE);
        BOOST_ASSERT(akq[x.getId()] != QUEEN);

        return complex.coincidenceIndex(x, y);
    }


	int INLINE_CRITICAL_CODE toKingCoeff(Cell x, Cell y, Cell y_star)
    {
        BOOST_ASSERT(akq[x.getId()] != QUEEN);
        BOOST_ASSERT(akq[y_star.getId()] == QUEEN);
        BOOST_ASSERT(akq[y.getId()] == KING);

        return -1 * complex.coincidenceIndex(x, y_star) / complex.coincidenceIndex(y, y_star);
    }

	void INLINE_CRITICAL_CODE followPath(Cell c)
    {
        std::stack<std::pair<Cell, int> > S; // no cycles!

        S.push(std::make_pair(c, 1));

        while (S.size())
        {
            Cell curr = S.top().first;
            int accumulatedWeight = S.top().second;
            S.pop();

            BOOST_ASSERT(akq[curr.getId()] != QUEEN);

            if (curr.getId() != c.getId() && akq[curr.getId()] == ACE)
            {
                // these are aces - small
                ++numPathsBetween[std::make_pair(c.getId(), curr.getId())];

                coeffs[std::make_pair(c.getId(), curr.getId())] += accumulatedWeight;

                continue;
            }

            int ourMorseValue = morse[curr.getId()];

            // case 1: to ace
            BOOST_FOREACH(Cell to, complex.iterators().bdCells(curr))
            {
                if (akq[to.getId()] == ACE && morse[to.getId()] < ourMorseValue)
                {
                    S.push(std::make_pair(to, accumulatedWeight * toAceCoeff(curr, to)));
                }
            }

            // case 2: to king
            BOOST_FOREACH(Cell bd, complex.iterators().bdCells(curr))
            {
                if (akq[bd.getId()] != QUEEN)
                    continue;

                Cell to = kerKing[bd.getId()];

                if (akq[to.getId()] != KING)
                {
                	continue;
                }

                BOOST_ASSERT(akq[to.getId()] == KING);

                if (morse[to.getId()] < ourMorseValue)
                {
                    S.push(std::make_pair(to, accumulatedWeight * toKingCoeff(curr, to, bd)));
                }
            }
        }
    }

    void reportPaths()
    {
        BOOST_FOREACH(Cell ace, aces)
        {
            markNullPaths(ace, ace);
        }

        BOOST_FOREACH(Cell ace, aces)
        {
            followPath(ace);
        }

        typedef std::pair<std::pair<int,int>,int> Triple;

        std::vector<OutputComplexType::Dim> dims(aces.size());

        std::map<size_t, size_t> from0;

        BOOST_FOREACH(Cell ace, aces)
        {
            size_t next = from0.size();
            dims[next] = ace.getDim();
            from0[ace.getId()] = next;
        }

        OutputComplexType::KappaMap kap;

        BOOST_FOREACH(Triple p, coeffs)
        {
            int coef = p.second;
            kap.push_back(boost::make_tuple(from0[p.first.first], from0[p.first.second], coef));
        }

        outputComplex = new OutputComplexType(3, dims, kap, 1);
    }

    template<typename ImplT1, typename ImplT2>
    static void reduce(const typename SComplex::template CellProxy<ImplT1>& a, const typename SComplex::template CellProxy<ImplT2>& b)
    {
        a.template setColor<2>();
        b.template setColor<2>();
    }

    template<typename ImplT>
    static void reduce(const typename SComplex::template CellProxy<ImplT>& cell)
    {
        cell.template setColor<2>();
    }

    typename Traits::Extract::result_type extract()
    {
		for (; extractDim <= max_d; )
		{
		  extractIt  = complex.iterators(1).dimCells(extractDim).begin();;
		  extractEnd = complex.iterators(1).dimCells(extractDim).end();;

			while(extractIt != extractEnd && extractIt->getColor() != 1)
			{
				++extractIt;
			}

			if (extractIt != extractEnd)
			{
				int v = 0;

				if (extractDim != 0)
					v = calcMorseValue(*extractIt);

				morse[extractIt->getId()] = v;
				aces.push_back(*extractIt);
				akq[extractIt->getId()] = ACE;

				// we increment the iterator here!
				return typename Traits::Extract::result_type::value_type(*extractIt);
			}

			if (++extractDim > max_d)
				break;

			extractIt = complex.iterators(1).dimCells(extractDim).begin();
			extractEnd = complex.iterators(1).dimCells(extractDim).end();
		}

        reportPaths();
        return typename Traits::Extract::result_type();
    }

    static typename Traits::ForceCoreduction::result_type forceCoreductionPair()
    {
        return typename Traits::ForceCoreduction::result_type();
    }

    OutputComplexType* getOutputComplex()
    {
        return outputComplex;
    }

protected:

	typename SComplex::ColoredIterators::Iterators::DimCells::iterator extractIt;
	typename SComplex::ColoredIterators::Iterators::DimCells::iterator extractEnd;

  int extractDim;
    std::vector<int> morse;
    std::vector<AKQType> akq;
    std::vector<Cell> kerKing;
    std::vector<Cell> aces;

    int max_d;
    std::map<std::pair<int,int>, int> coeffs;
    std::map<std::pair<int,int>, int> numPathsBetween; // only between aces - small

    OutputComplexType *outputComplex;
  using DefaultReduceStrategy<SComplexT>::complex;
};

#endif

