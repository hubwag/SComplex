#ifndef ALGS_AKQ_REDUCE_STRATEGY_HPP_
#define ALGS_AKQ_REDUCE_STRATEGY_HPP_

#include <boost/optional.hpp>
#include <boost/foreach.hpp>
#include <utility>
#include <stack>

#include <map>
#include <vector>

#include "redHom/complex/scomplex/SComplex.hpp"
#include "redHom/complex/scomplex/SComplexDefaultTraits.hpp"

#include <redHom/algorithm/Algorithms.hpp>
#include "redHom/algorithm/strategy/DefaultReduceStrategy.hpp"
#include "redHom/algorithm/Coreduction.hpp"

template<typename SComplexT>
class AKQReduceStrategy : public DefaultReduceStrategy<SComplexT>
{
protected:
  using DefaultReduceStrategyBase<SComplexT>::complex;

public:
    enum AKQType {UNSET, KING, QUEEN, ACE};

    typedef std::map<int,int> PathsInfo;

    std::vector<int> nullPathMemo;
    std::vector<PathsInfo> followMemoTable;

    typedef ::SComplex<SComplexDefaultTraits> OutputComplexType;
    typedef SComplexT SComplex;
    typedef DefaultReduceStrategyTraits<SComplex> Traits;
    typedef typename SComplex::Cell Cell;

    PathsInfo notSet;

    AKQReduceStrategy(SComplex& _complex): DefaultReduceStrategy<SComplex>(_complex)
    {
        maxExtractDim = _complex.getDim();
        extractDim = 0;

        int csize = _complex.size();

        morse.resize(csize, 0);
        nullPathMemo.resize(csize, -1);

        notSet[-1] = -1;
        followMemoTable.resize(csize, notSet);
        akq.resize(csize);
        kerKing.resize(csize, Cell(this->complex));
    }

    SComplex& getComplex() const
    {
        return complex;
    }

    template<typename CellT1, typename CellT2>
	void coreduce(CellT1 a, CellT2 b)
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

    typename Traits::Extract::result_type extract()
    {
        for (; extractDim <= maxExtractDim; extractDim++)
        {
            typename SComplex::ColoredIterators::Iterators::DimCells dimCells = complex.iterators(1).dimCells(extractDim);
            typename SComplex::ColoredIterators::Iterators::DimCells::iterator end = dimCells.end(),
                    it = dimCells.begin();

            if (it != end)
            {
                int v = 0;

                if (extractDim != 0)
                    v = calcMorseValue(*it);

                morse[it->getId()] = v;
                aces.push_back(*it);
                akq[it->getId()] = ACE;

                return typename Traits::Extract::result_type::value_type(*it);
            }
        }

        reportPaths();
        return typename Traits::Extract::result_type();
    }

    OutputComplexType* getOutputComplex()
    {
        return outputComplex;
    }

protected:

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
        akq[king.getId()] = KING;
        akq[queen.getId()] = QUEEN;
        int v = calcMorseValue(king);
        morse[king.getId()] = v;
        morse[queen.getId()] = v;
        kerKing[queen.getId()] = king;
    }

    int toAceCoeff(Cell x, Cell y)
    {
        BOOST_ASSERT(akq[y.getId()] == ACE);
        BOOST_ASSERT(akq[x.getId()] != QUEEN);

        return complex.coincidenceIndex(x, y);
    }

    int toKingCoeff(Cell x, Cell y, Cell y_star)
    {
        BOOST_ASSERT(akq[x.getId()] != QUEEN);
        BOOST_ASSERT(akq[y_star.getId()] == QUEEN);
        BOOST_ASSERT(akq[y.getId()] == KING);

        return -1 * complex.coincidenceIndex(x, y_star) / complex.coincidenceIndex(y, y_star);
    }

    void followPath(Cell c)
    {
        std::stack<std::pair<Cell, int> > S;

        S.push(std::make_pair(c, 1));

        while (S.size())
        {
            Cell curr = S.top().first;
            int accumulatedWeight = S.top().second;
            S.pop();

            // std::cout << curr.getId() << " d: " << curr.getDim() << " m: " << morse[curr.getId()] << std::endl;

            BOOST_ASSERT(akq[curr.getId()] != QUEEN);

            if (curr.getId() != c.getId() && akq[curr.getId()] == ACE)
            {
                coeffs[std::make_pair(c.getId(), curr.getId())] += accumulatedWeight;
                continue;
            }

            int ourMorseValue = morse[curr.getId()];

            // case 1: to ace
            BOOST_FOREACH(Cell to, complex.iterators().bdCells(curr))
            {
                if (akq[to.getId()] == ACE && morse[to.getId()] < ourMorseValue)
                {
                	//std::cout << curr.getId() << " d: " << curr.getDim() << " m: " << morse[curr.getId()] << "=>";
                	//std::cout << to.getId() << " d: " << to.getDim() << " m: " << morse[to.getId()] << std::endl;
                    S.push(std::make_pair(to, accumulatedWeight * toAceCoeff(curr, to)));
                }
            }

            // case 2: to king
            BOOST_FOREACH(Cell bd, complex.iterators().bdCells(curr))
            {
                if (akq[bd.getId()] != QUEEN)
                    continue;

                Cell to = kerKing[bd.getId()];

                if (akq[to.getId()] == KING && morse[to.getId()] < ourMorseValue)
                {
                	//std::cout << curr.getId() << " d: " << curr.getDim() << " m: " << morse[curr.getId()] << "=>";
                	//std::cout << to.getId() << " d: " << to.getDim() << " m: " << morse[to.getId()] << std::endl;
                    S.push(std::make_pair(to, accumulatedWeight * toKingCoeff(curr, to, bd)));
                }
            }
        }
    }

  const PathsInfo& followPathMemo(Cell curr, Cell from)
    {
    	PathsInfo &coeffs = followMemoTable[curr.getId()];

    	if (coeffs != notSet)
	  return coeffs;
	
    	if (akq[curr.getId()] == UNSET)
	  return emptyPathsInfo;

    	BOOST_ASSERT(akq[curr.getId()] != ACE || curr.getId() == from.getId());

		coeffs.clear();

        int ourMorseValue = morse[curr.getId()];

        // case 1: to ace
		BOOST_FOREACH(Cell to, this->complex.iterators().bdCells(curr))
		{
			if (akq[to.getId()] == ACE && morse[to.getId()] < ourMorseValue)
			{
				int newCoeff = toAceCoeff(curr, to);
				coeffs[to.getId()] += newCoeff;
			}
		}
		// case 2: to king
		BOOST_FOREACH(Cell bd, this->complex.iterators().bdCells(curr))
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
				int newCoeff = toKingCoeff(curr, to, bd);

				PathsInfo fromHere = followPathMemo(to, from);

				for (PathsInfo::const_iterator it = fromHere.begin(); it != fromHere.end(); ++it)
				{
					coeffs[it->first] += it->second*newCoeff;
				}
			}
		}

		return coeffs;
    }

    bool markNullPaths(Cell curr, Cell from)
    {
    	int& ok = nullPathMemo[curr.getId()];
    	if (ok != -1)
			return static_cast<bool>(ok); // true or false

    	if (akq[curr.getId()] == ACE && curr.getId() != from.getId())
			return true;

    	if (akq[curr.getId()] == UNSET)
    		return false;

        int ourMorseValue = morse[curr.getId()];

        ok = false;

        // case 1: to ace
		BOOST_FOREACH(Cell to, this->complex.iterators().bdCells(curr))
		{
			if (akq[to.getId()] == ACE && morse[to.getId()] < ourMorseValue)
			{
				ok = true;
				break;
			}
		}
		// case 2: to king
		BOOST_FOREACH(Cell bd, this->complex.iterators().bdCells(curr))
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

    void reportPaths()
    {
    	/*
        BOOST_FOREACH(Cell ace, aces)
        {
		    markNullPaths(ace, ace);
        }

        BOOST_FOREACH(Cell ace, aces)
        {
			followPath(ace);
        }
        */

        BOOST_FOREACH(Cell ace, aces)
        {
            PathsInfo cf = followPathMemo(ace, ace);
            for (PathsInfo::const_iterator it = cf.begin(); it != cf.end(); ++it)
            {
				coeffs[std::make_pair(ace.getId(), it->first)] += it->second;
            }
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

    std::vector<int> morse;
    std::vector<AKQType> akq;
    std::vector<Cell> kerKing;
    std::vector<Cell> aces;



	int extractDim;
    int maxExtractDim;
    std::map<std::pair<int,int>, int> coeffs;
    PathsInfo emptyPathsInfo;
    OutputComplexType *outputComplex;
};

#endif

