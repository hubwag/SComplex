#ifndef ALGS_AKQ_REDUCE_STRATEGY_HPP_
#define ALGS_AKQ_REDUCE_STRATEGY_HPP_

#include <boost/optional.hpp>
#include <boost/foreach.hpp>
#include <utility>
#include <stack>

#include "redHom/complex/scomplex/SComplex.hpp"
#include "redHom/complex/scomplex/SComplexDefaultTraits.hpp"

static const bool verboseAKQ = true;

template<typename SComplexT>
class AKQReduceStrategyTraits
{
public:

    template<typename ArgT>
struct GetReductionPair : public std::unary_function<const ArgT&,
                boost::optional<typename SComplexT::Cell> > {};

    template<typename ArgT>
struct GetCoreductionPair : public std::unary_function<const ArgT&,
                boost::optional<typename SComplexT::Cell> > {};

    struct ForceCoreduction
    {
        typedef boost::optional<std::pair<typename SComplexT::Cell,
        typename SComplexT::Cell> > result_type;
    };

    struct Extract
    {
        typedef boost::optional<typename SComplexT::Cell >  result_type;
    };
};

template<typename SComplexT>
class AKQReduceStrategyBase
{

public:
    enum AKQType {UNSET, KING, QUEEN, ACE};
    typedef SComplexT SComplex;
    typedef AKQReduceStrategyTraits<SComplex> Traits;
    typedef typename SComplex::Cell Cell;

    AKQReduceStrategyBase(SComplex& _complex): complex(_complex), dummyCell2(_complex),  dummyCell3(_complex)
    {
        std::cout << "CHECKING MAX DIM BRUTALLY! ADD getMaxDim to SCOMPLEX!!!!" << std::endl;
        std::cout << "SIZE, CARDINALITY OR MAX_ID?" << std::endl;
        max_d = getMaxDim();

        morse.resize(_complex.size());
        akq.resize(_complex.size());
        kerKing.resize(_complex.size(), *complex.iterators(1).allCells().begin());

        /*
        morse.resize(3000000);
        akq.resize(3000000);
        kerKing.resize(3000000);
        */
    }

    SComplex& getComplex() const
    {
        return complex;
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
        //BOOST_FOREACH(Cell el, complex.iterators().bdCells(c))
        //{
        //    val = max(val, morse[el.getId()]);
        //}

        for (typename SComplexT::Iterators::BdCells::iterator it = complex.iterators().bdCells(c).begin();
			it != complex.iterators().bdCells(c).end(); ++it)
        {
            val = max(val, morse[(*it).getId()]);
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

    int toAceCoeff(Cell x, Cell y)
    {
        BOOST_ASSERT(akq[y.getId()] == ACE);
        BOOST_ASSERT(akq[x.getId()] != QUEEN);

        return complex.coincidenceIndex(x, y);
    }

    //
    int toKingCoeff(Cell x, Cell y, Cell y_star)
    {
        BOOST_ASSERT(akq[x.getId()] != QUEEN);
        BOOST_ASSERT(akq[y_star.getId()] == QUEEN);
        BOOST_ASSERT(akq[y.getId()] == KING);

        return -1 * complex.coincidenceIndex(x, y_star) / complex.coincidenceIndex(y, y_star);
    }

    void followPath(Cell c)
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
                if (verboseAKQ)
                {
                    std::cout << "found path from: " << c.getId() << " to " << curr.getId() << std::endl;
                    std::cout << "between values: " << morse[c.getId()] << "and " << morse[curr.getId()] << " with coef product" << accumulatedWeight << std::endl;
                }
                coeffs[std::make_pair(c.getId(), curr.getId())] += accumulatedWeight;

                continue;
            }

            int ourMorseValue = morse[curr.getId()];

            // case 1: to ace
            /*BOOST_FOREACH(Cell to, complex.iterators().bdCells(curr))
            {
                if (akq[to.getId()] == ACE && morse[to.getId()] < ourMorseValue)
                {
                    S.push(std::make_pair(to, accumulatedWeight * toAceCoeff(curr, to)));
                }
            }*/
            for (typename SComplexT::Iterators::BdCells::iterator it = complex.iterators().bdCells(c).begin();
				 it != complex.iterators().bdCells(c).end(); ++it)
			{
				Cell to = *it;
				if (akq[to.getId()] == ACE && morse[to.getId()] < ourMorseValue)
                {
                    S.push(std::make_pair(to, accumulatedWeight * toAceCoeff(curr, to)));
                }
			}

            // case 2: to king
            //BOOST_FOREACH(Cell bd, complex.iterators().bdCells(curr))
            for (typename SComplexT::Iterators::BdCells::iterator it = complex.iterators().bdCells(c).begin();
				 it != complex.iterators().bdCells(c).end(); ++it)
			{
				Cell bd = *it;

                if (akq[bd.getId()] != QUEEN)
                    continue;

                Cell to = kerKing[bd.getId()];
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
        std::cout << " asow bylo: " << aces.size() << std::endl;

        if (verboseAKQ)
            for (size_t i = 0; i < aces.size(); i++)
            {
                std::cout << aces[i].getId() << " ";
            }
        std::cout << std::endl;

        BOOST_FOREACH(Cell ace, aces)
        {
            followPath(ace);
        }

        std::cout << "\n\n\n";

        typedef std::pair<std::pair<int,int>,int> Triple;

        BOOST_FOREACH(Triple p, numPathsBetween)
        {
            if (p.second == 1)
            {
                std::cout << "There was a single path!" << std::endl;
                break;
            }
        }

        if (verboseAKQ)
        {
            std::cout << "num paths between: \n";
            BOOST_FOREACH(Triple p, numPathsBetween)
            {
                std::cout << p.first.first << " " << p.first.second << " = " << p.second << std::endl;
            }

            std::cout << "coefficients: \n";
            BOOST_FOREACH(Triple p, coeffs)
            {
                std::cout << p.first.first << " " << p.first.second << " => " << p.second << std::endl;
            }
        }

        std::vector<size_t> dims(aces.size());

        std::map<size_t, size_t> from0;

        BOOST_FOREACH(Cell ace, aces)
        {
            size_t next = from0.size();
            dims[next] = ace.getDim();
            from0[ace.getId()] = next;
        }

        std::cout << "constructing general SComplex" << std::endl;

        OutputComplexType::KappaMap kap;

        BOOST_FOREACH(Triple p, coeffs)
        {
            int coef = p.second;

            kap.push_back(boost::make_tuple(from0[p.first.first], from0[p.first.second], coef));

            if (verboseAKQ)
                std::cout << from0[p.first.first] << "[d=" <<dims[from0[p.first.first]] << "]" << " : " <<  from0[p.first.second] << "[d="<<dims[from0[p.first.second]] << "]"  << " => " << coef << std::endl;
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
        // for (int d = 0; d <= complex.getMaxDim(); d++)
        for (int d = 0; d <= max_d; d++)
        {
            typename SComplex::ColoredIterators::Iterators::DimCells dimCells = complex.iterators(1).dimCells(d);
            typename SComplex::ColoredIterators::Iterators::DimCells::iterator end = dimCells.end(),
                    it = dimCells.begin();

            if (it != end)
            {
                int v = 0;

                if (d != 0)
                    v = calcMorseValue(*it);

                morse[it->getId()] = v;
                aces.push_back(*it);
                akq[it->getId()] = ACE;

                // cout << "AS!" << it->getId() << ' ' << it->getDim() << "with M = " << v << endl;

                return typename Traits::Extract::result_type::value_type(*it);
            }
        }

        reportPaths();
        return typename Traits::Extract::result_type();
    }

    static typename Traits::ForceCoreduction::result_type forceCoreductionPair()
    {
        return typename Traits::ForceCoreduction::result_type();
    }

    template<typename ArgT>
    typename Traits::template GetCoreductionPair<ArgT>::result_type
    getCoreductionPair(const ArgT& cell)
    {
        int times = 0;
        for (typename SComplexT::ColoredIterators::Iterators::BdCells::iterator it = complex.iterators(1).bdCells(cell).begin();
				 it != complex.iterators(1).bdCells(cell).end(); ++it)

    	//BOOST_FOREACH(typename SComplex::ColoredIterators::Iterators::BdCells::iterator::value_type v,
        //              complex.iterators(1).bdCells(cell))
		{
			// Cell v = *it;

            if (times == 0)
            {
                dummyCell3 = *it;
            }
            ++times;
            if (times == 2)
            {
                break;
            }
        }

        if (times == 1)
        {
            return typename Traits::template GetCoreductionPair<ArgT>::result_type(dummyCell3);
        }

        return typename Traits::template GetCoreductionPair<ArgT>::result_type();
    }

    template<typename ArgT>
    typename Traits::template GetReductionPair<ArgT>::result_type
    getReductionPair(const ArgT &cell)
    {
        int times = 0;

        //BOOST_FOREACH(typename SComplex::ColoredIterators::Iterators::CbdCells::iterator::value_type v,
        //              complex.iterators(1).cbdCells(cell))
		for (typename SComplexT::Iterators::CbdCells::iterator it = complex.iterators(1).cbdCells(cell).begin();
				 it != complex.iterators(1).cbdCells(cell).end(); ++it)
        {
        	Cell v = *it;

            if (times == 0)
            {
                dummyCell2 = v;
            }
            ++times;
            if (times == 2)
            {
                break;
            }
        }

        if (times == 1)
        {
            return typename Traits::template GetReductionPair<ArgT>::result_type(dummyCell2);
        }
        return typename Traits::template GetReductionPair<ArgT>::result_type();
    }


    size_t getMaxDim()
    {
        typename SComplex::Dim maxDim = 0;
        for (typename SComplex::ColoredIterators::Iterators::AllCells::iterator it = complex.template iterators<1>().allCells().begin(),
                end = complex.template iterators<1>().allCells().end();
                it != end; ++it)
        {
            maxDim = std::max(maxDim, (*it).getDim());
        }

        return maxDim;
    }

    typedef ::SComplex<SComplexDefaultTraits> OutputComplexType;

    OutputComplexType* getOutputComplex()
    {
        return outputComplex;
    }

protected:

    std::vector<int> morse;
    std::vector<AKQType> akq;
    std::vector<Cell> kerKing;
    std::vector<Cell> aces;

    int max_d;
    std::map<std::pair<int,int>, int> coeffs;
    std::map<std::pair<int,int>, int> numPathsBetween; // only between aces - small

    OutputComplexType *outputComplex;

    SComplex& complex;
    Cell dummyCell2, dummyCell3;
};

template<typename SComplexT>
class AKQReduceStrategy: public AKQReduceStrategyBase<SComplexT>
{
public:
    AKQReduceStrategy(SComplexT& _complex): AKQReduceStrategyBase<SComplexT>(_complex) {}
};

#endif

