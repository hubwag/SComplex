#ifndef SCOMPLEX_ALGS_DEFAULT_REDUCE_STRATEGY_HPP_
#define SCOMPLEX_ALGS_DEFAULT_REDUCE_STRATEGY_HPP_

#include <boost/optional.hpp>
#include <boost/foreach.hpp>
#include <utility>

#include "SComplex.hpp"
#include "SComplexDefaultTraits.hpp"
#include "SComplexAlgs.hpp"
//#include "CrHomS.hpp"

//#include <boost/tuples/tuple.hpp>

const bool verbose = false;

template<typename SComplexT>
class DefaultReduceStrategyTraits
{
public:

    // template<typename ImplT>
    // struct Proxy: pu	blic CellProxy<ImplT> {
    // 	 template<typename ImplT2>
    // 	 Proxy(const ImplT2& impl): CellProxy<ImplT>(impl) {}

    // 	 Proxy(const SComplexT& c): CellProxy<ImplT>(ImplT(c)) {}
    // };

    // template<typename ImplT>
    // struct Proxy<CellProxy<ImplT> >: public CellProxy<ImplT> {
    // 	 template<typename ImplT2>
    // 	 Proxy(const ImplT2& impl): CellProxy<ImplT>(impl) {}
    // };

    // template<typename ImplT>
    // static Proxy<ImplT*> makeProxy(const CellProxy<ImplT>& impl) {
    // 	 return Proxy<ImplT*>(impl.getImpl());
    // }

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
class DefaultReduceStrategyBase
{

public:
    typedef SComplexT SComplex;
    typedef DefaultReduceStrategyTraits<SComplex> Traits;
    typedef typename SComplex::Cell Cell;

    map<int, int> morse;
    map<Cell, char> akq;
    vector<Cell> aces;

    DefaultReduceStrategyBase(SComplex& _complex): complex(_complex), dummyCell2(_complex),  dummyCell3(_complex)
    {
    	cout << "CHECKING MAX DIM BRUTALLY! ADD getMaxDim to SCOMPLEX!!!!" << endl;
    	max_d = 0;
    	typename SComplex::Iterators::AllCells::iterator it = complex.iterators().allCells().begin();
    	typename SComplex::Iterators::AllCells::iterator end = complex.iterators().allCells().end();

    	for (; it != end; ++it)
			max_d = max<int>(max_d, it->getDim());

		/*
    	BOOST_FOREACH(typename SComplex::Iterators::AllCells::iterator::value_type v,
                      complex.iterators().allCells())
                      {
                      	max_d = max<int>(max_d, v.getDim());
                      }*/
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
    int calc_morse_value(const cellT1 &c)
    {
    	int val = 0;
    	BOOST_FOREACH(Cell el, complex.iterators().bdCells(c))
    	{
    		val = max(val, morse[el.getId()]);
    	}

    	return val + 1;
    }

    template<typename T1, typename T2>
    void king_gets_married(const T1 &king, const T2 &queen)
    {
    	akq[king] = 'k';
		akq[queen] = 'q';
        int v = calc_morse_value(king);
        morse[king.getId()] = v;
        morse[queen.getId()] = v;
        her_king.insert(make_pair(queen, king));
    }

    template<typename ImplT1, typename ImplT2>
    void coreduce(const typename SComplex::template CellProxy<ImplT1>& a, const typename SComplex::template CellProxy<ImplT2>& b)
    {
        BOOST_ASSERT(a.getDim() != b.getDim());

        if (a.getDim() > b.getDim()) {
            king_gets_married(a, b);
        } else {
            king_gets_married(b, a);
        }

        a.template setColor<2>();
        b.template setColor<2>();
    }

	int max_d;
    map<pair<int,int>, int> coeffs;
    map<Cell, Cell> her_king;
    map<pair<int,int>, int> num_paths_between;

    template<typename T, typename U>
    pair<T, U> make_sorted_pair(const T &a, const U &b)
    {
    	return make_pair(min(a,b), max(a,b));
    }

    int toAceCoeff(Cell x, Cell y)
    {
    	BOOST_ASSERT(akq[y] == 'a');
    	BOOST_ASSERT(akq[x] != 'q');

    	return complex.coincidenceIndex(x, y);
    }

	//
    int toKingCoeff(Cell x, Cell y, Cell y_star)
    {
    	BOOST_ASSERT(akq[x] != 'q');
    	BOOST_ASSERT(akq[y_star] == 'q');
    	BOOST_ASSERT(akq[y] == 'k');

    	return -1 * complex.coincidenceIndex(x, y_star) / complex.coincidenceIndex(y, y_star);
    }

    void follow_path(Cell c)
    {
    	stack<pair<Cell, int> > S; // no cycles!
    	stack<set<Cell> > Spath;
    	set<Cell> _s;
    	_s.insert(c);

    	S.push(make_pair(c, 1));

    	while(S.size())
    	{
    		Cell curr = S.top().first;
    		int accumulated_weight = S.top().second;
    		S.pop();

    		BOOST_ASSERT(akq[curr] != 'q');

    		if (curr.getId() != c.getId() && akq[curr] == 'a')
    		{
    			++num_paths_between[make_pair(c.getId(), curr.getId())];
    			if (verbose)
    			{
					cout << "found path from: " << c.getId() << " to " << curr.getId() << endl;
					cout << "between values: " << morse[c.getId()] << "and " << morse[curr.getId()] << " with coef product" << accumulated_weight << endl;
    			}
    			coeffs[make_pair(c.getId(), curr.getId())] += accumulated_weight;

    			continue;
    		}

    		int our_value = morse[curr.getId()];

			// case 1: to ace
    		BOOST_FOREACH(Cell to, complex.iterators().bdCells(curr))
    		{
    			if (akq[to] == 'a' && morse[to.getId()] < our_value)
    			{
					S.push(make_pair(to, accumulated_weight * toAceCoeff(curr, to)));
    			}
    		}

			// case 2: to king
    		BOOST_FOREACH(Cell bd, complex.iterators().bdCells(curr))
    		{
    			if (akq[bd] != 'q')
					continue;

    			// Cell to = her_king[bd];
    			Cell to = her_king.find(bd)->second;
    			BOOST_ASSERT(akq[to] == 'k');

				if (morse[to.getId()] < our_value)
    			{
    				{
						S.push(make_pair(to, accumulated_weight * toKingCoeff(curr, to, bd)));
    				}
    			}
    		}
    	}
    }

    typedef ::SComplex<SComplexDefaultTraits> OutputComplexType;
    OutputComplexType *outputComplex;

    void gender_reassignment(Cell &king, Cell &queen)
    {

    }

    void get_path(Cell &c, vector<Cell> &ordered_path_ret)
    {
    	stack<pair<Cell, int> > S; // no cycles!
    	stack<set<Cell> > Spath;
    	stack<vector<Cell> > Sordered;

    	vector<Cell> ordered_path;

    	set<Cell> _s;
    	_s.insert(c);
    	ordered_path.push_back(c);

    	Spath.push(_s);
    	Sordered.push(ordered_path);

    	S.push(make_pair(c, 1));

    	while(S.size())
    	{
    		Cell curr = S.top().first;
    		int accumulated_weight = S.top().second;
    		S.pop();

    		set<Cell> path = Spath.top();
    		vector<Cell> ordered_path = Sordered.top();
    		Spath.pop();
    		Sordered.pop();

    		BOOST_ASSERT(akq[curr] != 'q');

    		if (curr.getId() != c.getId() && akq[curr] == 'a')
    		{
    			if (verbose)
    			{
					cout << "found path from: " << c.getId() << " to " << curr.getId() << endl;
					cout << "between values: " << morse[c.getId()] << "and " << morse[curr.getId()] << " with coef product" << accumulated_weight << endl;
    			}


    			ordered_path_ret = ordered_path;

				break;
    		}

    		int our_value = morse[curr.getId()];

			// case 1: to ace
    		BOOST_FOREACH(Cell to, complex.iterators().bdCells(curr))
    		{
    			if (akq[to] == 'a')
    			{
    				BOOST_ASSERT(morse[to.getId()] < our_value);

    				if (path.count(to) == 0)
    				{
    					ordered_path.push_back(to);
    					path.insert(to);
						Spath.push(path);
						Sordered.push(ordered_path);
						S.push(make_pair(to, accumulated_weight * toAceCoeff(curr, to)));
    				}
    			}
    		}

			// case 2: to king
    		BOOST_FOREACH(Cell bd, complex.iterators().bdCells(curr))
    		{
    			if (akq[bd] != 'q')
					continue;

    			// Cell to = her_king[bd];
    			Cell to = her_king.find(bd)->second;
    			BOOST_ASSERT(akq[to] == 'k');

				if (path.count(to) == 0)
				{
					BOOST_ASSERT(morse[to.getId()] < our_value);
					path.insert(to);
					ordered_path.push_back(bd);
					ordered_path.push_back(to);
					Spath.push(path);
					Sordered.push(ordered_path);
					S.push(make_pair(to, accumulated_weight * toKingCoeff(curr, to, bd)));
				}

    		}
    	}
    }

    void report_paths()
    {
    	cout << " asow bylo: " << aces.size() << endl;

    	if (verbose)
    	for (size_t i = 0; i < aces.size(); i++)
    	{
			cout << aces[i].getId() << " ";
    	}
    	cout << endl;

    	BOOST_FOREACH(Cell ace, aces)
    	{
    		follow_path(ace);
    	}

    	typedef pair<pair<int,int>,int> Pair;

    	cout << "\n\n\n";

		/*
    	BOOST_FOREACH(Pair p, num_paths_between)
		{
			cout << "mamy jedyna sciezke pomiedzy : ";
			cout << p.first.first << " " << p.first.second << " = " << p.second << endl;


			vector<Cell> path;

			for (int i = 0; i < aces.size(); i++)
			{
				if (aces[i].getId() == p.first.first)
				{
					get_path(aces[i], path);
					break;
				}
			}

			cout << "sciezka: " << endl;
			for (int i = 0; i < path.size(); i++)
			{
				cout << path[i].getId() << " " << akq[path[i]] << " " << "[wart: " << morse[path[i].getId()] << "] ";
			}

			akq[path.back()] = 'q';
			akq[path[0]] = 'k';

			for (int i = (int)path.size() - 1; i >= 1; i-=2)
				her_king[path[i]] = path[i-1];
		}

		if (verbose)
		{
			cout << "num paths between: \n";
			BOOST_FOREACH(Pair p, num_paths_between)
			{
				cout << p.first.first << " " << p.first.second << " = " << p.second << endl;
			}

			cout << "coefficients: \n";
			BOOST_FOREACH(Pair p, coeffs)
			{
				cout << p.first.first << " " << p.first.second << " => " << p.second << endl;
			}
		}

		*/

	    vector<size_t> dims(aces.size());

	    map<size_t, size_t> from0;

	    BOOST_FOREACH(Cell ace, aces)
	    {
	    	size_t next = from0.size();
	    	dims[next] = ace.getDim();
	    	from0[ace.getId()] = next;
	    }

		OutputComplexType::KappaMap kap;

		time_t t = time(0);
		srand(t);

		cout << "constructing general SComplex" << endl;

		BOOST_FOREACH(Pair p, coeffs)
    	{
    		int coef = p.second;
    		// coef = rand()%2 == 0 ? 1 : -1;

    		kap.push_back(boost::make_tuple(from0[p.first.first], from0[p.first.second], coef));

    		if (verbose)
				cout << from0[p.first.first] << "[d=" <<dims[from0[p.first.first]] << "]" << " : " <<  from0[p.first.second] << "[d="<<dims[from0[p.first.second]] << "]"  << " => " << coef << endl;
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
					v = calc_morse_value(*it);

				morse[it->getId()] = v;
				aces.push_back(*it);
				akq[*it] = 'a';

				// cout << "AS!" << it->getId() << ' ' << it->getDim() << "with M = " << v << endl;

                return typename Traits::Extract::result_type::value_type(*it);
            }
        }
        cout << "no i tyle bylo asow..." << endl;
        report_paths();
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
        BOOST_FOREACH(typename SComplex::ColoredIterators::Iterators::BdCells::iterator::value_type v,
                      complex.iterators(1).bdCells(cell))
        {
            if (times == 0)
            {
                dummyCell3 = v;
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
        BOOST_FOREACH(typename SComplex::ColoredIterators::Iterators::CbdCells::iterator::value_type v,
                      complex.iterators(1).cbdCells(cell))
        {
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

protected:
    SComplex& complex;
    Cell dummyCell2, dummyCell3;
};

template<typename SComplexT>
class DefaultReduceStrategy: public DefaultReduceStrategyBase<SComplexT>
{

public:
    DefaultReduceStrategy(SComplexT& _complex): DefaultReduceStrategyBase<SComplexT>(_complex) {}
};


#endif

