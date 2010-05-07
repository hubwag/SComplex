#ifndef SCOMPLEX_ALGS_COAKQ_STRATEGY_HPP_
#define SCOMPLEX_ALGS_COAKQ_STRATEGY_HPP_

#include <boost/tuple/tuple.hpp>
#include <boost/optional.hpp>
#include <boost/foreach.hpp>
#include <utility>
#include <map>
#include <vector>
#include <stack>
#include <algorithm>


#include <iostream>


template<typename SComplexT>
class COAKQStrategyTraits {
public:
  
  template<typename ArgT>
  struct GetReductionPair: public std::unary_function<const ArgT&,
						      boost::optional<typename SComplexT::Cell> > {};

  struct ForceReduction {
    typedef boost::optional<std::pair<typename SComplexT::Cell,
				      typename SComplexT::Cell> > result_type;
  };

  struct Extract {
    typedef boost::optional<typename SComplexT::Cell >  result_type;
  };
};

template<typename SComplexT, typename OutputComplexT>
class COAKQStrategyBase {

public:
  typedef SComplexT SComplex;
  typedef COAKQStrategyTraits<SComplex> Traits;
  typedef typename SComplex::Cell Cell;


  COAKQStrategyBase(SComplex& _complex): complex(_complex), dummyCell2(_complex),  dummyCell3(_complex) {
    max_d = 0;
    BOOST_FOREACH(typename SComplex::Iterators::AllCells::iterator::value_type v,
		  complex.iterators().allCells())
      {
	max_d = std::max<int>(max_d, v.getDim());
      }
}
  
  SComplex& getComplex() const {
	 return complex;
  }

  OutputComplexT& getOutputComplex() {
    return *outputComplex;
  }
  
  template<typename ImplT>
  static bool reduced(const typename SComplex::template CellProxy<ImplT>& cell) {
  	 return cell.getColor() == 2;
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
        // her_king[queen] = king;
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
  template<typename ImplT1, typename ImplT2>	 
  void reduce(const typename SComplex::template CellProxy<ImplT1>& a, const typename SComplex::template CellProxy<ImplT2>& b)  {
        BOOST_ASSERT(a.getDim() != b.getDim());

        if (a.getDim() > b.getDim()) {
            king_gets_married(a, b);
        } else {
            king_gets_married(b, a);
        }

  	 a.template setColor<2>();
  	 b.template setColor<2>();
  }

  template<typename ImplT>
  static void reduce(const typename SComplex::template CellProxy<ImplT>& cell) {
	 cell.template setColor<2>();
  }
  
  typename Traits::Extract::result_type extract() {
      for (int d = max_d; d >= 0; d--)
        {
	  typename SComplex::ColoredIterators::Iterators::DimCells dimCells = complex.iterators(1).dimCells(d);
	  typename SComplex::ColoredIterators::Iterators::DimCells::iterator end = dimCells.end(),
	    it = dimCells.begin();

	  if (it != end)
            {
	      int v = 0;

	      if (d != max_d)
		v = calc_morse_value(*it);

	      morse[it->getId()] = v;
	      aces.push_back(*it);
	      akq[*it] = 'a';

	      return typename Traits::Extract::result_type::value_type(*it);
            }
        }

      report_paths();
      return typename Traits::Extract::result_type();
  }
  
  int toAceCoeff(Cell x, Cell y)
  {
    BOOST_ASSERT(akq[y] == 'a');
    BOOST_ASSERT(akq[x] != 'q');

    return complex.coincidenceIndex(x, y);
  }

  int toKingCoeff(Cell x, Cell y, Cell y_star)
  {
    BOOST_ASSERT(akq[x] != 'q');
    BOOST_ASSERT(akq[y_star] == 'q');
    BOOST_ASSERT(akq[y] == 'k');

    return -1 * complex.coincidenceIndex(x, y_star) / complex.coincidenceIndex(y, y_star);
  }

  void follow_path(Cell c)
  {
    std::stack<std::pair<Cell, int> > S; // no cycles!
    S.push(std::make_pair(c, 1));

    while(S.size())
      {
	Cell curr = S.top().first;
	int accumulated_weight = S.top().second;
	S.pop();

	BOOST_ASSERT(akq[curr] != 'q');

	if (curr.getId() != c.getId() && akq[curr] == 'a')
	  {
	    ++num_paths_between[make_pair(c.getId(), curr.getId())];

	    std::cout << "found path from: " << c.getId() << " to " << curr.getId() << std::endl;
	    std::cout << "between values: " << morse[c.getId()] << "and " << morse[curr.getId()] << " with coef product" << accumulated_weight << std::endl;
	      
	    coeffs[make_pair(c.getId(), curr.getId())] += accumulated_weight;
	    continue;
	  }

	int our_value = morse[curr.getId()];

	// case 1: to ace
	BOOST_FOREACH(Cell to, complex.iterators().cbdCells(curr))
	  {
	    if (akq[to] == 'a' && morse[to.getId()] < our_value)
	      {
		S.push(make_pair(to, accumulated_weight * toAceCoeff(curr, to)));
	      }
	  }

	// case 2: to king
	BOOST_FOREACH(Cell bd, complex.iterators().cbdCells(curr))
	  {
	    if (akq[bd] != 'q')
	      continue;

	    // Cell to = her_king[bd];
	    Cell to = her_king.find(bd)->second;
	    BOOST_ASSERT(akq[to] == 'k');
	    if (morse[to.getId()] < our_value)
	      {
		S.push(make_pair(to, accumulated_weight * toKingCoeff(curr, to, bd)));
	      }
	  }
      }
  }

  void report_paths()
  {
    std::cout << " asow bylo: " << aces.size() << std::endl;


    for (size_t i = 0; i < aces.size(); i++)
      {
	std::cout << aces[i].getId() << " ";
      }
    std::cout << std::endl;


    BOOST_FOREACH(Cell ace, aces)
      {
	follow_path(ace);
      }

    typedef std::pair<std::pair<int,int>,int> Pair;

    std::cout << "\n\n\n";


    std::cout << "num paths between: \n";
    BOOST_FOREACH(Pair p, num_paths_between)
      {
	std::cout << p.first.first << " " << p.first.second << " = " << p.second << std::endl;
      }

    std::cout << "coefficients: \n";
    BOOST_FOREACH(Pair p, coeffs)
      {
	std::cout << p.first.first << " " << p.first.second << " => " << p.second << std::endl;
      }
      

    std::vector<size_t> dims(aces.size());

    std::map<size_t, size_t> from0;

    BOOST_FOREACH(Cell ace, aces)
      {
	size_t next = from0.size();
	dims[next] = ace.getDim();
	from0[ace.getId()] = next;
      }

    typename OutputComplexT::KappaMap kap;

    time_t t = time(0);
    srand(t);

    std::cout << "constructing general SComplex" << std::endl;

    BOOST_FOREACH(Pair p, coeffs)
      {
	int coef = p.second;
	// coef = rand()%2 == 0 ? 1 : -1;

	kap.push_back(boost::make_tuple(from0[p.first.first], from0[p.first.second], coef));

	
	std::cout << from0[p.first.first] << "[d=" <<dims[from0[p.first.first]] << "]" << " : " << from0[p.first.second] << "[d="<<dims[from0[p.first.second]] << "]" << " => " << coef << std::endl;
      
      }

    outputComplex = new OutputComplexT(3, dims, kap, 1);
  }

     typename Traits::ForceReduction::result_type forceReductionPair() {
	  for (int d = max_d; d >= 0; d--)
	  {
	       typedef typename SComplex::ColoredIterators::Iterators::DimCells DimCells;
	       DimCells  dimCells = complex.iterators(1).dimCells(d);
	       typename DimCells::iterator end = dimCells.end(),
		    it = dimCells.begin();

	       while(it != end) {
		    typename Traits::template GetReductionPair<typename DimCells::iterator::value_type>::result_type
			 coface = getReductionPair(*it);
		    if (coface) {
			 return typename Traits::ForceReduction::result_type(std::make_pair(*it, *coface));
		    }
		    ++it;
	       }
	  }
  
	  return typename Traits::ForceReduction::result_type();
     }

  template<typename ArgT>
  typename Traits::template GetReductionPair<ArgT>::result_type
  getReductionPair(const ArgT& cell)
  {
  	 int times = 0;
  	 BOOST_FOREACH(typename SComplex::ColoredIterators::Iterators::CbdCells::iterator::value_type v,
  						complex.iterators(1).cbdCells(cell)) {
  		if (times == 0) {
  		  dummyCell3 = v;
  		}
  		++times;
  		if (times == 2) {
		  break;
  		}
  	 }

  	 if (times == 1) {
  		return typename Traits::template GetReductionPair<ArgT>::result_type(dummyCell3);
  	 }
  	 return typename Traits::template GetReductionPair<ArgT>::result_type();
  }


  size_t getMaxDim() {
	 typename SComplex::Dim maxDim = 0;
	 for (typename SComplex::ColoredIterators::Iterators::AllCells::iterator it = complex.template iterators<1>().allCells().begin(),
			  end = complex.template iterators<1>().allCells().end();
			it != end; ++it) {

		maxDim = std::max(maxDim, (*it).getDim());
	 }
	 
	 return maxDim;
  }

protected:
  SComplex& complex;
  Cell dummyCell2, dummyCell3;

  OutputComplexT *outputComplex;
  std::map<int, int> morse;
  std::map<Cell, char> akq;
  std::vector<Cell> aces;
  std::map<std::pair<int,int>, int> coeffs;
  std::map<Cell, Cell> her_king;
  std::map<std::pair<int,int>, int> num_paths_between;
  int max_d;
  
};

template<typename SComplexT, typename OutputComplexT>
class COAKQStrategy: public COAKQStrategyBase<SComplexT, OutputComplexT> {

public:
  COAKQStrategy(SComplexT& _complex): COAKQStrategyBase<SComplexT, OutputComplexT>(_complex) {}
};

#endif

