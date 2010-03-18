#include "CubSComplex.hpp"
#include <functional>
#include <boost/utility/enable_if.hpp>


template<>
class DefaultReduceStrategyTraits<CubSComplex> {
public:

  template<typename ImplT>
  struct Proxy: public CubSComplex::CubCellProxy<ImplT> {
	 template<typename ImplT2>
	 Proxy(const ImplT2& impl): CubSComplex::CubCellProxy<ImplT>(impl) {}
  };

  template<typename ImplT>
  struct Proxy<CubSComplex::CubCellProxy<ImplT> >: public CubSComplex::CubCellProxy<ImplT> {
	 template<typename ImplT2>
	 Proxy(const ImplT2& impl): CubSComplex::CubCellProxy<ImplT>(impl) {}
  };

  template<typename ImplT>
  static Proxy<ImplT*> makeProxy(const CubSComplex::CubCellProxy<ImplT>& impl) {
	 return Proxy<ImplT*>(impl.getImpl());
  }
  
  template<typename ArgT>
  struct GetReductionPair:  std::unary_function<const ArgT&,
																boost::optional<Proxy<CubSComplex::DynamicCell::Impl*> > > {};
  template<typename ArgT>
  struct GetCoreductionPair:  std::unary_function<const ArgT&,
																  boost::optional<Proxy<CubSComplex::DynamicCell::Impl*> > > {};
  struct ForceCoreduction {
	 typedef boost::optional<std::pair<Proxy<CubSComplex::DynamicCell::Impl*>,
												  Proxy<CubSComplex::DynamicCell::Impl*> > > result_type;
  };

  struct Extract {
	 typedef boost::optional<Proxy<CubSComplex::BitCoordCellImpl> >  result_type;
  };

};
