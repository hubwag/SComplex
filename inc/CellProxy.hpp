#ifndef CELL_PROXY_HPP
#define CELL_PROXY_HPP

#include <boost/ref.hpp>

template<typename CellImplT>
class CellProxy {
protected:
  CellImplT impl;
  //  boost::reference_wrapper<CellImplT> nonConstImpl;
  
public:
  typedef typename CellImplT::Color Color;
  typedef typename CellImplT::Dim Dim;

  typedef CellImplT Impl;
  
  CellProxy(const CellImplT& _impl): impl(_impl) {} //, nonConstImpl(impl) {}
  
  Color getColor() const{
	 return impl.getColor();
  }

  template<Color color>
  void setColor() {
	 impl.template setColor<color>();
  }

  void setColor(const Color& color) {
	 impl.setColor(color);
  }
	 
  Dim getDim() const {
	 return impl.getDim();
  }

  bool operator<(const CellProxy& b) const {
	 return impl < b.impl;
  }

  CellImplT* getImpl() const {
	 //return nonConstImpl.get_pointer();
	 return &(const_cast<CellProxy*>(this)->impl);
  }
};

template<typename CellImplT>
class CellProxy<CellImplT*> {
protected:
  CellImplT* impl;
  //  boost::reference_wrapper<CellImplT> nonConstImpl;
public:
  typedef typename CellImplT::Color Color;
  typedef typename CellImplT::Dim Dim;

  typedef CellImplT Impl;

  CellProxy(): impl(NULL) {} //, nonConstImpl(*impl) {}
  CellProxy(CellImplT* _impl): impl(_impl) {} //, nonConstImpl(*impl) {}

  CellProxy(const CellProxy<CellImplT>& other): impl(other.getImpl()) {}
																//nonConstImpl(*impl) {}

  Color getColor() const{
	 return impl->getColor();
  }

  template<Color color>
  void setColor() {
	 impl->template setColor<color>();
  }

  void setColor(const Color& color) {
	 impl->setColor(color);
  }
	 
  Dim getDim() const {
	 return impl->getDim();
  }

  bool operator<(const CellProxy& b) const {
	 return *impl < b->impl;
  }

  CellImplT* getImpl() const {
	 //return nonConstImpl.get_pointer();
	 return const_cast<CellProxy*>(this)->impl;
  }

};

#endif
