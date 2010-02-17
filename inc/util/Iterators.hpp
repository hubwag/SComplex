#ifndef __ITERATORS_HPP_
#define __ITERATORS_HPP_

#include <functional>

namespace Util {
  namespace Iterators {
	 
	 template<typename CollectionT, typename TransformF =
				 std::_Identity<typename CollectionT::iterator> > // a problem with std::identity, where it is?
	 class CollectionBeginEnd {
	 public:
		typedef typename TransformF::result_type iterator;
		
		explicit CollectionBeginEnd(CollectionT& _collection, const TransformF& _transform = TransformF()): collection(_collection), transform(_transform) {}
  
		iterator begin() {
		  return transform(collection.begin());
		}

		iterator end() {
		  return transform(collection.end());
		}

	 private:
		CollectionT& collection;
		TransformF transform;
	 };
  }
}

#endif
