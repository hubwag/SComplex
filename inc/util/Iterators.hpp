#ifndef __ITERATORS_HPP_
#define __ITERATORS_HPP_

#include <functional>

namespace Util {
  namespace Iterators {
	 
	 template<typename CollectionT>
	 class CollectionBeginEnd {
	 public:
		typedef typename CollectionT::iterator iterator;
		typedef typename CollectionT::const_iterator const_iterator;
		
		explicit CollectionBeginEnd(CollectionT& _collection): collection(_collection) {}
  
		iterator begin() {
		  return (collection.begin());
		}

		const_iterator begin() const {
		  return (collection.begin());
		}

		iterator end() {
		  return (collection.end());
		}

		const_iterator end() const {
		  return (collection.end());
		}

	 private:
		CollectionT& collection;
	 };
  }
}

#endif
