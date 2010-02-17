#ifndef __ITERATORS_HPP_
#define __ITERATORS_HPP_

namespace Util {
  namespace Iterators {
	 
	 template<typename CollectionT>
	 class CollectionBeginEnd {
	 public:
		explicit CollectionBeginEnd(CollectionT& _collection): collection(_collection) {}
  
		typename CollectionT::iterator begin() {
		  return collection.begin();
		}

		typename CollectionT::iterator end() {
		  return collection.end();
		}

	 private:
		CollectionT& collection;
	 };
  }
}

#endif
