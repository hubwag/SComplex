#ifndef __ITERATORS_HPP_
#define __ITERATORS_HPP_

#include <functional>
#include <boost/iterator/transform_iterator.hpp>

namespace Util {
  namespace Iterators {
	 
	 template<bool isConst, typename CollectionT, typename TransformT>
	 class CollectionBeginEnd {
	 public:
		typedef typename boost::mpl::if_c<isConst,
													 const CollectionT&,
													 CollectionT&
													 >::type Collection;
		

		typedef boost::transform_iterator<TransformT, typename boost::mpl::if_c<isConst, typename CollectionT::const_iterator, typename CollectionT::iterator>::type> iterator;
		
		typedef boost::transform_iterator<TransformT, typename CollectionT::const_iterator> const_iterator;
		
		explicit CollectionBeginEnd(Collection _collection): collection(_collection) {}

		iterator begin() const {
		  return iterator(collection.begin(), transform);
		}

		iterator end() const {
		  return iterator(collection.end(), transform);
		}

	 private:
		Collection collection;
		TransformT transform;
	 };

	 template<bool isConst, typename CollectionT, typename TransformT>
	 class BeginEnd {
		typedef typename boost::mpl::if_c<isConst, typename CollectionT::const_iterator, typename CollectionT::iterator>::type CollectionIterator;
		typedef typename boost::mpl::if_c<isConst,
													 const CollectionT&,
													 CollectionT&
													 >::type Collection;

	 public:
		typedef boost::transform_iterator<TransformT, CollectionIterator> iterator;
		
		typedef boost::transform_iterator<TransformT, typename CollectionT::const_iterator> const_iterator;
		
		explicit BeginEnd(Collection _collection, TransformT _transform): b(_collection.begin()), e(_collection.end()), transform(_transform) {}
		
		iterator begin() const {
		  return iterator(b, transform);
		}

		iterator end() const {
		  return iterator(e, transform);
		}

	 private:
		CollectionIterator b, e;
		TransformT transform;
	 };
	 
  }
}

#endif
