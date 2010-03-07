#ifndef __COLOR_LIST_NEIGHBOUR_MODEL_HPP_
#define __COLOR_LIST_NEIGHBOUR_MODEL_HPP_

#include "Iterators.hpp"
#include "ColorListModel.hpp"
#include <boost/range.hpp>

namespace Util {
  namespace Neighbours {


	 template<typename Object, typename Color>
	 class ColorListNeighbourModel {

		template<typename T>
		struct Dereference {
		  T& operator()(T* t) const {
			 return *t;
		  }
		};
		
	 public:
		struct NeighbourLink;
		typedef typename std::list<NeighbourLink*> NeighbourLinkPtrs;
		typedef typename std::vector<NeighbourLinkPtrs> NeighbourLinkPtrsByColor;  

		typedef typename NeighbourLinkPtrsByColor::value_type::iterator NeighbourLinkPtrsIterator;		
		struct NeighbourLink {
		  Object objectRef;
		  NeighbourLinkPtrsIterator neighbourLinkPtrsIterator; //an iterator to NeighbourLinkPtrs not in this, but in an instance of the class for a neighbour.

		  explicit NeighbourLink(Object _o) : objectRef(_o) {}
		};

		typedef typename std::vector<NeighbourLink> NeighbourLinks;		

		typedef typename boost::sub_range<NeighbourLinks> AllNeighbours;
		typedef Util::Iterators::RangeTransform<NeighbourLinkPtrs, Dereference<NeighbourLink> > NeighboursInColor;
		
		AllNeighbours allNeighbours() {
		  return AllNeighbours(neighbours);
		}

		NeighboursInColor neighboursInColor(const Color& color) {
		  return NeighboursInColor(neighboursByColor[color], Dereference<NeighbourLink>());
		}
		
		void init(size_t colors, size_t size) {
		  neighboursByColor.resize(colors);
		  neighbours.reserve(size);
		}
					
		NeighbourLinkPtrsIterator add(const Object& object, const Color& color) {
		  neighbours.push_back(NeighbourLink(object));
		  return neighboursByColor[color].insert(neighboursByColor[color].end(), &(*(neighbours.end() - 1)));	 
		}

		void changeColor(NeighbourLinkPtrsIterator it, const Color& oldColor, const Color& newColor) {
		  neighboursByColor[newColor].splice(neighboursByColor[newColor].end(), neighboursByColor[oldColor], it);
		}

	 private:
		NeighbourLinks neighbours;
		NeighbourLinkPtrsByColor neighboursByColor;
	 };
  }
}

#endif //  __COLOR_LIST_NEIGHBOUR_MODEL_HPP_
