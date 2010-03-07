#ifndef __COLOR_LIST_NEIGHBOUR_MODEL_HPP_
#define __COLOR_LIST_NEIGHBOUR_MODEL_HPP_

#include "Iterators.hpp"
#include "ColorListModel.hpp"
#include <boost/range.hpp>

namespace Util {
  namespace Neighbours {


	 template<typename Object, typename Color>
	 class ColorListNeighbourModel {
		
	 public:
		struct NeighbourLink;
		typedef ColorListModel<NeighbourLink, Color> ColoredNeighbours;

		typedef typename ColoredNeighbours::ObjectPtrsIterator NeighbourLinkPtrsIterator;
		
		struct NeighbourLink {
		  Object objectRef;
		  NeighbourLinkPtrsIterator neighbourLinkPtrsIterator; //an iterator to NeighbourLinkPtrs not in this, but in an instance of the class for a neighbour.

		  explicit NeighbourLink(Object _o) : objectRef(_o) {}
		};

		typedef typename ColoredNeighbours::AllObjects AllNeighbours;
	  
		typedef typename ColoredNeighbours::ObjectsInColor NeighboursInColor;
		
		AllNeighbours allNeighbours() {
		  return coloredNeighbours.allObjects();
		}

		NeighboursInColor neighboursInColor(const Color& color) {
		  return coloredNeighbours.objectsInColor(color);
		}
		
		void init(size_t colors, size_t size) {
		  coloredNeighbours.init(colors, size);
		}
					
		NeighbourLinkPtrsIterator add(const Object& object, const Color& color) {
		  return coloredNeighbours.add(NeighbourLink(object), color);
		}

		void changeColor(NeighbourLinkPtrsIterator it, const Color& oldColor, const Color& newColor) {
		  coloredNeighbours.changeColor(it, oldColor, newColor);
		}

	 private:
		ColoredNeighbours coloredNeighbours;
	 };
  }
}

#endif //  __COLOR_LIST_NEIGHBOUR_MODEL_HPP_
