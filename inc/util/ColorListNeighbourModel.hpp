#ifndef __COLOR_LIST_NEIGHBOUR_MODEL_HPP_
#define __COLOR_LIST_NEIGHBOUR_MODEL_HPP_

namespace Util {
  namespace Neighbours {
	 
	 template<typename ObjectT, typename Color>
	 class ColorListNeighbourModel {

	 public:
		
		typedef ObjectT Object;

		struct NeighbourLink;
		
		typedef typename std::list<NeighbourLink*> NeighbourLinkPtrs;
		typedef typename std::vector<NeighbourLinkPtrs> NeighbourLinkPtrsByColor;  
		typedef typename NeighbourLinkPtrsByColor::value_type::iterator NeighbourLinkPtrsIterator;
		
		struct NeighbourLink {
		  Object objectRef;
		  NeighbourLinkPtrsIterator neighbourLinkPtrsIterator; //an iterator to NeighbourLinkPtrs not in this, but in an instance of the class for a netighbour.

		  explicit NeighbourLink(Object _o) : objectRef(_o) {}
		};
		
		typedef typename std::vector<NeighbourLink> NeighbourLinks;		
		typedef typename NeighbourLinks::iterator iterator;
		typedef typename NeighbourLinks::const_iterator const_iterator;

  
		void init(size_t colors, size_t size) {
		  neighboursByColor.resize(colors);
		  neighbours.reserve(size);
		}

		iterator begin() {
		  return neighbours.begin();
		}

		iterator end() {
		  return neighbours.end();
		}

		const_iterator begin() const {
		  return neighbours.begin();
		}

		const_iterator end() const {
		  return neighbours.end();
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
