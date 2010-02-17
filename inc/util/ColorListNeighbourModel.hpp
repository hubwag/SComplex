#ifndef __COLOR_LIST_NEIGHBOUR_MODEL_HPP_
#define __COLOR_LIST_NEIGHBOUR_MODEL_HPP_

namespace Util {
  namespace Neighbours {
	 
	 template<typename ObjectT>
	 class ColorListNeighbourModel {
		typedef ObjectT Object;
		typedef boost::reference_wrapper<Object> ObjectRef;
		typedef typename std::list<ObjectRef> ObjectsList;
		typedef typename std::vector<ObjectsList> ObjectsByColor;  
		typedef typename std::vector<ObjectRef> ObjectsVect;
  
	 public:

		typedef typename ObjectsList::iterator NeighbourId;
		typedef typename ObjectsVect::iterator iterator;

  
		explicit ColorListNeighbourModel(size_t colors): neighboursByColor(colors) {}

		iterator begin() {
		  return neighbours.begin();
		}

		iterator end() {
		  return neighbours.end();
		}

		NeighbourId add(Object& object, const typename Object::Color& color) {
		  return neighboursByColor[color].insert(neighboursByColor[color].end(), object);	 
		}

		void changeColor(NeighbourId it, const typename Object::Color& oldColor, const typename Object::Color& newColor) {
		  neighboursByColor[newColor].splice(neighboursByColor[newColor].end(), neighboursByColor[oldColor], it);
		}

	 private:
		ObjectsVect neighbours;
		ObjectsByColor neighboursByColor;
	 };
  }
}

#endif //  __COLOR_LIST_NEIGHBOUR_MODEL_HPP_
