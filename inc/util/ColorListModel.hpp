#ifndef __COLOR_LIST_MODEL_HPP_
#define __COLOR_LIST_MODEL_HPP_

#include <vector>
#include <list>
#include <boost/range.hpp>

namespace Util {
  namespace Neighbours {


	 template<typename Object, typename Color>
	 class ColorListModel {

		template<typename T>
		struct Dereference {
		  T& operator()(T* t) const {
			 return *t;
		  }
		};

	 public:
		typedef std::vector<Object> Objects;	 
		typedef std::list<Object*> ObjectPtrs;
		typedef std::vector<ObjectPtrs> ObjectPtrsByColor;
		typedef typename ObjectPtrs::iterator ObjectPtrsIterator;

		typedef typename boost::sub_range<Objects> AllObjects;
		typedef Util::Iterators::RangeTransform<ObjectPtrs&, Dereference<Object> > ObjectsInColor;

		AllObjects allObjects() {
		  return AllObjects(objects);
		}

		ObjectsInColor objectsInColor(const Color& color) {
		  return ObjectsInColor(objectsInColor[color], Dereference<Object>());
		}
		
		void init(size_t colors, size_t size) {
		  objectsByColor.resize(colors);
		  objects.reserve(size);
		}
					
		ObjectPtrsIterator add(const Object& object, const Color& color) {
		  objects.push_back(object);
		  return objectsByColor[color].insert(objectsByColor[color].end(), &(*(objects.end() - 1)));	 
		}

		void changeColor(ObjectPtrsIterator it, const Color& oldColor, const Color& newColor) {
		  objectsByColor[newColor].splice(objectsByColor[newColor].end(), objectsByColor[oldColor], it);
		}

	 private:
		Objects objects;
		ObjectPtrsByColor objectsByColor;
	 };
  
  }
}

#endif
