#ifndef __SCOMPLEX_DEFUALT_TRAITS_HPP
#define __SCOMPLEX_DEFUALT_TRAITS_HPP

#include <cstring>
#include "util/ColoredObjectsListModel.hpp"


class SComplexDefaultTraits {

public:

  typedef size_t Size;
  typedef size_t Dim;
  typedef size_t Id;
  typedef size_t Color;

  template<typename ObjectT>
  class ColoredObjectsModel: public Util::Colors::ColoredObjectsListModel<ObjectT, Color> {

  };
  
};

#endif

