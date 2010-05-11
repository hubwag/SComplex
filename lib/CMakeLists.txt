
add_library(${PROJECT_NAME} STATIC 
complex/cubical/CubSComplex_Cell.cpp
complex/simplicial/SimplexCell.cpp
)

install(TARGETS ${PROJECT_NAME}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)


