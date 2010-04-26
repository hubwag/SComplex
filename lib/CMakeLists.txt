
add_library(${PROJECT_NAME} STATIC 
CubSComplex_Cell.cpp
SimplexCell.cpp
)

install(TARGETS ${PROJECT_NAME}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)


