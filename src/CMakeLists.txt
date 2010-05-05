
add_executable(CrHomS CrHomS.cpp)
target_link_libraries(CrHomS capd)

install(TARGETS CrHomS
  RUNTIME DESTINATION bin)
