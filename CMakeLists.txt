cmake_minimum_required(VERSION 2.6)

project(SComplex CXX)

#set(RUNTIME_OUTPUT_DIRECTORY build)

SET(Boost_ADDITIONAL_VERSIONS "1.42" "1.42.0")
find_package(Boost 1.42)


set(CODE ${PROJECT_SOURCE_DIR}/)
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/install)


#set(CMAKE_USE_RELATIVE_PATHS TRUE)
#set(CMAKE_SUPPRESS_REGENERATION TRUE)

enable_testing()

include_directories(${CODE}/inc)

include_directories(SYSTEM $ENV{CAPD_ROOT}/include)
link_directories($ENV{CAPD_ROOT}/lib)


include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
link_directories(SYSTEM ${Boost_INCLUDE_DIRS})


#add_library(/home/mj/workspace/capd/lib/capd STATIC IMPORTED)


add_subdirectory(${CODE}/src src)
add_subdirectory(${CODE}/test test)
add_subdirectory(${CODE}/lib lib)
add_subdirectory(${CODE}/inc inc)



include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE  
     "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set (CPACK_PACKAGE_VERSION_MAJOR 1)
set (CPACK_PACKAGE_VERSION_MINOR 0)

set(CPACK_SOURCE_IGNORE_FILES
"~$"
"^${PROJECT_SOURCE_DIR}/.#.*"
"^${PROJECT_SOURCE_DIR}/.git"
"^${PROJECT_SOURCE_DIR}/build/"
"^${PROJECT_SOURCE_DIR}/install/"
)

include (CPack)

