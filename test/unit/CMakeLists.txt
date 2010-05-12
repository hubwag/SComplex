
#SimplexSComplexTest.cpp

set(TEST_OUTPUT_FORMAT HRF)
set(TEST_LOG_LEVEL message)
set(TEST_REPORT_LEVEL short)
set(TEST_RUN *)

set(ARGS "${ARGS} -VV")

#add_definitions(-DBOOST_TEST_DYN_LINK)

add_executable(SComplexUnitTest
  SComplexUnitTestMain.cpp
  SComplexInitTest.cpp
  SComplexColorTest.cpp
  SComplexReductionTest.cpp
#  SimplexSComplexTest.cpp
  CubSComplexIteratorsTest.cpp
  CubSComplexReductionTest.cpp
)

add_executable(COAKQTest  COAKQTest.cpp   SComplexUnitTestMain.cpp)

target_link_libraries(SComplexUnitTest capd)
target_link_libraries(SComplexUnitTest ${Boost_LIBRARIES})

target_link_libraries(COAKQTest capd)
target_link_libraries(COAKQTest ${Boost_LIBRARIES})


if (WIN32)
   set (LD_VARNAME "PATH")
   set (LD_PATH "${Boost_LIBRARY_DIRS};$ENV{PATH}")

   #
   # IMPORTANT NOTE: The set_tests_properties(), below, internally
   # stores its name/value pairs with a semicolon delimiter.
   # because of this we must protect the semicolons in the path
   #
   string(REPLACE ";" "\\;" LD_PATH "${LD_PATH}")
else ()
   set (LD_VARNAME "LD_LIBRARY_PATH")
   set (LD_PATH "${Boost_LIBRARY_DIRS}:$ENV{LD_LIBRARY_PATH}")
endif ()

add_test(SComplexUnitTest SComplexUnitTest --output_format=${TEST_OUTPUT_FORMAT} --log_level=${TEST_LOG_LEVEL} --report_level=${TEST_REPORT_LEVEL} --run_test=${TEST_RUN})
set_tests_properties(SComplexUnitTest PROPERTIES ENVIRONMENT "${LD_VARNAME}=${LD_PATH}")

add_test(COAKQTest COAKQTest --output_format=${TEST_OUTPUT_FORMAT} --log_level=${TEST_LOG_LEVEL} --report_level=${TEST_REPORT_LEVEL} --run_test=${TEST_RUN})

set_tests_properties(COAKQTest PROPERTIES ENVIRONMENT "${LD_VARNAME}=${LD_PATH}")


