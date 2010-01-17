
BUILD_DIR=build
BINS_DIR=$(BUILD_DIR)/bin
LIBS_DIR=$(BUILD_DIR)/lib
OBJS_DIR=$(BUILD_DIR)/obj
RUN_DIR=$(BUILD_DIR)/run

INCS_DIR=inc
SRCS_DIR=src
TESTS_DIR=test


APP_NAME=CrHomS
LIB_NAME=libSComplex.a
TEST_RESULT_XML=$(RUN_DIR)/test_result.xml

TEST_OUTPUT_FORMAT=HRF
TEST_LOG_LEVEL=message
TEST_REPORT_LEVEL=short


BOOST_HOME=/home/juda/local/apps/boost-1-41-0
CAPD_HOME=/home/juda/workspace/capd

LOCAL_INC_PATHS= -I$(INCS_DIR)  
INC_PATHS = $(LOCAL_INC_PATHS) -isystem$(CAPD_HOME)/include -isystem$(BOOST_HOME)/include
LIB_PATHS = -L$(LIBS_DIR) -L$(CAPD_HOME)/lib -L$(BOOST_HOME)/lib


LIB_SRCS = $(SRCS_DIR)/CubSComplex_Cell.cpp
APP_SRCS = $(SRCS_DIR)/CrHomS.cpp

LIB_OBJS = $(LIB_SRCS:%.cpp=$(OBJS_DIR)/%.o)
APP_OBJS = $(APP_SRCS:%.cpp=$(OBJS_DIR)/%.o)

APP_LIBS = -lcapd -lSComplex


TEST_APP_NAME=CubSComplexTest
TEST_SRCS = $(TESTS_DIR)/CubSComplexReductionTest.cpp \
	$(TESTS_DIR)/CubSComplexTestMain.cpp \
	$(TESTS_DIR)/CubSComplexIteratorsTest.cpp \

TEST_OBJS = $(TEST_SRCS:%.cpp=$(OBJS_DIR)/%.o)

TEST_LIBS = $(APP_LIBS) -lboost_unit_test_framework


MKDIR=mkdir
CC=g++

DEBUG_ON= -g
DEBUG_OFF=-DBOOST_DISABLE_ASSERTS
DEBUG= $(DEBUG_ON) 
BOOST_FLAGS= -DBOOST_TEST_DYN_LINK
DEPEND_FLAGS=$(CFLAGS) $(LOCAL_INC_PATHS) $(BOOST_FLAGS)
WARNINGS= -Wall
# -Winline -Wdisabled-optimization
#OPT_FLAGS=-O2 -finline-functions -finline-limit=100000 --param large-function-growth=100000 --param inline-unit-growth=100000
OPT_FLAGS=-O2
PROJECT_FLAGS= $(WARNINGS) -pedantic $(OPT_FLAGS)
COMP_FLAGS=$(CFLAGS) $(DEBUG) $(PROJECT_FLAGS) $(INC_PATHS) $(BOOST_FLAGS) 



.PHONY: init depend clean libs apps all


$(BINS_DIR):
	$(MKDIR) $(BINS_DIR)

$(LIBS_DIR):
	$(MKDIR) $(LIBS_DIR)

$(OBJS_DIR):
	$(MKDIR) $(OBJS_DIR)
	$(MKDIR) $(OBJS_DIR)/$(SRCS_DIR)
	$(MKDIR) $(OBJS_DIR)/$(TESTS_DIR)

$(RUN_DIR):
	$(MKDIR) $(RUN_DIR)

all: init libs apps test

init: $(BINS_DIR) $(LIBS_DIR) $(OBJS_DIR) $(RUN_DIR)

libs: init $(LIBS_DIR)/$(LIB_NAME)

apps: init $(BINS_DIR)/$(APP_NAME)

test: init $(BINS_DIR)/$(TEST_APP_NAME)
	@time (LD_LIBRARY_PATH=$(BOOST_HOME)/lib $(PWD)/$(BINS_DIR)/$(TEST_APP_NAME) --output_format=$(TEST_OUTPUT_FORMAT) --log_level=$(TEST_LOG_LEVEL) --report_level=$(TEST_REPORT_LEVEL))

$(OBJS_DIR)/%.o: %.cpp
	@echo $(CC) $<
	@$(CC) $(COMP_FLAGS) -c $< -o $@


$(LIBS_DIR)/$(LIB_NAME): $(LIB_OBJS)
	@echo $(AR) $(LIB_NAME)
	@$(AR) rcs $(LIBS_DIR)/$(LIB_NAME) $(LIB_OBJS) 

$(BINS_DIR)/$(APP_NAME): init libs $(APP_OBJS)
	@echo $(APP_NAME)
	@$(CC) $(COMP_FLAGS) -o $(APP_NAME) $(APP_OBJS) $(LIB_PATHS) $(APP_LIBS)


$(BINS_DIR)/$(TEST_APP_NAME): init libs $(TEST_OBJS)
	@echo $(TEST_APP_NAME)
	@$(CC) -o $(BINS_DIR)/$(TEST_APP_NAME) $(LIB_PATHS) $(TEST_OBJS) $(TEST_LIBS) 


clean:
	$(RM) -r $(BUILD_DIR)/*

depend: $(APP_SRCS) $(LIB_SRCS)  $(TEST_SRCS)
	makedepend  -p$(OBJS_DIR)/ $(DEPEND_FLAGS) $^

#Bellow are dependencies generated by 'dpend' target (make depend)
# DO NOT DELETE

build/objs/src/CrHomS.o: inc/CubSComplex.hpp inc/CubSComplex_Cell.hpp
build/objs/src/CrHomS.o: inc/CubSComplex_IteratorProvider.hpp
build/objs/src/CrHomS.o: inc/CubSComplex_Iterators.hpp
build/objs/src/CrHomS.o: inc/CubSComplex_ColoredIterators.hpp
build/objs/src/CrHomS.o: inc/CubSComplex_Numerators.hpp inc/SComplexAlgs.hpp
build/objs/src/CrHomS.o: inc/SComplexAlgs_Coreduction.hpp
build/objs/src/CrHomS.o: inc/SComplexAlgs_DefaultReduceStrategy.hpp
build/objs/src/CrHomS.o: inc/SComplexAlgs_Shave.hpp
build/objs/src/CubSComplex_Cell.o: inc/CubSComplex.hpp
build/objs/src/CubSComplex_Cell.o: inc/CubSComplex_Cell.hpp
build/objs/src/CubSComplex_Cell.o: inc/CubSComplex_IteratorProvider.hpp
build/objs/src/CubSComplex_Cell.o: inc/CubSComplex_Iterators.hpp
build/objs/src/CubSComplex_Cell.o: inc/CubSComplex_ColoredIterators.hpp
build/objs/src/CubSComplex_Cell.o: inc/CubSComplex_Numerators.hpp
build/objs/test/CubSComplexReductionTest.o: inc/CubSComplex.hpp
build/objs/test/CubSComplexReductionTest.o: inc/CubSComplex_Cell.hpp
build/objs/test/CubSComplexReductionTest.o: inc/CubSComplex_IteratorProvider.hpp
build/objs/test/CubSComplexReductionTest.o: inc/CubSComplex_Iterators.hpp
build/objs/test/CubSComplexReductionTest.o: inc/CubSComplex_ColoredIterators.hpp
build/objs/test/CubSComplexReductionTest.o: inc/CubSComplex_Numerators.hpp
build/objs/test/CubSComplexReductionTest.o: inc/SComplexAlgs.hpp
build/objs/test/CubSComplexReductionTest.o: inc/SComplexAlgs_Coreduction.hpp
build/objs/test/CubSComplexReductionTest.o: inc/SComplexAlgs_DefaultReduceStrategy.hpp
build/objs/test/CubSComplexReductionTest.o: inc/SComplexAlgs_Shave.hpp
build/objs/test/CubSComplexIteratorsTest.o: inc/CubSComplex.hpp
build/objs/test/CubSComplexIteratorsTest.o: inc/CubSComplex_Cell.hpp
build/objs/test/CubSComplexIteratorsTest.o: inc/CubSComplex_IteratorProvider.hpp
build/objs/test/CubSComplexIteratorsTest.o: inc/CubSComplex_Iterators.hpp
build/objs/test/CubSComplexIteratorsTest.o: inc/CubSComplex_ColoredIterators.hpp
build/objs/test/CubSComplexIteratorsTest.o: inc/CubSComplex_Numerators.hpp
