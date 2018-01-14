
# solver library Makefile

include ../../../../Makefile.config

SOLVERLIB_BINARY = libsolver.so

SOLVERLIB_CXX_FLAGS	+= $(CXX_FLAGS) #-g #TODO
SOLVERLIB_SHARED_LINKER_FLAGS += $(SHARED_LINKER_FLAGS) -shared

SOLVERLIB_SRCS = $(shell find . -name \*.hpp -o -name \*.cpp)

SOLVERLIB_CXX_FLAGS += -I'$(BOOST_INCLUDE_DIR)' -I'$(LUA_INCLUDE_DIR)'
SOLVERLIB_SHARED_LINKER_FLAGS += $(BOOST_LIBRARIES)

#SOLVERLIB_SHARED_LINKER_FLAGS += -lgmp

SOLVERLIB_SHARED_LINKER_FLAGS += ../../../../CMakeFiles/star.dir/src/solver/base.cpp.o
SOLVERLIB_SHARED_LINKER_FLAGS += ../../../../CMakeFiles/star.dir/src/solver/external.cpp.o
SOLVERLIB_SHARED_LINKER_FLAGS += ../../../../CMakeFiles/star.dir/src/solver_loader/base.cpp.o
SOLVERLIB_SHARED_LINKER_FLAGS += ../../../../CMakeFiles/star.dir/src/solver_loader/parser/parser.cpp.o
SOLVERLIB_SHARED_LINKER_FLAGS += ../../../../CMakeFiles/star.dir/src/solver_loader/parser/sbml.cpp.o

ifeq ($(HAVE_DLIB),TRUE)
	SOLVERLIB_SHARED_LINKER_FLAGS += ../../../../CMakeFiles/star.dir/src/solver/lib/dlib/all/source.cpp.o
endif

ifeq ($(WEB_INTERFACE),ON)
	SOLVERLIB_CXX_FLAGS += -D'STAR_WEB_INTERFACE' -I'$(LIBMYSQLCLIENT_INCLUDE_DIR)'
	SOLVERLIB_SHARED_LINKER_FLAGS += $(LIBMYSQLCLIENT_LIBRARY)
endif

ifeq ($(HAVE_GNUPLOT),TRUE)
	SOLVERLIB_CXX_FLAGS += -D'HAVE_GNUPLOT'
endif

ifeq ($(HAVE_LIBSBML),TRUE)
	SOLVERLIB_CXX_FLAGS += -D'HAVE_LIBSBML' -I'$(LIBSBML_INCLUDE_DIR)'
	SOLVERLIB_SHARED_LINKER_FLAGS += $(LIBSBML_LIBRARY)
endif

ifeq ($(HAVE_MCR),TRUE)
	SOLVERLIB_CXX_FLAGS += -D'HAVE_MCR' -I'$(MCR_INCLUDE_DIR)'
	SOLVERLIB_SHARED_LINKER_FLAGS += $(MCR_LIBRARY) $(BINARY_DIR)/$(LIBMATLABSD_LIBRARY)
endif

ifeq ($(HAVE_NLOPT),TRUE)
	SOLVERLIB_CXX_FLAGS += -D'HAVE_NLOPT' -I'$(NLOPT_INCLUDE_DIR)'
	SOLVERLIB_SHARED_LINKER_FLAGS += $(NLOPT_LIBRARY)
endif

SOLVERLIB_OBJS = $(filter %.o, $(SOLVERLIB_SRCS:.cpp=.o))

$(SOLVERLIB_BINARY): $(SOLVERLIB_OBJS)
	@$(CXX_COMPILER) $(SOLVERLIB_OBJS) $(SOLVERLIB_SHARED_LINKER_FLAGS) -o $@

.c.o:
	@$(C_COMPILER) $(SOLVERLIB_C_FLAGS) -c -o $@ $<

.cpp.o:
	@$(CXX_COMPILER) $(SOLVERLIB_CXX_FLAGS) -c -o $@ $<

all: $(SOLVERLIB_BINARY)
