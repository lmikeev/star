
cmake_minimum_required(VERSION 2.8)
project(star C CXX)

set(star_VERSION_MAJOR 1)
set(star_VERSION_MINOR 0)

option(web-interface  "web-interface"    OFF)
option(codegen        "codegen"          ON)
option(no-sbml        "no-sbml"          ON)
option(no-mcr         "no-mcr"           ON)
option(no-nlopt       "no-nlopt"         OFF)
option(no-dlib        "no-dlib"          ON)

###############################################################################

set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/lib/x86_64-linux-gnu)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${star_SOURCE_DIR}/cmake)

#set(Boost_USE_STATIC_LIBS ON)

find_package(Boost COMPONENTS filesystem program_options system REQUIRED)
find_package(Lua REQUIRED)
find_package(Gnuplot QUIET)

if(NOT no-mcr)
  find_package(MCR QUIET)
endif()

if(NOT no-nlopt)
  find_package(NLopt QUIET)
endif()

if(NOT no-sbml)
  find_package(LibSBML REQUIRED)
endif()

if(web-interface)
  find_package(LibMySQLclient REQUIRED)
endif()

if(NOT codegen)
  find_package(LLVM REQUIRED)
  include_directories(${LLVM_INCLUDE_DIR})
  link_directories(${LLVM_LIBRARY_DIR})
endif()

include_directories(${Boost_INCLUDE_DIR})
include_directories(${LUA_INCLUDE_DIR})

###############################################################################

file(GLOB star_SOURCES *.cpp)
file(GLOB star_HEADERS *.hpp)

file(GLOB src_SOURCES src/*.cpp)
file(GLOB src_HEADERS src/*.hpp)

set(star_SOURCES ${star_SOURCES} ${src_SOURCES})
set(star_HEADERS ${star_HEADERS} ${src_HEADERS})


file(GLOB src_solver_loader_SOURCES src/solver_loader/*.cpp)
file(GLOB src_solver_loader_HEADERS src/solver_loader/*.hpp)

set(star_SOURCES ${star_SOURCES} ${src_solver_loader_SOURCES})
set(star_HEADERS ${star_HEADERS} ${src_solver_loader_HEADERS})

file(GLOB src_solver_loader_parser_SOURCES src/solver_loader/parser/*.cpp)
file(GLOB src_solver_loader_parser_HEADERS src/solver_loader/parser/*.hpp)

set(star_SOURCES ${star_SOURCES} ${src_solver_loader_parser_SOURCES})
set(star_HEADERS ${star_HEADERS} ${src_solver_loader_parser_HEADERS})

set(star_SOURCES ${star_SOURCES} src/solver/base.cpp)
set(star_HEADERS ${star_HEADERS} src/solver/base.hpp)

set(star_SOURCES ${star_SOURCES} src/solver/external.cpp)
set(star_HEADERS ${star_HEADERS} src/solver/external.hpp)

if(NOT no-dlib)
  set(star_SOURCES ${star_SOURCES} src/solver/lib/dlib/all/source.cpp)
  list(APPEND star_COMPILE_DEFINITIONS "DLIB_NO_GUI_SUPPORT")
endif()

###############################################################################

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")

add_definitions(-DNDEBUG)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

set(star_BINARY star)

add_executable(${star_BINARY} ${star_SOURCES} ${star_HEADERS})

target_link_libraries(${star_BINARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${LUA_LIBRARIES}
)

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  target_link_libraries(${star_BINARY} dl)
  if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
    set(PLATFORM_DIR "linux64")
  else("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
    set(PLATFORM_DIR "linux32")
  endif("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
  set(LIBMATLABSD_LIBRARY "libmatlabsd.so")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
    set(PLATFORM_DIR "win64")
  else("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
    set(PLATFORM_DIR "win32")
  endif("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
  set(LIBMATLABSD_LIBRARY "libmatlabsd.dll")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(PLATFORM_DIR "mac")
  set(LIBMATLABSD_LIBRARY "libmatlabsd.dylib")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

if(MCR_FOUND)
  if(NOT EXISTS "${star_BINARY_DIR}/${LIBMATLABSD_LIBRARY}")
    set(MCR_FOUND FALSE)
  endif()
endif()

if(web-interface)
  target_link_libraries(${star_BINARY} ${LIBMYSQLCLIENT_LIBRARY})
  list(APPEND star_COMPILE_DEFINITIONS "STAR_WEB_INTERFACE")
  include_directories(${LIBMYSQLCLIENT_INCLUDE_DIR})
endif()

if(codegen)
  list(APPEND star_COMPILE_DEFINITIONS "STAR_CODEGEN")
endif()

if(GNUPLOT_FOUND)
  list(APPEND star_COMPILE_DEFINITIONS "HAVE_GNUPLOT")
endif()

if(MCR_FOUND)
  list(APPEND star_COMPILE_DEFINITIONS "HAVE_MCR")
  include_directories(${MCR_INCLUDE_DIR})

  find_library(LIBMATLABSD matlabsd HINTS "${star_BINARY_DIR}")
  target_link_libraries(${star_BINARY} ${MCR_LIBRARY} ${LIBMATLABSD})
endif()

if(LLVM_FOUND)
  list(APPEND star_COMPILE_DEFINITIONS "LLVM_JIT")
  target_link_libraries(${star_BINARY} ${LLVM_MODULE_LIBS} z ncurses)
endif()

if(LIBSBML_FOUND)
  target_link_libraries(${star_BINARY} ${LIBSBML_LIBRARY})
  list(APPEND star_COMPILE_DEFINITIONS "HAVE_LIBSBML")
  include_directories(${LIBSBML_INCLUDE_DIR})
endif()

if(NLOPT_FOUND)
  list(APPEND star_COMPILE_DEFINITIONS "HAVE_NLOPT")
endif()

set_target_properties(${star_BINARY} PROPERTIES COMPILE_DEFINITIONS "${star_COMPILE_DEFINITIONS}")

#TODO
#set_target_properties(${star_BINARY} PROPERTIES COMPILE_FLAGS "-std=c++0x")

###############################################################################

add_subdirectory(src/solver/lib/lsoda)

###############################################################################

configure_file(
  "${star_SOURCE_DIR}/config.hpp.in"
  "${star_SOURCE_DIR}/config.hpp"
)

configure_file(
  "${star_SOURCE_DIR}/Makefile.config.in"
  "${star_SOURCE_DIR}/Makefile.config"
)

###############################################################################

if(LIBSBML_FOUND)
endif()


###############################################################################

#TODO: set default values for Unix Makefile
#param = {}
#init = {}
#tspan = linspace(0,1)
#ntraces = 100
#nmoments = 2

add_custom_target(
  run_src
  ${star_BINARY} -s "\"$(src)\""
)

add_custom_target(
  run_bg
  nohup ./${star_BINARY} -e "$(experiment)" > out.txt 2> err.txt < /dev/null &
)

add_custom_target(
  run_bg_src
  nohup ./$(star_BINARY) -s "\"$(src)\"" > out.txt 2> err.txt < /dev/null &
)

add_custom_target(
  run_transient
  ./${star_BINARY} -s "\"transient({model='$(model)',tspan=$(tspan),plot=true,stats=true});\""
)

add_custom_target(
  run_transient_hybrid
  ./${star_BINARY} -s "\"transient({model='$(model)',kinetics='hybrid',nmoments=$(nmoments),tspan=$(tspan),plot=true,stats=true});\""
)

add_custom_target(
  run_transient_det
  ./${star_BINARY} -s "\"transient({model='$(model)',kinetics='det',tspan=$(tspan),plot=true,stats=true});\""
)

add_custom_target(
  run_transient_ssa
  ./${star_BINARY} -s "\"transient({solver='ssa',model='$(model)',tspan=$(tspan),ntrajectories=$(ntrajectories),plot=true,stats=true});\""
)

add_custom_target(
  mfile_export_hybrid
  ./${star_BINARY} -s "\"export_mfile({model='$(model)',kinetics='hybrid',tspan=$(tspan),nmoments=$(nmoments)});\""
)

add_custom_target(
  mfile_export_det
  ./${star_BINARY} -s "\"export_mfile({model='$(model)',kinetics='det',tspan=$(tspan),nmoments=$(nmoments)});\""
)
