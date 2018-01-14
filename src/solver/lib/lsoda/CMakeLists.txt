
set(lsoda_BINARY lsoda_test)

file(GLOB lsoda_SOURCES *.c)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D'_LSODA_MAIN' -fPIC")
set(CMAKE_SHARED_LIBRARY_C_FLAGS "${CMAKE_SHARED_LIBRARY_C_FLAGS} -fPIC")

add_executable(${lsoda_BINARY} ${lsoda_SOURCES})
