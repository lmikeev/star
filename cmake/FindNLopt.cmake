
set(_POSSIBLE_NLOPT_INCLUDE include)
set(_POSSIBLE_NLOPT_LIBRARY nlopt)

find_path(NLOPT_INCLUDE_DIR nlopt.hpp
	HINTS
	$ENV{NLOPT_DIR}
	PATH_SUFFIXES ${_POSSIBLE_NLOPT_INCLUDE}
	PATHS
	~/Library/Frameworks
	/Library/Frameworks
	/usr/local
	/usr
	/sw # Fink
	/opt/local # DarwinPorts
	/opt/csw # Blastwave
	/opt
)

find_library(NLOPT_LIBRARY
	NAMES ${_POSSIBLE_NLOPT_LIBRARY}
	HINTS
	$ENV{NLOPT_DIR}
	PATH_SUFFIXES lib64 lib
	PATHS
	~/Library/Frameworks
	/Library/Frameworks
	/usr/local
	/usr
	/sw
	/opt/local
	/opt/csw
	/opt
)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set NLOPT_FOUND to TRUE if
# all listed variables are TRUE
find_package_handle_standard_args(NLopt
                                  REQUIRED_VARS NLOPT_LIBRARY NLOPT_INCLUDE_DIR)

mark_as_advanced(NLOPT_INCLUDE_DIR NLOPT_LIBRARY)
