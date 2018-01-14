
set(_POSSIBLE_LIBSBML_INCLUDE include include/sbml)
set(_POSSIBLE_LIBSBML_LIBRARY sbml)

find_path(LIBSBML_INCLUDE_DIR SBMLTypes.h
	HINTS
	$ENV{LIBSBML_DIR}
	PATH_SUFFIXES ${_POSSIBLE_LIBSBML_INCLUDE}
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

find_library(LIBSBML_LIBRARY
	NAMES ${_POSSIBLE_LIBSBML_LIBRARY}
	HINTS
	$ENV{LIBSBML_DIR}
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
# handle the QUIETLY and REQUIRED arguments and set LIBSBML_FOUND to TRUE if
# all listed variables are TRUE
find_package_handle_standard_args(LibSBML
                                  REQUIRED_VARS LIBSBML_LIBRARY LIBSBML_INCLUDE_DIR)

mark_as_advanced(LIBSBML_INCLUDE_DIR LIBSBML_LIBRARY)
