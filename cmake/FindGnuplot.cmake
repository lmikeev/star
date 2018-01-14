
find_program(GNUPLOT_EXECUTABLE gnuplot)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GNUPLOT_FOUND to TRUE if 
# all listed variables are TRUE
find_package_handle_standard_args(Gnuplot
                                  REQUIRED_VARS GNUPLOT_EXECUTABLE)

mark_as_advanced(GNUPLOT_EXECUTABLE)
