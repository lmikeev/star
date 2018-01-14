
REQUIREMENTS

- C++11-capable compiler
- CMake

- libLua
- libBoost

optional:

- Gnuplot (plotting)
- NLopt (optimization)
- libSBML (SBML import/export)

============================================================

COMPILATION

1. cmake CMakeLists.txt
2. make

============================================================

USAGE

  ./star -e <experiment_file_name>  

Please take a look at sample models in 'models' folder
and sample experiments in 'experiments' folder.

The results are saved into the 'output' folder.
The 'output' folder can be emptied by running
  ./star --cleanout

============================================================

DOCUMENTATION

  http://almacompute.mmci.uni-saarland.de/star/help

============================================================

CONTACT

  http://groups.google.com/group/star-tool
  mikeev@cs.uni-saarland.de
