/*
 *  fit.hpp
 *
 *  Created by Linar Mikeev (mikeev@cs.uni-saarland.de).
 *  Copyright (C) 2015 Saarland University. All rights reserved.
 *
 *
 *  This file is part of star.
 *
 *  star is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  star is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with star.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SOLVER_LOADER_FIT_HPP_
#define SOLVER_LOADER_FIT_HPP_

#include "optimize.hpp"

namespace solver_loader {

class fit : public optimize {
 protected:
  enum type get_type() const { return SL_FIT; }

  virtual char* sprint_solver_description(char* const buf) const {
    char* c = buf;

    c += sprint_description_solver_name(c, "fit");

    c += sprint_description_param(c, "time_series_data",
                                  time_series_data_src.c_str());

    return c;
  }

#ifdef STAR_CODEGEN
  bool write_subsolver_def(std::ostream& os) const {
    os << "#include \"../../../../src/solver/fit.hpp\"" << std::endl;
    os << std::endl;
    os << "namespace solver" << std::endl;
    os << "{" << std::endl;
    os << std::endl;
    base::write_solver_def(os, "fit", true);
    os << std::endl;
    os << "}" << std::endl;
    return true;
  }
#endif

 public:
  fit(experiment* const E) : optimize(E, true) {}
};
}

#endif
