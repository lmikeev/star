/*
 *  plot_distr.hpp
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

#ifndef SOLVER_LOADER_TASK_INFO_PLOT_DISTR_HPP_
#define SOLVER_LOADER_TASK_INFO_PLOT_DISTR_HPP_

#include "plot.hpp"

namespace solver_loader {
namespace task_info {

class plot_distr : public plot, public cmp {
 public:
  plot_distr(const plot_props& props, const std::vector<plot_expr>& exprs,
             const std::string& cnd_str = "", const std::string& cmp_path = "",
             const std::string& copy_data_path = "",
             const bool plot_stddev_ = false, const bool isdynamic_ = false)
      : plot(props, exprs, cnd_str, copy_data_path, plot_stddev_, isdynamic_),
        cmp(cmp_path) {}

  plot_distr(const plot_distr& ti)
      : plot_distr(ti.get_props(), ti.get_exprs(), ti.get_cond_str(),
                   ti.get_cmp_path()) {}

  virtual type get_type() const { return TASK_PLOT_DISTR; }

  virtual bool check(solver_loader::base const* const sl) {
    name = "plot_distr";

    if (!check_exprs(sl)) {
      return false;
    }

    return true;
  }

  base* clone() const { return new plot_distr(*this); }
};
}
}

#endif
