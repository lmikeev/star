/*
 *  plot_objf.hpp
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

#ifndef SOLVER_LOADER_TASK_INFO_PLOT_OBJF_HPP_
#define SOLVER_LOADER_TASK_INFO_PLOT_OBJF_HPP_

#include "plot_distr.hpp"

namespace solver_loader {
namespace task_info {

class plot_objf : public plot_distr {
 public:
  plot_objf(const plot_props& props, const std::vector<plot_expr>& exprs,
            const std::string& cnd_str = "", const std::string& cmp_path = "",
            const std::string& copy_data_path = "",
            const bool plot_stddev_ = false)
      : plot_distr(props, exprs, cnd_str, cmp_path, copy_data_path,
                   plot_stddev_) {}

  plot_objf(const plot_objf& ti) : plot_objf(ti.get_props(), ti.get_exprs()) {}

  type get_type() const { return TASK_PLOT_OBJF; }

  bool check(solver_loader::base const* const sl);
  bool post_check(solver_loader::base const* const sl);

  base* clone() const { return new plot_objf(*this); }
};
}
}

#endif
