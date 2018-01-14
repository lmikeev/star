/*
 *  plot_objf_2d.hpp
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

#ifndef SOLVER_LOADER_TASK_INFO_PLOT_OBJF_2D_HPP_
#define SOLVER_LOADER_TASK_INFO_PLOT_OBJF_2D_HPP_

#include "plot_distr_2d.hpp"

namespace solver_loader {
namespace task_info {

class plot_objf_2d : public plot_distr_2d {
 public:
  plot_objf_2d(const plot_props& props, const std::string& expr1,
               const std::string& expr2,
               const plot_surface_props& surface_props,
               const std::string& cnd_str = "",

               const std::string& dparam = "", const std::string& dparam2 = "")
      : plot_distr_2d(props, expr1, expr2, surface_props, cnd_str, dparam,
                      dparam2) {}

  plot_objf_2d(const plot_objf_2d& ti)
      : plot_objf_2d(ti.get_props(), ti.get_expr1(), ti.get_expr2(), ti) {}

  type get_type() const { return TASK_PLOT_OBJF_2D; }

  bool check(solver_loader::base const* const sl);
  bool post_check(solver_loader::base const* const sl);

  base* clone() const { return new plot_objf_2d(*this); }
};
}
}

#endif
