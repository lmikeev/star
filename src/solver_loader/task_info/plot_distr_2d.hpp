/*
 *  plot_distr_2d.hpp
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

#ifndef SOLVER_LOADER_TASK_INFO_PLOT_DISTR_2D_HPP_
#define SOLVER_LOADER_TASK_INFO_PLOT_DISTR_2D_HPP_

#include "base.hpp"

namespace solver_loader {
namespace task_info {

class plot_distr_2d : public plot_2d_base, public plot_surface_props {
 private:
 public:
  plot_distr_2d(const plot_props& props, const std::string& expr1,
                const std::string& expr2,
                const plot_surface_props& surface_props,
                const std::string& cnd_str = "", const std::string& dparam = "",
                const std::string& dparam2 = "", const bool isdynamic_ = false)
      : plot_2d_base(props, expr1, expr2, cnd_str, dparam, dparam2, isdynamic_),
        plot_surface_props(surface_props) {}

  plot_distr_2d(const plot_distr_2d& ti)
      : plot_distr_2d(ti.get_props(), ti.get_expr1(), ti.get_expr2(), ti,
                      ti.get_cond_str(), ti.get_dparam(), ti.get_dparam2()) {}

  virtual type get_type() const { return TASK_PLOT_DISTR_2D; }

  virtual bool check(solver_loader::base const* const sl);

  base* clone() const { return new plot_distr_2d(*this); }
};
}
}

#endif
