/*
 *  plot.hpp
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

#ifndef SOLVER_LOADER_TASK_INFO_PLOT_HPP_
#define SOLVER_LOADER_TASK_INFO_PLOT_HPP_

#include "base.hpp"

namespace solver_loader {
namespace task_info {

class plot_srch_pts;

class plot : public plot_base {
 protected:
  std::vector<plot_expr> exprs;

 private:
  const std::string copy_data_path;
  bool plot_stddev_;

 public:
  plot(const plot_props& props, const std::vector<plot_expr>& exprs,
       const std::string& cnd_str = "", const std::string& copy_data_path = "",
       const bool plot_stddev_ = false, const bool isdynamic_ = true)
      : plot_base(props, cnd_str, isdynamic_),
        exprs(exprs),
        copy_data_path(copy_data_path),
        plot_stddev_(plot_stddev_) {}

  plot(const plot& ti)
      : plot(ti.get_props(), ti.get_exprs(), ti.get_cond_str(),
             ti.get_copy_data_path(), ti.plot_stddev(), ti.is_dynamic()) {}

  virtual ~plot() {
    for (auto& e : exprs) {
      e.free();
    }
  }

  virtual type get_type() const { return TASK_PLOT; }

  const std::vector<plot_expr>& get_exprs() const { return exprs; }

  char const* get_copy_data_path() const { return copy_data_path.c_str(); }

  void do_not_plot_stddev() { plot_stddev_ = false; }

  bool plot_stddev() const { return plot_stddev_; }

  virtual bool check(solver_loader::base const* const sl) {
    name = "plot";

    if (!check_exprs(sl)) {
      return false;
    }
    return true;
  }

  bool check_exprs(solver_loader::base const* const sl);

  base* clone() const { return new plot(*this); }
};
}
}

#endif
