/*
 *  dump_moments.hpp
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

#ifndef SOLVER_LOADER_TASK_INFO_DUMP_MOMENTS_HPP_
#define SOLVER_LOADER_TASK_INFO_DUMP_MOMENTS_HPP_

#include "base.hpp"

namespace solver_loader {

namespace model_info {

class var;
}

namespace task_info {

class dump_moments : public base, public cond, public cmp {
 private:
  const std::vector<std::string> var_names;

 protected:
  std::vector<model_info::var const*> vars;
  int nmoments;
  const bool raw;

 public:
  dump_moments(const int nmoments = -1, const bool raw = false,
               const std::string& cnd = "", const std::string& cmp_path = "")
      : cond(cnd), cmp(cmp_path), nmoments(nmoments), raw(raw) {}

  dump_moments(const std::vector<std::string>& var_names,
               const int nmoments = -1, const bool raw = false,
               const std::string& cnd = "", const std::string& cmp_path = "")
      : cond(cnd),
        cmp(cmp_path),
        var_names(var_names),
        nmoments(nmoments),
        raw(raw) {}

  dump_moments(const dump_moments& ti)
      : dump_moments(ti.get_var_names(), ti.get_nmoments(), ti.dump_raw(),
                     ti.get_cond_str(), ti.get_cmp_path()) {}

  type get_type() const { return TASK_DUMP_MOMENTS; }

  const std::vector<std::string>& get_var_names() const { return var_names; }

  const std::vector<model_info::var const*>& get_vars() const { return vars; }

  int get_nmoments() const { return nmoments; }

  void adjust_nmoments(const int nmoments_, const bool force_ = false) {
    if (nmoments < 0) {
      nmoments = nmoments_;
    } else {
      if (force_ && nmoments > nmoments_) {
        nmoments = nmoments_;
      }
    }
  }

  bool dump_raw() const { return raw; }

  bool check(solver_loader::base const* const sl) {
    name = "moments";
    if (!check_var_names(sl, var_names, vars)) {
      return false;
    }
    if (!check_cond(sl, name)) {
      return false;
    }
    if (raw) {
      name += ", raw";
    } else {
      name += ", central";
    }
    return true;
  }

  base* clone() const { return new dump_moments(*this); }
};
}
}

#endif
