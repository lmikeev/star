/*
 *  scan.hpp
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

#ifndef SOLVER_LOADER_SCAN_HPP_
#define SOLVER_LOADER_SCAN_HPP_

#include "transient.hpp"
#include "parser/parser.hpp"

namespace solver_loader {

class scan : public transient {
 protected:
  enum type get_type() const { return SL_SCAN; }

  virtual char* sprint_solver_description(char* const buf) const {
    char* c = buf;

    c += sprint_description_solver_name(c, "scan");

    c += sprint_description_param(c, "objf", objf.c_str());

    return c;
  }

  bool check_objf() {
    if (objf == "") {
      last_err << "no objective function specified";
      return false;
    }

    std::size_t objf_len = 1;
    if (!sa_no_der) {
      objf_len += sa_params.size();
      if (!sa_no_2der) {
        objf_len += PAR2_LEN(sa_params.size());
      }
    }

    objf_et.resize(objf_len);
    objf_index.resize(objf_len);

    solver_loader::parser::state ps;
    if (!solver_loader::parser::parse_expr(objf.c_str(), this, objf_et[0],
                                           ps)) {
      last_error() << solver_loader::parser::get_error(ps);
      return false;
    }

    add_expr(objf_et[0], objf_index[0]);

    if (!sa_no_der) {
      for (std::size_t i = 0, ij = sa_params.size(); i < sa_params.size();
           i++) {
        objf_et[1 + i] = build_expr_der(objf_et[0], sa_params[i]);
        add_expr(objf_et[1 + i], objf_index[1 + i]);

        if (!sa_no_2der) {
          for (std::size_t j = 0; j <= i; j++, ij++) {
            objf_et[1 + ij] = build_expr_der(objf_et[1 + i], sa_params[j]);
            add_expr(objf_et[1 + ij], objf_index[1 + ij]);
          }
        }
      }

      set_stoch();
    }

    return true;
  }

  virtual bool afterload_check_options() {
    use_transient_subsolver = true;

    for (auto& p : iterate_params) {
      p.c = findcnst(p.param);
      if (p.c == nullptr) {
        last_err << "unknown parameter '" << p.param << "'";
        return false;
      }
      add_param(p.c);
    }
    for (auto& p : sample_params) {
      p.c = findcnst(p.param);
      if (p.c == nullptr) {
        last_err << "unknown parameter '" << p.param << "'";
        return false;
      }
      add_param(p.c);
    }

    for (auto& p : iterate_inits) {
      p.v = findvar(p.param);
      if (p.v == nullptr) {
        last_err << "unknown variable '" << p.param << "'";
        return false;
      }
      add_ivar(p.v);
    }
    for (auto& p : sample_inits) {
      p.v = findvar(p.param);
      if (p.v == nullptr) {
        last_err << "unknown variable '" << p.param << "'";
        return false;
      }
      add_ivar(p.v);
    }

    if (!check_objf()) {
      return false;
    }

    do_sa = false;
    sa_no_der = true;

    if (hsucc_has_rates_g) {
      hsucc_has_tr = true;
    }

    return true;
  }

  bool check_tasks() {
    for (auto const& t : tasks) {
      switch (t->get_type()) {
        case task_info::TASK_DUMP_STATS:
        case task_info::TASK_PLOT_OBJF:
        case task_info::TASK_PLOT_OBJF_2D:
          if (!t->check(this)) {
            return false;
          }
          break;

        default:
          break;
      }
    }

    return true;
  }

  bool post_check_tasks() {
    for (auto const& t : tasks) {
      switch (t->get_type()) {
        case task_info::TASK_DUMP_STATS:
        case task_info::TASK_PLOT_OBJF:
        case task_info::TASK_PLOT_OBJF_2D:
          if (!t->post_check(this)) {
            return false;
          }
          break;

        default:
          break;
      }
    }

    return true;
  }

#ifdef STAR_CODEGEN
  virtual bool write_subsolver_def(std::ostream& os) const {
    if (!transient::write_solver_def(os, true)) {
      return false;
    }
    return true;
  }

  bool write_solver_def(std::ostream& os, const bool = false) const {
    if (!write_subsolver_def(os)) {
      return false;
    }

    os << "#include \"../../../../src/solver/scan.hpp\"" << std::endl;
    os << std::endl;
    os << "namespace solver" << std::endl;
    os << "{" << std::endl;
    os << std::endl;
    base::write_solver_def(os, "scan<Subsolver>");
    os << std::endl;
    os << "} //solver" << std::endl;
    return true;
  }

  bool write_Makefile_defs(std::ostream& os) const {
    if (!transient::write_Makefile_defs(os)) {
      return false;
    }
    os << "SOLVERLIB_CXX_FLAGS  += -D'STAR_USE_TRANSIENT_SUBSOLVER'"
       << std::endl;

    return true;
  }
#endif

 public:
  scan(experiment* const E) : transient(E) {}
};
}

#endif
