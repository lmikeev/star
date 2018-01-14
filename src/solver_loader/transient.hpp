/*
 *  transient.hpp
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

#ifndef SOLVER_LOADER_TRANSIENT_HPP_
#define SOLVER_LOADER_TRANSIENT_HPP_

#include "base.hpp"

namespace solver_loader {

class transient : public base {
 protected:
  enum type get_type() const { return SL_TRANSIENT; }

  virtual char* sprint_solver_description(char* const buf) const {
    char* c = buf;

    c += sprint_description_solver_name(c, "transient");
    c += sprint_description_method_name(c);

    if (!is_stoch()) {
      c += sprint_description_param(c, "kinetics",
                                    is_hybrid() ? "hybrid" : "det");

      if (is_hybrid()) {
        c += sprint_description_vars(c, c_vars, "c_vars");
        c += sprint_description_vars(c, p_vars, "p_vars");
      }

      c += sprint_description_param(c, "nmoments", static_cast<int>(nmoments));
    }

    if (!method_stoch_only()) {
      c += sprint_description_param(c, "atol", abs_tol);
      c += sprint_description_param(c, "rtol", rel_tol);
    }

    return c;
  }

  bool method_stoch_only() const {
    return method_name == "ssa" || method_name == "ssa-fast" ||
           method_name == "fau" || method_name == "aau";
  }

  bool method_implicit() const {
    return method_name == "beuler" || method_name == "irk3" ||
           method_name == "ros23" || method_name == "heuler";
  }

  virtual bool preload_check_options() {
    if (method_stoch_only()) {
      set_stoch();
    }

    if (is_det()) {
      hsucc_has_tr = false;
      hsucc_has_rates_g = false;
    } else if (is_hybrid()) {
      hsucc_has_tr = true;

      if (!get("hsucc_has_rates_g", hsucc_has_rates_g)) {
        hsucc_has_rates_g = false;
      }
    } else {
      if (!get("hsucc_has_rates_g", hsucc_has_rates_g)) {
        hsucc_has_rates_g = true;
      }

      hsucc_has_tr = !hsucc_has_rates_g;
    }

    if (method_name == "") {
      if (is_det()) {
        method_name = "lsoda";
      } else {
        method_name = "rk45";
      }
    }

    bool export_mfile;
    if (get("export_mfile", export_mfile) && export_mfile) {
      hs_has_index = true;
      hsucc_has_tr = true;
    }

    if (method_name == "lsoda") {
      hs_has_index = true;
    } else if (method_name == "fau" || method_name == "aau" ||
               method_implicit()) {
      hs_has_exitrate = true;
      if (method_implicit()) {
        hs_has_flags = true;
      }
    }

    return true;
  }

  virtual bool afterload_check_options() {
    do_sa = !sa_params.empty();

    if (do_sa) {
      if (!get("sa_no_der", sa_no_der)) {
        sa_no_der = false;
      }

      if (!sa_no_der) {
        if (!get("sa_no_2der", sa_no_2der)) {
          sa_no_2der = false;
        }
      }
    }

    return true;
  }

  bool check_tasks() {
    for (auto& t : tasks) {
      switch (t->get_type()) {
        case task_info::TASK_DUMP_MOMENTS: {
          task_info::dump_moments* const ti =
              static_cast<task_info::dump_moments*>(t);
          if (is_stoch()) {
            ti->adjust_nmoments(2);

            if (ti->get_nmoments() > (int)nmoments) {
              nmoments = ti->get_nmoments();
            }
          } else {
            ti->adjust_nmoments(nmoments, true);

            if (is_det()) {
              ti->remove_cond();
            }
          }
          if (!t->check(this)) {
            return false;
          }
          break;
        }

        case task_info::TASK_PLOT: {
          if (is_det()) {
            task_info::cond* const ti = dynamic_cast<task_info::cond*>(t);
            ti->remove_cond();

            if (nmoments < 2) {
              task_info::plot* const ti_ = dynamic_cast<task_info::plot*>(t);
              ti_->do_not_plot_stddev();
            }
          }
          if (!t->check(this)) {
            return false;
          }
          break;
        }

        case task_info::TASK_PLOT_2D: {
          if (is_det()) {
            task_info::cond* const ti = dynamic_cast<task_info::cond*>(t);
            ti->remove_cond();
          }
          if (!t->check(this)) {
            return false;
          }
          break;
        }

        case task_info::TASK_DUMP_DISTR:
        case task_info::TASK_PLOT_DISTR:
        case task_info::TASK_PLOT_DISTR_2D: {
          if (!is_det()) {
            if (!t->check(this)) {
              return false;
            }
          }
          break;
        }

        case task_info::TASK_DUMP_STATS: {
          if (!t->check(this)) {
            return false;
          }
          break;
        }

        default:
          break;
      }
    }

    return true;
  }

#ifdef STAR_CODEGEN
  bool write_solver_def(std::ostream& os, const bool subsloader = false) const {
    bool methylation = false;
    get("methylation", methylation);
    if (methylation) {
      os << "#include \"../../../../src/solver/misc/methylation.hpp\""
         << std::endl;
      os << std::endl;
      os << "namespace solver" << std::endl;
      os << "{" << std::endl;
      os << std::endl;
      base::write_solver_def(os, "misc::methylation", subsloader);
      os << std::endl;
      os << "}" << std::endl;
      return true;
    }

    os << "#include \"../../../../src/solver/transient/all.hpp\"" << std::endl;
    os << std::endl;
    os << "namespace solver" << std::endl;
    os << "{" << std::endl;
    os << std::endl;

    static char solver_class[128];
    if (method_name == "euler") {
      sprintf(solver_class, "transient::EulerT<>");
    } else if (method_name == "rk23") {
      sprintf(solver_class, "transient::RK23_T<>");
    } else if (method_name == "rk45") {
      sprintf(solver_class, "transient::RK45_T<>");
    } else if (method_name == "lsoda") {
      sprintf(solver_class, "transient::LSODAwrapper");
    } else if (method_name == "matlab-ode45") {
      sprintf(solver_class,
              "transient::matlabODEwrapper<transient::MATLAB_ODE45>");
    } else if (method_name == "matlab-ode23") {
      sprintf(solver_class,
              "transient::matlabODEwrapper<transient::MATLAB_ODE23>");
    } else if (method_name == "matlab-ode113") {
      sprintf(solver_class,
              "transient::matlabODEwrapper<transient::MATLAB_ODE113>");
    } else if (method_name == "matlab-ode15s") {
      sprintf(solver_class,
              "transient::matlabODEwrapper<transient::MATLAB_ODE15S>");
    } else if (method_name == "matlab-ode23s") {
      sprintf(solver_class,
              "transient::matlabODEwrapper<transient::MATLAB_ODE23S>");
    } else if (method_name == "matlab-ode23t") {
      sprintf(solver_class,
              "transient::matlabODEwrapper<transient::MATLAB_ODE23T>");
    } else if (method_name == "matlab-ode23tb") {
      sprintf(solver_class,
              "transient::matlabODEwrapper<transient::MATLAB_ODE23TB>");
    } else if (method_name == "ssa") {
      sprintf(solver_class, "transient::SSA");
    } else if (method_name == "ssa-fast") {
      sprintf(solver_class, "transient::SSAfast");
    } else if (method_name == "fau") {
      sprintf(solver_class, "transient::FAU");
    } else if (method_name == "aau") {
      int uniformization_method = 0;
      get("uniformization_method", uniformization_method);
      switch (uniformization_method) {
        case 1:
          sprintf(solver_class, "transient::AAU<transient::AAU_EPOCH>");
          break;
        case 2:
          sprintf(solver_class, "transient::AAU<transient::AAU_APRIORI>");
          break;
        case 0:
        default:
          sprintf(solver_class, "transient::AAU<transient::AAU_STD>");
          break;
      }
    } else if (method_name == "beuler") {
      sprintf(solver_class, "transient::BEuler<>");
    } else if (method_name == "heuler") {
      sprintf(solver_class, "transient::HEuler<>");
    } else if (method_name == "bdf2") {
      sprintf(solver_class, "transient::BDF2<>");
    } else if (method_name == "irk3") {
      sprintf(solver_class, "transient::IRK3<>");
    } else if (method_name == "ros23") {
      sprintf(solver_class, "transient::Ros23<>");
    } else {
      last_err << "unknown method";
      return false;
    }

    base::write_solver_def(os, solver_class, subsloader);

    os << std::endl;
    os << "} //solver" << std::endl;
    return true;
  }

  bool write_Makefile_defs(std::ostream& os) const {
    if (method_name == "lsoda") {
      os << "SOLVERLIB_SHARED_LINKER_FLAGS  += $(shell find "
            "../../../../src/solver/lib/lsoda/CMakeFiles/lsoda_test.dir -name "
            "\\*.o)" << std::endl;
    }
    return true;
  }
#endif

 public:
  transient(experiment* const E) : base(E) {}
};
}

#endif
