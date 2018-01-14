/*
 *  optimize.hpp
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

#ifndef SOLVER_LOADER_OPTIMIZE_HPP_
#define SOLVER_LOADER_OPTIMIZE_HPP_

#include "scan.hpp"

namespace solver_loader {

class optimize : public scan {
 protected:
  enum type get_type() const { return SL_OPTIMIZE; }

  virtual char* sprint_solver_description(char* const buf) const {
    char* c = buf;

    c += sprint_description_solver_name(c, "optimize");

    c += sprint_description_param(c, "objf", objf.c_str());

    return c;
  }

  virtual bool preload_check_options() {
    if (!transient::preload_check_options()) {
      return false;
    }

    if (opt_method_name == "") {
#if HAVE_MCR
      opt_method_name = "matlab-fmincon";
#else
#if HAVE_NLOPT
      opt_method_name = "nlopt";
#else
#if HAVE_DLIB
      opt_method_name = "dlib-opt";
#else
      last_err << "no optimization solver available";
      return false;
#endif
#endif
#endif
    }

    if (opt_method_name == "dlib-opt") {
#if HAVE_DLIB
      sa_no_der = false;
#else
      last_err << "dlib is not available";
      return false;
#endif
    } else if (opt_method_name == "nlopt" ||
               opt_method_name == "nlopt-gn-direct" ||
               opt_method_name == "nlopt-gn-direct-l" ||
               opt_method_name == "nlopt-gn-direct-l-rand" ||
               opt_method_name == "nlopt-gn-direct-noscal" ||
               opt_method_name == "nlopt-gn-direct-l-noscal" ||
               opt_method_name == "nlopt-gn-direct-l-rand-noscal" ||
               opt_method_name == "nlopt-gn-mlsl") {
#if HAVE_NLOPT
      sa_no_der = true;
#else
      last_err << "nlopt is not available";
      return false;
#endif
    } else if (opt_method_name == "nlopt-gd-mlsl" ||
               opt_method_name == "nlopt-ld-mma") {
#if HAVE_NLOPT
      sa_no_der = false;
      sa_no_2der = true;
#else
      last_err << "nlopt is not available";
      return false;
#endif
    }

    else if (opt_method_name == "matlab-fmincon" ||
             opt_method_name == "matlab-globalsearch" ||
             opt_method_name == "matlab-multistart") {
#if HAVE_MCR
#else
      last_err << "mcr is not available";
      return false;
#endif
    } else if (opt_method_name == "matlab-patternsearch" ||
               opt_method_name == "matlab-ga" ||
               opt_method_name == "matlab-simulannealbnd") {
#if HAVE_MCR
      sa_no_der = true;
#else
      last_err << "mcr is not available";
      return false;
#endif
    } else {
      last_err << "unknown opt method '" << opt_method_name << "'";
      return false;
    }

    return true;
  }

  virtual bool afterload_check_options() {
    use_transient_subsolver = true;

    for (auto& p : estimate_params) {
      p.c = findcnst(p.param);
      if (p.c == nullptr) {
        last_err << "unknown parameter '" << p.param << "' [3]";
        return false;
      }
      add_param(p.c);
    }

    if (get_type() == SL_FIT) {
      for (auto& p : estimate_inits) {
        p.v = findvar(p.param);
        if (p.v == nullptr) {
          last_err << "unknown variable '" << p.param << "'";
          return false;
        }
        add_ivar(p.v);
      }

      for (auto& p : estimate_obs_errors) {
        p.v = findvar(p.param);
        if (p.v == nullptr) {
          last_err << "unknown variable '" << p.param << "'";
          return false;
        }
        add_evar(p.v);
      }
    }

    if (!no_objf && !check_objf()) {
      return false;
    }

    do_sa = true;

    if (hsucc_has_rates_g) {
      hsucc_has_tr = true;
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

    os << "#include \"../../../../src/solver/optimization/all.hpp\""
       << std::endl;
    os << std::endl;
    os << "namespace solver" << std::endl;
    os << "{" << std::endl;
    os << std::endl;

    static char solver_class[128];
    if (opt_method_name == "dlib-opt") {
      sprintf(solver_class, "optimization::dlibOptWrapper<Subsolver>");
    } else if (opt_method_name == "dlib-opt-box-constrained") {
      sprintf(solver_class,
              "optimization::dlibOptWrapper<Subsolver, "
              "optimization::DLIB_OPT_BOX_CONSTRAINED>");
    }
#if HAVE_NLOPT
    else if (opt_method_name == "nlopt") {
      sprintf(solver_class, "optimization::NLoptWrapper<Subsolver>");
    } else if (opt_method_name == "nlopt-gn-direct") {
      sprintf(solver_class,
              "optimization::NLoptWrapper<Subsolver, "
              "nlopt::algorithm::GN_DIRECT>");
    } else if (opt_method_name == "nlopt-gn-direct-l") {
      sprintf(solver_class,
              "optimization::NLoptWrapper<Subsolver, "
              "nlopt::algorithm::GN_DIRECT_L>");
    } else if (opt_method_name == "nlopt-gn-direct-l-rand") {
      sprintf(solver_class,
              "optimization::NLoptWrapper<Subsolver, "
              "nlopt::algorithm::GN_DIRECT_L_RAND>");
    } else if (opt_method_name == "nlopt-gn-direct-noscal") {
      sprintf(solver_class,
              "optimization::NLoptWrapper<Subsolver, "
              "nlopt::algorithm::GN_DIRECT_NOSCAL>");
    } else if (opt_method_name == "nlopt-gn-direct-l-noscal") {
      sprintf(solver_class,
              "optimization::NLoptWrapper<Subsolver, "
              "nlopt::algorithm::GN_DIRECT_L_NOSCAL>");
    } else if (opt_method_name == "nlopt-gn-direct-l-rand-noscal") {
      sprintf(solver_class,
              "optimization::NLoptWrapper<Subsolver, "
              "nlopt::algorithm::GN_DIRECT_L_RAND_NOSCAL>");
    } else if (opt_method_name == "nlopt-gn-mlsl") {
      sprintf(
          solver_class,
          "optimization::NLoptWrapper<Subsolver, nlopt::algorithm::GN_MLSL>");
    } else if (opt_method_name == "nlopt-gd-mlsl") {
      sprintf(
          solver_class,
          "optimization::NLoptWrapper<Subsolver, nlopt::algorithm::GD_MLSL>");
    } else if (opt_method_name == "nlopt-ld-mma") {
      sprintf(
          solver_class,
          "optimization::NLoptWrapper<Subsolver, nlopt::algorithm::LD_MMA>");
    }
#endif

#if HAVE_MCR
    else if (opt_method_name == "matlab-fmincon") {
      if (sa_no_der) {
        sprintf(solver_class,
                "optimization::matlabOptWrapper<Subsolver, "
                "optimization::MATLAB_FMINCON_F>");
      } else {
        if (sa_no_2der) {
          sprintf(solver_class,
                  "optimization::matlabOptWrapper<Subsolver, "
                  "optimization::MATLAB_FMINCON_FG>");
        } else {
          sprintf(solver_class,
                  "optimization::matlabOptWrapper<Subsolver, "
                  "optimization::MATLAB_FMINCON_FGH>");
        }
      }
    } else if (opt_method_name == "matlab-globalsearch") {
      if (sa_no_der) {
        sprintf(solver_class,
                "optimization::matlabOptWrapper<Subsolver, "
                "optimization::MATLAB_GLOBALSEARCH_F>");
      } else {
        if (sa_no_2der) {
          sprintf(solver_class,
                  "optimization::matlabOptWrapper<Subsolver, "
                  "optimization::MATLAB_GLOBALSEARCH_FG>");
        } else {
          sprintf(solver_class,
                  "optimization::matlabOptWrapper<Subsolver, "
                  "optimization::MATLAB_GLOBALSEARCH_FGH>");
        }
      }
    } else if (opt_method_name == "matlab-multistart") {
      sprintf(solver_class,
              "optimization::matlabOptWrapper<Subsolver, "
              "optimization::MATLAB_MULTISTART>");
    } else if (opt_method_name == "matlab-patternsearch") {
      sprintf(solver_class,
              "optimization::matlabOptWrapper<Subsolver, "
              "optimization::MATLAB_PATTERNSEARCH>");
    } else if (opt_method_name == "matlab-ga") {
      sprintf(solver_class,
              "optimization::matlabOptWrapper<Subsolver, "
              "optimization::MATLAB_GA>");
    } else if (opt_method_name == "matlab-simulannealbnd") {
      sprintf(solver_class,
              "optimization::matlabOptWrapper<Subsolver, "
              "optimization::MATLAB_SIMULANNEALBND>");
    }
#endif

    base::write_solver_def(os, solver_class);

    os << std::endl;
    os << "} //solver" << std::endl;
    return true;
  }

  bool write_Makefile_defs(std::ostream& os) const {
    if (!scan::write_Makefile_defs(os)) {
      return false;
    }
    return true;
  }
#endif

 private:
  bool no_objf;

 public:
  optimize(experiment* const E, const bool no_objf_ = false)
      : scan(E), no_objf(no_objf_) {
    bool methylation = false;
    if (get("methylation", methylation) && methylation) {
      no_objf = true;
    }
  }
};
}

#endif
