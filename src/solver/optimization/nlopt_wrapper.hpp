/*
 *  nlopt_wrapper.hpp
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

#ifndef SOLVER_OPTIMIZATION_NLOPT_WRAPPER_HPP_
#define SOLVER_OPTIMIZATION_NLOPT_WRAPPER_HPP_

#if HAVE_NLOPT
#include <nlopt.hpp>
#endif

#include "base.hpp"
#include "../external.hpp"

namespace solver {
namespace optimization {

template <class S
#if HAVE_NLOPT
          ,
          nlopt::algorithm alg_id = nlopt::algorithm::GD_MLSL
#endif
          >
class NLoptWrapper : public base<S> {
 protected:
  double nlopt_f(const std::vector<double>& x, std::vector<double>& grad,
                 void* data) {
    assert(x.size() == this->nlparams);

    srch_pt_s p;
    p.c.resize(this->nlparams);
    for (std::size_t i = 0; i < this->nlparams; i++) {
      p.c[i] = this->from_log_space(x[i]);
    }

    if (!this->run_objf(p)) {
      return false;
    }

    if (!grad.empty()) {
      for (std::size_t i = 0; i < this->nlparams; i++) {
        grad[i] =
            this->v_df_dc(p.f.data())[i] * p.c[i] * this->log_log_space_base;
      }
    }

    this->add_srch_pt(p);

    return p.f[0];
  }

  bool search(std::vector<srch_pt_s>& pts) {
#if HAVE_NLOPT
    nlopt::opt opt(alg_id, this->nlparams);

    std::vector<double> lb(this->nlparams), ub(this->nlparams);

    pts.resize(1);
    pts.front().c.resize(this->nlparams);
    pts.front().f.resize(this->f_len);

    std::size_t k = 0;
    for (std::size_t i = 0; i < this->sl->estimate_params.size(); i++, k++) {
      const est_param_info& ep = this->sl->estimate_params[i];
      lb[k] = this->to_log_space(ep.min);
      ub[k] = this->to_log_space(ep.max);
      pts.front().c[k] = this->to_log_space(ep.value0);
    }
    for (std::size_t i = 0; i < this->sl->estimate_inits.size(); i++, k++) {
      const est_param_info& ep = this->sl->estimate_inits[i];
      lb[k] = this->to_log_space(ep.min);
      ub[k] = this->to_log_space(ep.max);
      pts.front().c[k] = this->to_log_space(ep.value0);
    }
    for (std::size_t i = 0; i < this->sl->estimate_obs_errors.size();
         i++, k++) {
      const est_param_info& ep = this->sl->estimate_obs_errors[i];
      lb[k] = this->to_log_space(ep.min);
      ub[k] = this->to_log_space(ep.max);
      pts.front().c[k] = this->to_log_space(ep.value0);
    }

    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);

    if (this->reltol_x > 0.0 || this->abstol_x > 0.0 || this->reltol_f > 0.0 ||
        this->abstol_f > 0.0) {
      if (this->reltol_x > 0.0) {
        opt.set_xtol_rel(this->reltol_x);
      }

      if (this->abstol_x > 0.0) {
        opt.set_xtol_abs(this->abstol_x);
      }

      if (this->reltol_f > 0.0) {
        opt.set_ftol_rel(this->reltol_f);
      }

      if (this->abstol_f > 0.0) {
        opt.set_ftol_abs(this->abstol_f);
      }
    } else {
      this->reltol_x = 1e-3;
      opt.set_xtol_rel(this->reltol_x);
    }

    if (this->max_fevals > 0) {
      opt.set_maxeval(this->max_fevals);
    }

    if (this->max_runtime > 0.0) {
      opt.set_maxtime(this->max_runtime);
    }

    if (this->maximize) {
      opt.set_max_objective(::nlopt_f, nullptr);
    } else {
      opt.set_min_objective(::nlopt_f, nullptr);
    }

    nlopt::result result = opt.optimize(pts.front().c, pts.front().f.front());

    for (std::size_t i = 0; i < this->nlparams; i++) {
      pts.front().c[i] = this->from_log_space(pts.front().c[i]);
    }

    if (result < 0) {
      return false;
    }

    return true;
#endif

    this->sl->last_error() << "NLopt not available";
    return false;
  }

 public:
  NLoptWrapper(solver_loader::base const* const sl) : options(sl), base<S>(sl) {
    external_set_solver(this);
  }
};
}
}

#endif
