/*
 *  base.hpp
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

#ifndef SOLVER_TRANSIENT_BASE_HPP_
#define SOLVER_TRANSIENT_BASE_HPP_

#include "../base.hpp"

namespace solver {
namespace transient {

class base : public solver::base {
 public:
  base(solver_loader::base const* const sl)
      : options(sl),
        solver::base(sl),
        error(0.0),
        h_sum(0.0),
        niter(0),
        niter_failed(0) {
    printf("abs_tol = %lg\nrel_tol = %lg\n", sl->abs_tol, sl->rel_tol);
  }

  virtual ~base() {}

 protected:
  double hmin, hmax;

  double error;
  double h_sum;
  unsigned int niter, niter_failed;

  bool get_timepoint(const unsigned int timepoint_index,
                     double& timepoint) const {
    timepoint = sl->tspan[timepoint_index];
    return true;
  }

  double get_error() const { return error; }

  double get_avg_timestep() const { return h_sum / niter; };

  unsigned int get_niter() const { return niter; }

  unsigned int get_niter_failed() const { return niter_failed; }

  unsigned int ss_miniter, ss_maxiter;
  double ss_reltol, ss_abstol;

  bool update_progress(const double t) const {
#ifdef STAR_WEB_INTERFACE
    if (t < sl->tspan[sl->tspan.size() - 1])
#endif
    {
      const double progress =
          (t - sl->tspan.front()) / sl->tspan.back() * 100.0;
      return sl->update_progress(progress);
    }
    return true;
  }

  bool dump(const double t, const std::size_t timepoint_index) {
    int timepoint_id;
    if (!create_timepoint(t, timepoint_id)) {
      return false;
    }

    return dump(t, timepoint_index, timepoint_id);
  }

  virtual bool dump(const double, const std::size_t, const std::size_t) {
    return true;
  }

  virtual bool init_(solver_loader::model_info::ic const* const) {
    return true;
  }

  virtual bool init_() { return init_(this->sl->get_ics().front()); }

  virtual bool iterate0(const double, double&, const bool, double&, bool& redo,
                        bool&) {
    redo = true;
    return true;
  }

  virtual bool post_iterate() { return true; }

  virtual bool iterate(const double H, const double t0, double& h) {
    double t_ = 0.0;
    bool nofailed = true;
    bool done = false;
    while (!done) {
      hmin = 16.0 * std::numeric_limits<double>::epsilon() * t_;
      h = std::min(hmax, std::max(hmin, h));

      if (1.10 * h >= H - t_) {
        h = H - t_;
        done = true;
      }

      const double h_old = h;
      double lerr;
      bool redo = false;
      if (!iterate0(t0 + t_, h, done, lerr, redo, nofailed)) {
        return false;
      }

      if (redo) {
        niter_failed++;
        nofailed = false;
        done = false;
        continue;
      } else {
        nofailed = true;
      }

      h_sum += h_old;
      t_ += h_old;
      error = lerr;
      niter++;
    }

    post_iterate();

    return true;
  }

  virtual bool init() { return init_(); }

  void set_hmax() { hmax = 0.10 * (sl->tspan.back() - sl->tspan.front()); }

  virtual void comp_h0(double& h) {
    h = 0.0000010 * (sl->tspan.back() - sl->tspan.front());
  }

  virtual bool run_() {
    h_sum = 0.0;
    niter = 0;
    niter_failed = 0;

    set_hmax();

    auto ti = sl->tspan.begin();
    double t = *ti, h;
    comp_h0(h);

    if (!update_progress(t) || !dump(t, ti - sl->tspan.begin())) {
      return false;
    }

    ti++;

    while ((this->do_steady_state() || ti != sl->tspan.end()) &&
           (elapsed_time < max_runtime || max_runtime <= 0.0)) {
      if (!iterate(*ti - t, t, h)) {
        return false;
      }
      t = *ti;

      update_elapsed_time();

      if (!update_progress(t) || !dump(t, ti - sl->tspan.begin())) {
        return false;
      }

      ti++;
    }

    return true;
  }

  virtual bool run_objf(srch_pt_s& p) {
    update(p);

    if (!init_()) {
      return false;
    }

    if (!this->do_steady_state()) {
      set_hmax();
      double h;
      comp_h0(h);
      if (!iterate(sl->tspan.back() - sl->tspan.front(), sl->tspan.front(),
                   h)) {
        return false;
      }

      comp_objf(p.f.data());
    } else {
    }

    return true;
  }

  bool err_tol_not_met(const double t) {
    this->sl->last_error()
        << "Unable to meet integration tolerances without reducing the "
           "step size below the smallest value allowed (" << hmin
        << ") at time t = " << t;

    return false;
  }
};
}
}

#endif
