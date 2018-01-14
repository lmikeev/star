/*
 *  impl_base.hpp
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

#ifndef SOLVER_TRANSIENT_IMPL_HPP_
#define SOLVER_TRANSIENT_IMPL_HPP_

#include "rk_base.hpp"
#include "../linsys.hpp"

namespace solver {
namespace transient {

class impl_lstate_data_accessor : public rk_lstate_data_accessor {
 public:
  impl_lstate_data_accessor(options const* const o, const std::size_t n = 2)
      : rk_lstate_data_accessor(o, n) {}

  virtual void update(void* const d) const { memcpy(z(d), y(d), size_); }

  virtual void reset(void* const d) const { memcpy(y(d), z(d), size_); }
};

template <class lstate_data_accessor_t = impl_lstate_data_accessor,
          class hstate_data_accessor_t = hstate_data_accessor,
          class hstate_succ_data_accessor_t = hstate_succ_data_accessor>
class ImplBaseT : public RKbaseT<lstate_data_accessor_t, hstate_data_accessor_t,
                                 hstate_succ_data_accessor_t> {
 public:
  typedef RKbaseT<lstate_data_accessor_t, hstate_data_accessor_t,
                  hstate_succ_data_accessor_t> RKbase;

  typedef typename RKbase::state state;
  typedef typename RKbase::hstate hstate;
  typedef typename RKbase::hstate_succ hstate_succ;

  ImplBaseT(solver_loader::base const* const sl) : options(sl), RKbase(sl) {}

  virtual ~ImplBaseT() {}

  virtual std::size_t get_niter_total() const = 0;

 protected:
  using RKbase::comp_stages;
  virtual void comp_stages(const double, double&, double&, double&) {}

  bool iterate0(const double t, double& h, const bool done, double& lerr,
                bool& redo, bool& nofailed) {
    assert(this->get_states()->get_first_temp() == nullptr);

    double err, normy, normynew;

    comp_stages(h, err, normy, normynew);

    this->update_err(err, normy, normynew);

    redo = err > this->sl->rel_tol;

    const std::size_t nstates_before = this->get_states()->get_nactive();
    const double h_before = h;
    bool success = true;
    this->prob_sum = 0.0;

    if (!redo) {
      state* s = this->get_states()->get_first();
      while (s != nullptr) {
        state* const lnext = s->lnext;
        void* const ls_d = this->get_da()->get_lstate_data(s);
        const double p = this->get_da()->s_p(ls_d);

        if (!this->can_ignore_s(p)) {
          this->lstate_da->update(ls_d);
          this->prob_sum += p;
        } else {
          this->get_da()->hs_exitrate(s->hs) =
              std::fabs(this->get_da()->hs_exitrate(s->hs));
          this->get_states()->s_erase(s);
        }

        s = lnext;
      }

      s = this->get_states()->get_first_temp();
      while (s != nullptr) {
        state* const lnext = s->lnext;
        void* const ls_d = this->get_da()->get_lstate_data(s);
        this->lstate_da->zero(ls_d);
        this->get_da()->hs_unset_temp(s->hs);
        this->get_states()->erase_temp(s);
        s = lnext;
      }

      lerr = 1.0 - this->prob_sum;
    } else {
      state* s = this->get_states()->get_first();
      while (s != nullptr) {
        state* const lnext = s->lnext;

        this->get_states()->s_reset(s);

        void* const ls_d = this->get_da()->get_lstate_data(s);
        const double p = this->get_da()->s_p(ls_d);

        if (!this->can_ignore_s(p)) {
          this->prob_sum += p;
        } else {
          this->get_da()->hs_exitrate(s->hs) =
              std::fabs(this->get_da()->hs_exitrate(s->hs));
          this->get_states()->s_erase(s);
        }

        s = lnext;
      }

      s = this->get_states()->get_first_temp();
      while (s != nullptr) {
        state* const lnext = s->lnext;
        void* const ls_d = this->get_da()->get_lstate_data(s);
        this->get_da()->hs_unset_temp(s->hs);
        this->get_states()->restore_temp(s);
        this->lstate_da->reset(ls_d);

        const double p = this->get_da()->s_p(ls_d);

        if (!this->can_ignore_s(p)) {
          this->prob_sum += p;
        } else {
          this->get_da()->hs_exitrate(s->hs) =
              std::fabs(this->get_da()->hs_exitrate(s->hs));
          this->get_states()->s_erase(s);
        }

        s = lnext;
      }

      lerr = 999.0;

      if (h <= this->hmin) {
        success = false;
      }

      nofailed = false;
      h = std::max(this->hmin,
                   h * std::max(0.10, 0.80 * std::pow(this->sl->rel_tol / err,
                                                      this->get_pow())));
    }

    if (nofailed) {
      const double tmp = 1.25 * pow(err / this->sl->rel_tol, this->get_pow());
      if (tmp > 0.20) {
        h /= tmp;
      } else {
        h *= 5.0;
      }
    }

    if (done || !success || !((this->niter + this->niter_failed) & 0xfff)) {
      printf("%9u %5u\t%lu %le %le %lu(%lu,%lu) %le %d %le\n", this->niter,
             this->niter_failed, this->get_niter_total(), t, h_before,
             this->get_states()->get_nactive(), nstates_before,
             this->get_states()->get_ntotal(), err, redo, 1.0 - this->prob_sum);
    }

    if (!success) {
      return this->err_tol_not_met(t);
    }
    return true;
  }
};
}
}

#endif
