/*
 *  linsys.hpp
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

#ifndef SOLVER_LINSYS_HPP_
#define SOLVER_LINSYS_HPP_

#include <unordered_set>
#include "hash_list_based.hpp"

#define LINSYS_PRE                              \
  hstate* const hs = s->hs;                     \
  cstate const* const cs = this->da->hs_cs(hs); \
  double const* const x = nullptr;              \
  void* const ls_d = this->da->get_lstate_data(s);

#define LINSYS_LOOP                                                      \
  if (!neg) {                                                            \
    this->states->s_explore_unexplored(s);                               \
  }                                                                      \
  if (!neg || !this->da->hs_is_unexplored(hs)) {                         \
    hstate_succ* hsucc = hs->succ;                                       \
    double exitrate = 0.0;                                               \
    while (hsucc != nullptr) {                                           \
      void const* const hsucc_d = this->da->get_hstate_succ_data(hsucc); \
      double const* rates_g;                                             \
      double const* rates_h;                                             \
      this->hlbs->tr_rates(hsucc_d, cs, x, rates_g, rates_h);            \
      const bool succ_erased = this->da->hs_is_erased(hsucc->hs);        \
      const bool succ_temp = this->da->hs_is_temp(hsucc->hs);            \
      assert(!succ_erased || !succ_temp);                                \
      if (!neg || (!succ_erased && !succ_temp &&                         \
                   this->da->hs_exitrate(hsucc->hs) >= 0.0)) {           \
        if (succ_erased) {                                               \
          this->states->hs_restore(hsucc->hs);                           \
        } else if (succ_temp) {                                          \
          this->states->hs_restore_temp(hsucc->hs);                      \
        }                                                                \
        void* const lsucc_d = this->da->get_lstate_data(hsucc->hs->ls);

#define LINSYS_POST                            \
  }                                            \
  hsucc = hsucc->lnext;                        \
  exitrate += this->tr_rate(rates_g, rates_h); \
  }                                            \
  if (!neg) {                                  \
    this->da->hs_exitrate(hs) = exitrate;      \
  }                                            \
  }

namespace solver {

template <class lstate_data_accessor_t = lstate_data_accessor,
          class hstate_data_accessor_t = hstate_data_accessor,
          class hstate_succ_data_accessor_t = hstate_succ_data_accessor,
          typename lstate_t = lstate_h, typename hstate_t = hstate_h,
          typename hstate_succ_t = hstate_succ_h>
class linSysT {
 public:
  typedef hashListBasedT<lstate_data_accessor_t, hstate_data_accessor_t,
                         hstate_succ_data_accessor_t> hashListBased;

  typedef typename hashListBased::dataAccessor dataAccessor;
  typedef typename hashListBased::stateHashList stateHashList;

  typedef typename hashListBased::state state;
  typedef typename hashListBased::hstate hstate;
  typedef typename hashListBased::hstate_succ hstate_succ;

  linSysT(hashListBased* const hlbs)
      : hlbs(hlbs),
        da(hlbs->get_da()),
        states(hlbs->get_states()),
        omega(1.0),
        niter_total(0) {}

  virtual ~linSysT() {}

  virtual void propagate(state* const s, const std::size_t niter,
                         const bool neg = false) const = 0;

  virtual bool update(state* const s, const std::size_t niter, double const h,
                      double& jerr1, double& jnormy1, double& jnormynew1,
                      double& jerr2, double& jnormy2, double& jnormynew2,
                      double& jerr3, double& jnormy3, double& jnormynew3,
                      double& err, double& normy, double& normynew) const = 0;

  void solve(double const h, const double prob_sum, double& err, double& normy,
             double& normynew) {
    const double RTOL = this->hlbs->sl->rel_tol;

    double jerr1 = 1.0, jerr2 = 1.0, jerr3 = 1.0;
    std::size_t niter_ = 0;

    std::unordered_set<std::uintptr_t> chksum_;
    std::uintptr_t chksum = 0;

    while (jerr1 > RTOL || jerr2 > RTOL || jerr3 > RTOL ||
           !chksum_.count(chksum)) {
      chksum_.insert(chksum);

      std::size_t si = states->get_nactive();
      state* s = states->get_first();

      while (si--) {
        state* const lnext = s->lnext;

        propagate(s, niter_);

        s = lnext;
      }

      jerr1 = 0.0;
      double jnormy1 = 0.0;
      double jnormynew1 = 0.0;

      jerr2 = 0.0;
      double jnormy2 = 0.0;
      double jnormynew2 = 0.0;

      jerr3 = 0.0;
      double jnormy3 = 0.0;
      double jnormynew3 = 0.0;

      err = 0.0;
      normy = 0.0;
      normynew = 0.0;

      chksum = 0;

      s = states->get_first();
      si = states->get_nactive();
      while (si--) {
        state* const lnext = s->lnext;

        hlbs->s_comp_unexplored_exitrate(s);

        if (update(s, niter_, h, jerr1, jnormy1, jnormynew1, jerr2, jnormy2,
                   jnormynew2, jerr2, jnormy2, jnormynew2, err, normy,
                   normynew)) {
          chksum ^= reinterpret_cast<std::uintptr_t>(s);
        }

        s = lnext;
      }

      hlbs->update_err(jerr1, jnormy1, jnormynew1);
      hlbs->update_err(jerr2, jnormy2, jnormynew2);
      hlbs->update_err(jerr3, jnormy3, jnormynew3);

      niter_++;

      if (!(niter_ & 0xfff)) {
        printf("\t> %lu : %lu\t%le\t%le %le %le\t%lu\n", niter_,
               states->get_nactive(), 1.0 - prob_sum, jerr1, jerr2, jerr3,
               chksum);
      }
    }

    niter_total += niter_;
  }

  std::size_t get_niter_total() const { return niter_total; }

 protected:
  void s_make_temp(state* const s) const {
    da->hs_exitrate(s->hs) = std::fabs(da->hs_exitrate(s->hs));
    states->s_make_temp(s);
  }

  hashListBased* const hlbs;

  dataAccessor const* const da;
  stateHashList* const states;

  const double omega;

  std::size_t niter_total;

  lstate_data_accessor_t const* get_lstate_da() const {
    return hlbs->get_lstate_da();
  }

  double tr_rate(double const* const rates_g,
                 double const* const rates_h) const {
    double const* const x = nullptr;
    return hlbs->tr_rate(rates_g, rates_h, x);
  }

  double& v_p(double* const v) const { return hlbs->v_p(v); }

  double v_p(double const* const v) const { return hlbs->v_p(v); }

  bool can_ignore_tr(const double d, const double h) const {
    return hlbs->can_ignore_tr(d * h);
  }

  bool can_ignore_s(const double p) const { return hlbs->can_ignore_s(p); }
};
}

#endif
