/*
 *  beuler.hpp
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

#ifndef SOLVER_TRANSIENT_BEULER_HPP_
#define SOLVER_TRANSIENT_BEULER_HPP_

#include "impl_base.hpp"
#include "../linsys.hpp"

namespace solver {
namespace transient {

class beuler_lstate_data_accessor : public impl_lstate_data_accessor {
 protected:
  enum beuler_lstate_data_offset {
    do_u = do_z + 1,
    do_g,

    do_s_u,
    do_s_y,
    do_s_g
  };

 public:
  beuler_lstate_data_accessor(options const* const o,
                              const std::size_t n = do_s_g + 1)
      : impl_lstate_data_accessor(o, n) {}

  double* f0(void* const d) const { return s_g(d); }

  double const* f0(void const* const d) const { return s_g(d); }

  double* u(void* const d) const {
    return static_cast<double*>(d) + do_u * len;
  }

  double const* u(void const* const d) const {
    return static_cast<double const*>(d) + do_u * len;
  }

  double* g(void* const d) const {
    return static_cast<double*>(d) + do_g * len;
  }

  double const* g(void const* const d) const {
    return static_cast<double const*>(d) + do_g * len;
  }

  double* s_u(void* const d) const {
    return static_cast<double*>(d) + do_s_u * len;
  }

  double const* s_u(void const* const d) const {
    return static_cast<double const*>(d) + do_s_u * len;
  }

  double* s_y(void* const d) const {
    return static_cast<double*>(d) + do_s_y * len;
  }

  double const* s_y(void const* const d) const {
    return static_cast<double const*>(d) + do_s_y * len;
  }

  double* s_g(void* const d) const {
    return static_cast<double*>(d) + do_s_g * len;
  }

  double const* s_g(void const* const d) const {
    return static_cast<double const*>(d) + do_s_g * len;
  }

  virtual void init(void* const d) const {
    memcpy(z(d), y(d), size_);
    memcpy(u(d), y(d), size_);
    memcpy(g(d), y(d), size_);

    memset(s_y(d), 0, size_);
    memset(s_u(d), 0, size_);
    memset(s_g(d), 0, size_);
  }

  virtual void zero(void* const d) const {
    memset(z(d), 0, size_);
    memset(y(d), 0, size_);
    memset(u(d), 0, size_);
    memset(g(d), 0, size_);
  }

  virtual void update(void* const d) const {
    memcpy(z(d), y(d), size_);
    memcpy(u(d), y(d), size_);
    memcpy(g(d), y(d), size_);
  }

  virtual void reset(void* const d) const {
    memcpy(y(d), z(d), size_);
    memcpy(u(d), z(d), size_);
    memcpy(g(d), z(d), size_);
  }
};

template <class lstate_data_accessor_t = beuler_lstate_data_accessor,
          class hstate_data_accessor_t = hstate_data_accessor,
          class hstate_succ_data_accessor_t = hstate_succ_data_accessor>
class beuler_linSys
    : public linSysT<lstate_data_accessor_t, hstate_data_accessor_t,
                     hstate_succ_data_accessor_t> {
 public:
  typedef linSysT<lstate_data_accessor_t, hstate_data_accessor_t,
                  hstate_succ_data_accessor_t> linSys;
  typedef typename linSys::state state;
  typedef typename linSys::hstate hstate;
  typedef typename linSys::hstate_succ hstate_succ;
  typedef typename linSys::hashListBased hashListBased;

  beuler_linSys(hashListBased* const hlbs) : linSys(hlbs) {}

 protected:
  void propagate(state* const s, const std::size_t,
                 const bool neg = false) const {
    LINSYS_PRE

    double const* const u = this->get_lstate_da()->u(ls_d);
    double const* const y = this->get_lstate_da()->y(ls_d);
    double const* const g = this->get_lstate_da()->g(ls_d);

    LINSYS_LOOP

    double* const s_u = this->get_lstate_da()->s_u(lsucc_d);
    double* const s_y = this->get_lstate_da()->s_y(lsucc_d);
    double* const s_g = this->get_lstate_da()->s_g(lsucc_d);

    this->hlbs->s_comp_succ(rates_g, rates_h, u, s_u);
    this->hlbs->s_comp_succ(rates_g, rates_h, y, s_y);
    this->hlbs->s_comp_succ(rates_g, rates_h, g, s_g);

    LINSYS_POST
  }

  bool update(state* const s, const std::size_t niter, double const h,
              double& jerr1, double& jnormy1, double& jnormynew1, double& jerr2,
              double& jnormy2, double& jnormynew2, double& jerr3,
              double& jnormy3, double& jnormynew3, double& err, double& normy,
              double& normynew) const {
    hstate* const hs = s->hs;
    const double exitrate = this->da->hs_exitrate(hs);

    assert(exitrate >= 0.0);

    void* const ls_d = this->da->get_lstate_data(s);
    auto const sda = this->get_lstate_da();
    double const* const z = sda->z(ls_d);

    double* const u = sda->u(ls_d);
    double* const y = sda->y(ls_d);
    double* const g = sda->g(ls_d);

    double* const s_u = sda->s_u(ls_d);
    double* const s_y = sda->s_y(ls_d);
    double* const s_g = sda->s_g(ls_d);

    const double su = this->v_p(s_u);

    const double u_ =
        (this->v_p(z) + 0.50 * h * su) / (1.0 + 0.50 * h * exitrate);

    this->v_p(s_u) = 0.0;

    const double uold = this->v_p(u);
    const double unew = uold + this->omega * (u_ - uold);

    assert(unew >= 0.0);

    const double sy = this->v_p(s_y);

    const double y_ = (unew + 0.50 * h * sy) / (1.0 + 0.50 * h * exitrate);

    this->v_p(s_y) = 0.0;

    const double yold = this->v_p(y);
    const double ynew = yold + this->omega * (y_ - yold);

    assert(ynew >= 0.0);

    const double sg = this->v_p(s_g);

    const double g_ = (this->v_p(z) + h * sg) / (1.0 + h * exitrate);

    this->v_p(s_g) = 0.0;

    const double gold = this->v_p(g);
    const double gnew = gold + this->omega * (g_ - gold);

    assert(gnew >= 0.0);

    this->da->hs_exitrate(hs) = -exitrate;

    this->v_p(u) = unew - uold;
    this->v_p(y) = ynew - yold;
    this->v_p(g) = gnew - gold;

    propagate(s, niter, true);

    this->v_p(u) = unew;
    this->v_p(y) = ynew;
    this->v_p(g) = gnew;

    if (this->can_ignore_s(uold) && this->can_ignore_s(unew) &&
        this->can_ignore_s(yold) && this->can_ignore_s(ynew) &&
        this->can_ignore_s(gold) && this->can_ignore_s(gnew) &&
        this->can_ignore_s(unew * exitrate) &&
        this->can_ignore_s(ynew * exitrate) &&
        this->can_ignore_s(gnew * exitrate)) {
      this->s_make_temp(s);
      return false;
    }

    this->hlbs->update_err(jerr1, jnormy1, jnormynew1, unew - uold,
                           this->v_p(z), ynew);

    this->hlbs->update_err(jerr2, jnormy2, jnormynew2, ynew - yold, unew, ynew);

    this->hlbs->update_err(jerr3, jnormy3, jnormynew3, gnew - gold,
                           this->v_p(z), gnew);

    this->hlbs->update_err(err, normy, normynew, gnew - ynew, this->v_p(z),
                           ynew);

    return true;
  }
};

template <class lstate_data_accessor_t = beuler_lstate_data_accessor,
          class hstate_data_accessor_t = hstate_data_accessor,
          class hstate_succ_data_accessor_t = hstate_succ_data_accessor>
class BEuler : public ImplBaseT<lstate_data_accessor_t, hstate_data_accessor_t,
                                hstate_succ_data_accessor_t> {
 public:
  typedef ImplBaseT<lstate_data_accessor_t, hstate_data_accessor_t,
                    hstate_succ_data_accessor_t> ImplBase;

  typedef typename ImplBase::state state;
  typedef typename ImplBase::hstate hstate;
  typedef typename ImplBase::hstate_succ hstate_succ;

  BEuler(solver_loader::base const* const sl)
      : options(sl), ImplBase(sl), ls(new beuler_linSys<>(this)) {}

  std::size_t get_niter_total() const { return ls->get_niter_total(); }

 private:
  beuler_linSys<>* const ls;

 protected:
  virtual double get_pow() const { return 1.0 / (1.0 + 1.0); }

  void comp_stages(const double h, double& err, double& normy,
                   double& normynew) {
    ls->solve(h, this->prob_sum, err, normy, normynew);
  }
};
}
}

#endif
