/*
 *  euler.hpp
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

#ifndef SOLVER_TRANSIENT_EULER_HPP_
#define SOLVER_TRANSIENT_EULER_HPP_

#include "rk_base.hpp"

namespace solver {
namespace transient {

class euler_lstate_data_accessor : public rk_lstate_data_accessor {
 protected:
  enum euler_lstate_data_offset { do_k1 = do_z + 1, do_k2 };

 public:
  euler_lstate_data_accessor(options const* const o,
                             const std::size_t n = do_k2 + 1)
      : rk_lstate_data_accessor(o, n) {}

  double* f0(void* const d) const { return k1(d); }

  double const* f0(void const* const d) const { return k1(d); }

  double* k1(void* const d) const {
    return static_cast<double*>(d) + do_k1 * len;
  }

  double const* k1(void const* const d) const {
    return static_cast<double const*>(d) + do_k1 * len;
  }

  double* k2(void* const d) const {
    return static_cast<double*>(d) + do_k2 * len;
  }

  double const* k2(void const* const d) const {
    return static_cast<double const*>(d) + do_k2 * len;
  }

  virtual void zero_tmp(void* const d) const { memset(k1(d), 0, 2 * size_); }
};

template <class lstate_data_accessor_t = euler_lstate_data_accessor,
          class hstate_data_accessor_t = hstate_data_accessor,
          class hstate_succ_data_accessor_t = hstate_succ_data_accessor>
class EulerT : public RKbaseT<lstate_data_accessor_t, hstate_data_accessor_t,
                              hstate_succ_data_accessor_t> {
 public:
  typedef RKbaseT<lstate_data_accessor_t, hstate_data_accessor_t,
                  hstate_succ_data_accessor_t> RKbase;
  typedef typename RKbase::state state;

  EulerT(solver_loader::base const* const sl) : options(sl), RKbase(sl) {}
  virtual ~EulerT() {}

 protected:
  virtual double get_pow() const { return 1.0 / (1.0 + 1.0); }

  template <unsigned int k>
  struct k2t {
    typedef typename RKbaseT<lstate_data_accessor_t, hstate_data_accessor_t,
                             hstate_succ_data_accessor_t>::template k2t<k> type;
  };

  virtual void s_comp_k_set_y(void* const, const double,
                              typename k2t<1>::type) const {}

  virtual void s_comp_k_set_y(void* const ls_d, const double h,
                              typename k2t<2>::type) const {
    double* const y = this->lstate_da->y(ls_d);
    double const* const z = this->lstate_da->z(ls_d);
    double* const k1 = this->lstate_da->k1(ls_d);

    this->update_nonneg(y, k1);

    for (std::size_t i = 0; i < this->s_y_len; i++) {
      y[i] = z[i] + 0.50 * h * k1[i];
    }
  }

  virtual void final_stage_update_nonneg(void* const ls_d) const {
    double const* const y = this->lstate_da->y(ls_d);
    double* const k2 = this->lstate_da->k2(ls_d);
    this->update_nonneg(y, k2);
  }

  virtual void s_ynew_d(double* const z, void* const ls_d, const double h,
                        const std::size_t i, const bool mbnn, double& ynew,
                        double& d) const {
    double* const y = this->lstate_da->y(ls_d);
    double const* const k1 = this->lstate_da->k1(ls_d);
    double const* const k2 = this->lstate_da->k2(ls_d);

    ynew = z[i] + h * k1[i];

    d = 0.50 * (k1[i] - k2[i]);

    this->update_nonneg(ynew, mbnn, y[i]);
  }

  virtual void comp_stages(const double h) {
    this->template comp_k<1>(h);
    this->template comp_k<2>(h);
  }

  void s_comp_k_succ(cstate const* const cs, double const* const x,
                     double const* const y, void const* const hsucc_d,
                     void* const ls_d, void* const lsucc_d,
                     typename k2t<1>::type) const {
    this->s_comp_succ(cs, x, y, hsucc_d, this->lstate_da->k1(ls_d),
                      this->lstate_da->k1(lsucc_d));
  }

  void s_comp_k_succ(cstate const* const cs, double const* const x,
                     double const* const y, void const* const hsucc_d,
                     void* const ls_d, void* const lsucc_d,
                     typename k2t<2>::type) const {
    this->s_comp_succ(cs, x, y, hsucc_d, this->lstate_da->k2(ls_d),
                      this->lstate_da->k2(lsucc_d));
  }

  void s_comp_k_c(cstate const* const cs, double const* const x,
                  double const* const y, void* const ls_d,
                  typename k2t<1>::type) const {
    this->s_comp_c(cs, x, y, this->lstate_da->k1(ls_d));
  }

  void s_comp_k_c(cstate const* const cs, double const* const x,
                  double const* const y, void* const ls_d,
                  typename k2t<2>::type) const {
    this->s_comp_c(cs, x, y, this->lstate_da->k2(ls_d));
  }

  void s_comp_k_c_raw(cstate const* const cs, double const* const x,
                      double const* const y, void* const ls_d,
                      typename k2t<1>::type) const {
    this->s_comp_c_raw(cs, x, y, this->lstate_da->k1(ls_d));
  }

  void s_comp_k_c_raw(cstate const* const cs, double const* const x,
                      double const* const y, void* const ls_d,
                      typename k2t<2>::type) const {
    this->s_comp_c_raw(cs, x, y, this->lstate_da->k2(ls_d));
  }

  void s_comp_k_succ_cov(cstate const* const cs, double const* const x,
                         double const* const y, void const* const hsucc_d,
                         double const* const dx, void* const ls_d,
                         void* const lsucc_d, typename k2t<1>::type) const {
    this->s_comp_succ_cov(cs, x, y, hsucc_d, dx, this->lstate_da->k1(ls_d),
                          this->lstate_da->k1(lsucc_d));
  }

  void s_comp_k_succ_cov(cstate const* const cs, double const* const x,
                         double const* const y, void const* const hsucc_d,
                         double const* const dx, void* const ls_d,
                         void* const lsucc_d, typename k2t<2>::type) const {
    this->s_comp_succ_cov(cs, x, y, hsucc_d, dx, this->lstate_da->k2(ls_d),
                          this->lstate_da->k2(lsucc_d));
  }

  void s_comp_k_c_cov(cstate const* const cs, double const* const x,
                      double const* const y, void* const ls_d,
                      typename k2t<1>::type) const {
    this->s_comp_c_cov(cs, x, y, this->lstate_da->k1(ls_d));
  }

  void s_comp_k_c_cov(cstate const* const cs, double const* const x,
                      double const* const y, void* const ls_d,
                      typename k2t<2>::type) const {
    this->s_comp_c_cov(cs, x, y, this->lstate_da->k2(ls_d));
  }
};
}
}

#endif
