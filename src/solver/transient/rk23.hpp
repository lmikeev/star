/*
 *  rk23.hpp
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

#ifndef SOLVER_TRANSIENT_RK23_HPP_
#define SOLVER_TRANSIENT_RK23_HPP_

#include "euler.hpp"

namespace solver {
namespace transient {

class rk23_lstate_data_accessor : public euler_lstate_data_accessor {
 protected:
  enum rk23_lstate_data_offset { do_k3 = do_k2 + 1, do_k4 };

 public:
  rk23_lstate_data_accessor(options const* const o,
                            const std::size_t n = do_k4 + 1)
      : euler_lstate_data_accessor(o, n) {}

  double* k3(void* const d) const {
    return static_cast<double*>(d) + do_k3 * len;
  }

  double const* k3(void const* const d) const {
    return static_cast<double const*>(d) + do_k3 * len;
  }

  double* k4(void* const d) const {
    return static_cast<double*>(d) + do_k4 * len;
  }

  double const* k4(void const* const d) const {
    return static_cast<double const*>(d) + do_k4 * len;
  }

  virtual void zero_tmp(void* const d) const { memset(k1(d), 0, 4 * size_); }
};

template <class lstate_data_accessor_t = rk23_lstate_data_accessor,
          class hstate_data_accessor_t = hstate_data_accessor,
          class hstate_succ_data_accessor_t = hstate_succ_data_accessor>
class RK23_T : public EulerT<lstate_data_accessor_t, hstate_data_accessor_t,
                             hstate_succ_data_accessor_t> {
 public:
  typedef EulerT<lstate_data_accessor_t, hstate_data_accessor_t,
                 hstate_succ_data_accessor_t> Euler;

  RK23_T(solver_loader::base const* const sl) : options(sl), Euler(sl) {}
  virtual ~RK23_T() {}

 protected:
  virtual double get_pow() const { return 1.0 / (1.0 + 2.0); }

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
      static const double k2_c1 = 1.0 / 2.0;

      y[i] = z[i] + h * k2_c1 * k1[i];
    }
  }

  virtual void s_comp_k_set_y(void* const ls_d, const double h,
                              typename k2t<3>::type) const {
    double* const y = this->lstate_da->y(ls_d);
    double const* const z = this->lstate_da->z(ls_d);
    double* const k2 = this->lstate_da->k2(ls_d);

    this->update_nonneg(y, k2);

    for (std::size_t i = 0; i < this->s_y_len; i++) {
      static const double k3_c2 = 3.0 / 4.0;

      y[i] = z[i] + h * k3_c2 * k2[i];
    }
  }

  virtual void s_comp_k_set_y(void* const ls_d, const double h,
                              typename k2t<4>::type) const {
    double* const y = this->lstate_da->y(ls_d);
    double const* const z = this->lstate_da->z(ls_d);
    double const* const k1 = this->lstate_da->k1(ls_d);
    double const* const k2 = this->lstate_da->k2(ls_d);
    double* const k3 = this->lstate_da->k3(ls_d);

    this->update_nonneg(y, k3);

    for (std::size_t i = 0; i < this->s_y_len; i++) {
      static const double k4_c1 = 2.0 / 9.0;
      static const double k4_c2 = 1.0 / 3.0;
      static const double k4_c3 = 4.0 / 9.0;

      y[i] = z[i] + h * (k4_c1 * k1[i] + k4_c2 * k2[i] + k4_c3 * k3[i]);
    }
  }

  virtual void final_stage_update_nonneg(void* const ls_d) const {
    double const* const y = this->lstate_da->y(ls_d);
    double* const k4 = this->lstate_da->k4(ls_d);
    this->update_nonneg(y, k4);
  }

  virtual void s_ynew_d(double* const z, void* const ls_d, const double h,
                        const std::size_t i, const bool mbnn, double& ynew,
                        double& d) const {
    double* const y = this->lstate_da->y(ls_d);
    double const* const k1 = this->lstate_da->k1(ls_d);
    double const* const k2 = this->lstate_da->k2(ls_d);
    double const* const k3 = this->lstate_da->k3(ls_d);
    double const* const k4 = this->lstate_da->k4(ls_d);

    static const double d_c1 = -15.0 / 216.0;
    static const double d_c2 = 1.0 / 12.0;
    static const double d_c3 = 1.0 / 9.0;
    static const double d_c4 = -1.0 / 8.0;

    static const double y_c1 = 2.0 / 9.0;
    static const double y_c2 = 1.0 / 3.0;
    static const double y_c3 = 4.0 / 9.0;

    d = d_c1 * k1[i] + d_c2 * k2[i] + d_c3 * k3[i] + d_c4 * k4[i];

    ynew = z[i] + h * (y_c1 * k1[i] + y_c2 * k2[i] + y_c3 * k3[i]);

    this->update_nonneg(ynew, mbnn, y[i]);
  }

  virtual void comp_stages(const double h) {
    this->template comp_k<1>(h);
    this->template comp_k<2>(h);
    this->template comp_k<3>(h);
    this->template comp_k<4>(h);
  }

  using Euler::s_comp_k_succ;
  using Euler::s_comp_k_c;
  using Euler::s_comp_k_c_raw;
  using Euler::s_comp_k_succ_cov;
  using Euler::s_comp_k_c_cov;

  void s_comp_k_succ(cstate const* const cs, double const* const x,
                     double const* const y, void const* const hsucc_d,
                     void* const ls_d, void* const lsucc_d,
                     typename k2t<3>::type) const {
    this->s_comp_succ(cs, x, y, hsucc_d, this->lstate_da->k3(ls_d),
                      this->lstate_da->k3(lsucc_d));
  }

  void s_comp_k_succ(cstate const* const cs, double const* const x,
                     double const* const y, void const* const hsucc_d,
                     void* const ls_d, void* const lsucc_d,
                     typename k2t<4>::type) const {
    this->s_comp_succ(cs, x, y, hsucc_d, this->lstate_da->k4(ls_d),
                      this->lstate_da->k4(lsucc_d));
  }

  void s_comp_k_c(cstate const* const cs, double const* const x,
                  double const* const y, void* const ls_d,
                  typename k2t<3>::type) const {
    this->s_comp_c(cs, x, y, this->lstate_da->k3(ls_d));
  }

  void s_comp_k_c(cstate const* const cs, double const* const x,
                  double const* const y, void* const ls_d,
                  typename k2t<4>::type) const {
    this->s_comp_c(cs, x, y, this->lstate_da->k4(ls_d));
  }

  void s_comp_k_c_raw(cstate const* const cs, double const* const x,
                      double const* const y, void* const ls_d,
                      typename k2t<3>::type) const {
    this->s_comp_c_raw(cs, x, y, this->lstate_da->k3(ls_d));
  }

  void s_comp_k_c_raw(cstate const* const cs, double const* const x,
                      double const* const y, void* const ls_d,
                      typename k2t<4>::type) const {
    this->s_comp_c_raw(cs, x, y, this->lstate_da->k4(ls_d));
  }

  void s_comp_k_succ_cov(cstate const* const cs, double const* const x,
                         double const* const y, void const* const hsucc_d,
                         double const* const dx, void* const ls_d,
                         void* const lsucc_d, typename k2t<3>::type) const {
    this->s_comp_succ_cov(cs, x, y, hsucc_d, dx, this->lstate_da->k3(ls_d),
                          this->lstate_da->k3(lsucc_d));
  }

  void s_comp_k_succ_cov(cstate const* const cs, double const* const x,
                         double const* const y, void const* const hsucc_d,
                         double const* const dx, void* const ls_d,
                         void* const lsucc_d, typename k2t<4>::type) const {
    this->s_comp_succ_cov(cs, x, y, hsucc_d, dx, this->lstate_da->k4(ls_d),
                          this->lstate_da->k4(lsucc_d));
  }

  void s_comp_k_c_cov(cstate const* const cs, double const* const x,
                      double const* const y, void* const ls_d,
                      typename k2t<3>::type) const {
    this->s_comp_c_cov(cs, x, y, this->lstate_da->k3(ls_d));
  }

  void s_comp_k_c_cov(cstate const* const cs, double const* const x,
                      double const* const y, void* const ls_d,
                      typename k2t<4>::type) const {
    this->s_comp_c_cov(cs, x, y, this->lstate_da->k4(ls_d));
  }
};
}
}

#endif
