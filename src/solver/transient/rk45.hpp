/*
 *  rk45.hpp
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

#ifndef SOLVER_TRANSIENT_RK45_HPP_
#define SOLVER_TRANSIENT_RK45_HPP_

#include "rk23.hpp"

namespace solver {
namespace transient {

class rk45_lstate_data_accessor : public rk23_lstate_data_accessor {
 protected:
  enum rk45_lstate_data_offset { do_k5 = do_k4 + 1, do_k6, do_k7 };

 public:
  rk45_lstate_data_accessor(options const* const o,
                            const std::size_t n = do_k7 + 1)
      : rk23_lstate_data_accessor(o, n) {}

  double* k5(void* const d) const {
    return static_cast<double*>(d) + do_k5 * len;
  }

  double const* k5(void const* const d) const {
    return static_cast<double const*>(d) + do_k5 * len;
  }

  double* k6(void* const d) const {
    return static_cast<double*>(d) + do_k6 * len;
  }

  double const* k6(void const* const d) const {
    return static_cast<double const*>(d) + do_k6 * len;
  }

  double* k7(void* const d) const {
    return static_cast<double*>(d) + do_k7 * len;
  }

  double const* k7(void const* const d) const {
    return static_cast<double const*>(d) + do_k7 * len;
  }

  virtual void zero_tmp(void* const d) const { memset(k1(d), 0, 7 * size_); }
};

template <class lstate_data_accessor_t = rk45_lstate_data_accessor,
          class hstate_data_accessor_t = hstate_data_accessor,
          class hstate_succ_data_accessor_t = hstate_succ_data_accessor>
class RK45_T : public RK23_T<lstate_data_accessor_t, hstate_data_accessor_t,
                             hstate_succ_data_accessor_t> {
 public:
  typedef RK23_T<lstate_data_accessor_t, hstate_data_accessor_t,
                 hstate_succ_data_accessor_t> RK23;

  RK45_T(solver_loader::base const* const sl) : options(sl), RK23(sl) {}
  virtual ~RK45_T() {}

 protected:
  virtual double get_pow() const { return 1.0 / (1.0 + 4.0); }

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
      static const double k2_c1 = 1.0 / 5.0;

      y[i] = z[i] + h * k2_c1 * k1[i];
    }
  }

  virtual void s_comp_k_set_y(void* const ls_d, const double h,
                              typename k2t<3>::type) const {
    double* const y = this->lstate_da->y(ls_d);
    double const* const z = this->lstate_da->z(ls_d);
    double const* const k1 = this->lstate_da->k1(ls_d);
    double* const k2 = this->lstate_da->k2(ls_d);

    this->update_nonneg(y, k2);

    for (std::size_t i = 0; i < this->s_y_len; i++) {
      static const double k3_c1 = 3.0 / 40.0;
      static const double k3_c2 = 9.0 / 40.0;

      y[i] = z[i] + h * (k3_c1 * k1[i] + k3_c2 * k2[i]);
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
      static const double k4_c1 = 44.0 / 45.0;
      static const double k4_c2 = -56.0 / 15.0;
      static const double k4_c3 = 32.0 / 9.0;

      y[i] = z[i] + h * (k4_c1 * k1[i] + k4_c2 * k2[i] + k4_c3 * k3[i]);
    }
  }

  virtual void s_comp_k_set_y(void* const ls_d, const double h,
                              typename k2t<5>::type) const {
    double* const y = this->lstate_da->y(ls_d);
    double const* const z = this->lstate_da->z(ls_d);
    double const* const k1 = this->lstate_da->k1(ls_d);
    double const* const k2 = this->lstate_da->k2(ls_d);
    double const* const k3 = this->lstate_da->k3(ls_d);
    double* const k4 = this->lstate_da->k4(ls_d);

    this->update_nonneg(y, k4);

    for (std::size_t i = 0; i < this->s_y_len; i++) {
      static const double k5_c1 = 19372.0 / 6561.0;
      static const double k5_c2 = -25360.0 / 2187.0;
      static const double k5_c3 = 64448.0 / 6561.0;
      static const double k5_c4 = -212.0 / 729.0;

      y[i] =
          z[i] +
          h * (k5_c1 * k1[i] + k5_c2 * k2[i] + k5_c3 * k3[i] + k5_c4 * k4[i]);
    }
  }

  virtual void s_comp_k_set_y(void* const ls_d, const double h,
                              typename k2t<6>::type) const {
    double* const y = this->lstate_da->y(ls_d);
    double const* const z = this->lstate_da->z(ls_d);
    double const* const k1 = this->lstate_da->k1(ls_d);
    double const* const k2 = this->lstate_da->k2(ls_d);
    double const* const k3 = this->lstate_da->k3(ls_d);
    double const* const k4 = this->lstate_da->k4(ls_d);
    double* const k5 = this->lstate_da->k5(ls_d);

    this->update_nonneg(y, k5);

    for (std::size_t i = 0; i < this->s_y_len; i++) {
      static const double k6_c1 = 9017.0 / 3168.0;
      static const double k6_c2 = -355.0 / 33.0;
      static const double k6_c3 = 46732.0 / 5247.0;
      static const double k6_c4 = 49.0 / 176.0;
      static const double k6_c5 = -5103.0 / 18656.0;

      y[i] = z[i] +
             h * (k6_c1 * k1[i] + k6_c2 * k2[i] + k6_c3 * k3[i] +
                  k6_c4 * k4[i] + k6_c5 * k5[i]);
    }
  }

  virtual void s_comp_k_set_y(void* const ls_d, const double h,
                              typename k2t<7>::type) const {
    double* const y = this->lstate_da->y(ls_d);
    double const* const z = this->lstate_da->z(ls_d);
    double const* const k1 = this->lstate_da->k1(ls_d);
    double const* const k3 = this->lstate_da->k3(ls_d);
    double const* const k4 = this->lstate_da->k4(ls_d);
    double const* const k5 = this->lstate_da->k5(ls_d);
    double* const k6 = this->lstate_da->k6(ls_d);

    this->update_nonneg(y, k6);

    for (std::size_t i = 0; i < this->s_y_len; i++) {
      static const double k7_c1 = 35.0 / 384.0;
      static const double k7_c3 = 500.0 / 1113.0;
      static const double k7_c4 = 125.0 / 192.0;
      static const double k7_c5 = -2187.0 / 6784.0;
      static const double k7_c6 = 11.0 / 84.0;

      y[i] = z[i] +
             h * (k7_c1 * k1[i] + k7_c3 * k3[i] + k7_c4 * k4[i] +
                  k7_c5 * k5[i] + k7_c6 * k6[i]);
    }
  }

  virtual void final_stage_update_nonneg(void* const ls_d) const {
    double const* const y = this->lstate_da->y(ls_d);
    double* const k7 = this->lstate_da->k7(ls_d);
    this->update_nonneg(y, k7);
  }

  virtual void s_ynew_d(double* const z, void* const ls_d, const double h,
                        const std::size_t i, const bool mbnn, double& ynew,
                        double& d) const {
    double* const y = this->lstate_da->y(ls_d);
    double const* const k1 = this->lstate_da->k1(ls_d);
    double const* const k3 = this->lstate_da->k3(ls_d);
    double const* const k4 = this->lstate_da->k4(ls_d);
    double const* const k5 = this->lstate_da->k5(ls_d);
    double const* const k6 = this->lstate_da->k6(ls_d);
    double const* const k7 = this->lstate_da->k7(ls_d);

    static const double d_c1 = 71.0 / 57600.0;
    static const double d_c3 = -71.0 / 16695.0;
    static const double d_c4 = 71.0 / 1920.0;
    static const double d_c5 = -17253.0 / 339200.0;
    static const double d_c6 = 22.0 / 525.0;
    static const double d_c7 = -1.0 / 40.0;

    static const double y_c1 = 35.0 / 384.0;
    static const double y_c3 = 500.0 / 1113.0;
    static const double y_c4 = 125.0 / 192.0;
    static const double y_c5 = -2187.0 / 6784.0;
    static const double y_c6 = 11.0 / 84.0;

    d = d_c1 * k1[i] + d_c3 * k3[i] + d_c4 * k4[i] + d_c5 * k5[i] +
        d_c6 * k6[i] + d_c7 * k7[i];

    ynew = z[i] +
           h * (y_c1 * k1[i] + y_c3 * k3[i] + y_c4 * k4[i] + y_c5 * k5[i] +
                y_c6 * k6[i]);

    this->update_nonneg(ynew, mbnn, y[i]);
  }

  virtual void comp_stages(const double h) {
    this->template comp_k<1>(h);
    this->template comp_k<2>(h);
    this->template comp_k<3>(h);
    this->template comp_k<4>(h);
    this->template comp_k<5>(h);
    this->template comp_k<6>(h);
    this->template comp_k<7>(h);
  }

  using RK23::s_comp_k_succ;
  using RK23::s_comp_k_c;
  using RK23::s_comp_k_c_raw;
  using RK23::s_comp_k_succ_cov;
  using RK23::s_comp_k_c_cov;

  void s_comp_k_succ(cstate const* const cs, double const* const x,
                     double const* const y, void const* const hsucc_d,
                     void* const ls_d, void* const lsucc_d,
                     typename k2t<5>::type) const {
    this->s_comp_succ(cs, x, y, hsucc_d, this->lstate_da->k5(ls_d),
                      this->lstate_da->k5(lsucc_d));
  }

  void s_comp_k_succ(cstate const* const cs, double const* const x,
                     double const* const y, void const* const hsucc_d,
                     void* const ls_d, void* const lsucc_d,
                     typename k2t<6>::type) const {
    this->s_comp_succ(cs, x, y, hsucc_d, this->lstate_da->k6(ls_d),
                      this->lstate_da->k6(lsucc_d));
  }

  void s_comp_k_succ(cstate const* const cs, double const* const x,
                     double const* const y, void const* const hsucc_d,
                     void* const ls_d, void* const lsucc_d,
                     typename k2t<7>::type) const {
    this->s_comp_succ(cs, x, y, hsucc_d, this->lstate_da->k7(ls_d),
                      this->lstate_da->k7(lsucc_d));
  }

  void s_comp_k_c(cstate const* const cs, double const* const x,
                  double const* const y, void* const ls_d,
                  typename k2t<5>::type) const {
    this->s_comp_c(cs, x, y, this->lstate_da->k5(ls_d));
  }

  void s_comp_k_c(cstate const* const cs, double const* const x,
                  double const* const y, void* const ls_d,
                  typename k2t<6>::type) const {
    this->s_comp_c(cs, x, y, this->lstate_da->k6(ls_d));
  }

  void s_comp_k_c(cstate const* const cs, double const* const x,
                  double const* const y, void* const ls_d,
                  typename k2t<7>::type) const {
    this->s_comp_c(cs, x, y, this->lstate_da->k7(ls_d));
  }

  void s_comp_k_c_raw(cstate const* const cs, double const* const x,
                      double const* const y, void* const ls_d,
                      typename k2t<5>::type) const {
    this->s_comp_c_raw(cs, x, y, this->lstate_da->k5(ls_d));
  }

  void s_comp_k_c_raw(cstate const* const cs, double const* const x,
                      double const* const y, void* const ls_d,
                      typename k2t<6>::type) const {
    this->s_comp_c_raw(cs, x, y, this->lstate_da->k6(ls_d));
  }

  void s_comp_k_c_raw(cstate const* const cs, double const* const x,
                      double const* const y, void* const ls_d,
                      typename k2t<7>::type) const {
    this->s_comp_c_raw(cs, x, y, this->lstate_da->k7(ls_d));
  }

  void s_comp_k_succ_cov(cstate const* const cs, double const* const x,
                         double const* const y, void const* const hsucc_d,
                         double const* const dx, void* const ls_d,
                         void* const lsucc_d, typename k2t<5>::type) const {
    this->s_comp_succ_cov(cs, x, y, hsucc_d, dx, this->lstate_da->k5(ls_d),
                          this->lstate_da->k5(lsucc_d));
  }

  void s_comp_k_succ_cov(cstate const* const cs, double const* const x,
                         double const* const y, void const* const hsucc_d,
                         double const* const dx, void* const ls_d,
                         void* const lsucc_d, typename k2t<6>::type) const {
    this->s_comp_succ_cov(cs, x, y, hsucc_d, dx, this->lstate_da->k6(ls_d),
                          this->lstate_da->k6(lsucc_d));
  }

  void s_comp_k_succ_cov(cstate const* const cs, double const* const x,
                         double const* const y, void const* const hsucc_d,
                         double const* const dx, void* const ls_d,
                         void* const lsucc_d, typename k2t<7>::type) const {
    this->s_comp_succ_cov(cs, x, y, hsucc_d, dx, this->lstate_da->k7(ls_d),
                          this->lstate_da->k7(lsucc_d));
  }

  void s_comp_k_c_cov(cstate const* const cs, double const* const x,
                      double const* const y, void* const ls_d,
                      typename k2t<5>::type) const {
    this->s_comp_c_cov(cs, x, y, this->lstate_da->k5(ls_d));
  }

  void s_comp_k_c_cov(cstate const* const cs, double const* const x,
                      double const* const y, void* const ls_d,
                      typename k2t<6>::type) const {
    this->s_comp_c_cov(cs, x, y, this->lstate_da->k6(ls_d));
  }

  void s_comp_k_c_cov(cstate const* const cs, double const* const x,
                      double const* const y, void* const ls_d,
                      typename k2t<7>::type) const {
    this->s_comp_c_cov(cs, x, y, this->lstate_da->k7(ls_d));
  }
};
}
}

#endif
