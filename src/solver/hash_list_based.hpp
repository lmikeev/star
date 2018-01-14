/*
 *  hash_list_based.hpp
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

#ifndef SOLVER_HASH_LIST_BASED_HPP_
#define SOLVER_HASH_LIST_BASED_HPP_

#include <iostream>
#include <fstream>
#include <queue>
#include "state_hash_list.hpp"

namespace solver {

template <class lstate_data_accessor_t = lstate_data_accessor,
          class hstate_data_accessor_t = hstate_data_accessor,
          class hstate_succ_data_accessor_t = hstate_succ_data_accessor,
          typename lstate_t = lstate_h, typename hstate_t = hstate_h,
          typename hstate_succ_t = hstate_succ_h>
class hashListBasedT : virtual public options {
 public:
  typedef dataAccessorT<lstate_data_accessor_t, hstate_data_accessor_t,
                        hstate_succ_data_accessor_t, lstate_t, hstate_t,
                        hstate_succ_t> dataAccessor;

  typedef stateHashListT<lstate_data_accessor_t, hstate_data_accessor_t,
                         hstate_succ_data_accessor_t> stateHashList;

  typedef lstate_t state;
  typedef hstate_t hstate;
  typedef hstate_succ_t hstate_succ;

  hashListBasedT(solver_loader::base const* const sl)
      : options(sl),
        lstate_da(new lstate_data_accessor_t(this)),
        hstate_da(new hstate_data_accessor_t(this)),
        hstate_succ_da(new hstate_succ_data_accessor_t(this)),
        lstate_allocator(new container::allocator<lstate_t>(lstate_da)),
        hstate_allocator(new container::allocator<hstate_t>(hstate_da)),
        hstate_succ_allocator(
            new container::allocator<hstate_succ_t>(hstate_succ_da)),
        da(new dataAccessor(this, lstate_da, hstate_da, hstate_succ_da,
                            lstate_allocator, hstate_allocator,
                            hstate_succ_allocator)),
        states(new stateHashList(this, lstate_da, hstate_da, hstate_succ_da,
                                 lstate_allocator, hstate_allocator,
                                 hstate_succ_allocator, da)),
        err_threshold(sl->abs_tol / sl->rel_tol),
        prob_threshold(std::min(sl->abs_tol, 1e-10)),
        tr_prob_threshold(0.010 * prob_threshold),
        cov_(xlen(nvars, nmoments)),
        dx_(npvars),
        ez_(mlt_u_len),
        ecovz_(mlt_u_len),
        rate_(mlt_u_len),
        p_rate_g_(mlt_u_len),
        sI(nmoments),
        sIk(nmoments),
        sJ(nmoments),
        sJk(nmoments),
        sK(nmoments),
        sI_(nvars),
        sIk_(nvars),
        sJ_(nvars),
        sK_(nvars) {
    printf("prob_threshold = %lg\n", prob_threshold);
    printf("tr_prob_threshold = %lg\n", tr_prob_threshold);
    printf("p_vars_len = %lu\n", sl->p_vars.size());
    printf("c_vars_len = %lu\n", sl->c_vars.size());
    printf("ls_data_size = %lu (%lu)\n", lstate_da->get_size(),
           sizeof(lstate_t));
    printf("cs_size = %lu\n", hstate_da->get_cs_size());
    printf("hs_data_size = %lu (%lu)\n", hstate_da->get_size(),
           sizeof(hstate_t));
    printf("hs_has_index = %d\n", hs_has_index());
    printf("hs_has_exitrate = %d\n", hs_has_exitrate());
    printf("hsucc_data_size = %lu (%lu)\n", hstate_succ_da->get_size(),
           sizeof(hstate_succ_t));
    printf("x_len = %lu\n", s_x_len);
    printf("y_len = %lu\n", s_y_len);

    printf("do_sa = %d\n", do_sa());
    printf("sa_no_der = %d\n", sa_no_der());
    printf("sa_no_2der = %d\n", sa_no_2der());

    printf("nmoments = %lu\n", nmoments);
    printf("det_centered = %d %d\n", is_det_centered(), sl->det_centered);

    printf("mlt_u = %lu\n", mlt_v);
    printf("mlt_du_dc = %lu\n", mlt_du_dc);
    printf("mlt_d2u_dc2 = %lu\n", mlt_d2u_dc2);
    printf("mlt_u_len = %lu\n", mlt_u_len);

    printf("f_len = %lu\n", f_len);

    printf("mlt_v = %lu\n", mlt_v);
    printf("mlt_dv_dc = %lu\n", mlt_dv_dc);
    printf("mlt_dv_di = %lu\n", mlt_dv_di);
    printf("mlt_dv_de = %lu\n", mlt_dv_de);
    printf("mlt_d2v_dc2 = %lu\n", mlt_d2v_dc2);
    printf("mlt_v_len = %lu\n", mlt_v_len);

    printf("v_rate_h_dx_len = %lu\n", v_rate_h_dx_len);
    printf("v_rate_h_offset = %lu\n", v_rate_h_offset);
    printf("v_drate_h_dx_offset = %lu\n", v_drate_h_dx_offset);
    printf("v_rate_h_len = %lu\n", v_rate_h_len);

    printf("v_rate_h_cnt = %lu\n", v_rate_h_cnt);
    printf("v_rates_h_len = %lu\n", v_rates_h_len);

    printf("tr_rates_g_len = %lu\n", tr_rates_g_len);
    printf("tr_rates_h_len = %lu\n", tr_rates_h_len);
  }

  stateHashList* create_states() const {
    return new stateHashList(this, lstate_da, hstate_da, hstate_succ_da,
                             lstate_allocator, hstate_allocator,
                             hstate_succ_allocator, da);
  }

  lstate_data_accessor_t const* get_lstate_da() const { return lstate_da; }

  dataAccessor const* get_da() const { return da; }

  stateHashList* get_states() const { return states; }

  bool can_ignore_s(const double p) const { return p < prob_threshold; }

  bool can_ignore_tr(const double p) const { return p < tr_prob_threshold; }

  double comp_err(const double d, const double z, const double ynew) const {
    return std::fabs(d) /
           std::max(std::max(std::fabs(z), std::fabs(ynew)), err_threshold);
  }

  void update_err(double& err, double& normy, double& normynew, const double d,
                  const double z, const double ynew) const {
    if (sl->norm_control) {
      err += d * d;
      normy += z * z;
      normynew += ynew * ynew;
    } else {
      err = std::max(
          err, std::fabs(d) / std::max(std::max(std::fabs(z), std::fabs(ynew)),
                                       err_threshold));
    }
  }

  void update_err(double& err, double& normy, double& normynew) const {
    if (sl->norm_control) {
      err = std::sqrt(err);
      normy = std::sqrt(normy);
      normynew = std::sqrt(normynew);

      err /= std::max(std::max(normy, normynew), err_threshold);
    }
  }

  void update_err_nn(double& err, double& normy, double& normynew,
                     bool& is_nonneg, double& nn_err, const bool mustbe_nonneg,
                     const double d, const double z, const double ynew) const {
    update_err(err, normy, normynew, d, z, ynew);

    if (mustbe_nonneg && ynew < 0.0) {
      is_nonneg = true;
      if (sl->norm_control) {
        nn_err += ynew * ynew;
      } else {
        nn_err = std::max(nn_err, -ynew / err_threshold);
      }
    }
  }

  void update_err_nn(double& err, double& normy, double& normynew,
                     double& nn_err) const {
    update_err(err, normy, normynew);

    if (sl->norm_control) {
      nn_err = std::sqrt(nn_err);

      nn_err /= std::max(std::max(normy, normynew), err_threshold);
    }
  }

  void add_init_state(const solver_loader::model_info::ic_s& s0,
                      int const* const iv) {
    states->add_init_state(s0, iv);
  }

  void load_init(solver_loader::model_info::ic const* const init) {
    states->load_init(init);
  }

  void tr_rates(transition const* const tr, cstate const* const cs,
                double const* const x, double const*& rates_g,
                double const*& rates_h) const {
    da->tr_rates(tr, cs, x, rates_g, rates_h);
  }

  void tr_rates(void const* const hsucc_d, cstate const* const cs,
                double const* const x, double const*& rates_g,
                double const*& rates_h) const {
    da->tr_rates(hsucc_d, cs, x, rates_g, rates_h);
  }

  double tr_rate(double const* const rates_g, double const* const rates_h,
                 double const* const x) const {
    double rate = v_rate_g(rates_g);
    if (!is_stoch()) {
      rate *= s_eh(x, rates_h);
    }
    return rate;
  }

  double tr_rate(void const* const hsucc_d, cstate const* const cs,
                 double const* const x) const {
    double const* rates_g;
    double const* rates_h;
    tr_rates(hsucc_d, cs, x, rates_g, rates_h);
    return tr_rate(rates_g, rates_h, x);
  }

  double tr_rate(transition const* const tr, cstate const* const cs,
                 double const* const x) const {
    double const* rates_g;
    double const* rates_h;
    tr_rates(tr, cs, x, rates_g, rates_h);
    return tr_rate(rates_g, rates_h, x);
  }

  void s_comp_succ_i(cstate const* const cs, double const* const x,
                     double const* const y, double const* const y_i,
                     void const* const hsucc_d, const double c_r,
                     const double c_i, double* const succ_dy,
                     double* const succ_dy_i) const {
    double const* rates_g;
    double const* rates_h;

    tr_rates(cs, x, hsucc_d, rates_g, rates_h);

    s_comp_succ_i(rates_g, rates_h, y, y_i, c_r, c_i, succ_dy, succ_dy_i);
  }

  void s_comp_succ_i(double const* const rates_g, double const* const,
                     double const* const y, double const* const y_i,
                     const double c_r, const double c_i, double* const succ_dy,
                     double* const succ_dy_i) const {
    const double p = v_p(y);
    const double p_i = v_p(y_i);
    const double rate_g = v_rate_g(rates_g);
    double rate = rate_g;

    const double d_p = p * rate;
    const double d_p_i = p_i * rate;

    v_p(succ_dy) += c_r * d_p + c_i * d_p_i;
    v_p(succ_dy_i) += c_r * d_p_i - c_i * d_p;
  }

  void s_comp_succ(cstate const* const cs, double const* const x,
                   double const* const y, void const* const hsucc_d,
                   double* const succ_dy) const {
    double const* rates_g;
    double const* rates_h;

    tr_rates(cs, x, hsucc_d, rates_g, rates_h);

    s_comp_succ(rates_g, rates_h, y, succ_dy);
  }

  void s_comp_succ(double const* const rates_g, double const* const rates_h,
                   double const* const y, double* const succ_dy) const {
    const double p = v_p(y);
    const double rate_g = v_rate_g(rates_g);
    double rate = rate_g;

    if (is_hybrid()) {
      rate *= v_rate_h(rates_h);
    }

    const double d_p = p * rate;

    v_p(succ_dy) += d_p;

    if (do_sa() && !sa_no_der()) {
      std::size_t kl = 0;

      for (std::size_t i = 0, ij = 0; i < nparams; i++) {
        const double drate_g_dci = v_drate_g_dc(rates_g, i);
        const double dp_dci = v_dp_dc(y, i);
        const double d_dp_dci = rate * dp_dci + drate_g_dci * p;

        v_dp_dc(succ_dy, i) += d_dp_dci;

        if (!sa_no_2der()) {
          for (std::size_t j = 0; j <= i; j++, ij++, kl++) {
            const double d_d2p_dci_dcj =
                v_drate_g_dc(rates_g, j) * dp_dci + rate * v_d2p_dc2(y, ij) +
                v_d2rate_g_dc2(rates_g, ij) * p + drate_g_dci * v_dp_dc(y, j);

            v_d2p_dc2(succ_dy, ij) += d_d2p_dci_dcj;
          }
        }
      }

      for (std::size_t i = 0; i < nivars; i++) {
        const double dp_dii = v_dp_di(y, i);
        const double d_dp_dii = rate * dp_dii;

        v_dp_di(succ_dy, i) += d_dp_dii;

        if (!sa_no_2der()) {
          for (std::size_t j = 0; j < nparams; j++, kl++) {
            const double d_d2p_dii_dcj =
                rate * v_d2p_dc2(y, kl) + v_drate_g_dc(rates_g, j) * dp_dii;

            v_d2p_dc2(succ_dy, kl) += d_d2p_dii_dcj;
          }
          for (std::size_t j = 0; j <= i; j++, kl++) {
            const double d_d2p_dii_dij = rate * v_d2p_dc2(y, kl);

            v_d2p_dc2(succ_dy, kl) += d_d2p_dii_dij;
          }
        }
      }

      for (std::size_t i = 0; i < nevars; i++) {
        const double dp_dei = v_dp_de(y, i);
        const double d_dp_dei = rate * dp_dei;

        v_dp_de(succ_dy, i) += d_dp_dei;

        if (!sa_no_2der()) {
          for (std::size_t j = 0; j < nparams; j++, kl++) {
            const double d_d2p_dei_dcj =
                rate * v_d2p_dc2(y, kl) + v_drate_g_dc(rates_g, j) * dp_dei;

            v_d2p_dc2(succ_dy, kl) += d_d2p_dei_dcj;
          }
          for (std::size_t j = 0; j < nivars; j++, kl++) {
            const double d_d2p_dei_dij = rate * v_d2p_dc2(y, kl);

            v_d2p_dc2(succ_dy, kl) += d_d2p_dei_dij;
          }
          for (std::size_t j = 0; j <= i; j++, kl++) {
            const double d_d2p_dei_dej = rate * v_d2p_dc2(y, kl);

            v_d2p_dc2(succ_dy, kl) += d_d2p_dei_dej;
          }
        }
      }
    }
  }

  void s_comp_unexplored_exitrate(state* const s) {
    hstate* const hs = s->hs;
    if (da->hs_is_unexplored(hs)) {
      cstate const* const cs = da->hs_cs(hs);
      void const* const ls_d = da->get_lstate_data(s);
      double const* const y = lstate_da->y(ls_d);
      double const* const x = da->s_x_cov(y);

      double exitrate = 0.0;
      for (auto const& tr : trs_stoch) {
        if (tr->is_enabled(cs)) {
          exitrate += tr_rate(tr, cs, x);
        }
      }
      da->hs_exitrate(hs) = exitrate;
    }
  }

 protected:
  lstate_data_accessor_t const* const lstate_da;
  hstate_data_accessor_t const* const hstate_da;
  hstate_succ_data_accessor_t const* const hstate_succ_da;

  container::allocator<lstate_t>* const lstate_allocator;
  container::allocator<hstate_t>* const hstate_allocator;
  container::allocator<hstate_succ_t>* const hstate_succ_allocator;

  dataAccessor const* const da;

  stateHashList* const states;

  const double err_threshold;
  const double prob_threshold;
  const double tr_prob_threshold;

  mutable std::vector<double> cov_;
  mutable std::vector<double> dx_;
  mutable std::vector<double> ez_;
  mutable std::vector<double> ecovz_;
  mutable std::vector<double> rate_;
  mutable std::vector<double> p_rate_g_;

  mutable std::vector<std::size_t> sI, sIk, sJ, sJk, sK;
  mutable std::vector<std::size_t> sI_, sIk_, sJ_, sK_;

  double s_ez_(double const* const x, double const* const rates_h,
               std::size_t* const K_, std::size_t* const K,
               std::size_t* const Jk_, const std::size_t d,
               const std::size_t hci0, const std::size_t ci0) const {
    double ez = 0.0;

    std::size_t Jk_len;
    cov_I_I(npvars, sJk.data(), Jk_len, Jk_);

    if (Jk_len != 1) {
      double c = 1.0;
      for (std::size_t i = 0; i < npvars; i++) {
        c *= _1_fac[K_[i]];
      }

      const double t_ = c * rates_h[hci0 + cov_indexl(K, d)];
      if (!Jk_len) {
        ez += t_;
      } else if (Jk_len > 1) {
        ez += t_ * x[ci0 + cov_indexl(sJk.data(), Jk_len)];
      }
    }

    if (Jk_len + 1 <= nmoments) {
      const std::size_t hcov_len_ = covlen(npvars, d);
      const std::size_t cov_len_ = Jk_len ? covlen(npvars, Jk_len) : 0;

      for (std::size_t i = 0; i <= K[d - 1]; i++) {
        K[d] = i;
        K_[i]++;

        Jk_[i]++;

        ez += s_ez_(x, rates_h, K_, K, Jk_, d + 1, hci0 + hcov_len_,
                    ci0 + cov_len_);

        Jk_[i]--;

        K_[i]--;
      }
    }

    return ez;
  }

  double s_ez(double const* const x, double const* const rates_h,
              std::size_t* const J_, std::size_t* const J,
              const std::size_t J_len) const {
    double ez = 0.0;
    std::size_t ci0 = J_len ? cov_index0(npvars, J_len) : 0;
    std::size_t hci0 = 0;

    if (!J_len) {
      ez += v_rate_h(rates_h);
    } else if (J_len > 1) {
      ez += v_rate_h(rates_h) * x[ci0 + cov_indexl(J, J_len)];
    }

    if (J_len + 1 <= nmoments) {
      ci0 += (J_len ? covlen(npvars, J_len) : 0);
      hci0 += 1;

      for (std::size_t i = 0; i < npvars; i++) {
        sK[0] = i;
        sK_[i]++;

        assert(sK_[i] == 1);

        J_[i]++;

        ez += s_ez_(x, rates_h, sK_.data(), sK.data(), J_, 1, hci0, ci0);

        J_[i]--;

        sK_[i]--;

        assert(sK_[i] == 0);
      }
    }

    return ez;
  }

  double s_ez(double const* const x, double const* const rates_h,
              std::size_t* const J_) const {
    std::size_t J_len;
    cov_I_I(npvars, sJ.data(), J_len, J_);
    return s_ez(x, rates_h, J_, sJ.data(), J_len);
  }

  double s_eh(double const* const x, double const* const rates_h) const {
    for (std::size_t i = 0; i < npvars; i++) {
      assert(sJ_[i] == 0);
    }
    return s_ez(x, rates_h, sJ_.data());
  }

  double s_cov_z(double const* const x, double const* const rates_h,
                 double const* const ch, std::size_t const* const I_,
                 std::size_t* const J_, const std::size_t d,
                 double const* const dx = nullptr,
                 const bool fmx = true) const {
    double z = 0.0;

    assert(J_[d] == 0);

    const std::size_t I_mx = I_[d] + ((!fmx || d + 1 < npvars) ? 1 : 0);

    for (std::size_t i = 0; i < I_mx; i++) {
      J_[d] = i;

      if (d + 1 < npvars) {
        z += s_cov_z(x, rates_h, ch, I_, J_, d + 1, dx, fmx && i == I_[d]);
      } else {
        double k = 1.0;
        for (std::size_t j = 0; j < npvars; j++) {
          k *= cmb[I_[j]][J_[j]];
        }

        if (dx == nullptr) {
          for (std::size_t j = 0; j < npvars; j++) {
            const double ch_ = ch[j];
            for (std::size_t l = 0; l < I_[j] - J_[j]; l++) {
              k *= ch_;
            }
          }
        }

        if (std::fabs(k) > 0.0) {
          if (dx != nullptr) {
            for (std::size_t j = 0; j < npvars; j++) {
              const double ch_ = ch[j];
              for (std::size_t l = 0; l < I_[j] - J_[j]; l++) {
                k *= dx[j] + ch_;
              }
            }
          }
          const double ez = s_ez(x, rates_h, J_);
          z += k * ez;
        }
      }
    }
    J_[d] = 0;

    return z;
  }

  void s_ez_raw_(double const* const cov, double const* const rates_h,
                 std::size_t* const K_, std::size_t* const K,
                 const std::size_t d, const std::size_t hci0,
                 const std::size_t ci0, double* const ez) const {
    if (d > 1) {
      double c = 1.0;
      for (std::size_t i = 0; i < npvars; i++) {
        c *= _1_fac[K_[i]];
      }
      const std::size_t ci_ = cov_indexl(K, d);

      v_f(ez) += c * rates_h[hci0 + ci_] * cov[ci0 + ci_];

      if (do_sa() && !sa_no_der()) {
        for (std::size_t j = 0, jk = 0; j < nparams; j++) {
          v_df_dc(ez)[j] +=
              c * (v_drate_h_dc(rates_h, j)[hci0 + ci_] * cov[ci0 + ci_] +
                   rates_h[hci0 + ci_] *
                       cov[(mlt_du_dc + j) * s_x_len + ci0 + ci_]);

          if (!sa_no_2der()) {
            for (std::size_t k = 0; k <= j; k++, jk++) {
              v_d2f_dc2(ez)[jk] +=
                  c *
                  (v_drate_h_dc(rates_h, j)[hci0 + ci_] *
                       cov[(mlt_du_dc + k) * s_x_len + ci0 + ci_] +
                   v_d2rate_h_dc2(rates_h, jk)[hci0 + ci_] * cov[ci0 + ci_] +
                   v_drate_h_dc(rates_h, k)[hci0 + ci_] *
                       cov[(mlt_du_dc + j) * s_x_len + ci0 + ci_] +
                   rates_h[hci0 + ci_] *
                       cov[(mlt_d2u_dc2 + jk) * s_x_len + ci0 + ci_]);
            }
          }
        }
      }
    }

    if (d + 1 <= nmoments) {
      const std::size_t cov_len_ = covlen(npvars, d);

      for (std::size_t i = 0; i <= K[d - 1]; i++) {
        K[d] = i;
        K_[i]++;

        s_ez_raw_(cov, rates_h, K_, K, d + 1, hci0 + cov_len_, ci0 + cov_len_,
                  ez);

        K_[i]--;
      }
    }
  }

  double* s_ez_raw(double const* const cov, double const* const rates_h,
                   std::size_t* const, std::size_t* const J,
                   const std::size_t J_len) const {
    std::size_t hci0 =
        (J_len ? 1 + cov_indexg(npvars, J, J_len) : 0) * v_rate_h_len;

    v_f(ez_.data()) = rates_h[hci0];

    if (do_sa() && !sa_no_der()) {
      for (std::size_t j = 0, jk = 0; j < nparams; j++) {
        v_df_dc(ez_.data())[j] = v_drate_h_dc(rates_h, j)[hci0];

        if (!sa_no_2der()) {
          for (std::size_t k = 0; k <= j; k++, jk++) {
            v_d2f_dc2(ez_.data())[jk] = v_d2rate_h_dc2(rates_h, jk)[hci0];
          }
        }
      }
    }

    if (nmoments > 1) {
      hci0 += 1;

      for (std::size_t i = 0; i < npvars; i++) {
        sK[0] = i;
        sK_[i]++;

        assert(sK_[i] == 1);

        s_ez_raw_(cov, rates_h, sK_.data(), sK.data(), 1, hci0, 0, ez_.data());

        sK_[i]--;

        assert(sK_[i] == 0);
      }
    }

    return ez_.data();
  }

  double* s_ez_raw(double const* const cov, double const* const rates_h,
                   std::size_t* const J_) const {
    std::size_t J_len;
    cov_I_I(npvars, sJ.data(), J_len, J_);
    return s_ez_raw(cov, rates_h, J_, sJ.data(), J_len);
  }

  double* s_eh_raw(double const* const cov, double const* const rates_h) const {
    for (std::size_t i = 0; i < npvars; i++) {
      assert(sJ_[i] == 0);
    }
    return s_ez_raw(cov, rates_h, sJ_.data());
  }

  void s_cov_z_raw(double const* const cov, double const* const rates_h,
                   double const* const ch, std::size_t const* const I_,
                   std::size_t* const J_, const std::size_t d, const bool fmx,
                   double* const z) const {
    assert(J_[d] == 0);

    const std::size_t I_mx = I_[d] + ((!fmx || d + 1 < npvars) ? 1 : 0);

    for (std::size_t i = 0; i < I_mx; i++) {
      J_[d] = i;

      if (d + 1 < npvars) {
        s_cov_z_raw(cov, rates_h, ch, I_, J_, d + 1, fmx && i == I_[d], z);
      } else {
        double k = 1.0;
        for (std::size_t j = 0; j < npvars; j++) {
          k *= cmb[I_[j]][J_[j]];
          const double ch_ = ch[j];
          for (std::size_t l = 0; l < I_[j] - J_[j]; l++) {
            k *= ch_;
          }
        }
        if (std::fabs(k) > 0.0) {
          double const* const ez = s_ez_raw(cov, rates_h, J_);
          for (std::size_t j = 0; j < mlt_u_len; j++) {
            z[j] += k * ez[j];
          }
        }
      }
    }
    J_[d] = 0;
  }

  double* s_cov_z_raw(double const* const cov, double const* const rates_h,
                      double const* const ch, std::size_t const* const I_,
                      std::size_t* const J_, const std::size_t d,
                      const bool fmx = true) const {
    std::fill(ecovz_.begin(), ecovz_.end(), 0.0);
    s_cov_z_raw(cov, rates_h, ch, I_, J_, d, fmx, ecovz_.data());
    return ecovz_.data();
  }

  void s_comp_succ(cstate const* const cs, double const* const x,
                   double const* const y, void const* const hsucc_d,
                   double* const s_dy, double* const succ_dy) const {
    double const* rates_g;
    double const* rates_h;
    tr_rates(hsucc_d, cs, x, rates_g, rates_h);
    s_comp_succ(rates_g, rates_h, x, y, hsucc_d, s_dy, succ_dy);
  }

  void s_comp_succ(double const* const rates_g, double const* const rates_h,
                   double const* const x, double const* const y,
                   void const* const hsucc_d, double* const s_dy,
                   double* const succ_dy) const {
    const double p = v_p(y);

    v_f(rate_.data()) = v_rate_g(rates_g);

    v_f(p_rate_g_.data()) = p * v_rate_g(rates_g);

    if (do_sa() && !sa_no_der()) {
      for (std::size_t j = 0, jk = 0; j < nparams; j++) {
        v_df_dc(rate_.data())[j] = v_drate_g_dc(rates_g, j);

        v_df_dc(p_rate_g_.data())[j] =
            p * v_drate_g_dc(rates_g, j) + v_dp_dc(y, j) * v_rate_g(rates_g);

        if (!sa_no_2der()) {
          for (std::size_t k = 0; k <= j; k++, jk++) {
            v_d2f_dc2(rate_.data())[jk] = v_d2rate_g_dc2(rates_g, jk);

            v_d2f_dc2(p_rate_g_.data())[jk] =
                p * v_d2rate_g_dc2(rates_g, jk) +
                v_dp_dc(y, k) * v_drate_g_dc(rates_g, j) +
                v_dp_dc(y, j) * v_drate_g_dc(rates_g, k) +
                v_d2p_dc2(y, jk) * v_rate_g(rates_g);
          }
        }
      }
    }

    double eh;
    if (is_hybrid()) {
      eh = s_eh(x, rates_h);

      v_f(rate_.data()) *= eh;
    }

    const double d_p = p * v_f(rate_.data());

    v_p(s_dy) -= d_p;
    v_p(succ_dy) += d_p;

    if (do_rare_event()) {
      const double p2 = v_p2(y);
      const double rate2 = v_rate2_g(rates_g);

      const double d_p2 = p2 * rate2;

      v_p2(s_dy) -= d_p2;
      v_p2(succ_dy) += d_p2;
    }

    if (do_sa() && !sa_no_der()) {
      std::size_t kl = 0;

      for (std::size_t i = 0, ij = 0; i < nparams; i++) {
        const double dp_dci = v_dp_dc(y, i);
        const double drate_dci = v_df_dc(rate_.data())[i];
        const double d_dp_dci = v_f(rate_.data()) * dp_dci + drate_dci * p;

        v_dp_dc(s_dy, i) -= d_dp_dci;
        v_dp_dc(succ_dy, i) += d_dp_dci;

        if (!sa_no_2der()) {
          for (std::size_t j = 0; j <= i; j++, ij++, kl++) {
            const double d_d2p_dci_dcj = v_df_dc(rate_.data())[j] * dp_dci +
                                         v_f(rate_.data()) * v_d2p_dc2(y, ij) +
                                         v_d2rate_g_dc2(rates_g, ij) * p +
                                         drate_dci * v_dp_dc(y, j);

            v_d2p_dc2(s_dy, kl) -= d_d2p_dci_dcj;
            v_d2p_dc2(succ_dy, kl) += d_d2p_dci_dcj;
          }
        }
      }

      for (std::size_t i = 0; i < nivars; i++) {
        const double dp_dii = v_dp_di(y, i);
        const double d_dp_dii = v_f(rate_.data()) * dp_dii;

        v_dp_di(s_dy, i) -= d_dp_dii;
        v_dp_di(succ_dy, i) += d_dp_dii;

        if (!sa_no_2der()) {
          for (std::size_t j = 0; j < nparams; j++, kl++) {
            const double d_d2p_dii_dcj = v_f(rate_.data()) * v_d2p_dc2(y, kl) +
                                         v_df_dc(rate_.data())[j] * dp_dii;

            v_d2p_dc2(s_dy, kl) -= d_d2p_dii_dcj;
            v_d2p_dc2(succ_dy, kl) += d_d2p_dii_dcj;
          }
          for (std::size_t j = 0; j <= i; j++, kl++) {
            const double d_d2p_dii_dij = v_f(rate_.data()) * v_d2p_dc2(y, kl);

            v_d2p_dc2(s_dy, kl) -= d_d2p_dii_dij;
            v_d2p_dc2(succ_dy, kl) += d_d2p_dii_dij;
          }
        }
      }

      for (std::size_t i = 0; i < nevars; i++) {
        const double dp_dei = v_dp_de(y, i);
        const double d_dp_dei = v_f(rate_.data()) * dp_dei;

        v_dp_de(s_dy, i) -= d_dp_dei;
        v_dp_de(succ_dy, i) += d_dp_dei;

        if (!sa_no_2der()) {
          for (std::size_t j = 0; j < nparams; j++, kl++) {
            const double d_d2p_dei_dcj = v_f(rate_.data()) * v_d2p_dc2(y, kl) +
                                         v_df_dc(rate_.data())[j] * dp_dei;

            v_d2p_dc2(s_dy, kl) -= d_d2p_dei_dcj;
            v_d2p_dc2(succ_dy, kl) += d_d2p_dei_dcj;
          }
          for (std::size_t j = 0; j < nivars; j++, kl++) {
            const double d_d2p_dei_dij = v_f(rate_.data()) * v_d2p_dc2(y, kl);

            v_d2p_dc2(s_dy, kl) -= d_d2p_dei_dij;
            v_d2p_dc2(succ_dy, kl) += d_d2p_dei_dij;
          }
          for (std::size_t j = 0; j <= i; j++, kl++) {
            const double d_d2p_dei_dej = v_f(rate_.data()) * v_d2p_dc2(y, kl);

            v_d2p_dc2(s_dy, kl) -= d_d2p_dei_dej;
            v_d2p_dc2(succ_dy, kl) += d_d2p_dei_dej;
          }
        }
      }
    }

    if (is_hybrid()) {
      transition const* const tr =
          sl->get_model()->get_transitions()[da->hsucc_tri(hsucc_d)];

      if (this->is_det_centered()) {
        for (std::size_t i = 0; i < npvars; i++) {
          sI[0] = i;
          sI_[i]++;

          assert(sI_[i] == 1);

          const double ez = s_ez(x, rates_h, sI_.data(), sI.data(), 1);

          const double dpx_ = v_f(p_rate_g_.data()) * (ez + x[i] * eh);

          const double ch = tr->change()[i];

          v_px(s_dy)[i] -= dpx_;
          v_px(succ_dy)[i] += dpx_ + ch * v_f(p_rate_g_.data()) * eh;

          sI_[i]--;

          assert(sI_[i] == 0);
        }
      } else {
        double const* const ch = tr->change();

        for (std::size_t i = 0; i < npvars; i++) {
          sI[0] = i;
          sI_[i]++;

          assert(sI_[i] == 1);

          double const* const ez =
              s_ez_raw(x, rates_h, sI_.data(), sI.data(), 1);

          v_px(s_dy)[i] -= v_f(p_rate_g_.data()) * v_f(ez);
          v_px(succ_dy)[i] += v_f(p_rate_g_.data()) * (v_f(ez) + ch[i] * eh);

          if (nmoments > 1) {
            s_dcov_raw(s_dy, succ_dy, x, p_rate_g_.data(), rates_h, ch,
                       sI.data(), sI_.data(), 1, v_px_offset + npvars);
          }

          sI_[i]--;

          assert(sI_[i] == 0);
        }
      }

      if (do_sa() && !sa_no_der()) {
      }
    }
  }

  void s_comp_c(cstate const* const cs, double const* const x,
                double const* const y, double* const s_dy,
                const double t = 0.0) const {
    const double p = v_p(y);
    for (auto const& tr : trs_det) {
      if (tr->is_enabled(cs, t)) {
        double const* const ch = tr->change();
        const double dpx_ = p * tr_rate(tr, cs, x);
        for (std::size_t i = 0; i < npvars; i++) {
          v_px(s_dy)[i] += ch[i] * dpx_;
        }
      }
    }
  }

  void s_comp_c_raw(cstate const* const cs, double const* const x,
                    double const* const y, double* const s_dy,
                    const double t = 0.0) const {
    const double p = v_p(y);

    for (auto const& tr : trs_det) {
      if (tr->is_enabled(cs, t)) {
        double const* rates_g;
        double const* rates_h;
        tr_rates(tr, cs, x, rates_g, rates_h);
        const double p_rate_g = p * v_rate_g(rates_g);
        const double dpx_ = p * tr_rate(rates_g, rates_h, x);
        double const* const ch = tr->change();

        for (std::size_t i = 0; i < npvars; i++) {
          sI[0] = i;
          sI_[i]++;

          assert(sI_[i] == 1);

          v_px(s_dy)[i] += ch[i] * dpx_;

          if (nmoments > 1) {
            s_dcov_c_raw(s_dy, x, p_rate_g, rates_h, ch, sI.data(), sI_.data(),
                         1, v_px_offset + npvars);
          }

          sI_[i]--;

          assert(sI_[i] == 0);
        }
      }
    }
  }

  void s_dcov_c_raw(double* const s_dy, double const* const x,
                    const double p_rate_g, double const* const rates_h,
                    double const* const ch, std::size_t* const I,
                    std::size_t* const I_, const std::size_t d,
                    const std::size_t ci0) const {
    const std::size_t cov_len_ = covlen(npvars, d + 1);

    for (std::size_t i = 0; i <= I[d - 1]; i++) {
      I[d] = i;
      I_[i]++;

      const std::size_t ci = ci0 + cov_indexl(I, d + 1);

      double const* const ecovz =
          s_cov_z_raw(x, rates_h, ch, I_, sJ_.data(), 0);

      s_dy[ci] += p_rate_g * v_f(ecovz);

      if (d + 1 < nmoments) {
        s_dcov_c_raw(s_dy, x, p_rate_g, rates_h, ch, I, I_, d + 1,
                     ci0 + cov_len_);
      }

      I_[i]--;
    }
  }

  void s_dcov_raw(double* const s_dy, double* const succ_dy,
                  double const* const x, double const* const p_rate_g,
                  double const* const rates_h, double const* const ch,
                  std::size_t* const I, std::size_t* const I_,
                  const std::size_t d, const std::size_t ci0) const {
    const std::size_t cov_len_ = covlen(npvars, d + 1);

    for (std::size_t i = 0; i <= I[d - 1]; i++) {
      I[d] = i;
      I_[i]++;

      const std::size_t ci = cov_indexl(I, d + 1);

      double const* const ez = s_ez_raw(x, rates_h, I_, I, d + 1);

      s_dy[ci0 + ci] -= v_f(p_rate_g) * v_f(ez);

      double const* const ecovz =
          s_cov_z_raw(x, rates_h, ch, I_, sJ_.data(), 0, false);

      succ_dy[ci0 + ci] += v_f(p_rate_g) * v_f(ecovz);

      if (d + 1 < nmoments) {
        s_dcov_raw(s_dy, succ_dy, x, p_rate_g, rates_h, ch, I, I_, d + 1,
                   ci0 + cov_len_);
      }

      I_[i]--;
    }
  }

  void s_comp_succ_cov(cstate const* const cs, double const* const x,
                       double const* const y, void const* const hsucc_d,
                       double const* const dx, double* const s_dy,
                       double* const succ_dy) const {
    double const* rates_g;
    double const* rates_h;

    tr_rates(hsucc_d, cs, x, rates_g, rates_h);

    const double p_rate_g = v_p(y) * v_rate_g(rates_g);

    double const* const ch = da->tr_change(hsucc_d);

    for (std::size_t i = 0; i < npvars; i++) {
      sI[0] = i;
      sI_[i]++;

      assert(sI_[i] == 1);

      s_dcov(s_dy, succ_dy, x, dx, p_rate_g, rates_h, ch, sI.data(), sI_.data(),
             1, v_px_offset + npvars);

      sI_[i]--;

      assert(sI_[i] == 0);
    }
  }

  void s_dcov(double* const s_dy, double* const succ_dy, double const* const x,
              double const* const dx, const double p_rate_g,
              double const* const rates_h, double const* const ch,
              std::size_t* const I, std::size_t* const I_, const std::size_t d,
              const std::size_t ci0) const {
    const std::size_t cov_len_ = covlen(npvars, d + 1);

    for (std::size_t i = 0; i <= I[d - 1]; i++) {
      I[d] = i;
      I_[i]++;

      const std::size_t ci = cov_indexl(I, d + 1);

      s_dy[ci0 + ci] -= p_rate_g * s_ez(x, rates_h, I_, I, d + 1);

      succ_dy[ci0 + ci] +=
          p_rate_g * s_cov_z(x, rates_h, ch, I_, sJ_.data(), 0, dx, false);

      if (d + 1 < nmoments) {
        s_dcov(s_dy, succ_dy, x, dx, p_rate_g, rates_h, ch, I, I_, d + 1,
               ci0 + cov_len_);
      }

      I_[i]--;
    }
  }

  void s_dcov_2(double* s_dy, double const* const x, const std::size_t vi,
                std::size_t* const I, std::size_t* const I_,
                const std::size_t d, const std::size_t ci0) const {
    const std::size_t cov_len_ = covlen(npvars, d + 1);

    for (std::size_t i = 0; i <= I[d - 1]; i++) {
      I[d] = i;
      I_[i]++;

      if (I_[vi]) {
        I_[vi]--;

        std::size_t Ik_len;
        cov_I_I(npvars, sIk.data(), Ik_len, I_);

        I_[vi]++;

        if (Ik_len != 1) {
          double dy_ = I_[vi] * (v_px(s_dy)[vi] - x[vi] * v_p(s_dy));

          if (Ik_len > 1) {
            dy_ *= x[cov_indexg(npvars, sIk.data(), Ik_len)];
          }

          s_dy[ci0 + cov_indexl(I, d + 1)] -= dy_;
        }
      }

      if (d + 1 < nmoments) {
        s_dcov_2(s_dy, x, vi, I, I_, d + 1, ci0 + cov_len_);
      }

      I_[i]--;
    }
  }

  void s_dcov_c(double* const s_dy, double const* const x,
                const double p_rate_g, double const* const rates_h,
                double const* const ch, std::size_t* const I,
                std::size_t* const I_, const std::size_t d,
                const std::size_t ci0) const {
    const std::size_t cov_len_ = covlen(npvars, d + 1);

    for (std::size_t i = 0; i <= I[d - 1]; i++) {
      I[d] = i;
      I_[i]++;

      s_dy[ci0 + cov_indexl(I, d + 1)] +=
          p_rate_g * s_cov_z(x, rates_h, ch, I_, sJ_.data(), 0);

      if (d + 1 < nmoments) {
        s_dcov_c(s_dy, x, p_rate_g, rates_h, ch, I, I_, d + 1, ci0 + cov_len_);
      }

      I_[i]--;
    }
  }

  void s_comp_c_cov(cstate const* const cs, double const* const x,
                    double const* const y, double* const s_dy,
                    const double t = 0.0) const {
    const double p = v_p(y);

    for (auto const& tr : trs_det) {
      if (tr->is_enabled(cs, t)) {
        double const* rates_g;
        double const* rates_h;
        tr_rates(tr, cs, x, rates_g, rates_h);

        const double p_rate_g = p * v_rate_g(rates_g);
        double const* const ch = tr->change();

        for (std::size_t i = 0; i < npvars; i++) {
          sI[0] = i;
          sI_[i]++;

          assert(sI_[i] == 1);

          s_dcov_c(s_dy, x, p_rate_g, rates_h, ch, sI.data(), sI_.data(), 1,
                   v_px_offset + npvars);

          sI_[i]--;

          assert(sI_[i] == 0);
        }
      }
    }

    for (std::size_t i = 0; i < npvars; i++) {
      sI[0] = i;
      sI_[i]++;

      assert(sI_[i] == 1);

      for (std::size_t j = 0; j < npvars; j++) {
        s_dcov_2(s_dy, x, j, sI.data(), sI_.data(), 1, v_px_offset + npvars);
      }

      sI_[i]--;

      assert(sI_[i] == 0);
    }
  }

  void center_moments(double const* const x, double* const cov) const {
    for (std::size_t i = 0; i < npvars; i++) {
      cov[i] = x[i];

      if (nmoments > 1) {
        sI[0] = i;
        sI_[i]++;

        assert(sI_[i] == 1);

        states->center_moments_(x, cov, nmoments, npvars, sI.data(), sI_.data(),
                                1, npvars);

        sI_[i]--;

        assert(sI_[i] == 0);
      }
    }
  }

  void read_states_y(const double* const y) const {
    double const* z = y;
    for (state* s = states->get_first(); s != nullptr; s = s->lnext) {
      double* const s_y = lstate_da->y(da->get_lstate_data(s));

      memcpy(s_y, z, s_y_len * sizeof(double));
      z += s_y_len;
    }
  }

  void write_states_y(double* const y) const {
    double* z = y;
    for (state* s = states->get_first(); s != nullptr; s = s->lnext) {
      double const* const s_y = lstate_da->y(da->get_lstate_data(s));

      memcpy(z, s_y, s_y_len * sizeof(double));
      z += s_y_len;
    }
  }

  void write_states_dy(const double t, double const* const y,
                       double* const dy) const {
    read_states_y(y);

    const std::size_t nmodes = states->get_nactive();
    const std::size_t neqs = nmodes * s_y_len;
    memset(dy, 0, neqs * sizeof(double));

    for (state* s = states->get_first(); s != nullptr; s = s->lnext) {
      hstate const* const hs = s->hs;
      void const* const hs_d = da->get_hstate_data(hs);
      cstate const* const cs = da->hs_d_cs(hs_d);
      const std::size_t s_i = da->hs_d_index(hs_d);
      double const* const y = lstate_da->y(da->get_lstate_data(s));
      double const* const x = da->s_x_cov(y);

      if (!is_stoch() && !is_det_centered()) {
        this->center_moments(x, cov_.data());
      }

      double* const s_dy = dy + s_i * s_y_len;

      if (!da->hs_is_unexplored(hs)) {
        hstate_succ const* hsucc = hs->succ;
        while (hsucc != nullptr) {
          const std::size_t succ_i = da->hs_index(hsucc->hs);
          void const* const hsucc_d = da->get_hstate_succ_data(hsucc);

          double* const succ_dy = dy + succ_i * s_y_len;

          s_comp_succ(cs, is_det_centered() ? x : cov_.data(), y, hsucc_d, s_dy,
                      succ_dy);

          hsucc = hsucc->lnext;
        }
      }

      if (!is_stoch()) {
        if (is_det_centered()) {
          s_comp_c(cs, x, y, s_dy, t);
        } else {
          s_comp_c_raw(cs, cov_.data(), y, s_dy, t);
        }
      }
    }

    if (!is_stoch() && is_det_centered() && nmoments > 1) {
      for (state* s = states->get_first(); s != nullptr; s = s->lnext) {
        hstate const* const hs = s->hs;
        void const* const hs_d = da->get_hstate_data(hs);
        cstate const* const cs = da->hs_d_cs(hs_d);
        const std::size_t s_i = da->hs_d_index(hs_d);
        double const* const y = lstate_da->y(da->get_lstate_data(s));
        double const* const x = da->s_x_cov(y);

        double* const s_dy = dy + s_i * s_y_len;

        if (!da->hs_is_unexplored(hs)) {
          hstate_succ const* hsucc = hs->succ;
          while (hsucc != nullptr) {
            const std::size_t succ_i = da->hs_index(hsucc->hs);
            void const* const hsucc_d = da->get_hstate_succ_data(hsucc);
            void const* const lsucc_d = da->get_lstate_data(hsucc->hs->ls);
            double const* const y_succ = lstate_da->y(lsucc_d);
            const double p_succ = v_p(y_succ);

            double const* dx;
            if (p_succ > 0.0) {
              double const* const px_succ = v_px(y_succ);
              for (std::size_t i = 0; i < npvars; i++) {
                dx_[i] = x[i] - px_succ[i] / p_succ;
              }
              dx = dx_.data();
            } else {
              dx = x;
            }

            double* const succ_dy = dy + succ_i * s_y_len;

            s_comp_succ_cov(cs, x, y, hsucc_d, dx, s_dy, succ_dy);

            hsucc = hsucc->lnext;
          }
        }

        s_comp_c_cov(cs, x, y, s_dy, t);
      }
    }
  }

  void write_states_d2y(const double, double const* const y,
                        double* const d2y) const {
    read_states_y(y);

    const std::size_t neqs = states->get_nactive() * s_y_len;
    memset(d2y, 0, neqs * neqs * sizeof(double));

    assert(false);
  }
};
}

#endif
