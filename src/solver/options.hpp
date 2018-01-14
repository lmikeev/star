/*
 *  options.hpp
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

#ifndef SOLVER_OPTIONS_HPP_
#define SOLVER_OPTIONS_HPP_

#include <cstddef>
#include <functional>
#include "model.hpp"
#include "../solver_loader/base.hpp"

namespace solver {

class options {
 public:
  options(solver_loader::base const* const sl)
      : sl(sl),

        nvars(sl->get_vars().size()),
        npvars(sl->p_vars.size()),
        ncvars(sl->c_vars.size()),

        ntransition_classes(sl->get_transitions().size()),

        nparams(sl->sa_params.size()),
        nivars(sl->sa_ivars.size()),
        nevars(sl->sa_evars.size()),

        nlparams(nparams + nivars + nevars),

        param_vals(std::vector<double>(nparams)),
        init_vals(std::vector<double>(nivars)),
        obserr_vals(std::vector<double>(nevars)),

        nmoments(sl->nmoments),

        s_x_len(xlen_c(npvars, nmoments)),

        mlt_v(0),
        mlt_dv_dc(mlt_v + 1),
        mlt_dv_di(mlt_dv_dc + ((do_sa() && !sa_no_der()) ? nparams : 0)),
        mlt_dv_de(mlt_dv_di + ((do_sa() && !sa_no_der()) ? nivars : 0)),
        mlt_d2v_dc2(mlt_dv_de + ((do_sa() && !sa_no_der()) ? nevars : 0)),
        mlt_v_len(mlt_d2v_dc2 +
                  ((do_sa() && !sa_no_2der()) ? PAR2_LEN(nlparams) : 0)),

        v_p_offset(0),
        v_p2_offset(v_p_offset + 1),
        v_px_offset(v_p2_offset + (do_rare_event() ? 1 : 0)),
        v_p_len(v_px_offset + s_x_len),

        s_y_len(v_p_len * mlt_v_len),

        mlt_u(0),
        mlt_du_dc(mlt_u + 1),
        mlt_d2u_dc2(mlt_du_dc + ((do_sa() && !sa_no_der()) ? nparams : 0)),
        mlt_u_len(mlt_d2u_dc2 +
                  ((do_sa() && !sa_no_2der()) ? PAR2_LEN(nparams) : 0)),

        v_rate_g_offset(0),
        v_rate2_g_offset(v_rate_g_offset + 1),
        v_rate_g_len(v_rate2_g_offset + (do_rare_event() ? 1 : 0)),

        tr_rates_g_len(v_rate_g_len * mlt_u_len),

        v_rate_h_offset(0),
        v_drate_h_dx_offset(v_rate_h_offset + 1),

        v_rate_h_dx_len((nmoments > 1) ? s_x_len : 0),

        v_rate_h_len(v_drate_h_dx_offset + v_rate_h_dx_len),

        v_rate_h_cnt(1 + ((!is_det_centered()) ? s_x_len : 0)),

        v_rates_h_len(v_rate_h_len * v_rate_h_cnt),

        tr_rates_h_len(v_rates_h_len * mlt_u_len),

        f_len(mlt_v_len),

        cmb_len(std::max(nvars + nmoments, nmoments + 1)),
        cmb(std::vector<std::vector<std::size_t>>(
            cmb_len, std::vector<std::size_t>(cmb_len, 0))),

        fac(std::vector<uint64_t>(nmoments + 1)),
        _1_fac(std::vector<double>(nmoments + 1)),

        exprs(std::vector<std::function<double(
            cstate const* const, double const* const, double const* const)>>(
            sl->get_exprs().size())),

        csvSep(',') {
    for (std::size_t i = 0; i < cmb_len; i++) {
      cmb[i][0] = 1;
    }
    for (std::size_t i = 1; i < cmb_len; i++) {
      for (std::size_t j = 1; j <= i; j++) {
        cmb[i][j] = cmb[i - 1][j - 1] + cmb[i - 1][j];
      }
    }

    fac[0] = 1;
    fac[1] = 1;
    _1_fac[0] = 1.0;
    _1_fac[1] = 1.0;
    for (std::size_t i = 2; i <= nmoments; i++) {
      fac[i] = fac[i - 1] * i;
      _1_fac[i] = _1_fac[i - 1] / (double)i;
    }
  }

  int CI(const std::size_t i, const std::size_t j) const {
    assert(i >= j);
    return (((i + 1) * i) >> 1) + j;
  }

  int CIu(const std::size_t i, const std::size_t j) const {
    return (i >= j) ? CI(i, j) : CI(j, i);
  }

  double& v_p(double* const v) const { return v[v_p_offset]; }

  double v_p(double const* const v) const { return v[v_p_offset]; }

  double& v_p2(double* const v) const { return v[v_p2_offset]; }

  double v_p2(double const* const v) const { return v[v_p2_offset]; }

  double* v_px(double* const v) const { return &v[v_px_offset]; }

  double const* v_px(double const* const v) const { return &v[v_px_offset]; }

  double& v_dp_dc(double* const v, const std::size_t ci) const {
    return v[(mlt_dv_dc + ci) * v_p_len + v_p_offset];
  }

  double v_dp_dc(double const* const v, const std::size_t ci) const {
    return v[(mlt_dv_dc + ci) * v_p_len + v_p_offset];
  }

  double& v_dp_di(double* const v, const std::size_t ii) const {
    return v[(mlt_dv_di + ii) * v_p_len + v_p_offset];
  }

  double v_dp_di(double const* const v, const std::size_t ii) const {
    return v[(mlt_dv_di + ii) * v_p_len + v_p_offset];
  }

  double& v_dp_de(double* const v, const std::size_t ei) const {
    return v[(mlt_dv_de + ei) * v_p_len + v_p_offset];
  }

  double v_dp_de(double const* const v, const std::size_t ei) const {
    return v[(mlt_dv_de + ei) * v_p_len + v_p_offset];
  }

  double& v_d2p_dc2(double* const v, const std::size_t cij) const {
    return v[(mlt_d2v_dc2 + cij) * v_p_len + v_p_offset];
  }

  double v_d2p_dc2(double const* const v, const std::size_t cij) const {
    return v[(mlt_d2v_dc2 + cij) * v_p_len + v_p_offset];
  }

  double* v_dpx_dc(double* const v, const std::size_t ci) const {
    return v + (mlt_dv_dc + ci) * v_p_len + v_px_offset;
  }

  double const* v_dpx_dc(double const* const v, const std::size_t ci) const {
    return v + (mlt_dv_dc + ci) * v_p_len + v_px_offset;
  }

  double* v_d2px_dc2(double* const v, const std::size_t cij) const {
    return v + (mlt_d2v_dc2 + cij) * v_p_len + v_px_offset;
  }

  double const* v_d2px_dc2(double const* const v, const std::size_t cij) const {
    return v + (mlt_d2v_dc2 + cij) * v_p_len + v_px_offset;
  }

  double& v_rate_g(double* const v) const { return v[v_rate_g_offset]; }

  double v_rate_g(double const* const v) const { return v[v_rate_g_offset]; }

  double& v_rate2_g(double* const v) const { return v[v_rate2_g_offset]; }

  double v_rate2_g(double const* const v) const { return v[v_rate2_g_offset]; }

  double& v_drate_g_dc(double* const v, const std::size_t ci) const {
    return v[(mlt_du_dc + ci) * v_rate_g_len + v_rate_g_offset];
  }

  double v_drate_g_dc(double const* const v, const std::size_t ci) const {
    return v[(mlt_du_dc + ci) * v_rate_g_len + v_rate_g_offset];
  }

  double& v_d2rate_g_dc2(double* const v, const std::size_t cij) const {
    return v[(mlt_d2u_dc2 + cij) * v_rate_g_len + v_rate_g_offset];
  }

  double v_d2rate_g_dc2(double const* const v, const std::size_t cij) const {
    return v[(mlt_d2u_dc2 + cij) * v_rate_g_len + v_rate_g_offset];
  }

  double& v_rate_h(double* const v) const { return v[v_rate_h_offset]; }

  double v_rate_h(double const* const v) const { return v[v_rate_h_offset]; }

  double* v_drate_h_dx(double* const v) const {
    return v + v_drate_h_dx_offset;
  }

  double const* v_drate_h_dx(double const* const v) const {
    return v + v_drate_h_dx_offset;
  }

  double* v_drate_h_dc(double* const v, const std::size_t ci) const {
    return v + (mlt_du_dc + ci) * v_rates_h_len + v_rate_h_offset;
  }

  double const* v_drate_h_dc(double const* const v,
                             const std::size_t ci) const {
    return v + (mlt_du_dc + ci) * v_rates_h_len + v_rate_h_offset;
  }

  double* v_d2rate_h_dc2(double* const v, const std::size_t cij) const {
    return v + (mlt_d2u_dc2 + cij) * v_rates_h_len + v_rate_h_offset;
  }

  double const* v_d2rate_h_dc2(double const* const v,
                               const std::size_t cij) const {
    return v + (mlt_d2u_dc2 + cij) * v_rates_h_len + v_rate_h_offset;
  }

  double& v_f(double* const v) const { return v[mlt_u]; }

  double v_f(double const* const v) const { return v[mlt_u]; }

  double* v_df_dc(double* const v) const { return &v[mlt_du_dc]; }

  double const* v_df_dc(double const* const v) const { return &v[mlt_du_dc]; }

  double* v_d2f_dc2(double* const v) const { return &v[mlt_d2u_dc2]; }

  double const* v_d2f_dc2(double const* const v) const {
    return &v[mlt_d2u_dc2];
  }

  std::size_t get_s_x_len() const { return s_x_len; }

  std::size_t get_s_y_len() const { return s_y_len; }

  std::size_t get_tr_rates_g_len() const { return tr_rates_g_len; }

  std::size_t get_tr_rates_h_len() const { return tr_rates_h_len; }

  bool is_det() const {
    return
#ifdef STAR_CODEGEN
#ifdef STAR_DET
        true
#else
        false
#endif
#else
        sl->is_det()
#endif
        ;
  }

  bool is_hybrid() const {
    return
#ifdef STAR_CODEGEN
#ifdef STAR_HYBRID
        true
#else
        false
#endif
#else
        sl->is_hybrid()
#endif
        ;
  }

  bool is_stoch() const { return !is_det() && !is_hybrid(); }

  bool is_det_centered() const {
    return
#ifdef STAR_CODEGEN
#ifdef STAR_DET_CENTERED
        true
#else
        false
#endif
#else
        sl->det_centered
#endif
        ;
  }

  bool hs_has_index() const {
    return
#ifdef STAR_CODEGEN
#ifdef STAR_HS_HAS_INDEX
        true
#else
        false
#endif

#else
        sl->hs_has_index
#endif
        ;
  }

  bool hs_has_flags() const {
    return
#ifdef STAR_CODEGEN
#ifdef STAR_HS_HAS_FLAGS
        true
#else
        false
#endif

#else
        sl->hs_has_flags
#endif
        ;
  }

  bool hs_has_exitrate() const {
    return
#ifdef STAR_CODEGEN
#ifdef STAR_HS_HAS_EXITRATE
        true
#else
        false
#endif
#else
        sl->hs_has_exitrate
#endif
        ;
  }

  bool hsucc_has_tr() const {
    return
#ifdef STAR_CODEGEN
#ifdef STAR_HSUCC_HAS_TR
        true
#else
        false
#endif
#else
        sl->hsucc_has_tr
#endif
        ;
  }

  bool hsucc_has_rates_g() const {
    return
#ifdef STAR_CODEGEN
#ifdef STAR_HSUCC_HAS_RATES_G
        true
#else
        false
#endif
#else
        sl->hsucc_has_rates_g
#endif
        ;
  }

  bool do_steady_state() const {
    return
#ifdef STAR_CODEGEN

#ifdef STAR_DO_STEADY_STATE
        true
#else
        false
#endif
#else
        sl->do_steady_state
#endif
        ;
  }

  bool do_rare_event() const {
    return
#ifdef STAR_CODEGEN

#ifdef STAR_DO_RARE_EVENT
        true
#else
        false
#endif
#else
        sl->do_rare_event
#endif
        ;
  }

  bool do_sa() const {
    return
#ifdef STAR_CODEGEN
#ifdef STAR_DO_SA
        true
#else
        false
#endif
#else
        sl->do_sa
#endif
        ;
  }

  bool sa_no_der() const {
    return
#ifdef STAR_CODEGEN
#ifdef STAR_SA_NO_DER
        true
#else
        false
#endif
#else
        sl->sa_no_der
#endif
        ;
  }

  bool sa_no_2der() const {
    return
#ifdef STAR_CODEGEN
#ifdef STAR_SA_NO_DER
        true
#else
#ifdef STAR_SA_NO_2DER
        true
#else
        false
#endif
#endif
#else
        sl->sa_no_2der
#endif
        ;
  }

  bool is_transient_subsolver() const {
    return
#ifdef STAR_CODEGEN
#ifdef STAR_USE_TRANSIENT_SUBSOLVER
        true
#else
        false
#endif
#else
        sl->use_transient_subsolver
#endif
        ;
  }

  std::vector<std::function<
      double(cstate const* const, double const* const, double const* const)>>&
  get_exprs() {
    return exprs;
  }

  std::size_t covlen(const std::size_t nvars_,
                     const std::size_t nmoments_) const {
    return cmb[nvars_ + nmoments_ - 1][nmoments_];
  }

  std::size_t xlen(const int nvars_, const std::size_t nmoments_) const {
    std::size_t l = 0;
    for (std::size_t i = 1; i <= nmoments_; i++) {
      l += covlen(nvars_, i);
    }
    return l;
  }

  std::size_t cov_indexl(std::size_t const* const I,
                         const std::size_t nmoments_) const {
    return sl->cov_indexl(I, nmoments_);
  }

  std::size_t cov_index_(
      std::size_t const* const I, const std::size_t d,
      const std::size_t nmoments_,
      const std::vector<solver_loader::model_info::var const*>& vars_) const {
    uint64_t k = 1;
    for (std::size_t i = 0; i < nmoments_ - d; i++) {
      k *= vars_[I[d]]->get_index() + i;
    }
    for (std::size_t i = 2; i <= nmoments_ - d; i++) {
      k /= i;
    }
    if (d + 1 < nmoments_) {
      k += cov_index_(I, d + 1, nmoments_, vars_);
    }
    return static_cast<std::size_t>(k);
  }

  std::size_t cov_indexl(
      std::size_t const* const I, const std::size_t nmoments_,
      const std::vector<solver_loader::model_info::var const*>& vars_) const {
    return cov_index_(I, 0, nmoments_, vars_);
  }

  std::size_t cov_index_(std::size_t const* const I, const std::size_t d,
                         const std::size_t nmoments_,
                         const std::vector<std::size_t>& vars_) const {
    uint64_t k = 1;
    for (std::size_t i = 0; i < nmoments_ - d; i++) {
      k *= vars_[I[d]] + i;
    }
    for (std::size_t i = 2; i <= nmoments_ - d; i++) {
      k /= i;
    }
    if (d + 1 < nmoments_) {
      k += cov_index_(I, d + 1, nmoments_, vars_);
    }
    return static_cast<std::size_t>(k);
  }

  std::size_t cov_indexl(std::size_t const* const I,
                         const std::size_t nmoments_,
                         const std::vector<std::size_t>& vars_) const {
    return cov_index_(I, 0, nmoments_, vars_);
  }

  std::size_t cov_index0(const std::size_t nvars_,
                         const std::size_t nmoments_) const {
    return xlen(nvars_, nmoments_ - 1);
  }

  std::size_t cov_indexg(const std::size_t nvars_, std::size_t const* const I,
                         const std::size_t nmoments_) const {
    return cov_index0(nvars_, nmoments_) + cov_indexl(I, nmoments_);
  }

  std::size_t cov_indexg(const std::size_t nvars_, std::size_t const* const I,
                         const std::size_t nmoments_,
                         const std::vector<std::size_t>& vars_) const {
    return cov_index0(nvars_, nmoments_) + cov_indexl(I, nmoments_, vars_);
  }

  std::size_t cov_indexg(
      const std::size_t nvars_, std::size_t const* const I,
      const std::size_t nmoments_,
      const std::vector<solver_loader::model_info::var const*>& vars_) const {
    return cov_index0(nvars_, nmoments_) + cov_indexl(I, nmoments_, vars_);
  }

  void cov_I_I(const std::size_t nvars_, std::size_t* const I,
               std::size_t& nmoments_, std::size_t const* const I_) const {
    nmoments_ = 0;
    for (int i = nvars_ - 1; i >= 0; i--) {
      for (std::size_t j = 0; j < I_[i]; j++) {
        I[nmoments_] = i;
        nmoments_++;
      }
    }
  }

  void update_trs() const {
    trs_stoch.clear();
    trs_det.clear();
    for (std::size_t i = 0; i < ntransition_classes; i++) {
      auto const& tri = std::find_if(
          sl->trs_stoch.begin(), sl->trs_stoch.end(),
          [&](solver_loader::model_info::transition const* const t) {
            return t == sl->get_transitions()[i];
          });

      if (tri != sl->trs_stoch.end()) {
        trs_stoch.push_back(sl->get_model()->get_transitions()[i]);
      } else {
        trs_det.push_back(sl->get_model()->get_transitions()[i]);
      }
    }

    if (!is_transient_subsolver()) {
      printf("trs_stoch:");
      for (auto const& tr : trs_stoch) {
        printf(" %lu", tr->get_index());
      }
      printf("\n");

      printf("trs_det:");
      for (auto const& tr : trs_det) {
        printf(" %lu", tr->get_index());
      }
      printf("\n");
    }
  }

  solver_loader::base const* const sl;

  const std::size_t nvars;
  const std::size_t npvars;
  const std::size_t ncvars;

  const std::size_t ntransition_classes;

  const std::size_t nparams;
  const std::size_t nivars;
  const std::size_t nevars;

  const std::size_t nlparams;

  std::vector<double> param_vals;
  std::vector<double> init_vals;
  std::vector<double> obserr_vals;

  const std::size_t nmoments;

  const std::size_t s_x_len;

  const std::size_t mlt_v;
  const std::size_t mlt_dv_dc;
  const std::size_t mlt_dv_di;
  const std::size_t mlt_dv_de;
  const std::size_t mlt_d2v_dc2;
  const std::size_t mlt_v_len;

  const std::size_t v_p_offset;
  const std::size_t v_p2_offset;
  const std::size_t v_px_offset;
  const std::size_t v_p_len;

  const std::size_t s_y_len;

  const std::size_t mlt_u;
  const std::size_t mlt_du_dc;
  const std::size_t mlt_d2u_dc2;
  const std::size_t mlt_u_len;

  const std::size_t v_rate_g_offset;
  const std::size_t v_rate2_g_offset;
  const std::size_t v_rate_g_len;

  const std::size_t tr_rates_g_len;

  const std::size_t v_rate_h_offset;
  const std::size_t v_drate_h_dx_offset;
  const std::size_t v_rate_h_dx_len;

  const std::size_t v_rate_h_len;
  const std::size_t v_rate_h_cnt;

  const std::size_t v_rates_h_len;

  const std::size_t tr_rates_h_len;

  const std::size_t f_len;

  const std::size_t cmb_len;
  std::vector<std::vector<std::size_t>> cmb;

  std::vector<uint64_t> fac;
  std::vector<double> _1_fac;

  std::vector<std::function<double(cstate const* const, double const* const,
                                   double const* const)>> exprs;

  const char csvSep;

  mutable std::vector<transition const*> trs_stoch, trs_det;
};
}

#endif
