/*
 *  rk_base.hpp
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

#ifndef SOLVER_TRANSIENT_RK_BASE_HPP_
#define SOLVER_TRANSIENT_RK_BASE_HPP_

#include "../hash_list_based.hpp"
#include "../mfile_exporter.hpp"
#include "base.hpp"

namespace solver {
namespace transient {

class rk_lstate_data_accessor : public lstate_data_accessor {
 protected:
  enum rk_lstate_data_offset { do_z = do_y + 1 };

 public:
  rk_lstate_data_accessor(options const* const o,
                          const std::size_t n = do_z + 1)
      : lstate_data_accessor(o, n) {}

  double* z(void* const d) const {
    return static_cast<double*>(d) + do_z * len;
  }

  double const* z(void const* const d) const {
    return static_cast<double const*>(d) + do_z * len;
  }

  virtual void init(void* const d) const {
    memcpy(z(d), y(d), size_);
    zero_tmp(d);
  }

  virtual void reset(void* const d) const { memcpy(y(d), z(d), size_); }

  virtual void zero(void* const d) const {
    memset(y(d), 0, size_);
    memset(z(d), 0, size_);
  }

  virtual void zero_tmp(void* const) const {}
};

template <class lstate_data_accessor_t = rk_lstate_data_accessor,
          class hstate_data_accessor_t = hstate_data_accessor,
          class hstate_succ_data_accessor_t = hstate_succ_data_accessor>
class RKbaseT
    : public base,
      public hashListBasedT<lstate_data_accessor_t, hstate_data_accessor_t,
                            hstate_succ_data_accessor_t> {
 public:
  typedef hashListBasedT<lstate_data_accessor_t, hstate_data_accessor_t,
                         hstate_succ_data_accessor_t> hashListBased;

  typedef mfileExporterT<lstate_data_accessor_t, hstate_data_accessor_t,
                         hstate_succ_data_accessor_t> mfileExporter;

  typedef typename hashListBased::state state;
  typedef typename hashListBased::hstate hstate;
  typedef typename hashListBased::hstate_succ hstate_succ;

  template <unsigned int k>
  struct k2t {};

  RKbaseT(solver_loader::base const* const sl)
      : options(sl),
        base(sl),
        hashListBased(sl),
        mfile_exporter(
            new mfileExporter(this, this->get_da(), this->get_states())),
        X_(nvars * mlt_u_len),
        F_(exprs.size() * mlt_u_len),
        X2_(nvars * mlt_u_len),
        F2_(exprs.size() * mlt_u_len),
        dx_(npvars),
        cov_(s_x_len),
        prob_sum(0.0) {}

  virtual ~RKbaseT() {}

 private:
  mfileExporter const* const mfile_exporter;

  std::vector<double> X_;
  std::vector<double> F_;
  std::vector<double> X2_;
  std::vector<double> F2_;

 protected:
  std::vector<double> dx_;
  std::vector<double> cov_;

  double prob_sum;

  virtual bool init_(solver_loader::model_info::ic const* const idistr) {
    hashListBased::load_init(idistr);

    this->get_states()->init_states();

    return true;
  }

  bool add_tasks() {
    update_trs();

    bool export_mfile;
    if (sl->get("export_mfile", export_mfile) && export_mfile) {
      namespace fs = boost::filesystem;
      fs::path dump_mfile_path = sl->get_path() / "output.m";
      this->get_states()->bfs();
      if (!mfile_exporter->write_mfile(
              sl->get_experiment()->convws(dump_mfile_path.wstring()))) {
        return false;
      }
    }

    for (auto& t : sl->tasks) {
      switch (t->get_type()) {
        case solver_loader::task_info::TASK_DUMP_STATS: {
          if (!add_task_dump_stats(
                  static_cast<solver_loader::task_info::dump_stats*>(t))) {
            return false;
          }
          break;
        }

        case solver_loader::task_info::TASK_DUMP_MOMENTS: {
          if (!add_task_dump_moments(
                  static_cast<solver_loader::task_info::dump_moments*>(t))) {
            return false;
          }
          break;
        }

        case solver_loader::task_info::TASK_DUMP_DISTR: {
          if (!is_det()) {
            if (!add_task_dump_distr(
                    static_cast<solver_loader::task_info::dump_distr*>(t))) {
              return false;
            }
          }
          break;
        }

        case solver_loader::task_info::TASK_PLOT: {
          if (!add_task_plot(static_cast<solver_loader::task_info::plot*>(t))) {
            return false;
          }
          break;
        }

        case solver_loader::task_info::TASK_PLOT_2D: {
          if (!add_task_plot_2d(
                  static_cast<solver_loader::task_info::plot_2d*>(t))) {
            return false;
          }
          break;
        }

        case solver_loader::task_info::TASK_PLOT_DISTR: {
          if (!is_det()) {
            if (!add_task_plot_distr(
                    static_cast<solver_loader::task_info::plot_distr*>(t))) {
              return false;
            }
          }
          break;
        }

        case solver_loader::task_info::TASK_PLOT_DISTR_2D: {
          if (!is_det()) {
            if (!add_task_plot_distr_2d(
                    static_cast<solver_loader::task_info::plot_distr_2d*>(t))) {
              return false;
            }
          }
          break;
        }

        default:
          break;
      }
    }
    return true;
  }

  bool dump(const double t, const std::size_t timepoint_index,
            const std::size_t timepoint_id) {
    comp_exprs(X_.data(), F_.data(), X2_.data(), F2_.data());

    for (auto& tsk : tasks) {
      switch (tsk->get_info()->get_type()) {
        case solver_loader::task_info::TASK_DUMP_STATS: {
          if (!dump_task_dump_stats(tsk, t)) {
            return false;
          }
          break;
        }

        case solver_loader::task_info::TASK_DUMP_MOMENTS: {
          if (!dump_task_dump_moments(tsk, timepoint_index, t)) {
            return false;
          }
          break;
        }

        case solver_loader::task_info::TASK_DUMP_DISTR: {
          if (!is_det()) {
            if (!dump_task_dump_distr(tsk, timepoint_index, timepoint_id)) {
              return false;
            }
          }
          break;
        }

        case solver_loader::task_info::TASK_PLOT: {
          if (!dump_task_plot(tsk, t, X_.data(), F_.data(), X2_.data(),
                              F2_.data())) {
            return false;
          }
          break;
        }

        case solver_loader::task_info::TASK_PLOT_2D: {
          if (!dump_task_plot_2d(tsk, t, X_.data(), F_.data())) {
            return false;
          }
          break;
        }

        case solver_loader::task_info::TASK_PLOT_DISTR: {
          if (!dump_task_plot_distr(tsk, timepoint_index, timepoint_id)) {
            return false;
          }
          break;
        }

        case solver_loader::task_info::TASK_PLOT_DISTR_2D: {
          if (!dump_task_plot_distr_2d(tsk, timepoint_index, timepoint_id)) {
            return false;
          }
          break;
        }

        default:
          break;
      }
    }
    return true;
  }

  void comp_exprs(double* const X, double* const F, double* const Xstd,
                  double* const Fstd) const {
    this->get_states()->comp_exprs(X, F, Xstd, Fstd);
  }

  void comp_moments(
      double* const X,
      solver_loader::task_info::dump_moments const* const ti) const {
    this->get_states()->comp_moments(X, ti);
  }

  void comp_objf(double* const f) const { this->get_states()->comp_objf(f); }

  void update_rates() { this->get_states()->update_hsucc_data(); }

  void dump_distr(std::ofstream& os,
                  solver_loader::task_info::dump_distr const* const ti) const {
    this->get_states()->dump_distr(os, ti);
  }

  void pre_dump_stats(std::ofstream& os) const {
    os << "time";
    os << csvSep << "run-time";
    if (!is_det()) {
      os << csvSep << "nactive";
      os << csvSep << "ntotal";
      os << csvSep << "niter";
      os << csvSep << "niter-failed";
      os << csvSep << "avg-timestep";
    }
    os << std::endl;
  }

  void dump_stats(std::ofstream& os, const double t = 0.0) const {
    os << t;
    os << csvSep << get_elapsed_time();
    if (!is_det()) {
      os << csvSep << this->get_states()->get_nactive();
      os << csvSep << this->get_states()->get_ntotal();
      os << csvSep << get_niter();
      os << csvSep << get_niter_failed();
      os << csvSep << get_avg_timestep();
    }
    os << std::endl;
  }

  bool dump_plot_distr_data(task const* const tsk) const {
    return this->get_states()->dump_plot_distr_data(
        static_cast<solver_loader::task_info::plot_distr const*>(
            tsk->get_info()),
        tsk->get_dumps().back());
  }

  bool dump_plot_distr_2d_data(task const* const tsk) const {
    return this->get_states()->dump_plot_distr_2d_data(

        static_cast<solver_loader::task_info::plot_distr_2d const*>(
            tsk->get_info()),
        tsk->get_dumps().back());
  }

  void distr_err_precomp() { this->get_states()->distr_err_precomp(); }

  void distr_err_update(int const* const v, const double p, double& d,
                        double& d2, double& d_max, double& r, double& r2,
                        double& r_max, double& d_r_max) {
    this->get_states()->distr_err_update(v, p, d, d2, d_max, r, r2, r_max,
                                         d_r_max);
  }

  void distr_err_postcomp(double& d, double& d2, double& d_max, double& r,
                          double& r2, double& r_max, double& d_r_max) {
    this->get_states()->distr_err_postcomp(d, d2, d_max, r, r2, r_max, d_r_max);
  }

  virtual void comp_stages(const double) {}

  virtual double get_pow() const { return 1.0 / (1.0 + 1.0); }

  void update_nonneg(double const* const y, double* const f) const {
    std::size_t i = this->v_p_offset;
    for (; i < this->v_px_offset + this->npvars; i++) {
      if (y[i] < 0.0) {
        f[i] = std::max(f[i], 0.0);
      }
    }

    if (!is_det_centered()) {
      for (; i < this->v_p_len; i++) {
        if (y[i] < 0.0) {
          f[i] = std::max(f[i], 0.0);
        }
      }
    } else {
    }
  }

  void update_nonneg(const double ynew, const bool mbnn, double& y) const {
    y = ynew;
    if (mbnn) {
      y = std::max(y, 0.0);
    }
  }

  virtual void final_stage_update_nonneg(void* const) const {}

  virtual void s_ynew_d(double* const, void* const, const double,
                        const std::size_t, const bool, double&, double&) const {
  }

  virtual void s_comp(void* const ls_d, const double h, double& err,
                      double& normy, double& normynew, bool& is_nonneg,
                      double& nn_err) const {
    final_stage_update_nonneg(ls_d);

    double* const z = this->lstate_da->z(ls_d);
    double ynew, d;

    std::size_t i = this->v_p_offset;
    bool mbnn = true;
    for (; i < this->v_px_offset + this->npvars; i++) {
      s_ynew_d(z, ls_d, h, i, mbnn, ynew, d);
      this->update_err_nn(err, normy, normynew, is_nonneg, nn_err, mbnn, d,
                          z[i], ynew);
    }

    if (!is_det_centered()) {
      for (; i < this->v_p_len; i++) {
        s_ynew_d(z, ls_d, h, i, mbnn, ynew, d);
        this->update_err_nn(err, normy, normynew, is_nonneg, nn_err, mbnn, d,
                            z[i], ynew);
      }
    }

    mbnn = false;
    for (; i < this->s_y_len; i++) {
      s_ynew_d(z, ls_d, h, i, mbnn, ynew, d);
      this->update_err_nn(err, normy, normynew, is_nonneg, nn_err, mbnn, d,
                          z[i], ynew);
    }

    this->get_states()->s_zero_tmp(ls_d);
  }

  void comp_f0() {
    state* s = this->get_states()->get_first();
    const std::size_t nactive = this->get_states()->get_nactive();
    std::size_t si = nactive;
    while (si--) {
      hstate* const hs = s->hs;
      cstate const* const cs = this->get_da()->hs_cs(hs);
      void* const ls_d = this->get_da()->get_lstate_data(s);
      double const* const y = this->get_lstate_da()->y(ls_d);
      const double p = this->v_p(y);
      double exitrate = 0.0;

      if (p > 0.0) {
        double const* const x = this->get_da()->s_x_cov(y);
        if (!this->is_stoch() && !this->is_det_centered()) {
          this->center_moments(x, cov_.data());
        }

        double* const f0 = this->get_lstate_da()->f0(ls_d);

        this->get_states()->s_explore_unexplored(s);

        hstate_succ* hsucc = hs->succ;
        while (hsucc != nullptr) {
          void const* const hsucc_d =
              this->get_da()->get_hstate_succ_data(hsucc);
          double const* rates_g;
          double const* rates_h;
          this->tr_rates(hsucc_d, cs, x, rates_g, rates_h);

          const double rate = this->tr_rate(rates_g, rates_h, x);
          const double d_y = p * rate;

          if (!this->get_da()->hs_is_erased(hs) || !this->can_ignore_tr(d_y)) {
            this->get_states()->hs_restore_erased(hsucc->hs);

            void* const lsucc_d =
                this->get_da()->get_lstate_data(hsucc->hs->ls);
            double* const succ_f0 = this->get_lstate_da()->f0(lsucc_d);

            this->s_comp_succ(rates_g, rates_h,
                              this->is_det_centered() ? x : cov_.data(), y,
                              hsucc_d, f0, succ_f0);
          }

          hsucc = hsucc->lnext;
          exitrate += rate;
        }

        if (!this->is_stoch()) {
          if (this->is_det_centered()) {
            this->s_comp_c(cs, x, y, f0);
          } else {
            this->s_comp_c_raw(cs, cov_.data(), y, f0);
          }
        }
      }

      if (hs_has_exitrate()) {
        this->get_da()->hs_exitrate(hs) = exitrate;
      }

      s = s->lnext;
    }

    if (!this->is_stoch() && this->is_det_centered() && nmoments > 1) {
      s = this->get_states()->get_first();
      std::size_t si = nactive;
      while (si--) {
        hstate const* const hs = s->hs;
        cstate const* const cs = this->get_da()->hs_cs(hs);
        void* const ls_d = this->get_da()->get_lstate_data(s);
        double const* const y = this->lstate_da->y(ls_d);
        const double p = this->v_p(y);

        if (p > 0.0) {
          double const* const x = this->get_da()->s_x_cov(y);
          double* const f0 = this->get_lstate_da()->f0(ls_d);

          hstate_succ const* hsucc = hs->succ;
          while (hsucc != nullptr) {
            if (!this->get_da()->hs_is_erased(hsucc->hs)) {
              void const* const hsucc_d =
                  this->get_da()->get_hstate_succ_data(hsucc);
              void* const lsucc_d =
                  this->get_da()->get_lstate_data(hsucc->hs->ls);
              double* const succ_f0 = this->get_lstate_da()->f0(lsucc_d);
              double const* const y_succ = this->lstate_da->y(lsucc_d);
              const double p_succ = this->v_p(y_succ);

              double const* dx;
              if (p_succ > 0.0) {
                double const* const px_succ = this->v_px(y_succ);
                for (std::size_t i = 0; i < npvars; i++) {
                  dx_[i] = x[i] - px_succ[i] / p_succ;
                }
                dx = dx_.data();
              } else {
                dx = x;
              }

              this->s_comp_succ_cov(cs, x, y, hsucc_d, dx, f0, succ_f0);
            }

            hsucc = hsucc->lnext;
          }

          this->s_comp_c_cov(cs, x, y, f0);
        }

        s = s->lnext;
      }
    }
  }

  void comp_h0(double& h) {
    if (this->get_states()->get_first() == nullptr ||
        this->lstate_da->f0(this->get_da()->get_lstate_data(
            this->get_states()->get_first())) == nullptr) {
      return;
    }
    comp_f0();

    double rh = 0.0;
    prob_sum = 0.0;
    if (this->sl->norm_control) {
      double normy = 0.0;
      double normf0 = 0.0;
      state* s = this->get_states()->get_first();
      std::size_t si = this->get_states()->get_nactive();
      while (si--) {
        state* const lnext = s->lnext;
        void* const ls_d = this->get_da()->get_lstate_data(s);
        double const* const y = this->lstate_da->y(ls_d);
        double* const f0 = this->lstate_da->f0(ls_d);
        for (std::size_t i = 0; i < this->s_y_len; i++) {
          normy += y[i] * y[i];
          normf0 += f0[i] * f0[i];
          f0[i] = 0.0;
        }
        rh =
            std::sqrt(normf0) / std::max(std::sqrt(normy), this->err_threshold);
        if (!this->can_ignore_s(v_p(y))) {
          prob_sum += v_p(y);
        } else {
          this->get_states()->s_erase(s);
        }
        s = lnext;
      }
    } else {
      state* s = this->get_states()->get_first();
      std::size_t si = this->get_states()->get_nactive();
      while (si--) {
        state* const lnext = s->lnext;
        void* const ls_d = this->get_da()->get_lstate_data(s);
        double const* const y = this->lstate_da->y(ls_d);
        double* const f0 = this->lstate_da->f0(ls_d);
        for (std::size_t i = 0; i < this->s_y_len; i++) {
          rh = std::max(
              rh, std::fabs(f0[i]) / std::max(fabs(y[i]), this->err_threshold));
          f0[i] = 0.0;
        }
        if (!this->can_ignore_s(v_p(y))) {
          prob_sum += v_p(y);
        } else {
          this->get_states()->s_erase(s);
        }
        s = lnext;
      }
    }

    rh /= 0.80 * std::pow(this->sl->rel_tol, get_pow());

    h = hmax;
    if (h * rh > 1.0) {
      h = 1.0 / rh;
    }

    h = std::max(h, hmin);
  }

  virtual bool iterate0(const double t, double& h, const bool done,
                        double& lerr, bool& redo, bool& nofailed) {
    comp_stages(h);

    double err = 0.0;
    double normy = 0.0;
    double normynew = 0.0;
    bool is_nonneg = true;
    double nn_err = 0.0;
    for (state* s = this->get_states()->get_first(); s != nullptr;
         s = s->lnext) {
      s_comp(this->get_da()->get_lstate_data(s), h, err, normy, normynew,
             is_nonneg, nn_err);
    }

    this->update_err_nn(err, normy, normynew, nn_err);

    err *= h;

    redo = err > this->sl->rel_tol;

    bool nn_rejectstep = false;
    if (!redo && !is_nonneg && nn_err > this->sl->rel_tol) {
      err = nn_err;
      nn_rejectstep = true;
      redo = true;
    }

    const std::size_t nstates_before = this->get_states()->get_nactive();
    const double h_before = h;
    bool success = true;
    prob_sum = 0.0;

    if (!redo) {
      state* s = this->get_states()->get_first();
      std::size_t si = this->get_states()->get_nactive();
      while (si--) {
        state* const lnext = s->lnext;
        void* const ls_d = this->get_da()->get_lstate_data(s);
        double const* const y = this->lstate_da->y(ls_d);
        double* const z = this->lstate_da->z(ls_d);
        double ymx = 0.0;
        for (std::size_t i = 0; i < this->s_y_len; i++) {
          ymx = std::max(ymx, std::abs(y[i]));
          z[i] = y[i];
        }

        if (!this->can_ignore_s(ymx)) {
          prob_sum += this->v_p(y);
        } else {
          this->get_states()->s_erase(s);
        }

        s = lnext;
      }
      lerr = 1.0 - prob_sum;
    } else {
      state* s = this->get_states()->get_first();
      std::size_t si = this->get_states()->get_nactive();
      while (si) {
        state* const lnext = s->lnext;

        this->get_states()->s_reset(s);

        void* const ls_d = this->get_da()->get_lstate_data(s);
        const double p = this->get_da()->s_p(ls_d);

        if (!this->can_ignore_s(p)) {
          prob_sum += p;
        } else {
          this->get_states()->s_erase(s);
        }

        s = lnext;
        si--;
      }
      lerr = 999.0;

      if (h <= hmin) {
        success = false;
      }

      if (nofailed) {
        nofailed = false;
        h = std::max(
            hmin,
            h * ((nn_rejectstep)
                     ? 0.50
                     : std::max(0.10, 0.80 * std::pow(this->sl->rel_tol / err,
                                                      get_pow()))));
      } else {
        h = std::max(hmin, 0.50 * h);
      }
    }

    if (nofailed) {
      const double tmp = 1.25 * pow(err / this->sl->rel_tol, get_pow());
      if (tmp > 0.20) {
        h /= tmp;
      } else {
        h *= 5.0;
      }
    }

    if (!is_transient_subsolver()) {
      if (done || !success || !((this->niter + this->niter_failed) & 0xfff)) {
        printf("%9u %5u\t%le %le %lu(%lu,%lu) %le %le %d %le\n", this->niter,
               this->niter_failed, t, h_before,
               this->get_states()->get_nactive(), nstates_before,
               this->get_states()->get_ntotal(), err, nn_err, redo, lerr);
      }
    }

    if (!success) {
      return err_tol_not_met(t);
    }
    return true;
  }

  template <unsigned int k>
  void comp_k(const double h) {
    state* s = this->get_states()->get_first();
    const std::size_t nactive = this->get_states()->get_nactive();
    std::size_t si = nactive;
    while (si--) {
      hstate* const hs = s->hs;
      cstate const* const cs = this->get_da()->hs_cs(hs);
      void* const ls_d = this->get_da()->get_lstate_data(s);

      s_comp_k_set_y<k>(ls_d, h);

      double const* const y = this->lstate_da->y(ls_d);
      const double p = this->v_p(y);

      if (p > 0.0) {
        this->get_states()->s_explore_unexplored(s);

        double const* const x = this->get_da()->s_x_cov(y);
        if (!this->is_stoch() && !this->is_det_centered()) {
          this->center_moments(x, cov_.data());
        }

        hstate_succ* hsucc = hs->succ;
        while (hsucc != nullptr) {
          void const* const hsucc_d =
              this->get_da()->get_hstate_succ_data(hsucc);
          double const* rates_g;
          double const* rates_h;
          this->get_da()->tr_rates(hsucc_d, cs, x, rates_g, rates_h);

          const double rate = this->tr_rate(rates_g, rates_h, x);
          const double d_p = p * rate;

          if (!this->get_da()->hs_is_erased(hs) ||
              !this->can_ignore_tr(d_p * h)) {
            this->get_states()->hs_restore_erased(hsucc->hs);

            void* const lsucc_d =
                this->get_da()->get_lstate_data(hsucc->hs->ls);

            s_comp_k_succ<k>(cs, this->is_det_centered() ? x : cov_.data(), y,
                             hsucc_d, ls_d, lsucc_d);
          }

          hsucc = hsucc->lnext;
        }

        if (!this->is_stoch()) {
          if (this->is_det_centered()) {
            s_comp_k_c<k>(cs, x, y, ls_d);
          } else {
            s_comp_k_c_raw<k>(cs, cov_.data(), y, ls_d);
          }
        }
      }

      s = s->lnext;
    }

    if (!this->is_stoch() && this->is_det_centered() && nmoments > 1) {
      s = this->get_states()->get_first();
      std::size_t si = nactive;
      while (si--) {
        hstate const* const hs = s->hs;
        cstate const* const cs = this->get_da()->hs_cs(hs);
        void* const ls_d = this->get_da()->get_lstate_data(s);
        double const* const y = this->lstate_da->y(ls_d);
        const double p = this->v_p(y);

        if (p > 0.0) {
          double const* const x = this->get_da()->s_x_cov(y);

          hstate_succ const* hsucc = hs->succ;
          while (hsucc != nullptr) {
            if (!this->get_da()->hs_is_erased(hsucc->hs)) {
              void const* const hsucc_d =
                  this->get_da()->get_hstate_succ_data(hsucc);
              void* const lsucc_d =
                  this->get_da()->get_lstate_data(hsucc->hs->ls);
              double const* const y_succ = this->lstate_da->y(lsucc_d);
              const double p_succ = this->v_p(y_succ);

              double const* dx;
              if (p_succ > 0.0) {
                double const* const px_succ = this->v_px(y_succ);
                for (std::size_t i = 0; i < npvars; i++) {
                  dx_[i] = x[i] - px_succ[i] / p_succ;
                }
                dx = dx_.data();
              } else {
                dx = x;
              }

              s_comp_k_succ_cov<k>(cs, x, y, hsucc_d, dx, ls_d, lsucc_d);
            }

            hsucc = hsucc->lnext;
          }

          s_comp_k_c_cov<k>(cs, x, y, ls_d);
        }

        s = s->lnext;
      }
    }
  }

  template <unsigned int k>
  void s_comp_k_set_y(void* const ls_d, const double h) const {
    s_comp_k_set_y(ls_d, h, k2t<k>());
  }

  virtual void s_comp_k_set_y(void* const, const double, k2t<1>) const {}

  virtual void s_comp_k_set_y(void* const, const double, k2t<2>) const {}

  virtual void s_comp_k_set_y(void* const, const double, k2t<3>) const {}

  virtual void s_comp_k_set_y(void* const, const double, k2t<4>) const {}

  virtual void s_comp_k_set_y(void* const, const double, k2t<5>) const {}

  virtual void s_comp_k_set_y(void* const, const double, k2t<6>) const {}

  virtual void s_comp_k_set_y(void* const, const double, k2t<7>) const {}

  virtual void s_comp_k_set_y(void* const, const double, k2t<8>) const {}

  virtual void s_comp_k_set_y(void* const, const double, k2t<9>) const {}

  virtual void s_comp_k_set_y(void* const, const double, k2t<10>) const {}

  virtual void s_comp_k_set_y(void* const, const double, k2t<11>) const {}

  virtual void s_comp_k_set_y(void* const, const double, k2t<12>) const {}

  virtual void s_comp_k_set_y(void* const, const double, k2t<13>) const {}

  template <unsigned int k>
  void s_comp_k_succ(cstate const* const cs, double const* const x,
                     double const* const y, void const* const hsucc_d,
                     void* const ls_d, void* const lsucc_d) const {
    s_comp_k_succ(cs, x, y, hsucc_d, ls_d, lsucc_d, k2t<k>());
  }

  virtual void s_comp_k_succ(cstate const* const, double const* const,
                             double const* const, void const* const,
                             void* const, void* const, k2t<1>) const {}

  virtual void s_comp_k_succ(cstate const* const, double const* const,
                             double const* const, void const* const,
                             void* const, void* const, k2t<2>) const {}

  virtual void s_comp_k_succ(cstate const* const, double const* const,
                             double const* const, void const* const,
                             void* const, void* const, k2t<3>) const {}

  virtual void s_comp_k_succ(cstate const* const, double const* const,
                             double const* const, void const* const,
                             void* const, void* const, k2t<4>) const {}

  virtual void s_comp_k_succ(cstate const* const, double const* const,
                             double const* const, void const* const,
                             void* const, void* const, k2t<5>) const {}

  virtual void s_comp_k_succ(cstate const* const, double const* const,
                             double const* const, void const* const,
                             void* const, void* const, k2t<6>) const {}

  virtual void s_comp_k_succ(cstate const* const, double const* const,
                             double const* const, void const* const,
                             void* const, void* const, k2t<7>) const {}

  virtual void s_comp_k_succ(cstate const* const, double const* const,
                             double const* const, void const* const,
                             void* const, void* const, k2t<8>) const {}

  virtual void s_comp_k_succ(cstate const* const, double const* const,
                             double const* const, void const* const,
                             void* const, void* const, k2t<9>) const {}

  virtual void s_comp_k_succ(cstate const* const, double const* const,
                             double const* const, void const* const,
                             void* const, void* const, k2t<10>) const {}

  virtual void s_comp_k_succ(cstate const* const, double const* const,
                             double const* const, void const* const,
                             void* const, void* const, k2t<11>) const {}

  virtual void s_comp_k_succ(cstate const* const, double const* const,
                             double const* const, void const* const,
                             void* const, void* const, k2t<12>) const {}

  virtual void s_comp_k_succ(cstate const* const, double const* const,
                             double const* const, void const* const,
                             void* const, void* const, k2t<13>) const {}

  template <unsigned int k>
  void s_comp_k_c(cstate const* const cs, double const* const x,
                  double const* const y, void* const ls_d) const {
    s_comp_k_c(cs, x, y, ls_d, k2t<k>());
  }

  virtual void s_comp_k_c(cstate const* const, double const* const,
                          double const* const, void* const, k2t<1>) const {}

  virtual void s_comp_k_c(cstate const* const, double const* const,
                          double const* const, void* const, k2t<2>) const {}

  virtual void s_comp_k_c(cstate const* const, double const* const,
                          double const* const, void* const, k2t<3>) const {}

  virtual void s_comp_k_c(cstate const* const, double const* const,
                          double const* const, void* const, k2t<4>) const {}

  virtual void s_comp_k_c(cstate const* const, double const* const,
                          double const* const, void* const, k2t<5>) const {}

  virtual void s_comp_k_c(cstate const* const, double const* const,
                          double const* const, void* const, k2t<6>) const {}

  virtual void s_comp_k_c(cstate const* const, double const* const,
                          double const* const, void* const, k2t<7>) const {}

  virtual void s_comp_k_c(cstate const* const, double const* const,
                          double const* const, void* const, k2t<8>) const {}

  virtual void s_comp_k_c(cstate const* const, double const* const,
                          double const* const, void* const, k2t<9>) const {}

  virtual void s_comp_k_c(cstate const* const, double const* const,
                          double const* const, void* const, k2t<10>) const {}

  virtual void s_comp_k_c(cstate const* const, double const* const,
                          double const* const, void* const, k2t<11>) const {}

  virtual void s_comp_k_c(cstate const* const, double const* const,
                          double const* const, void* const, k2t<12>) const {}

  virtual void s_comp_k_c(cstate const* const, double const* const,
                          double const* const, void* const, k2t<13>) const {}

  template <unsigned int k>
  void s_comp_k_c_raw(cstate const* const cs, double const* const x,
                      double const* const y, void* const ls_d) const {
    s_comp_k_c_raw(cs, x, y, ls_d, k2t<k>());
  }

  virtual void s_comp_k_c_raw(cstate const* const, double const* const,
                              double const* const, void* const, k2t<1>) const {}

  virtual void s_comp_k_c_raw(cstate const* const, double const* const,
                              double const* const, void* const, k2t<2>) const {}

  virtual void s_comp_k_c_raw(cstate const* const, double const* const,
                              double const* const, void* const, k2t<3>) const {}

  virtual void s_comp_k_c_raw(cstate const* const, double const* const,
                              double const* const, void* const, k2t<4>) const {}

  virtual void s_comp_k_c_raw(cstate const* const, double const* const,
                              double const* const, void* const, k2t<5>) const {}

  virtual void s_comp_k_c_raw(cstate const* const, double const* const,
                              double const* const, void* const, k2t<6>) const {}

  virtual void s_comp_k_c_raw(cstate const* const, double const* const,
                              double const* const, void* const, k2t<7>) const {}

  virtual void s_comp_k_c_raw(cstate const* const, double const* const,
                              double const* const, void* const, k2t<8>) const {}

  virtual void s_comp_k_c_raw(cstate const* const, double const* const,
                              double const* const, void* const, k2t<9>) const {}

  virtual void s_comp_k_c_raw(cstate const* const, double const* const,
                              double const* const, void* const, k2t<10>) const {
  }

  virtual void s_comp_k_c_raw(cstate const* const, double const* const,
                              double const* const, void* const, k2t<11>) const {
  }

  virtual void s_comp_k_c_raw(cstate const* const, double const* const,
                              double const* const, void* const, k2t<12>) const {
  }

  virtual void s_comp_k_c_raw(cstate const* const, double const* const,
                              double const* const, void* const, k2t<13>) const {
  }

  template <unsigned int k>
  void s_comp_k_succ_cov(cstate const* const cs, double const* const x,
                         double const* const y, void const* const hsucc_d,
                         double const* const dx, void* const ls_d,
                         void* const lsucc_d) const {
    s_comp_k_succ_cov(cs, x, y, hsucc_d, dx, ls_d, lsucc_d, k2t<k>());
  }

  virtual void s_comp_k_succ_cov(cstate const* const, double const* const,
                                 double const* const, void const* const,
                                 double const* const, void* const, void* const,
                                 k2t<1>) const {}

  virtual void s_comp_k_succ_cov(cstate const* const, double const* const,
                                 double const* const, void const* const,
                                 double const* const, void* const, void* const,
                                 k2t<2>) const {}

  virtual void s_comp_k_succ_cov(cstate const* const, double const* const,
                                 double const* const, void const* const,
                                 double const* const, void* const, void* const,
                                 k2t<3>) const {}

  virtual void s_comp_k_succ_cov(cstate const* const, double const* const,
                                 double const* const, void const* const,
                                 double const* const, void* const, void* const,
                                 k2t<4>) const {}

  virtual void s_comp_k_succ_cov(cstate const* const, double const* const,
                                 double const* const, void const* const,
                                 double const* const, void* const, void* const,
                                 k2t<5>) const {}

  virtual void s_comp_k_succ_cov(cstate const* const, double const* const,
                                 double const* const, void const* const,
                                 double const* const, void* const, void* const,
                                 k2t<6>) const {}

  virtual void s_comp_k_succ_cov(cstate const* const, double const* const,
                                 double const* const, void const* const,
                                 double const* const, void* const, void* const,
                                 k2t<7>) const {}

  virtual void s_comp_k_succ_cov(cstate const* const, double const* const,
                                 double const* const, void const* const,
                                 double const* const, void* const, void* const,
                                 k2t<8>) const {}

  virtual void s_comp_k_succ_cov(cstate const* const, double const* const,
                                 double const* const, void const* const,
                                 double const* const, void* const, void* const,
                                 k2t<9>) const {}

  virtual void s_comp_k_succ_cov(cstate const* const, double const* const,
                                 double const* const, void const* const,
                                 double const* const, void* const, void* const,
                                 k2t<10>) const {}

  virtual void s_comp_k_succ_cov(cstate const* const, double const* const,
                                 double const* const, void const* const,
                                 double const* const, void* const, void* const,
                                 k2t<11>) const {}

  virtual void s_comp_k_succ_cov(cstate const* const, double const* const,
                                 double const* const, void const* const,
                                 double const* const, void* const, void* const,
                                 k2t<12>) const {}

  virtual void s_comp_k_succ_cov(cstate const* const, double const* const,
                                 double const* const, void const* const,
                                 double const* const, void* const, void* const,
                                 k2t<13>) const {}

  template <unsigned int k>
  void s_comp_k_c_cov(cstate const* const cs, double const* const x,
                      double const* const y, void* const ls_d) const {
    s_comp_k_c_cov(cs, x, y, ls_d, k2t<k>());
  }

  virtual void s_comp_k_c_cov(cstate const* const, double const* const,
                              double const* const, void* const, k2t<1>) const {}

  virtual void s_comp_k_c_cov(cstate const* const, double const* const,
                              double const* const, void* const, k2t<2>) const {}

  virtual void s_comp_k_c_cov(cstate const* const, double const* const,
                              double const* const, void* const, k2t<3>) const {}

  virtual void s_comp_k_c_cov(cstate const* const, double const* const,
                              double const* const, void* const, k2t<4>) const {}

  virtual void s_comp_k_c_cov(cstate const* const, double const* const,
                              double const* const, void* const, k2t<5>) const {}

  virtual void s_comp_k_c_cov(cstate const* const, double const* const,
                              double const* const, void* const, k2t<6>) const {}

  virtual void s_comp_k_c_cov(cstate const* const, double const* const,
                              double const* const, void* const, k2t<7>) const {}

  virtual void s_comp_k_c_cov(cstate const* const, double const* const,
                              double const* const, void* const, k2t<8>) const {}

  virtual void s_comp_k_c_cov(cstate const* const, double const* const,
                              double const* const, void* const, k2t<9>) const {}

  virtual void s_comp_k_c_cov(cstate const* const, double const* const,
                              double const* const, void* const, k2t<10>) const {
  }

  virtual void s_comp_k_c_cov(cstate const* const, double const* const,
                              double const* const, void* const, k2t<11>) const {
  }

  virtual void s_comp_k_c_cov(cstate const* const, double const* const,
                              double const* const, void* const, k2t<12>) const {
  }

  virtual void s_comp_k_c_cov(cstate const* const, double const* const,
                              double const* const, void* const, k2t<13>) const {
  }

  bool post_iterate() { return true; }
};
}
}

#endif
