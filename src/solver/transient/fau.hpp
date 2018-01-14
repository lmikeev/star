/*
 *  fau.hpp
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

#ifndef SOLVER_TRANSIENT_FAU_HPP_
#define SOLVER_TRANSIENT_FAU_HPP_

#include "rk_base.hpp"

namespace solver {
namespace transient {

class fau_lstate_data_accessor : public lstate_data_accessor {
 private:
  enum fau_lstate_data_offset { do_z = do_y + 1, do_k };

 public:
  fau_lstate_data_accessor(options const* const o)
      : lstate_data_accessor(o, do_k + 1) {}

  double* z(void* const d) const {
    return static_cast<double*>(d) + do_z * len;
  }

  double const* z(void const* const d) const {
    return static_cast<double const*>(d) + do_z * len;
  }

  double* k(void* const d) const {
    return static_cast<double*>(d) + do_k * len;
  }

  double const* k(void const* const d) const {
    return static_cast<double const*>(d) + do_k * len;
  }

  virtual void init(void* const d) const {
    memset(z(d), 0, size_);
    memset(k(d), 0, size_);
  }

  virtual void zero(void* const d) const {
    memset(y(d), 0, size_);
    memset(z(d), 0, size_);
  }

  virtual void prep(void* const d) const {
    memcpy(z(d), y(d), size_);
    memset(y(d), 0, size_);
  }

  virtual void zero_tmp(void* const) const {}
  virtual void reset(void* const) const {}
};

class FAU : public RKbaseT<fau_lstate_data_accessor> {
 public:
  typedef RKbaseT<fau_lstate_data_accessor> hashListBased;

  typedef typename hashListBased::state state;
  typedef typename hashListBased::hstate hstate;
  typedef typename hashListBased::hstate_succ hstate_succ;

  FAU(solver_loader::base const* const sl)
      : options(sl),
        hashListBased(sl),
        uniformization_rate(10000.0),
        fox_glynn_epsilon(1e-20),
        bdp_epsilon(1e-20),
        fg_last_h_(std::numeric_limits<double>::max()) {}

  virtual ~FAU() {}

 private:
  double uniformization_rate, fox_glynn_epsilon, bdp_epsilon;

  double fg_last_h_;
  std::size_t fg_L, fg_R;
  std::vector<double> fg_w;
  std::vector<double> bw;

  void update_fox_glynn(const double H, std::size_t& L, std::size_t& R,
                        std::vector<double>& w) {
    if (std::fabs(fg_last_h_ - H) < std::numeric_limits<double>::epsilon()) {
      return;
    }

    const double eps = fox_glynn_epsilon;

    double l = uniformization_rate * H;
    std::size_t m = (std::size_t)l;

    if (l < 25) {
      L = 0;
    } else {
      double bl = (1.0 + 1.0 / l) * std::exp(1.0 / (8.0 * l));
      int k = 4;

      while (eps / 2.0 <
             ((bl * std::exp(-(k * k) / 2.0)) / (k * sqrt(2.0 * M_PI)))) {
        k++;
      }
      L = (std::size_t)std::floor(m - k * sqrt(l) - 3.0 / 2.0);
    }

    if (l < 400) {
      l = 400;
    }

    double al = (1.0 + 1.0 / l) * std::exp(1.0 / 16.0) * std::sqrt(2.0);
    int k = 4;
    double dkl = std::exp(-(2.0 / 9.0) * (k * std::sqrt(2.0 * l) + 3.0 / 2.0));
    dkl = 1.0 / (1.0 - dkl);

    while (eps / 2.0 < ((al * dkl * std::exp(-(k * k) / 2.0)) /
                        (k * std::sqrt(2.0 * M_PI)))) {
      k++;
      dkl = std::exp(-(2.0 / 9.0) * (k * std::sqrt(2.0 * l) + 3.0 / 2.0));
      dkl = 1.0 / (1.0 - dkl);
    }
    R = (std::size_t)std::ceil(m + k * std::sqrt(2.0 * l) + 3.0 / 2.0);

    w.resize(R - L + 1);

    const double wm = (1.7e308 / 1.e10) / (R - L);

    std::size_t j = m;
    w[m - L] = wm;
    while (j > L) {
      w[j - 1 - L] = (j / l) * w[j - L];
      j--;
    }
    j = m;
    while (j < R) {
      w[j + 1 - L] = (l / (j + 1)) * w[j - L];
      j++;
    }

    double W = 0;
    std::size_t _s = L, _t = R;
    while (_s < _t) {
      if (w[_s - L] <= w[_t - L]) {
        W += w[_s - L];
        _s++;
      } else {
        W += w[_t - L];
        _t--;
      }
    }
    W += w[_s - L];

    for (std::size_t i = 0; i < R - L + 1; i++) {
      w[i] /= W;
    }
    fg_last_h_ = H;
  }

 protected:
  bool iterate(const double H, const double, double&) {
    state* s = get_states()->get_first();
    std::size_t si = get_states()->get_nactive();
    while (si--) {
      state* const lnext = s->lnext;
      void* const ls_d = get_da()->get_lstate_data(s);
      if (this->can_ignore_s(v_p(this->lstate_da->y(ls_d)))) {
        get_states()->s_erase(s);
      } else {
        this->lstate_da->prep(ls_d);
      }
      s = lnext;
    }

    update_fox_glynn(H, fg_L, fg_R, fg_w);

    std::size_t bk = 0;
    double bpc = 1.0;
    double ba, bb = 1.0;

    const int bsize = 1024;
    std::size_t bL = 0, bR = 0, bN = 1, bN_ = bsize;

    bw.resize(bsize);
    std::fill(bw.begin(), bw.end(), 0.0);

    bool success = true;
    while (success) {
      double max_exitrate = 0.0;

      state* s = get_states()->get_first();
      std::size_t si = get_states()->get_nactive();
      while (si--) {
        state* const lnext = s->lnext;
        void* const ls_d = get_da()->get_lstate_data(s);

        v_p(this->lstate_da->z(ls_d)) += v_p(this->lstate_da->k(ls_d));
        v_p(this->lstate_da->k(ls_d)) = 0.0;

        this->s_comp_unexplored_exitrate(s);
        max_exitrate =
            std::max(max_exitrate, this->get_da()->hs_exitrate(s->hs));

        s = lnext;
      }

      if (std::fabs(max_exitrate) < std::numeric_limits<double>::epsilon()) {
        this->sl->last_error() << "max_exitrate = 0";
        success = false;
        break;
      }

      bpc = 0.0;
      double bp = (bk) ? 0.0 : 1.0;

      ba = 1.0 - bb;
      bb = 1.0 - max_exitrate / uniformization_rate;

      if (ba < 0.0 || bb < 0.0) {
        this->sl->last_error() << "Please choose higher uniformization rate";
        success = false;
        break;
      }

      std::size_t j = 0;
      std::size_t j_ = 0;

      const double eps = bdp_epsilon;

      while (j < bN) {
        const double a_bw = ba * bw[j];
        bw[j] = bp;

        if (j_ == j && bw[j] < eps) {
          j_ = j + 1;
        }

        bp = a_bw + bb * bw[j];

        const std::size_t bj = j + bL;
        if (bj >= fg_L && bj <= fg_R) {
          bpc += fg_w[bj - fg_L] * bw[j];
        }

        if (j == bN - 1 && bw[j] >= eps) {
          if (j_ > 0) {
            for (std::size_t i = j_; i < bN; i++) {
              bw[i - j_] = bw[i];
              bw[i] = 0.0;
            }
            bN -= j_;
            bL += j_;
            j -= j_;
            j_ = 0;
          }

          bR++;
          bN++;
          if (bN > bN_) {
            bN_ += bsize;
            bw.resize(bN_);
          }
        }
        j++;
      }

      if (j_ > 0) {
        for (std::size_t i = j_; i < bN; i++) {
          bw[i - j_] = bw[i];
          bw[i] = 0.0;
        }
        bN -= j_;
        bL += j_;
        j_ = 0;
      }

      double psum = 0.0;
      for (state* s = get_states()->get_first(); s != nullptr; s = s->lnext) {
        void* const ls_d = get_da()->get_lstate_data(s);
        v_p(this->lstate_da->y(ls_d)) += bpc * v_p(this->lstate_da->z(ls_d));
        psum += v_p(this->lstate_da->y(ls_d));
      }

      printf(" psum = %.16lf, bpc = %.16lf, %lf\n", psum, bpc, max_exitrate);

      if (bpc < std::numeric_limits<double>::epsilon()) {
        break;
      }

      s = get_states()->get_first();
      si = get_states()->get_nactive();
      while (si--) {
        hstate* const hs = s->hs;
        cstate const* const cs = get_da()->hs_cs(hs);
        void* const ls_d = get_da()->get_lstate_data(s);
        double const* const x = nullptr;

        get_states()->s_explore_unexplored(s);

        double exitrate = 0.0;

        hstate_succ* hsucc = hs->succ;
        while (hsucc != nullptr) {
          void const* const hsucc_d = get_da()->get_hstate_succ_data(hsucc);

          const double rate = tr_rate(hsucc_d, cs, x);
          const double d_p =
              rate / max_exitrate * v_p(this->lstate_da->z(ls_d));

          if (!this->can_ignore_tr(d_p)) {
            get_states()->hs_restore_erased(hsucc->hs);

            void* const lsucc_d =
                this->get_da()->get_lstate_data(hsucc->hs->ls);

            v_p(this->lstate_da->k(ls_d)) -= d_p;
            v_p(this->lstate_da->k(lsucc_d)) += d_p;
          }

          exitrate += rate;
          hsucc = hsucc->lnext;
        }

        get_da()->hs_exitrate(s->hs) = exitrate;

        s = s->lnext;
      }

      bk++;
    }

    return success;
  }
};
}
}

#endif
