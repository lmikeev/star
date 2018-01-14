/*
 *  fit.hpp
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

#ifndef SOLVER_FIT_HPP_
#define SOLVER_FIT_HPP_

#include <dirent.h>
#include "transient/rk45.hpp"

namespace solver {

struct obs_pt_s {
  double t;
  std::vector<int> x;
};

struct time_series_s {
  std::vector<std::size_t> observables;
  std::vector<obs_pt_s> observations;
};

class fit : public transient::RK45_T<> {
 private:
  std::vector<time_series_s> time_series_data;

  solver_loader::model_info::ic_s init_cond_s0;
  std::vector<std::vector<double>> normpdf0, dnormpdf0, d2normpdf0;

  mutable std::vector<double> dw_de_, d2w_de2_;

 protected:
  bool weight(int const* const v, double* const y, const obs_pt_s& o,
              const time_series_s& ts) const {
    double w = 1.0;
    std::fill(dw_de_.begin(), dw_de_.end(), 0.0);
    std::fill(d2w_de2_.begin(), d2w_de2_.end(), 0.0);

    for (std::size_t i = 0; i < ts.observables.size(); i++) {
      const int ei = sl->get_vars()[ts.observables[i]]->get_index_e();
      if (ei < 0) {
        if (v[ts.observables[i]] != o.x[i]) {
          return false;
        }
      } else {
        const int d = std::abs(v[ts.observables[i]] - o.x[i]);
        if (d < normpdf0[ei].size()) {
          if (!sa_no_der()) {
            if (!sa_no_2der()) {
              for (std::size_t ej = 0; ej <= ei; ej++) {
                const std::size_t eij = ei * (ei + 1) / 2 + ej;

                d2w_de2_[eij] = d2w_de2_[eij] * normpdf0[ei][d] +
                                dw_de_[ei] * dnormpdf0[ej][d] +
                                dw_de_[ej] * dnormpdf0[ei][d] +
                                w * d2normpdf0[eij][d];
              }
            }

            dw_de_[ei] = dw_de_[ei] * normpdf0[ei][d] + w * dnormpdf0[ei][d];
          }

          w = w * normpdf0[ei][d];
        } else {
          return false;
        }
      }
    }

    if (!sa_no_der()) {
      std::size_t kl = 0;

      for (std::size_t i = 0; i < nparams; i++) {
        if (!sa_no_2der()) {
          for (std::size_t j = 0; j <= i; kl++) {
            v_d2p_dc2(y, kl) = v_d2p_dc2(y, kl) * w;
          }
        }

        v_dp_dc(y, i) = v_dp_dc(y, i) * w;
      }

      for (std::size_t i = 0; i < nivars; i++) {
        if (!sa_no_2der()) {
          for (std::size_t j = 0; j < nparams; j++, kl++) {
            v_d2p_dc2(y, kl) = v_d2p_dc2(y, kl) * w;
          }
          for (std::size_t j = 0; j <= i; j++, kl++) {
            v_d2p_dc2(y, kl) = v_d2p_dc2(y, kl) * w;
          }
        }

        v_dp_di(y, i) = v_dp_di(y, i) * w;
      }

      for (std::size_t i = 0, ij = 0; i < nevars; i++) {
        if (!sa_no_2der()) {
          for (std::size_t j = 0; j < nparams; j++, kl++) {
            v_d2p_dc2(y, kl) = v_d2p_dc2(y, kl) * w + v_dp_dc(y, j) * dw_de_[i];
          }
          for (std::size_t j = 0; j < nivars; j++, kl++) {
            v_d2p_dc2(y, kl) = v_d2p_dc2(y, kl) * w + v_dp_di(y, j) * dw_de_[i];
          }
          for (std::size_t j = 0; j <= i; j++, kl++, ij++) {
            v_d2p_dc2(y, kl) =
                v_d2p_dc2(y, kl) * w + v_dp_de(y, i) * dw_de_[j] +
                v_dp_de(y, j) * dw_de_[i] + v_p(y) * d2w_de2_[ij];
          }
        }

        v_dp_de(y, i) = v_dp_de(y, i) * w + v_p(y) * dw_de_[i];
      }
    }

    v_p(y) = v_p(y) * w;

    return true;
  }

  bool run_objf(srch_pt_s& p) {
    update(p);

    for (std::size_t i = 0; i < nevars; i++) {
      if (this->obserr_vals[i] > std::numeric_limits<double>::epsilon()) {
        normpdf0[i].clear();

        if (!sa_no_der()) {
          dnormpdf0[i].clear();

          if (!sa_no_2der()) {
            d2normpdf0[i].clear();
          }
        }

        const double _1_2pi = 1.0 / std::sqrt(2.0 * M_PI);
        const double _1_err = 1.0 / this->obserr_vals[i];
        const double _1_err2 = _1_err * _1_err;

        int d = 0;
        while (1) {
          const double d2 = d * d;
          const double f =
              _1_2pi * std::sqrt(_1_err) * std::exp(-0.50 * d2 * _1_err);
          if (f > 0.0) {
            normpdf0[i].push_back(f);

            if (!sa_no_der()) {
              const double z = 0.50 * _1_err * (d2 * _1_err - 1.0);
              const double df = f * z;
              dnormpdf0[i].push_back(df);

              if (!sa_no_2der()) {
                const double dz = 0.50 * _1_err2 * (1.0 - 2.0 * d2 * _1_err);
                const double d2f = df * z + f * dz;
                d2normpdf0[i].push_back(d2f);
              }
            }
          } else {
            break;
          }

          d++;
        }
      } else {
        normpdf0[i].assign(1, 1.0);

        if (!sa_no_der()) {
          dnormpdf0[i].assign(1, 0.0);

          if (!sa_no_2der()) {
            d2normpdf0[i].assign(1, 0.0);
          }
        }
      }
    }

    set_hmax();

    for (auto const& ts : time_series_data) {
      std::vector<double> Lcoeff;
      Lcoeff.push_back(1.0);

      if (nivars) {
        this->get_states()->load_gaussian_init(init_cond_s0);
      } else {
        if (!this->init_(this->sl->get_ics().front())) {
          return false;
        }
      }

      auto oi = ts.observations.begin();
      double t = oi->t, h;
      comp_h0(h);

      oi++;
      while (oi != ts.observations.end()) {
        if (!this->iterate(oi->t - t, t, h)) {
          return false;
        }

        state* s = this->get_states()->get_first();
        std::size_t si = this->get_states()->get_nactive();
        double psum = 0.0;

        const std::size_t tmp = si;

        while (si--) {
          state* const lnext = s->lnext;

          void* const ls_d = this->get_da()->get_lstate_data(s);
          double* const y = this->lstate_da->y(ls_d);

          if (weight(this->get_da()->cs_get(this->get_da()->hs_cs(s->hs)), y,
                     *oi, ts)) {
            const double p = this->v_p(y);

            psum += p;
          } else {
            this->get_states()->s_erase(s);
          }

          s = lnext;
        }

        if (!this->get_states()->empty()) {
          const double _1_psum = 1.0 / psum;

          for (state* s = this->get_states()->get_first(); s != nullptr;
               s = s->lnext) {
            void* const ls_d = this->get_da()->get_lstate_data(s);
            double* const y = this->lstate_da->y(ls_d);
            for (std::size_t i = 0; i < this->s_y_len; i++) {
              y[i] *= _1_psum;
            }

            this->lstate_da->init(ls_d);
          }
        } else {
          Lcoeff.clear();
        }

        Lcoeff.back() *= psum;
        if (Lcoeff.back() > 1e20 || Lcoeff.back() < 1e-20) {
          Lcoeff.push_back(1.0);
        }

        t = oi->t;
        oi++;
      }

      if (!Lcoeff.empty()) {
        for (auto const lc : Lcoeff) {
          this->v_f(p.f.data()) -= std::log(lc);
        }

        if (!sa_no_der()) {
          std::vector<double> dsum(this->nlparams, 0.0);
          std::vector<double> d2sum(PAR2_LEN(this->nlparams), 0.0);

          for (state* s = this->get_states()->get_first(); s != nullptr;
               s = s->lnext) {
            double const* const y =
                this->lstate_da->y(this->get_da()->get_lstate_data(s));

            for (std::size_t i = 0, ij = 0; i < this->nlparams; i++) {
              dsum[i] += v_dp_dc(y, i);

              if (!sa_no_2der()) {
                for (std::size_t j = 0; j <= i; j++, ij++) {
                  d2sum[ij] += v_d2p_dc2(y, ij);
                }
              }
            }
          }

          for (std::size_t i = 0, ij = 0; i < this->nlparams; i++) {
            this->v_df_dc(p.f.data())[i] -= dsum[i];

            if (!sa_no_2der()) {
              for (std::size_t j = 0; j <= i; j++, ij++) {
                this->v_d2f_dc2(p.f.data())[ij] +=
                    dsum[i] * dsum[j] - d2sum[ij];
              }
            }
          }
        }
      } else {
        this->v_f(p.f.data()) += 1e6;
      }
    }

    for (std::size_t i = 0; i < this->f_len; i++) {
      printf(" %lf", p.f[i]);
    }
    puts("");

    return true;
  }

  bool load_data(const int repeat_index) {
    static char data_dir[1024];

#ifdef STAR_WEB_INTERFACE
    char* dpath = nullptr;
    if (!sl->get_experiment()->get_dbconnector()->get_data_path(
            dpath, sl->get_experiment()->get_id(),
            this->sl->time_series_data_src.c_str())) {
      this->sl->last_error() << "error loading data";
      return false;
    }
    sprintf(data_dir, "/var/www/star/app/webroot/d/%s/", dpath);
#else
    if (repeat_index >= 1) {
      sprintf(data_dir, "%s/%d/", this->sl->time_series_data_src.c_str(),
              repeat_index);
    } else {
      sprintf(data_dir, "%s/", this->sl->time_series_data_src.c_str());
    }
#endif

    printf("Loading '%s' ...\n", data_dir);

    DIR* dp = opendir(data_dir);
    if (dp == nullptr) {
      this->sl->last_error() << "error(" << errno << ") loading data";
      return false;
    }

    dirent* dirp;

    int ntraces = -1;
    sl->get("ntraces", ntraces);

    bool success = true;

    std::size_t tracei = 0;
    while (ntraces <= 0 || tracei < ntraces) {
      time_series_s ts;

      dirp = readdir(dp);
      if (dirp == nullptr) {
        break;
      }

      if (*dirp->d_name == '.') {
        continue;
      }

      static char buf[1024];
      sprintf(buf, "%s/%s", data_dir, dirp->d_name);

      puts(buf);

      std::ifstream in(buf);
      if (!in.is_open()) {
        this->sl->last_error() << "couldn't open observations file '" << buf
                               << "'";
        success = false;
        break;
      }

      in.getline(buf, sizeof(buf));
      char* c = buf;
      while (*c != ',' && *c != '\0') {
        c++;
      }

      char c_ = *c;
      *c = '\0';
      if (strcmp(buf, "t") != 0 && strcmp(buf, "time") != 0) {
        this->sl->last_error() << "invalid header format: 't'/'time' expected";
        success = false;
        break;
      }

      std::vector<bool> vo;

      while (c_ != '\0') {
        char* const oname = ++c;
        while (*c != ',' && *c != '\0') {
          c++;
        }
        c_ = *c;
        *c = '\0';

        solver_loader::model_info::var const* const v =
            this->sl->findvar(static_cast<std::string>(oname));
        if (v == nullptr) {
          this->sl->last_error() << "invalid header format: unknown variable '"
                                 << oname << "'";
          success = false;
          break;
        }

        if (this->sl->is_var_observable(v)) {
          vo.push_back(true);
          ts.observables.push_back(v->get_index());
        } else {
          vo.push_back(false);
        }
      }

      if (!success) {
        break;
      }

      if (vo.empty()) {
        this->sl->last_error() << "invalid header format: no observables";
        success = false;
        break;
      }

      if (ts.observables.empty()) {
        puts("skipping unobserved data ...\n");
        continue;
      }

      const std::size_t nobs = vo.size();
      const std::size_t nobservables = ts.observables.size();

      while (!in.eof()) {
        obs_pt_s o;
        o.x.resize(nobservables);

        in.getline(buf, sizeof(buf));
        char* c = buf;
        while (*c != ',' && *c != '\0') {
          c++;
        }

        char c_ = *c;
        *c = '\0';
        o.t = atof(buf);

        std::size_t oi = 0;
        std::size_t oj = 0;

        while (c_ != '\0') {
          if (oi >= nobs) {
            this->sl->last_error() << "invalid data: too many observations";
            success = false;
            break;
          }

          char* const oval = ++c;
          while (*c != ',' && *c != '\0') {
            c++;
          }
          c_ = *c;
          *c = '\0';

          if (vo[oi++]) {
            o.x[oj++] = atoi(oval);
          }
        }

        if (!oi) {
          continue;
        }

        if (!success) {
          break;
        }

        if (oi < nobs) {
          this->sl->last_error() << "invalid data: too few observations";
          success = false;
          break;
        }

        if (!ts.observations.empty() && o.t <= ts.observations.back().t) {
          this->sl->last_error() << "invalid data: non-increasing time";
          success = false;
          break;
        }

        ts.observations.push_back(o);
      }

      if (!success) {
        break;
      }

      if (ts.observations.size() < 2) {
        this->sl->last_error() << "invalid data: too few observations";
        success = false;
        break;
      }

      time_series_data.push_back(ts);
      tracei++;
    }

    closedir(dp);

    std::cout << tracei << " trajectory(ies) loaded" << std::endl;

    if (!tracei) {
      this->sl->last_error() << "no data loaded";
      success = false;
    }

    return success;
  }

 public:
  fit(solver_loader::base const* const sl)
      : transient::RK45_T<>(sl),
        options(sl),
        normpdf0(nevars, std::vector<double>()),
        dnormpdf0(nevars, std::vector<double>()),
        d2normpdf0(nevars, std::vector<double>()),
        dw_de_(nevars),
        d2w_de2_(PAR2_LEN(nevars)) {
    if (!sa_no_der()) {
      init_cond_s0.dp_di.resize(nivars, 0.0);

      if (!sa_no_2der()) {
        init_cond_s0.d2p_di2.resize(PAR2_LEN(nivars), 0.0);
      }
    }

    for (auto const& s : this->sl->get_ics().front()->get_states().front().li) {
      auto const ii =
          std::find_if(init_cond_s0.li.begin(), init_cond_s0.li.end(),
                       [&](const solver_loader::model_info::ic_s_i& i) {
            return i.v == s.v;
          });

      if (ii == init_cond_s0.li.end()) {
        init_cond_s0.li.push_back(s);
      }
    }
  }
};
}

#endif
