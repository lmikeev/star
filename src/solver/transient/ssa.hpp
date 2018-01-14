/*
 *  ssa.hpp
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

#ifndef SOLVER_TRANSIENT_SSA_HPP_
#define SOLVER_TRANSIENT_SSA_HPP_

#include "base.hpp"
#include "../hash_list_based.hpp"

namespace solver {
namespace transient {

class SSA : public base, public hashListBasedT<> {
 public:
  typedef hashListBasedT<> hashListBased;
  typedef hashListBased::stateHashList stateHashList;

  struct dump_s {
    stateHashList* states;
    int timepoint_id;
    std::vector<double> X, F, Xstd, Fstd;
  };

  SSA(solver_loader::base const* const sl)
      : options(sl),
        base(sl),
        hashListBased(sl),
        vv_(nvars),
        rates_(ntransition_classes),
        rates_g_(get_tr_rates_g_len()),
        dumps(sl->tspan.size()) {
    printf("Initializing ... \n");

    const std::size_t nexprs = exprs.size();

    for (std::size_t i = 0; i < dumps.size(); i++) {
      dumps[i].states = i ? this->create_states() : this->get_states();
      dumps[i].X.resize(nvars);
      dumps[i].F.resize(nexprs);
      dumps[i].Xstd.resize(nvars);
      dumps[i].Fstd.resize(nexprs);
    }

    printf("Starting ... \n");
  }

  typedef hashListBasedT<>::state state;
  typedef hashListBasedT<>::hstate hstate;

 private:
  std::vector<int> vv_;
  std::vector<double> rates_;
  std::vector<double> rates_g_;

  char tmp[128];

 protected:
  std::vector<dump_s> dumps;

  int vv(const int val, const int iv) const {
    if (this->sl->o_err[iv] < std::numeric_limits<double>::epsilon()) {
      return val;
    }

    static std::default_random_engine generator;
    static std::normal_distribution<double> N01(0.0, 1.0);

    solver_loader::model_info::var const* const v = this->sl->get_vars()[iv];
    const double r = val + this->sl->o_err[iv] * N01(generator);
    return std::min(this->sl->get_var_maxvalue(v),
                    std::max(this->sl->get_var_minvalue(v), (int)(r + 0.50)));
  }

  bool run_() {
    bool success = true;

    unsigned int nrepeat = 1;
#ifndef STAR_WEB_INTERFACE
    sl->get("nrepeat", nrepeat);
#endif

    unsigned int ntrajectories = 1;
    sl->get("ntrajectories", ntrajectories);

    bool dump_trajectories = false;
#ifndef STAR_WEB_INTERFACE
    sl->get("dump_trajectories", dump_trajectories);
#endif

    bool dump_trajectories_by_time = false;
#ifndef STAR_WEB_INTERFACE
    sl->get("dump_trajectories_by_time", dump_trajectories_by_time);
#endif

    namespace fs = boost::filesystem;
    fs::path dump_tr_path = sl->get_path() / "tr";
    if (dump_trajectories) {
      try {
        fs::create_directory(dump_tr_path);
      } catch (const fs::filesystem_error& e) {
        sl->last_error() << e.what();
        return false;
      }
    }

    const std::size_t upmod = std::max((int)(0.10 * ntrajectories), 1);

    std::uniform_real_distribution<double> U1(0.0, 1.0), U2(0.0, 1.0);
    std::normal_distribution<double> N01(0.0, 1.0);

    auto& d0 = dumps[0];
    auto states0 = d0.states;

    hstate* hs = states0->get_hstate_allocator()->get_last();
    states0->get_hstate_allocator()->malloc();

    cstate* cs = get_da()->hs_cs(hs);

    for (std::size_t repeat_index = 1; repeat_index <= nrepeat;
         repeat_index++) {
      sprintf(tmp, "%02lu", repeat_index);

      fs::path dump_tr_path_r = dump_tr_path / tmp;
      if (dump_trajectories || dump_trajectories_by_time) {
        try {
          fs::create_directory(dump_tr_path_r);
        } catch (const fs::filesystem_error& e) {
          sl->last_error() << e.what();
          return false;
        }

        if (dump_trajectories_by_time) {
          for (std::size_t ti = 0; ti < sl->tspan.size(); ti++) {
            sprintf(tmp, "t%03lu.csv", ti);

            fs::path dump_tr_path_rti = dump_tr_path_r / tmp;
            FILE* hf_dump_by_time = fopen(
                sl->get_experiment()->convws(dump_tr_path_rti.wstring()), "w");
            fprintf(hf_dump_by_time, "trajectory");
            if (this->sl->o_vars.empty()) {
              for (std::size_t i = 0; i < nvars; i++) {
                fprintf(hf_dump_by_time, "%c%s", csvSep,
                        sl->get_vars()[i]->get_name().c_str());
              }
            } else {
              for (int o : this->sl->o_vars) {
                fprintf(hf_dump_by_time, "%c%s", csvSep,
                        sl->get_vars()[o]->get_name().c_str());
              }
            }
            fprintf(hf_dump_by_time, "\n");
            fclose(hf_dump_by_time);
          }
        }
      }

      for (std::size_t trajectory_index = 1; trajectory_index <= ntrajectories;
           trajectory_index++) {
        auto ti = sl->tspan.begin();
        double t = *ti;
        std::size_t dumpi = 0;

        const double u = U1(rand_engine);

        state const* s0 = states0->get_first();
        double ps = get_da()->s_p(s0);
        while (u > ps) {
          s0 = s0->lnext;

          assert(s0 != nullptr);

          ps += get_da()->s_p(s0);
        }

        ti++;

        memcpy(cs, get_da()->hs_cs(s0->hs), hstate_da->get_cs_size());

        FILE* hf_dump = nullptr;
        if (dump_trajectories || dump_trajectories_by_time) {
          int const* const v = get_da()->cs_get(cs);

          if (dump_trajectories) {
            sprintf(tmp, "%06lu.csv", trajectory_index);

            fs::path dump_tr_path_rt = dump_tr_path_r / tmp;

            hf_dump = fopen(
                sl->get_experiment()->convws(dump_tr_path_rt.wstring()), "w");
            fprintf(hf_dump, "time");
            if (this->sl->o_vars.empty()) {
              for (std::size_t i = 0; i < nvars; i++) {
                fprintf(hf_dump, "%c%s", csvSep,
                        sl->get_vars()[i]->get_name().c_str());
              }
            } else {
              for (int o : this->sl->o_vars) {
                fprintf(hf_dump, "%c%s", csvSep,
                        sl->get_vars()[o]->get_name().c_str());
              }
            }

            fprintf(hf_dump, "\n%lf", t);
            if (this->sl->o_vars.empty()) {
              for (std::size_t i = 0; i < nvars; i++) {
                fprintf(hf_dump, "%c%d", csvSep, v[i]);
              }
            } else {
              for (int o : this->sl->o_vars) {
                fprintf(hf_dump, "%c%d", csvSep, v[o]);
              }
            }
            fprintf(hf_dump, "\n");
          }

          if (dump_trajectories_by_time) {
            sprintf(tmp, "t%03u.csv", 0u);

            fs::path dump_tr_path_rti = dump_tr_path_r / tmp;
            FILE* hf_dump_by_time = fopen(
                sl->get_experiment()->convws(dump_tr_path_rti.wstring()), "a");

            fprintf(hf_dump_by_time, "%lu", trajectory_index);
            if (this->sl->o_vars.empty()) {
              for (std::size_t i = 0; i < nvars; i++) {
                fprintf(hf_dump_by_time, "%c%d", csvSep, v[i]);
              }
            } else {
              for (int o : this->sl->o_vars) {
                fprintf(hf_dump_by_time, "%c%d", csvSep, v[o]);
              }
            }
            fprintf(hf_dump_by_time, "\n");
            fclose(hf_dump_by_time);
          }
        }

        while (success && ti != sl->tspan.end() &&
               (elapsed_time < max_runtime || max_runtime <= 0.0)) {
          double exitrate = 0.0;
          for (std::size_t i = 0; i < ntransition_classes; i++) {
            transition const* const tr = sl->get_model()->get_transitions()[i];
            if (tr->is_enabled(cs)) {
              tr->rates_g(cs, rates_g_.data());
              rates_[i] = v_rate_g(rates_g_.data());
            } else {
              rates_[i] = 0.0;
            }
            exitrate += rates_[i];
          }

          transition const* tr_next;
          double h;

          if (exitrate > std::numeric_limits<double>::epsilon()) {
            const double u = exitrate * U1(rand_engine);
            double rs = rates_[0];
            std::size_t i = 0;
            while (u > rs) {
              i++;

              assert(i < ntransition_classes);

              rs += rates_[i];
            }

            tr_next = sl->get_model()->get_transitions()[i];

            h = -1.0 / exitrate * std::log(U2(rand_engine));
          } else {
            tr_next = nullptr;
            h = sl->tspan[sl->tspan.size() - 1] - t;
          }

          const double t_next = t + h;

          if (hf_dump != nullptr || dump_trajectories_by_time) {
            int const* const v = get_da()->cs_get(cs);
            if (this->sl->o_vars.empty()) {
              for (std::size_t i = 0; i < nvars; i++) {
                vv_[i] = vv(v[i], i);
              }
            } else {
              for (int o : this->sl->o_vars) {
                vv_[o] = vv(v[o], o);
              }
            }
          }

          while (ti != sl->tspan.end() && t_next >= *ti) {
            t = *ti;

            if (hf_dump != nullptr) {
              fprintf(hf_dump, "%lf", t);
              if (this->sl->o_vars.empty()) {
                for (std::size_t i = 0; i < nvars; i++) {
                  fprintf(hf_dump, "%c%d", csvSep, vv_[i]);
                }
              } else {
                for (int o : this->sl->o_vars) {
                  fprintf(hf_dump, "%c%d", csvSep, vv_[o]);
                }
              }
              fprintf(hf_dump, "\n");
            }

            if (dump_trajectories_by_time) {
              sprintf(tmp, "t%03lu.csv", ti - sl->tspan.begin());

              fs::path dump_tr_path_rti = dump_tr_path_r / tmp;
              FILE* hf_dump_by_time = fopen(
                  sl->get_experiment()->convws(dump_tr_path_rti.wstring()),
                  "a");

              fprintf(hf_dump_by_time, "%lu", trajectory_index);
              if (this->sl->o_vars.empty()) {
                for (std::size_t i = 0; i < nvars; i++) {
                  fprintf(hf_dump_by_time, "%c%d", csvSep, vv_[i]);
                }
              } else {
                for (int o : this->sl->o_vars) {
                  fprintf(hf_dump_by_time, "%c%d", csvSep, vv_[o]);
                }
              }
              fprintf(hf_dump_by_time, "\n");
              fclose(hf_dump_by_time);
            }

            dumpi++;

            auto& di = dumps[dumpi];
            auto statesi = di.states;

            hstate* hsucc = statesi->get_hstate_allocator()->get_last();
            cstate* csucc = get_da()->hs_cs(hsucc);

            memcpy(csucc, cs, hstate_da->get_cs_size());

            get_da()->hs_update_hash(hsucc);

            bool fnew;
            state* const s = statesi->add(hsucc, fnew);

            if (fnew) {
              get_da()->s_p(s) = 1.0;
              get_da()->hs_init(s->hs);
              statesi->get_hstate_allocator()->malloc();
            } else {
              get_da()->s_p(s)++;
            }

            ti++;
          }

          if (tr_next != nullptr) {
            tr_next->update(cs, t);
          }

          t = t_next;
        }

        if (trajectory_index % upmod == 0) {
          if (!update_dumps(trajectory_index) ||
              !sl->update_progress(100.0 * trajectory_index /
                                   (double)ntrajectories)) {
            return false;
          }
        }

        if (hf_dump != nullptr) {
          fclose(hf_dump);
        }
      }

      if (ntrajectories % upmod != 0) {
        if (!update_dumps(ntrajectories) || !sl->update_progress(100.0)) {
          return false;
        }
      }
    }

    return success;
  }

  bool init() {
    this->load_init(this->sl->get_ics().front());

    return true;
  }

  bool add_tasks() {
    for (std::size_t i = 0; i < dumps.size(); i++) {
      if (!create_timepoint(sl->tspan[i], dumps[i].timepoint_id)) {
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
          if (!add_task_dump_distr(
                  static_cast<solver_loader::task_info::dump_distr*>(t))) {
            return false;
          }

          task* const tsk = tasks.back();

          for (std::size_t i = 0; i < dumps.size(); i++) {
            if (!add_dump(tsk, i, dumps[i].timepoint_id, "distr", "d", "csv")) {
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
          if (!add_task_plot_distr(
                  static_cast<solver_loader::task_info::plot_distr*>(t))) {
            return false;
          }

          task* const tsk = tasks.back();

          for (std::size_t i = 0; i < dumps.size(); i++) {
            if (!add_dump(tsk, i, dumps[i].timepoint_id, "plot_distr data",
                          "plot_distr_data", "csv")) {
              return false;
            }
            if (!add_plot_distr(tsk, i, dumps[i].timepoint_id,
                                tsk->get_dumps().back(), "plot_distr",
                                "plot_distr", "png")) {
              return false;
            }
          }
          break;
        }

        case solver_loader::task_info::TASK_PLOT_DISTR_2D: {
          if (!add_task_plot_distr_2d(
                  static_cast<solver_loader::task_info::plot_distr_2d*>(t))) {
            return false;
          }

          task* const tsk = tasks.back();

          for (std::size_t i = 0; i < dumps.size(); i++) {
            if (!add_dump(tsk, i, dumps[i].timepoint_id, "plot_distr_2d data",
                          "plot_distr_2d_data", "csv")) {
              return false;
            }
            if (!add_plot_distr_2d(tsk, i, dumps[i].timepoint_id,
                                   tsk->get_dumps().back(), "plot_distr_2d",
                                   "plot_distr_2d", "png")) {
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

  bool update_dumps(const std::size_t trajectory_index) {
    update_elapsed_time();

    const double p_scale = 1.0 / (double)trajectory_index;

    for (std::size_t k = 0; k < dumps.size(); k++) {
      auto& d = dumps[k];
      d.states->comp_exprs(d.X.data(), d.F.data(), d.Xstd.data(), d.Fstd.data(),
                           k ? p_scale : 1.0);
    }

    for (auto& tsk : tasks) {
      switch (tsk->get_info()->get_type()) {
        case solver_loader::task_info::TASK_DUMP_STATS: {
          if (!dump_task_dump_stats(tsk, trajectory_index)) {
            return false;
          }
          break;
        }

        case solver_loader::task_info::TASK_DUMP_MOMENTS: {
          if (!update_dump_moments(tsk, trajectory_index)) {
            return false;
          }
          break;
        }

        case solver_loader::task_info::TASK_DUMP_DISTR: {
          if (!update_dump_distr(tsk, trajectory_index)) {
            return false;
          }
          break;
        }

        case solver_loader::task_info::TASK_PLOT: {
          if (!update_plot(tsk, trajectory_index)) {
            return false;
          }
          break;
        }

        case solver_loader::task_info::TASK_PLOT_2D: {
          if (!update_plot_2d(tsk, trajectory_index)) {
            return false;
          }
          break;
        }

        case solver_loader::task_info::TASK_PLOT_DISTR: {
          if (!update_plot_distr(tsk, trajectory_index)) {
            return false;
          }
          break;
        }

        case solver_loader::task_info::TASK_PLOT_DISTR_2D: {
          if (!update_plot_distr_2d(tsk, trajectory_index)) {
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

  bool update_dump_moments(task const* const tsk,
                           const std::size_t trajectory_index) const {
    solver_loader::task_info::dump_moments const* const ti =
        static_cast<solver_loader::task_info::dump_moments const*>(
            tsk->get_info());

    const std::size_t nvars_ = ti->get_vars().size();
    const std::size_t nmoments_ = ti->get_nmoments();

    std::vector<std::size_t> I(nmoments_, 0);
    std::vector<std::size_t> I_(nvars_, 0);

    for (std::size_t d = 0; d < nmoments_; d++) {
      std::ofstream os(tsk->get_dump(d)->get_local_path());
      os.precision(16);

      os << "time" << csvSep;

      for (std::size_t i = 0; i < nvars_; i++) {
        I[0] = i;
        I_[i]++;

        pre_dump_moments(os, ti, d + 1, I.data(), I_.data(), 1);

        I_[i]--;
      }

      os << std::endl;
      os.close();
    }

    const std::size_t x_len_ = xlen(nvars_, nmoments_);

    std::vector<double> X(x_len_);

    const double p_scale = 1.0 / (double)trajectory_index;

    for (std::size_t k = 0; k < dumps.size(); k++) {
      std::fill(X.begin(), X.end(), 0.0);

      dumps[k].states->comp_moments(X.data(), ti, k ? p_scale : 1.0);

      std::size_t i0 = 0;
      for (std::size_t d = 0; d < nmoments_; d++) {
        std::ofstream os(tsk->get_dump(d)->get_local_path(),
                         std::ios::out | std::ios::app);
        os.precision(16);

        os << sl->tspan[k] << csvSep;

        const std::size_t l = covlen(nvars_, d + 1);

        for (std::size_t i = 0; i < l; i++) {
          os << X[i0 + i];
          if (i < l - 1) {
            os << csvSep;
          }
        }

        os << std::endl;
        os.close();

        i0 += l;
      }
    }

    if (!tsk->update_web_copies()) {
      return false;
    }

    return true;
  }

  bool update_dump_distr(task const* const tsk,
                         const std::size_t trajectory_index) const {
    solver_loader::task_info::dump_distr const* const ti =
        static_cast<solver_loader::task_info::dump_distr const*>(
            tsk->get_info());

    const double p_scale = 1.0 / (double)trajectory_index;

    for (std::size_t k = 0; k < dumps.size(); k++) {
      std::ofstream os(tsk->get_dumps()[k]->get_local_path());
      os.precision(16);
      write_distr_header(ti, os);
      dumps[k].states->dump_distr(os, ti, k ? p_scale : 1.0);
      os.close();
    }

    if (!tsk->update_web_copies()) {
      return false;
    }

    return true;
  }

  bool update_plot(task const* const tsk, const std::size_t) const {
    solver_loader::task_info::plot const* const ti =
        static_cast<solver_loader::task_info::plot const*>(tsk->get_info());

    solver::dump const* const dmp = tsk->get_dumps().back();

    std::ofstream os(dmp->get_local_path());
    write_plot_data_header(ti, os);

    for (std::size_t k = 0; k < dumps.size(); k++) {
      auto const& d = dumps[k];

      os << sl->tspan[k];
      if (ti->get_exprs().empty()) {
        for (std::size_t i = 0; i < nvars; i++) {
          os << csvSep << d.X[i];
        }
        if (ti->plot_stddev()) {
          for (std::size_t i = 0; i < nvars; i++) {
            os << csvSep << d.Xstd[i];
          }
        }
      } else {
        for (std::size_t i = 0; i < ti->get_exprs().size(); i++) {
          auto const& e = ti->get_exprs()[i];
          os << csvSep << d.F[e.expr_index];
        }
        if (ti->plot_stddev()) {
          for (std::size_t i = 0; i < ti->get_exprs().size(); i++) {
            auto const& e = ti->get_exprs()[i];
            os << csvSep << d.Fstd[e.expr_index];
          }
        }
      }
      os << std::endl;
    }
    os.close();

#ifdef STAR_WEB_INTERFACE
    if (!ti->is_dynamic())
#endif
    {
      if (!plot(ti, dmp, tsk->get_plots().back())) {
        return false;
      }
    }
#ifdef STAR_WEB_INTERFACE
    else {
      if (!sl->get_experiment()->get_dbconnector()->update_dyn_plot_stats(
              tsk->get_dyn_plot_id())) {
        return false;
      }
    }
#endif

    if (!tsk->update_web_copies()) {
      return false;
    }

    return true;
  }

  bool update_plot_2d(task const* const tsk, const std::size_t) const {
    solver_loader::task_info::plot_2d const* const ti =
        static_cast<solver_loader::task_info::plot_2d const*>(tsk->get_info());

    solver::dump const* const dmp = tsk->get_dumps().back();

    std::ofstream os(dmp->get_local_path());
    write_plot_2d_data_header(ti, os);

    for (std::size_t k = 0; k < dumps.size(); k++) {
      auto const& d = dumps[k];
      os << sl->tspan[k] << csvSep << d.F[ti->get_expr1_index()] << csvSep
         << d.F[ti->get_expr2_index()] << std::endl;
    }
    os.close();

#ifdef STAR_WEB_INTERFACE
    if (!ti->is_dynamic())
#endif
    {
      if (!plot_2d(ti, dmp, tsk->get_plots().back())) {
        return false;
      }
    }
#ifdef STAR_WEB_INTERFACE
    else {
      if (!sl->get_experiment()->get_dbconnector()->update_dyn_2d_plot_stats(
              tsk->get_dyn_plot_id())) {
        return false;
      }
    }
#endif

    if (!tsk->update_web_copies()) {
      return false;
    }

    return true;
  }

  bool update_plot_distr(task const* const tsk,
                         const std::size_t trajectory_index) const {
    solver_loader::task_info::plot_distr const* const ti =
        static_cast<solver_loader::task_info::plot_distr const*>(
            tsk->get_info());

    const double p_scale = 1.0 / (double)trajectory_index;

    for (std::size_t k = 0; k < dumps.size(); k++) {
      auto const& d = dumps[k];
      if (!d.states->dump_plot_distr_data(ti, tsk->get_dumps()[k],
                                          k ? p_scale : 1.0)) {
        return false;
      }
      if (!plot_distr(ti, tsk->get_dumps()[k], tsk->get_plots()[k])) {
        return false;
      }
    }

    if (!tsk->update_web_copies()) {
      return false;
    }

    return true;
  }

  bool update_plot_distr_2d(task const* const tsk,
                            const std::size_t trajectory_index) const {
    solver_loader::task_info::plot_distr_2d const* const ti =
        static_cast<solver_loader::task_info::plot_distr_2d const*>(
            tsk->get_info());

    const double p_scale = 1.0 / (double)trajectory_index;

    for (std::size_t k = 0; k < dumps.size(); k++) {
      auto const& d = dumps[k];
      if (!d.states->dump_plot_distr_2d_data(ti, tsk->get_dumps()[k],
                                             k ? p_scale : 1.0)) {
        return false;
      }
      if (!plot_distr_2d(ti, tsk->get_dumps()[k], tsk->get_plots()[k])) {
        return false;
      }
    }

    if (!tsk->update_web_copies()) {
      return false;
    }

    return true;
  }

  void pre_dump_stats(std::ofstream& os) const {
    os << "trajectory" << csvSep << "run-time" << std::endl;
  }

  void dump_stats(std::ofstream& os, const double t = 0.0) const {
    os << (int)t << csvSep << get_elapsed_time() << std::endl;
  }
};
}
}

#endif
