/*
 *  scan.hpp
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

#ifndef SOLVER_SCAN_HPP_
#define SOLVER_SCAN_HPP_

#include "base.hpp"

namespace solver {

template <class S>
class scan : public S {
 private:
  uint64_t nsrch_pts;
  unsigned int nrepeat;

  solver_loader::model_info::ic* init_cond;

 protected:
  std::vector<srch_pt_s> srch_pts;
  unsigned int update_dumps_mod, last_dump_index;

  int max_fevals;

  bool iterate_init(const std::size_t k, srch_pt_s& p) {
    if (k < this->sl->iterate_inits.size()) {
      auto const& ip = this->sl->iterate_inits[k];
      for (std::size_t i = 0; i < ip.values.size(); i++) {
        p.c[this->nparams + k] = ip.values[i];

        if (!iterate_init(k + 1, p)) {
          return false;
        }
      }
    } else {
      iterate_params(0, p);
    }

    return true;
  }

  bool iterate_params(const std::size_t k, srch_pt_s& p) {
    if (k < this->sl->iterate_params.size()) {
      auto const& ip = this->sl->iterate_params[k];
      for (std::size_t i = 0; i < ip.values.size(); i++) {
        p.c[k] = ip.values[i];

        if (!iterate_params(k + 1, p)) {
          return false;
        }
      }
    } else {
      for (std::size_t i = 0; i < nrepeat; i++) {
        if (!sample_init(0, p)) {
          return false;
        }
      }
    }

    return true;
  }

  bool sample_init(const std::size_t k, srch_pt_s& p) {
    if (k < this->sl->sample_inits.size()) {
      auto const& sp = this->sl->sample_inits[k];

      p.c[this->nparams + this->sl->iterate_inits.size() + k] =
          sp.f->gen_sample();

      if (!sample_init(k + 1, p)) {
        return false;
      }
    } else {
      sample_params(0, p);
    }

    return true;
  }

  bool sample_params(const std::size_t k, srch_pt_s& p) {
    if (k < this->sl->sample_params.size()) {
      auto const& sp = this->sl->sample_params[k];

      p.c[this->sl->iterate_params.size() + k] = sp.f->gen_sample();

      if (!sample_params(k + 1, p)) {
        return false;
      }
    } else {
      if (!this->run_objf(p)) {
        return false;
      }

      if (!add_srch_pt(p)) {
        return false;
      }

      if (max_fevals > 0 && (int)srch_pts.size() >= max_fevals) {
        return true;
      }
    }

    return true;
  }

  virtual bool init() {
    nsrch_pts = nrepeat;
    for (std::size_t i = 0; i < this->sl->iterate_params.size(); i++) {
      nsrch_pts *= this->sl->iterate_params[i].values.size();
    }
    for (std::size_t i = 0; i < this->sl->iterate_inits.size(); i++) {
      nsrch_pts *= this->sl->iterate_inits[i].values.size();
    }

    if (max_fevals > 0 && nsrch_pts > max_fevals) {
      nsrch_pts = max_fevals;
    }

    srch_pts.clear();
    update_dumps_mod = std::max((int)(0.10 * nsrch_pts), 1);
    last_dump_index = 0;

    return true;
  }

  virtual bool init_() {
    if (init_cond != nullptr) {
      for (auto& i : init_cond->get_states().front().li) {
        i.value->set(this->init_vals[i.v->get_index_i()]);
      }
      return S::init_(init_cond);
    }
    return S::init_(this->sl->get_ics().front());
  }

  bool run_() {
    srch_pt_s p;
    p.c.resize(this->nlparams);

    if (!iterate_init(0, p)) {
      return false;
    }

    if (srch_pts.size() % update_dumps_mod != 0) {
      if (!update_dumps() || !this->sl->update_progress(100.0)) {
        return false;
      }
    }

    return true;
  }

  bool add_srch_pt(const srch_pt_s& p) {
    for (auto const c : p.c) {
      printf(" %lf", c);
    }
    printf("\t");
    for (std::size_t i = 0; i < this->f_len; i++) {
      printf(" %lf", p.f[i]);
    }
    printf("\n");

    srch_pts.push_back(p);

    if (srch_pts.size() % update_dumps_mod == 0) {
      if (!update_dumps() || !update_progress()) {
        return false;
      }
    }

    return true;
  }

  virtual bool update_progress() {
    return this->sl->update_progress(100.0 * srch_pts.size() /
                                     (double)nsrch_pts);
  }

  bool dump_plot_objf_2d_data(task const* const tsk) const {
    solver_loader::task_info::plot_objf_2d const* const ti =
        static_cast<solver_loader::task_info::plot_objf_2d const*>(
            tsk->get_info());

    const std::size_t ci1 =
        (ti->get_et1()->is_cnst())
            ? static_cast<solver_loader::parser::cnstAST const* const>(
                  ti->get_et1())
                  ->get_cnst()
                  ->get_index_p()
            : this->nparams +
                  static_cast<solver_loader::parser::varAST const* const>(
                      ti->get_et1())
                      ->get_var()
                      ->get_index_i();

    const std::size_t ci2 =
        (ti->get_et2()->is_cnst())
            ? static_cast<solver_loader::parser::cnstAST const* const>(
                  ti->get_et2())
                  ->get_cnst()
                  ->get_index_p()
            : this->nparams +
                  static_cast<solver_loader::parser::varAST const* const>(
                      ti->get_et2())
                      ->get_var()
                      ->get_index_i();

    dump const* const dmp = tsk->get_dumps().back();
    std::ofstream os(dmp->get_local_path(), std::ios::out | std::ios::app);
    for (std::size_t i = last_dump_index; i < srch_pts.size(); i++) {
      auto const& pt = srch_pts[i];
      os << pt.c[ci1] << this->csvSep << pt.c[ci2] << this->csvSep
         << this->v_f(pt.f.data());
      if (this->do_sa() && !this->sa_no_der()) {
        const double df_dc1 = this->v_df_dc(pt.f.data())[ci1];
        const double df_dc2 = this->v_df_dc(pt.f.data())[ci2];
        os << this->csvSep << df_dc1 << this->csvSep << df_dc2;
      }
      os << std::endl;
    }
    os.close();

    return true;
  }

  bool add_tasks() {
    for (auto& t : this->sl->tasks) {
      switch (t->get_type()) {
        case solver_loader::task_info::TASK_DUMP_STATS: {
          if (!this->add_task_dump_stats(
                  static_cast<solver_loader::task_info::dump_stats*>(t))) {
            return false;
          }
          break;
        }

        case solver_loader::task_info::TASK_PLOT_OBJF: {
          if (!this->add_task_plot_objf(
                  static_cast<solver_loader::task_info::plot_objf*>(t))) {
            return false;
          }
          break;
        }

        case solver_loader::task_info::TASK_PLOT_OBJF_2D: {
          if (!this->add_task_plot_objf_2d(
                  static_cast<solver_loader::task_info::plot_objf_2d*>(t))) {
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

  bool update_dumps() {
    this->update_elapsed_time();

    for (auto& tsk : this->tasks) {
      switch (tsk->get_info()->get_type()) {
        case solver_loader::task_info::TASK_DUMP_STATS: {
          if (!this->dump_task_dump_stats(tsk, srch_pts.size())) {
            return false;
          }
          break;
        }

        case solver_loader::task_info::TASK_PLOT_OBJF: {
          if (!this->dump_task_plot_objf(tsk, srch_pts)) {
            return false;
          }
          break;
        }

        case solver_loader::task_info::TASK_PLOT_OBJF_2D: {
          if (!this->dump_task_plot_objf_2d(tsk)) {
            return false;
          }
          break;
        }

        default:
          break;
      }
    }
    last_dump_index = srch_pts.size();
    return true;
  }

  void pre_dump_stats(std::ofstream& os) const {
    os << "search-point" << this->csvSep << "run-time" << std::endl;
  }

  void dump_stats(std::ofstream& os, const double t = 0.0) const {
    os << (int)t << this->csvSep << this->get_elapsed_time() << std::endl;
  }

 public:
  scan(solver_loader::base const* const sl)
      : options(sl), S(sl), nrepeat(1), max_fevals(-1) {
    if (!this->sl->sample_params.empty()) {
      sl->get("nrepeat", nrepeat);
    }

    this->sl->get("max_fevals", max_fevals);

    if (!this->sl->iterate_inits.empty() || !this->sl->sample_inits.empty()) {
      std::vector<solver_loader::model_info::ic_s> ls(1);
      ls.front().p = 1.0;
      for (auto& p : this->sl->iterate_inits) {
        solver_loader::model_info::ic_s_i i;
        i.v = p.v;
        i.value = new solver_loader::model_info::val(nullptr);
        ls.front().li.push_back(i);
      }
      for (auto& p : this->sl->sample_inits) {
        solver_loader::model_info::ic_s_i i;
        i.v = p.v;
        i.value = new solver_loader::model_info::val(nullptr);
        ls.front().li.push_back(i);
      }
      init_cond = new solver_loader::model_info::ic("", this->sl, ls);
    } else {
      init_cond = nullptr;
    }
  }
};
}

#endif
