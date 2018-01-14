/*
 *  base.hpp
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

#ifndef SOLVER_BASE_HPP_
#define SOLVER_BASE_HPP_

#include <sstream>
#include <inttypes.h>
#include <iostream>
#include <algorithm>
#include <map>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include "task.hpp"
#include "options.hpp"
#include "../solver_loader/base.hpp"
#include "../experiment.hpp"

#if HAVE_DLIB
#include "lib/dlib/optimization.h"

typedef dlib::matrix<double, 0, 1> column_vector;
typedef dlib::matrix<double> general_matrix;
#endif

#if HAVE_MCR
#include "../solver/matlab/libmatlabsd.hpp"
#endif

namespace solver {

struct srch_pt_s {
  std::vector<double> c;
  std::vector<double> f;
};

class base : virtual public options {
 public:
  base(solver_loader::base const* const sl)
      : options(sl),
        rand_engine(rand_device()),
        subsolver(nullptr),
        elapsed_time(0.0),
        max_runtime(-1.0) {
    this->sl->get("max_runtime", max_runtime);

    for (std::size_t i = 0; i < nparams; i++) {
      param_vals[i] = sl->sa_params[i]->get_value()->get();
    }
  }

  virtual ~base() {
    for (auto& t : tasks) {
    }
  }

  solver_loader::base const* get_solver_loader() const { return sl; }

  const std::vector<task*>& get_tasks() const { return tasks; }

  virtual bool get_timepoint(const unsigned int, double&) const {
    return false;
  }

  virtual bool init() { return true; }

  virtual bool init_() { return true; }

  virtual bool add_tasks() { return true; }

  virtual bool run_() = 0;

  bool run() {
#if HAVE_MCR

#endif

    start_timer();

    if (!init()) {
      return false;
    }

    if (!add_tasks()) {
      return false;
    }

    return run_();
  }

  base* get_subsolver() { return subsolver; }

  virtual void lsoda_f(double, double*, double*, void*) {}

  virtual double nlopt_f(const std::vector<double>&, std::vector<double>&,
                         void*) {
    return 0.0;
  }

#if HAVE_DLIB
  virtual double dlib_opt_f(const column_vector&) { return 0.0; }

  virtual column_vector dlib_opt_g(const column_vector&) {
    return column_vector();
  }

  virtual general_matrix dlib_opt_h(const column_vector&) {
    return general_matrix();
  }
#endif

#if HAVE_MCR
  virtual bool MW_CALL_CONV matlab_ft(int, mxArray* [], int, mxArray* []) {
    return false;
  }

  virtual bool MW_CALL_CONV matlab_fgh(int, mxArray* [], int, mxArray* []) {
    return false;
  }
#endif

 private:
  std::chrono::high_resolution_clock::time_point start_time;

  std::random_device rand_device;

  mutable char tmp[1 << 20];

 protected:
  std::mt19937 rand_engine;

  std::vector<task*> tasks;

  base* subsolver;

  double elapsed_time;

  double max_runtime;

  void start_timer() { start_time = std::chrono::high_resolution_clock::now(); }

  void update_elapsed_time() {
    using namespace std::chrono;
    const high_resolution_clock::time_point t = high_resolution_clock::now();
    elapsed_time = duration_cast<duration<double>>(t - start_time).count();
  }

  double get_elapsed_time() const { return elapsed_time; }

  void foo(const int) const {}

  task* get_task(const unsigned int i) const { return tasks[i]; }

  bool add_dump(task* const tsk, char const* const dump_name,
                char const* const fname, char const* const fext = "csv",
                const int index = -1) {
    tsk->add_dump(dump_name, fname, fext, index);

#ifdef STAR_WEB_INTERFACE
    if (!sl->get_experiment()->get_dbconnector()->create_dump(
            tsk->get_id(), dump_name,
            tsk->get_dumps().back()->get_web_fname())) {
      return false;
    }
#endif

    return true;
  }

  bool add_dump(task* const tsk, const std::size_t timepoint_index,
                const std::size_t timepoint_id, char const* const dump_name,
                char const* const fname, char const* const fext = "csv") {
    tsk->add_dump(timepoint_index, dump_name, fname, fext);

#ifdef STAR_WEB_INTERFACE
    if (!sl->get_experiment()->get_dbconnector()->create_dump(
            tsk->get_id(), timepoint_id, dump_name,
            tsk->get_dumps().back()->get_web_fname())) {
      return false;
    }
#else

    foo(timepoint_id);
#endif

    return true;
  }

  bool add_plot(task* const tsk, dump const* const dmp,
                char const* const plot_name, char const* const fname,
                char const* const fext = "png", const int index = -1) {
#ifdef STAR_WEB_INTERFACE
    solver_loader::task_info::plot const* const ti =
        static_cast<solver_loader::task_info::plot const*>(tsk->get_info());
    if (ti->is_dynamic()) {
      const unsigned int dyn_plot_id =
          sl->get_experiment()->get_dbconnector()->create_dyn_plot(
              tsk->get_id(), plot_name, dmp->get_web_fname());
      if (!dyn_plot_id) {
        return false;
      }
      tsk->set_dyn_plot_id(dyn_plot_id);
    } else
#else

    foo(*dmp->get_name());
#endif
    {
      tsk->add_plot(plot_name, fname, fext, index);
    }
    return true;
  }

  bool add_plot_2d(task* const tsk, dump const* const dmp,
                   char const* const plot_name, char const* const fname,
                   char const* const fext = "png", const int index = -1) {
#ifdef STAR_WEB_INTERFACE
    solver_loader::task_info::plot_2d const* const ti =
        static_cast<solver_loader::task_info::plot_2d const*>(tsk->get_info());
    if (ti->is_dynamic()) {
      const unsigned int dyn_2d_plot_id =
          sl->get_experiment()->get_dbconnector()->create_dyn_2d_plot(
              tsk->get_id(), plot_name, dmp->get_web_fname());
      if (!dyn_2d_plot_id) {
        return false;
      }
      tsk->set_dyn_plot_id(dyn_2d_plot_id);
    } else
#else

    foo(*dmp->get_name());
#endif
    {
      tsk->add_plot(plot_name, fname, fext, index);
    }
    return true;
  }

  bool add_plot_distr(task* const tsk, const std::size_t timepoint_index,
                      const std::size_t timepoint_id, dump const* const dmp,
                      char const* const plot_name, char const* const fname,
                      char const* const fext = "png") {
    tsk->add_plot(timepoint_index, plot_name, fname, fext);

#ifdef STAR_WEB_INTERFACE
    if (!sl->get_experiment()->get_dbconnector()->create_plot(
            tsk->get_id(), timepoint_id, plot_name,
            tsk->get_plots().back()->get_web_fname())) {
      return false;
    }
#else

    foo(timepoint_id);
    foo(*dmp->get_name());
#endif

    return true;
  }

  bool add_plot_distr_2d(task* const tsk, const std::size_t timepoint_index,
                         const std::size_t timepoint_id, dump const* const dmp,
                         char const* const plot_name, char const* const fname,
                         char const* const fext = "png") {
    tsk->add_plot(timepoint_index, plot_name, fname, fext);

#ifdef STAR_WEB_INTERFACE
    if (!sl->get_experiment()->get_dbconnector()->create_plot(
            tsk->get_id(), timepoint_id, plot_name,
            tsk->get_plots().back()->get_web_fname())) {
      return false;
    }
#else

    foo(timepoint_id);
    foo(*dmp->get_name());
#endif

    return true;
  }

  bool add_plot_objf(task* const tsk, char const* const plot_name,
                     char const* const fname, char const* const fext = "png") {
    tsk->add_plot(plot_name, fname, fext);

#ifdef STAR_WEB_INTERFACE
    if (!sl->get_experiment()->get_dbconnector()->create_plot(
            tsk->get_id(), plot_name,
            tsk->get_plots().back()->get_web_fname())) {
      return false;
    }
#endif

    return true;
  }

  bool add_plot_objf_2d(task* const tsk, char const* const plot_name,
                        char const* const fname,
                        char const* const fext = "png") {
    tsk->add_plot(plot_name, fname, fext);

#ifdef STAR_WEB_INTERFACE
    if (!sl->get_experiment()->get_dbconnector()->create_plot(
            tsk->get_id(), plot_name,
            tsk->get_plots().back()->get_web_fname())) {
      return false;
    }
#endif

    return true;
  }

  bool create_timepoint(const double t, int& timepoint_id) const {
#ifdef STAR_WEB_INTERFACE
    timepoint_id = sl->get_experiment()->get_dbconnector()->create_timepoint(
        sl->get_id(), t);
    if (timepoint_id == 0) {
      return false;
    }
#else
    static int timepoint_id_;
    timepoint_id = ++timepoint_id_;

    foo(static_cast<unsigned int>(t));
#endif

    return true;
  }

  bool create_task(char const* const task_name, int& task_id) const {
#ifdef STAR_WEB_INTERFACE
    task_id = sl->get_experiment()->get_dbconnector()->create_task(sl->get_id(),
                                                                   task_name);
    if (task_id == 0) {
      return false;
    }
#else
    task_id = tasks.size() + 1;

    foo(*task_name);
#endif

    return true;
  }

  bool create_task_dir(task const* const tsk) const {
    boost::filesystem::path p = sl->get_path() / tsk->get_dir();
    sprintf(tmp, "mkdir %s", p.c_str());
    if (system(tmp)) {
      return false;
    }
    return true;
  }

  bool add_task(solver_loader::task_info::base* const ti) {
    int tsk_id;
    if (!create_task(ti->get_name(), tsk_id)) {
      return false;
    }

    tasks.push_back(new task(this, tsk_id, ti, sl->get_path()));

    if (!create_task_dir(tasks.back())) {
      return false;
    }

    return true;
  }

  void pre_dump_moments(std::ofstream& os,
                        solver_loader::task_info::dump_moments const* const ti,
                        const std::size_t nmoments_, std::size_t* const I,
                        std::size_t* const I_, const std::size_t d) const {
    if (d < nmoments_) {
      for (std::size_t i = 0; i <= I[d - 1]; i++) {
        I[d] = i;
        I_[i]++;

        pre_dump_moments(os, ti, nmoments_, I, I_, d + 1);

        I_[i]--;
      }
    } else {
      bool f_ = false;
      for (std::size_t i = 0; i < d; i++) {
        if (i) {
          os << "^";
        }
        os << ti->get_vars()[I[i]]->get_name();

        if (I[i] < ti->get_vars().size() - 1) {
          f_ = true;
        }
      }

      if (f_) {
        os << csvSep;
      }
    }
  }

  void create_nrm_file(task const* const tsk) const {
    sprintf(tmp, "%s/nrm.csv", tsk->get_local_dir());
    std::ofstream os_nrm(tmp);
    if (os_nrm.is_open()) {
      os_nrm << "timepoint" << csvSep;
      sl->write_nrm_header(os_nrm) << std::endl;
      os_nrm.close();
    } else {
      std::cerr << "Couldn't open file '" << tmp << "'" << std::endl;
    }
  }

  bool add_task_dump_moments(solver_loader::task_info::dump_moments* const ti,
                             const bool dump_time = true) {
    if (!add_task(ti)) {
      return false;
    }

    task* const tsk = tasks.back();

    for (int i = 1; i <= ti->get_nmoments(); i++) {
      static char buf[16];
      sprintf(buf, "m_%d", i);
      if (!add_dump(tsk, buf, "m", "csv", i)) {
        return false;
      }
    }

    if (ti->do_cmp()) {
      sprintf(tmp, "%s/nrm.csv", tsk->get_local_dir());
      std::ofstream os_nrm(tmp);
      if (dump_time) {
        os_nrm << "time" << csvSep;
      }
      os_nrm << "moment" << csvSep;
      sl->write_nrm_header(os_nrm) << std::endl;
      os_nrm.close();
    }

    const std::size_t nvars_ = ti->get_vars().size();
    const std::size_t nmoments_ = ti->get_nmoments();

    std::vector<std::size_t> I(nmoments_, 0);
    std::vector<std::size_t> I_(nvars_, 0);

    for (std::size_t d = 0; d < nmoments_; d++) {
      std::ofstream os(tsk->get_dump(d)->get_local_path(),
                       std::ios::out | std::ios::app);
      os.precision(16);

      if (dump_time) {
        os << "time" << csvSep;
      }

      for (std::size_t i = 0; i < nvars_; i++) {
        I[0] = i;
        I_[i]++;

        pre_dump_moments(os, ti, d + 1, I.data(), I_.data(), 1);

        I_[i]--;
      }

      os << std::endl;
      os.close();
    }

    if (!tsk->update_web_copies()) {
      return false;
    }

    return true;
  }

  bool dump_task_dump_moments(task const* const tsk,
                              const std::size_t timepoint_index,
                              double const t = 0.0,
                              const bool dump_time = true) {
    solver_loader::task_info::dump_moments const* const ti =
        static_cast<solver_loader::task_info::dump_moments const*>(
            tsk->get_info());

    const std::size_t nvars_ = ti->get_vars().size();
    const std::size_t nmoments_ = ti->get_nmoments();

    std::vector<double> X(xlen(nvars_, nmoments_), 0.0);

    comp_moments(X.data(), ti);

    sprintf(tmp, "%s/mlist_%lu.csv", tsk->get_local_dir(), timepoint_index);
    std::ofstream os2(tmp);
    os2.precision(16);

    std::ifstream is_cmp;
    std::ofstream os_nrm;
    if (ti->do_cmp()) {
      sprintf(tmp, "%s/mlist_%lu.csv", ti->get_cmp_path().c_str(),
              timepoint_index);
      is_cmp.open(tmp);
      if (is_cmp.is_open()) {
        sprintf(tmp, "%s/nrm.csv", tsk->get_local_dir());
        os_nrm.open(tmp, std::ios::out | std::ios::app);
        if (!os_nrm.is_open()) {
          std::cerr << "Couldn't create file '" << tmp << "'" << std::endl;
        }
      } else {
        std::cerr << "Couldn't open file '" << tmp << "'" << std::endl;
      }
    }

    if (dump_time && os_nrm.is_open()) {
      os_nrm << t << std::endl;
    }

    std::size_t i0 = 0;
    for (std::size_t mi = 0; mi < nmoments_; mi++) {
      std::ofstream os(tsk->get_dump(mi)->get_local_path(),
                       std::ios::out | std::ios::app);
      os.precision(16);

      if (dump_time) {
        os << t << csvSep;
      }

      const std::size_t l = covlen(nvars_, mi + 1);

      double d = 0.0, d2 = 0.0, d_max = 0.0;
      double r = 0.0, r2 = 0.0, r_max = 0.0;
      double d_r_max = 0.0;

      for (std::size_t i = 0; i < l; i++) {
        os << X[i0 + i];
        if (i < l - 1) {
          os << csvSep;
        }

        os2 << X[i0 + i] << std::endl;

        if (os_nrm.is_open()) {
          double x;
          is_cmp >> x;
          sl->update_err(X[i0 + i], x, d, d2, d_max, r, r2, r_max, d_r_max);
        }
      }

      os << std::endl;
      os.close();

      os2 << std::endl;

      if (os_nrm.is_open()) {
        os_nrm << csvSep << mi + 1 << csvSep;
        sl->write_nrm(os_nrm, d, d2, d_max, r, r2, r_max, d_r_max) << std::endl;
      }

      i0 += l;
    }

    os2.close();

    if (os_nrm.is_open()) {
      is_cmp.close();

      os_nrm << std::endl;
      os_nrm.close();
    }

#ifndef STAR_WEB_INTERFACE
    std::string copy_mean_path;
    if (sl->get("copy_mean_path", copy_mean_path)) {
      sprintf(tmp, "cp %s %s", tsk->get_dump(0)->get_local_path(),
              copy_mean_path.c_str());
      if (system(tmp)) {
        return false;
      }
    }
#endif

    if (!tsk->update_web_copies()) {
      return false;
    }

    return true;
  }

  bool add_task_dump_stats(solver_loader::task_info::dump_stats* const ti) {
    if (!add_task(ti)) {
      return false;
    }

    task* const tsk = tasks.back();

    if (!add_dump(tsk, "stats", "stats", "csv")) {
      return false;
    }

    std::ofstream os(tsk->get_dump(0)->get_local_path(),
                     std::ios::out | std::ios::app);
    os.precision(16);
    pre_dump_stats(os);
    os.close();

    if (!tsk->update_web_copies()) {
      return false;
    }

    return true;
  }

  bool dump_task_dump_stats(task const* const tsk, double const t = 0.0) {
    std::ofstream os(tsk->get_dump(0)->get_local_path(),
                     std::ios::out | std::ios::app);
    dump_stats(os, t);
    os.close();

    if (!tsk->update_web_copies()) {
      return false;
    }

    return true;
  }

  virtual void distr_err_precomp() {}

  virtual void distr_err_update(int const* const, const double, double&,
                                double&, double&, double&, double&, double&,
                                double&) {}

  virtual void distr_err_postcomp(double&, double&, double&, double&, double&,
                                  double&, double&) {}

  bool add_task_dump_distr(solver_loader::task_info::dump_distr* const ti) {
    if (!add_task(ti)) {
      return false;
    }

    if (ti->do_cmp()) {
      create_nrm_file(tasks.back());
    }

    return true;
  }

  void write_distr_header(solver_loader::task_info::dump_distr const* const,
                          std::ofstream& os) const {
    for (auto const v : sl->p_vars) {
      os << v->get_name() << csvSep;
    }
    for (auto const v : sl->c_vars) {
      os << v->get_name() << csvSep;
    }
    os << "probability";
    if (do_sa() && !sa_no_der()) {
      for (std::size_t i = 0; i < nparams; i++) {
        os << csvSep
           << sl->get_cnsts()[sl->sa_params[i]->get_index_p()]->get_name();
      }
    }
    os << std::endl;
  }

  bool dump_task_dump_distr(task* const tsk, const std::size_t timepoint_index,
                            const std::size_t timepoint_id) {
    if (!add_dump(tsk, timepoint_index, timepoint_id, "distr", "d", "csv")) {
      return false;
    }

    solver_loader::task_info::dump_distr const* const ti =
        static_cast<solver_loader::task_info::dump_distr const*>(
            tsk->get_info());

    dump const* const dmp = tsk->get_dumps().back();

    std::ofstream os(dmp->get_local_path());
    os.precision(16);
    write_distr_header(ti, os);
    dump_distr(os, ti);
    os.close();

    dmp->update_web_copy();

    if (ti->do_cmp()) {
      sprintf(tmp, "%s/d_%lu.csv", ti->get_cmp_path().c_str(), timepoint_index);
      std::ifstream is_cmp(tmp);
      if (is_cmp.is_open()) {
        sprintf(tmp, "%s/nrm.csv", tsk->get_local_dir());
        std::ofstream os_nrm(tmp, std::ios::out | std::ios::app);
        if (os_nrm.is_open()) {
          double d = 0.0, d2 = 0.0, d_max = 0.0;
          double r = 0.0, r2 = 0.0, r_max = 0.0;
          double d_r_max = 0.0;

          distr_err_precomp();

          std::string line;

          getline(is_cmp, line);
          while (getline(is_cmp, line)) {
            static std::vector<int> v(1024);
            double p;
            std::istringstream iss(line.c_str());
            std::string s;
            std::size_t vi = 0;
            while (getline(iss, s, ',')) {
              if (vi < nvars) {
                v[vi++] = atoi(s.c_str());
              } else {
                p = atof(s.c_str());
              }
            }

            distr_err_update(v.data(), p, d, d2, d_max, r, r2, r_max, d_r_max);
          }
          is_cmp.close();

          distr_err_postcomp(d, d2, d_max, r, r2, r_max, d_r_max);

          os_nrm << timepoint_index << csvSep;
          sl->write_nrm(os_nrm, d, d2, d_max, r, r2, r_max, d_r_max)
              << std::endl;
          os_nrm.close();
        } else {
          std::cerr << "Couldn't create file '" << tmp << "'" << std::endl;
        }
      } else {
        std::cerr << "Couldn't open file '" << tmp << "'" << std::endl;
      }
    }

    return true;
  }

  void write_plot_data_header(solver_loader::task_info::plot const* const ti,
                              std::ofstream& os) const {
    os << "time";
    if (ti->get_exprs().empty()) {
      for (std::size_t i = 0; i < nvars; i++) {
        os << csvSep << sl->get_vars()[i]->get_name();
      }
    } else {
      std::string tmp_;
      for (std::size_t i = 0; i < ti->get_exprs().size(); i++) {
        auto const& e = ti->get_exprs()[i];
        if (e.line_props.get_title(tmp_)) {
          os << csvSep << tmp_;
        } else {
          os << csvSep << e.expr;
        }
      }
    }
    os << std::endl;
  }

  bool add_task_plot(solver_loader::task_info::plot* const ti) {
    if (!add_task(ti)) {
      return false;
    }

    task* const tsk = tasks.back();

    if (!add_dump(tsk, "plot data", "plot_data", "csv")) {
      return false;
    }

    dump const* const dmp = tsk->get_dumps().back();
    std::ofstream os(dmp->get_local_path());
    write_plot_data_header(ti, os);
    os.close();

    if (!dmp->update_web_copy()) {
      return false;
    }

    if (!add_plot(tsk, dmp, "plot", "plot", "png")) {
      return false;
    }

    return true;
  }

  std::size_t var_expr_index_(
      const solver_loader::task_info::plot_expr& e) const {
    assert(e.et->is_var());

    const std::size_t vi =
        static_cast<solver_loader::parser::varAST const* const>(e.et)
            ->get_var()
            ->get_index();

    if (e.dparami == nullptr) {
      return vi;
    }

    const std::size_t ci = e.dparami->get_index_p();

    if (e.dparamj == nullptr) {
      return (mlt_du_dc + ci) * nvars + vi;
    }

    const std::size_t cij = CI(ci, e.dparamj->get_index_p());

    return (mlt_d2u_dc2 + cij) * nvars + vi;
  }

  bool dump_task_plot(task const* const tsk, double const t,
                      double const* const X, double const* const F,
                      double const* const Xstd, double const* const Fstd) {
    solver_loader::task_info::plot const* const ti =
        static_cast<solver_loader::task_info::plot const*>(tsk->get_info());

    dump const* const dmp = tsk->get_dumps().back();

    std::ofstream os(dmp->get_local_path(), std::ios::out | std::ios::app);
    os.precision(16);
    os << t;
    if (ti->get_exprs().empty()) {
      for (std::size_t i = 0; i < nvars; i++) {
        os << csvSep << X[sl->get_vars()[i]->get_index()];
      }
      if (ti->plot_stddev()) {
        for (std::size_t i = 0; i < nvars; i++) {
          os << csvSep << Xstd[sl->get_vars()[i]->get_index()];
        }
      }
    } else {
      for (auto const& e : ti->get_exprs()) {
        if (e.et->is_var()) {
          os << csvSep << X[var_expr_index_(e)];
        } else {
          os << csvSep << F[e.expr_index];
        }
      }
      if (ti->plot_stddev()) {
        for (auto const& e : ti->get_exprs()) {
          if (e.et->is_var()) {
            os << csvSep << Xstd[var_expr_index_(e)];
          } else {
            os << csvSep << Fstd[e.expr_index];
          }
        }
      }
    }
    os << std::endl;
    os.close();

#ifdef STAR_WEB_INTERFACE
    if (!ti->is_dynamic())
#endif
    {
      if (!plot(ti, tsk->get_dumps().back(), tsk->get_plots().back())) {
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

#ifndef STAR_WEB_INTERFACE
    if (ti->get_copy_data_path()[0] != '\0') {
      sprintf(tmp, "cp %s %s", dmp->get_local_path(), ti->get_copy_data_path());
      if (system(tmp)) {
        return false;
      }
    }
#endif

    if (!tsk->update_web_copies()) {
      return false;
    }

    return true;
  }

  void write_plot_2d_data_header(
      solver_loader::task_info::plot_2d const* const ti,
      std::ofstream& os) const {
    os << "time" << csvSep << ti->get_expr1() << csvSep << ti->get_expr2()
       << std::endl;
  }

  bool add_task_plot_2d(solver_loader::task_info::plot_2d* const ti) {
    if (!add_task(ti)) {
      return false;
    }

    task* const tsk = tasks.back();

    if (!add_dump(tsk, "plot_2d data", "plot_2d_data", "csv")) {
      return false;
    }

    dump const* const dmp = tsk->get_dumps().back();
    std::ofstream os(dmp->get_local_path());
    write_plot_2d_data_header(ti, os);
    os.close();

    if (!dmp->update_web_copy()) {
      return false;
    }

    if (!add_plot_2d(tsk, dmp, "plot_2d", "plot_2d", "png")) {
      return false;
    }

    return true;
  }

  bool dump_task_plot_2d(task const* const tsk, double const t,
                         double const* const X, double const* const F) {
    solver_loader::task_info::plot_2d const* const ti =
        static_cast<solver_loader::task_info::plot_2d const*>(tsk->get_info());

    dump const* const dmp = tsk->get_dumps().back();

    std::ofstream os(dmp->get_local_path(), std::ios::out | std::ios::app);
    os << t << csvSep << F[ti->get_expr1_index()] << csvSep
       << F[ti->get_expr2_index()] << std::endl;
    os.close();

    foo(static_cast<int>(*X));

#ifdef STAR_WEB_INTERFACE
    if (!ti->is_dynamic())
#endif
    {
      if (!plot_2d(ti, tsk->get_dumps().back(), tsk->get_plots().back())) {
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

  bool add_task_plot_distr(solver_loader::task_info::plot_distr* const ti) {
    if (!add_task(ti)) {
      return false;
    }

    if (ti->do_cmp()) {
      create_nrm_file(tasks.back());
    }

    return true;
  }

  bool dump_task_plot_distr(task* const tsk, const std::size_t timepoint_index,
                            const std::size_t timepoint_id) {
    solver_loader::task_info::plot_distr const* const ti =
        static_cast<solver_loader::task_info::plot_distr const*>(
            tsk->get_info());

    if (!add_dump(tsk, timepoint_index, timepoint_id, "plot_distr data",
                  "plot_distr_data", "csv")) {
      return false;
    }

    dump const* const dmp = tsk->get_dumps().back();

    if (!dump_plot_distr_data(tsk)) {
      return false;
    }

    if (!dmp->update_web_copy()) {
      return false;
    }

    if (!add_plot_distr(tsk, timepoint_index, timepoint_id, dmp, "plot_distr",
                        "plot_distr", "png")) {
      return false;
    }

    if (!plot_distr(ti, tsk->get_dumps().back(), tsk->get_plots().back())) {
      return false;
    }

    if (!tsk->get_plots().back()->update_web_copy()) {
      return false;
    }

    return true;
  }

  bool add_task_plot_distr_2d(
      solver_loader::task_info::plot_distr_2d* const ti) {
    if (!add_task(ti)) {
      return false;
    }

    return true;
  }

  bool dump_task_plot_distr_2d(task* const tsk,
                               const std::size_t timepoint_index,
                               const std::size_t timepoint_id) {
    solver_loader::task_info::plot_distr_2d const* const ti =
        static_cast<solver_loader::task_info::plot_distr_2d const*>(
            tsk->get_info());

    if (!add_dump(tsk, timepoint_index, timepoint_id, "plot_distr_2d data",
                  "plot_distr_2d_data", "csv")) {
      return false;
    }

    dump const* const dmp = tsk->get_dumps().back();

    if (!dump_plot_distr_2d_data(tsk)) {
      return false;
    }

    if (!dmp->update_web_copy()) {
      return false;
    }

    if (!add_plot_distr_2d(tsk, timepoint_index, timepoint_id, dmp,
                           "plot_distr_2d", "plot_distr_2d", "png")) {
      return false;
    }

    if (!plot_distr_2d(ti, tsk->get_dumps().back(), tsk->get_plots().back())) {
      return false;
    }

    if (!tsk->get_plots().back()->update_web_copy()) {
      return false;
    }

    return true;
  }

  bool add_task_plot_objf(solver_loader::task_info::plot_objf* const ti) {
    if (!add_task(ti)) {
      return false;
    }

    task* const tsk = tasks.back();

    if (ti->get_exprs().empty()) {
      for (std::size_t i = 0; i < nparams; i++) {
        if (!add_dump(tsk, "plot_objf data", "plot_objf_data", "csv", i)) {
          return false;
        }

        dump const* const dmp = tsk->get_dumps().back();
        std::ofstream os(dmp->get_local_path());
        os << sl->sa_params[i]->get_name() << csvSep << sl->objf << std::endl;
        os.close();
      }
      for (std::size_t i = 0; i < nivars; i++) {
        if (!add_dump(tsk, "plot_objf data", "plot_objf_data", "csv",
                      nparams + i)) {
          return false;
        }

        dump const* const dmp = tsk->get_dumps().back();
        std::ofstream os(dmp->get_local_path());
        os << sl->sa_ivars[i]->get_name() << csvSep << sl->objf << std::endl;
        os.close();
      }
    } else {
      for (std::size_t i = 0; i < ti->get_exprs().size(); i++) {
        if (!add_dump(tsk, "plot_objf data", "plot_objf_data", "csv", i)) {
          return false;
        }

        dump const* const dmp = tsk->get_dumps().back();
        std::ofstream os(dmp->get_local_path());
        os << ti->get_exprs()[i].expr << csvSep << sl->objf << std::endl;
        os.close();
      }
    }

    if (!tsk->update_web_copies()) {
      return false;
    }

    if (!add_plot_objf(tsk, "plot_objf", "plot_objf", "png")) {
      return false;
    }

    return true;
  }

  bool dump_task_plot_objf(task const* const tsk,
                           const std::vector<srch_pt_s>& srch_pts) {
    solver_loader::task_info::plot_objf const* const ti =
        static_cast<solver_loader::task_info::plot_objf const*>(
            tsk->get_info());

    std::vector<std::pair<double, double>> v;

    if (ti->get_exprs().empty()) {
      for (std::size_t i = 0; i < nparams; i++) {
        for (auto const& pt : srch_pts) {
          v.push_back(std::make_pair(pt.c[i], v_f(pt.f.data())));
        }
        std::sort(v.begin(), v.end());
        std::ofstream os(tsk->get_dumps()[i]->get_local_path());
        os << sl->sa_params[i]->get_name() << csvSep << sl->objf << std::endl;
        for (std::size_t j = 0; j < v.size(); j++) {
          os << v[j].first << csvSep << v[j].second << std::endl;
        }
        os.close();
        v.clear();
      }
      for (std::size_t i = 0; i < nivars; i++) {
        for (auto const& pt : srch_pts) {
          v.push_back(std::make_pair(pt.c[nparams + i], v_f(pt.f.data())));
        }
        std::sort(v.begin(), v.end());
        std::ofstream os(tsk->get_dumps()[nparams + i]->get_local_path());
        os << sl->sa_ivars[i]->get_name() << csvSep << sl->objf << std::endl;
        for (std::size_t j = 0; j < v.size(); j++) {
          os << v[j].first << csvSep << v[j].second << std::endl;
        }
        os.close();
        v.clear();
      }
    } else {
      for (std::size_t i = 0; i < ti->get_exprs().size(); i++) {
        for (auto const& pt : srch_pts) {
          const std::size_t ci =
              (ti->get_exprs()[i].et->is_cnst())
                  ? static_cast<solver_loader::parser::cnstAST const* const>(
                        ti->get_exprs()[i].et)
                        ->get_cnst()
                        ->get_index_p()
                  : nparams +
                        static_cast<solver_loader::parser::varAST const* const>(
                            ti->get_exprs()[i].et)
                            ->get_var()
                            ->get_index_i();

          v.push_back(std::make_pair(pt.c[ci], v_f(pt.f.data())));
        }
        std::sort(v.begin(), v.end());
        std::ofstream os(tsk->get_dumps()[i]->get_local_path());
        os << ti->get_exprs()[i].expr << csvSep << sl->objf << std::endl;
        for (std::size_t j = 0; j < v.size(); j++) {
          os << v[j].first << csvSep << v[j].second << std::endl;
        }
        os.close();
        v.clear();
      }
    }

    if (!plot_objf(tsk)) {
      return false;
    }

    if (!tsk->update_web_copies()) {
      return false;
    }

    return true;
  }

  bool add_task_plot_objf_2d(solver_loader::task_info::plot_objf_2d* const ti) {
    if (!add_task(ti)) {
      return false;
    }

    task* const tsk = tasks.back();

    if (!add_dump(tsk, "plot_objf_2d data", "plot_objf_2d_data", "csv")) {
      return false;
    }

    dump const* const dmp = tsk->get_dumps().back();
    std::ofstream os(dmp->get_local_path());
    os << ti->get_expr1() << csvSep << ti->get_expr2() << csvSep << sl->objf
       << std::endl;
    os.close();

    if (!dmp->update_web_copy()) {
      return false;
    }

    if (!add_plot_objf_2d(tsk, "plot_objf_2d", "plot_objf_2d", "png")) {
      return false;
    }

    return true;
  }

  bool dump_task_plot_objf_2d(task const* const tsk) {
    if (!dump_plot_objf_2d_data(tsk)) {
      return false;
    }

    if (!plot_objf_2d(tsk)) {
      return false;
    }

    if (!tsk->update_web_copies()) {
      return false;
    }

    return true;
  }

  bool plot(solver_loader::task_info::plot const* const ti,
            dump const* const dmp, dump const* const plt) const;
  bool plot_2d(solver_loader::task_info::plot_2d const* const ti,
               dump const* const dmp, dump const* const plt) const;
  bool plot_distr(solver_loader::task_info::plot_distr const* const ti,
                  dump const* const dmp, dump const* const plt) const;
  bool plot_distr_2d(solver_loader::task_info::plot_distr_2d const* const ti,
                     dump const* const dmp, dump const* const plt) const;
  bool plot_objf(task const* const tsk) const;
  bool plot_objf_2d(task const* const tsk) const;

  virtual void comp_exprs(double* const, double* const, double* const,
                          double* const) const {}

  virtual void comp_moments(
      double* const,
      solver_loader::task_info::dump_moments const* const) const {}

  virtual void pre_dump_stats(std::ofstream&) const {}

  virtual void dump_stats(std::ofstream&, double const = 0.0) const {}

  virtual void dump_distr(
      std::ofstream&, solver_loader::task_info::dump_distr const* const) const {
  }

  virtual bool dump_plot_distr_data(task const* const) const { return false; }

  virtual bool dump_plot_distr_2d_data(task const* const) const {
    return false;
  }

  virtual bool dump_plot_objf_2d_data(task const* const) const { return false; }

  virtual void comp_objf(double* const) const {}

  virtual bool run_objf(srch_pt_s&) { return false; }

  void update(srch_pt_s& p) {
    p.f.assign(this->f_len, 0.0);

    this->param_vals.assign(p.c.begin(), p.c.begin() + nparams);
    this->init_vals.assign(p.c.begin() + nparams,
                           p.c.begin() + nparams + nivars);
    this->obserr_vals.assign(p.c.begin() + nparams + nivars,
                             p.c.begin() + nparams + nivars + nevars);

    this->update_rates();
  }

  virtual void update_rates() {}

  virtual bool load_data(const int) { return true; }
};
}

#endif
