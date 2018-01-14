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

#ifndef SOLVER_OPTIMIZATION_BASE_HPP_
#define SOLVER_OPTIMIZATION_BASE_HPP_

#include "../scan.hpp"

namespace solver {
namespace optimization {

template <class S>
class base : public scan<S> {
 protected:
  double to_log_space(const double x) const {
    return std::log(x) / std::log(log_space_base);
  }

  double from_log_space(const double x) const {
    return std::pow(log_space_base, x);
  }

  virtual bool search(std::vector<srch_pt_s>& pts) = 0;

  bool update_progress() {
    static char tmp[128];
    sprintf(tmp, "%lu points", this->srch_pts.size());
    return this->sl->update_progress(tmp);
  }

  virtual bool init() {
#if HAVE_MCR
    if (!libmatlabsdInitialize()) {
      this->sl->last_error() << "libmatlabsdInitialize = false";
      return false;
    }
#endif

    this->srch_pts.clear();
    this->update_dumps_mod = 100;
    this->last_dump_index = 0;

    if (!this->add_task(
            new solver_loader::task_info::base("parameter estimates"))) {
      return false;
    }

    tsk_pest = this->tasks.back();

    if (!this->add_dump(tsk_pest, "estimates", "estimates", "csv")) {
      return false;
    }

    std::ofstream os(tsk_pest->get_dumps().front()->get_local_path());
    for (std::size_t i = 0; i < this->nparams; i++) {
      os << this->sl->sa_params[i]->get_name() << this->csvSep;
    }
    for (std::size_t i = 0; i < this->nivars; i++) {
      os << this->sl->sa_ivars[i]->get_name() << "(0)" << this->csvSep;
    }
    for (std::size_t i = 0; i < this->nevars; i++) {
      os << "err[" << this->sl->sa_evars[i]->get_name() << "]" << this->csvSep;
    }
    os << "objf" << std::endl;
    os.close();

    return true;
  }

  bool run_() {
    unsigned int nrepeat = 1;
#ifndef STAR_WEB_INTERFACE
    this->sl->get("nrepeat", nrepeat);
#endif

    for (std::size_t repeat_index = 1; repeat_index <= nrepeat;
         repeat_index++) {
      if (!this->load_data(nrepeat > 1 ? repeat_index : -1)) {
        return false;
      }

      std::vector<srch_pt_s> pts;

      if (!search(pts)) {
        return false;
      }

      std::ofstream os(tsk_pest->get_dumps().front()->get_local_path(),
                       std::ios::out | std::ios::app);
      os.precision(16);
      for (auto const& pt : pts) {
        for (std::size_t i = 0; i < this->nlparams; i++) {
          os << pt.c[i] << this->csvSep;
        }
        os << this->v_f(pt.f.data()) << std::endl;
      }
      os.close();

      if (this->srch_pts.size() % this->update_dumps_mod != 0) {
        if (!this->update_dumps() || !this->sl->update_progress(100.0)) {
          return false;
        }
      }
    }

    return true;
  }

  bool maximize;

  double reltol_x;
  double abstol_x;
  double reltol_f;
  double abstol_f;

  task* tsk_pest;

  const double log_space_base;
  const double log_log_space_base;

 public:
  base(solver_loader::base const* const sl)
      : options(sl),
        scan<S>(sl),
        maximize(false),
        reltol_x(-1.0),
        abstol_x(-1.0),
        reltol_f(-1.0),
        abstol_f(-1.0),
        tsk_pest(nullptr),
        log_space_base(10.0),
        log_log_space_base(std::log(log_space_base)) {
    this->sl->get("maximize", maximize);

    this->sl->get("reltol_x", reltol_x);
    this->sl->get("abstol_x", abstol_x);
    this->sl->get("reltol_f", reltol_f);
    this->sl->get("abstol_f", abstol_f);
  }
};
}
}

#endif
