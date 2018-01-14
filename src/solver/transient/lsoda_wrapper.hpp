/*
 *  lsoda_wrapper.hpp
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

#ifndef SOLVER_TRANSIENT_LSODA_WRAPPER_HPP_
#define SOLVER_TRANSIENT_LSODA_WRAPPER_HPP_

#include <vector>
#include <algorithm>
#include "rk_base.hpp"
#include "../external.hpp"
#include "../lib/lsoda/lsoda.hpp"

namespace solver {
namespace transient {

class LSODAwrapper : public RKbaseT<> {
 public:
  typedef RKbaseT<> hashListBased;

  typedef typename hashListBased::state state;
  typedef typename hashListBased::hstate hstate;
  typedef typename hashListBased::hstate_succ hstate_succ;

  LSODAwrapper(solver_loader::base const* const sl)
      : options(sl), hashListBased(sl), neqs(0) {
    external_set_solver(this);
  }

  ~LSODAwrapper() { y_.clear(); }

 private:
  std::size_t neqs;
  std::vector<double> y_;

 protected:
  bool init_(solver_loader::model_info::ic const* const idistr) {
    hashListBased::load_init(idistr);

    get_states()->bfs();

    neqs = this->get_states()->get_nactive() * this->s_y_len;

    y_.resize(neqs);

    printf("%lu state(s), %lu equation(s)\n", this->get_states()->get_nactive(),
           neqs);

    return true;
  }

  bool iterate(const double H, const double t0, double&) {
    this->write_states_y(y_.data());

    int neq = neqs;
    double t = t0;
    double tout = t0 + H;
    int itol = 1;
    int itask = 1;
    int istate = 1;
    int iopt = 0;

    int iwork1, iwork2, iwork5, iwork6, iwork7, iwork8, iwork9;
    double rwork1, rwork5, rwork6, rwork7;
    iwork1 = iwork2 = iwork5 = iwork6 = iwork7 = iwork8 = iwork9 = 0;
    rwork1 = rwork5 = rwork6 = rwork7 = 0.0;

    double abs_tol = std::max(1e-6, sl->abs_tol);
    double rel_tol = std::max(1e-3, sl->rel_tol);

    int jt = 2;

    std::vector<double> yy(neqs + 1);
    for (std::size_t i = 0; i < neqs; ++i) {
      yy[i + 1] = y_[i];
    }

    lsoda(::lsoda_f, neq, yy.data(), &t, tout, itol, &rel_tol, &abs_tol, itask,
          &istate, iopt, jt, iwork1, iwork2, iwork5, iwork6, iwork7, iwork8,
          iwork9, rwork1, rwork5, rwork6, rwork7, 0);

    this->read_states_y(yy.data() + 1);

    return true;
  }

 public:
  void lsoda_f(double t, double* y, double* ydot, void*) {
    this->write_states_dy(t, y, ydot);
  }
};
}
}

#endif
