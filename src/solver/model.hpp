/*
 *  model.hpp
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

#ifndef SOLVER_MODEL_HPP_
#define SOLVER_MODEL_HPP_

#include <vector>
#include <cmath>

struct cstate;

namespace solver {

class options;

class transition {
 protected:
  options const* const o;
  const std::size_t index;

 public:
  transition(options const* const o, const std::size_t index)
      : o(o), index(index) {}
  virtual ~transition() {}

  std::size_t get_index() const { return index; }

  virtual bool is_enabled(cstate const* const, const double = 0.0) const {
    return true;
  }

  virtual void update(cstate* const, const double = 0.0) const {}

  virtual double const* change() const { return nullptr; }

  virtual void rates_g(cstate const* const, double* const,
                       const double = 0.0) const {}

  virtual void rates_h(cstate const* const, double const* const, double* const,
                       const double = 0.0) const {}
};

class model {
 private:
  std::vector<transition const*> transitions;

 public:
  model(options const* const) {}

  virtual ~model() {}

  void add_transition(transition const* const tr) { transitions.push_back(tr); }

  const std::vector<transition const*>& get_transitions() const {
    return transitions;
  }

  virtual void cs_set(cstate* const, int const* const) const {}

  virtual void cs_get(int* const, cstate const* const) const {}
};

constexpr std::size_t cmb_c(const int n, const int k) {
  return (k <= n) ? ((k > 0) ? cmb_c(n - 1, k - 1) + cmb_c(n - 1, k) : 1) : 0;
}

constexpr std::size_t covlen_c(const std::size_t npvars,
                               const std::size_t nmoments) {
  return cmb_c(npvars + nmoments - 1, nmoments);
}

constexpr std::size_t xlen_c(const int i, const int npvars,
                             const std::size_t nmoments) {
  return i >= 1 ? covlen_c(npvars, i) + xlen_c(i - 1, npvars, nmoments) : 0;
}

constexpr std::size_t xlen_c(const int npvars, const int nmoments) {
  return xlen_c(nmoments, npvars, nmoments);
}
}

#endif
