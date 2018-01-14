/*
 *  transition.hpp
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

#ifndef SOLVER_LOADER_MODEL_INFO_TRANSITION_HPP_
#define SOLVER_LOADER_MODEL_INFO_TRANSITION_HPP_

#include <vector>
#include "obj.hpp"

namespace solver_loader {

namespace parser {

class AST;
}

namespace model_info {

class var;

struct transition_update_item {
  var const* v;
  parser::AST const* u;
};

class transition : public obj {
 private:
  parser::AST const* guard;
  parser::AST const* rate;
  const std::vector<transition_update_item> updates;

 public:
  transition(const std::string& name, trsys const* const ts,
             parser::AST const* const guard, parser::AST const* const rate,
             const std::vector<transition_update_item>& updates)
      : obj(name, ts), guard(guard), rate(rate), updates(updates) {}
  virtual ~transition() {}

  virtual bool is_chemreaction() const { return false; }

  parser::AST const* get_guard() const { return guard; }

  parser::AST const* get_rate() const { return rate; }

  const std::vector<transition_update_item>& get_updates() const {
    return updates;
  }
};

struct chemreaction_item {
  var const* v;
  double c;
};

class chemreaction : public transition {
 private:
  parser::AST const* opt_guard;

  const std::vector<chemreaction_item> reactants;
  const std::vector<chemreaction_item> products;
  const std::vector<chemreaction_item> stoichiometry;

  const bool isreversible;
  chemreaction const* const reverse;

 public:
  chemreaction(const std::string& name, trsys const* const ts,
               parser::AST const* const guard, parser::AST const* const rate,
               const std::vector<transition_update_item>& updates,
               const std::vector<chemreaction_item>& reactants,
               const std::vector<chemreaction_item>& products,
               const std::vector<chemreaction_item>& stoichiometry,
               parser::AST const* const opt_guard = nullptr,
               const bool isreversible = false,
               chemreaction const* const reverse = nullptr)
      : transition(name, ts, guard, rate, updates),
        opt_guard(opt_guard),
        reactants(reactants),
        products(products),
        stoichiometry(stoichiometry),
        isreversible(isreversible),
        reverse(reverse) {}

  bool is_chemreaction() const { return true; }

  parser::AST const* get_opt_guard() const { return opt_guard; }

  const std::vector<chemreaction_item>& get_reactants() const {
    return reactants;
  }

  const std::vector<chemreaction_item>& get_products() const {
    return products;
  }

  const std::vector<chemreaction_item>& get_stoichiometry() const {
    return stoichiometry;
  }

  bool is_reversible() const { return isreversible; }

  chemreaction const* get_reverse() const { return reverse; }
};
}
}

#endif
