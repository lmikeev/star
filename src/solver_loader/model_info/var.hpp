/*
 *  var.hpp
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

#ifndef SOLVER_LOADER_MODEL_INFO_VAR_HPP_
#define SOLVER_LOADER_MODEL_INFO_VAR_HPP_

#include "obj.hpp"

namespace solver_loader {

namespace parser {

class AST;
}

namespace model_info {

class cnst;

class var : public obj {
 private:
  type const* const typ;

  cnst const* const compartment;

  mutable parser::AST const* init;
  mutable int index;
  mutable int index_i, index_e;

 public:
  var(const std::string& name, trsys const* const ts, type const* const typ,
      cnst const* const compartment = nullptr)
      : obj(name, ts),
        typ(typ),
        compartment(compartment),
        init(nullptr),
        index(-1),
        index_i(-1),
        index_e(-1) {}
  virtual ~var() {}

  type const* get_type() const { return typ; }

  cnst const* get_compartment() const { return compartment; }

  void set_init(parser::AST const* const et) const { init = et; }

  parser::AST const* get_init() const { return init; }

  void set_index(const int i) const { index = i; }

  int get_index() const { return index; }

  void set_index_i(const int i) const { index_i = i; }

  int get_index_i() const { return index_i; }

  void set_index_e(const int i) const { index_e = i; }

  int get_index_e() const { return index_e; }
};
}
}

#endif
