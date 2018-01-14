/*
 *  val.hpp
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

#ifndef SOLVER_LOADER_MODEL_INFO_VAL_HPP_
#define SOLVER_LOADER_MODEL_INFO_VAL_HPP_

namespace solver_loader {

namespace parser {

class AST;
}

namespace model_info {

typedef double val_t;

class val {
 private:
  parser::AST const* const expr;

  val_t value;
  bool isset;

 public:
  val(parser::AST const* const expr) : expr(expr), value(0.0), isset(false) {}

  val(val_t v) : expr(nullptr), value(v), isset(true) {}

  parser::AST const* get_expr() const { return expr; }

  bool is_set() const { return isset; }

  val_t get() const { return value; }

  void set(const val_t v) {
    value = v;
    isset = true;
  }
};
}
}

#endif
