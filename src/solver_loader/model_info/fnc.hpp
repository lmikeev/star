/*
 *  fnc.hpp
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

#ifndef SOLVER_LOADER_MODEL_INFO_FNC_HPP_
#define SOLVER_LOADER_MODEL_INFO_FNC_HPP_

#include <vector>
#include "obj.hpp"

namespace solver_loader {

namespace parser {

class AST;
}

namespace model_info {

class type;

struct fnc_arg {
  type const* typ;
  std::string name;
};

class fnc : public obj {
 private:
  type const* const rtyp;
  std::vector<fnc_arg> const* const args;
  parser::AST const* const rexpr;

 public:
  fnc(const std::string& name, trsys const* const ts, type const* const rtyp)
      : obj(name, ts), rtyp(rtyp), args(nullptr), rexpr(nullptr) {}

  fnc(const std::string& name, trsys const* const ts, type const* const rtyp,
      std::vector<fnc_arg> const* const args, parser::AST const* const rexpr)
      : obj(name, ts), rtyp(rtyp), args(args), rexpr(rexpr) {}

  virtual ~fnc() {}

  type const* get_rtype() const { return rtyp; }

  std::vector<fnc_arg> const* get_args() const { return args; }

  parser::AST const* get_rexpr() const { return rexpr; }

  std::size_t get_narg() const { return args->size(); }
};
}
}

#endif
