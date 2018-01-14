/*
 *  ic.hpp
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

#ifndef SOLVER_LOADER_MODEL_INFO_IC_HPP_
#define SOLVER_LOADER_MODEL_INFO_IC_HPP_

#include <vector>
#include "obj.hpp"
#include "val.hpp"

namespace solver_loader {

namespace parser {

class AST;
}

namespace model_info {

class var;

struct ic_s_i {
  var const* v;
  val* value;

  var const* get_var() { return v; }

  val* get_value() { return value; }
};

struct ic_s {
  double p;
  std::vector<double> dp_di;
  std::vector<double> d2p_di2;

  std::vector<ic_s_i> li;
};

class ic : public obj {
 private:
  std::vector<ic_s> ls;

 public:
  ic(const std::string& name, trsys const* const ts)
      : obj(name, ts), ls(std::vector<ic_s>()) {}

  ic(const std::string& name, trsys const* const ts,
     const std::vector<ic_s>& ls)
      : obj(name, ts), ls(ls) {}

  virtual ~ic() {}

  const std::vector<ic_s>& get_states() const { return ls; }

  std::vector<ic_s>& get_states() { return ls; }
};
}
}

#endif
