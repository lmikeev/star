/*
 *  cnst.hpp
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

#ifndef SOLVER_LOADER_MODEL_INFO_CNST_HPP_
#define SOLVER_LOADER_MODEL_INFO_CNST_HPP_

#include "obj.hpp"
#include "val.hpp"

namespace solver_loader {

namespace model_info {

class type;

class cnst : public obj {
 private:
  type const* const typ;
  val* const value;

  mutable int index;
  mutable int index_p;

 public:
  cnst(const std::string& name, trsys const* const ts, type const* const typ,
       val* const value)
      : obj(name, ts), typ(typ), value(value), index(-1), index_p(-1) {}

  type const* get_type() const { return typ; }

  val* get_value() const { return value; }

  void set_index(const int i) const { index = i; }

  int get_index() const { return index; }

  void set_index_p(const int i) const { index_p = i; }

  int get_index_p() const { return index_p; }
};
}
}

#endif
