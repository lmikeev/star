/*
 *  obj.hpp
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

#ifndef SOLVER_LOADER_MODEL_INFO_OBJ_HPP_
#define SOLVER_LOADER_MODEL_INFO_OBJ_HPP_

#include <string>

namespace solver_loader {
namespace model_info {

class trsys;

class obj {
 private:
  const std::string name;

 protected:
  trsys const* const ts;

 public:
  obj(const std::string& name, trsys const* const ts) : name(name), ts(ts) {}

  const std::string& get_name() const { return name; }

  bool is_unnamed() const { return name == ""; }

  trsys const* get_trsys() const { return ts; }
};
}
}

#endif
