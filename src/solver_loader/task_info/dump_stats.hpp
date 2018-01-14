/*
 *  dump_stats.hpp
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

#ifndef SOLVER_LOADER_TASK_INFO_DUMP_STATS_HPP_
#define SOLVER_LOADER_TASK_INFO_DUMP_STATS_HPP_

#include "base.hpp"

namespace solver_loader {
namespace task_info {

class dump_stats : public base {
 private:
  const std::vector<std::string> statnames;

 public:
  dump_stats() {}

  dump_stats(const std::vector<std::string>& statnames)
      : statnames(statnames) {}

  dump_stats(const dump_stats& ti) : dump_stats(ti.get_statnames()) {}

  type get_type() const { return TASK_DUMP_STATS; }

  const std::vector<std::string>& get_statnames() const { return statnames; }

  bool check(solver_loader::base const* const) {
    name = "stats";
    return true;
  }

  base* clone() const { return new dump_stats(*this); }
};
}
}

#endif
