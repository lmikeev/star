/*
 *  dump_samples.hpp
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

#ifndef SOLVER_LOADER_TASK_INFO_DUMP_SAMPLES_HPP_
#define SOLVER_LOADER_TASK_INFO_DUMP_SAMPLES_HPP_

#include "base.hpp"

namespace solver_loader {
namespace task_info {

class dump_samples : public base {
 private:
  const bool bytime;

 public:
  dump_samples(const bool bytime = false) : bytime(bytime) {}

  dump_samples(const dump_samples& ti) : dump_samples(ti.sortbytime()) {}

  type get_type() const { return TASK_DUMP_SAMPLES; }

  bool sortbytime() const { return bytime; }

  bool check(solver_loader::base const* const) {
    if (bytime) {
      name = "samples[time]";
    } else {
      name = "samples[trace]";
    }
    return true;
  }

  base* clone() const { return new dump_samples(*this); }
};
}
}

#endif
