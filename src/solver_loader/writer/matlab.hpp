/*
 *  matlab.hpp
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

#ifndef SOLVER_LOADER_WRITER_MATLAB_HPP_
#define SOLVER_LOADER_WRITER_MATLAB_HPP_

#include "base.hpp"

namespace solver_loader {
namespace writer {

class matlab : public base {
 protected:
  std::ostream& write_cnst_(std::ostream& os, model_info::cnst const* const c,
                            solver_loader::base const* const = nullptr) const {
    os << "c" << c->get_index();
    return os;
  }

  std::ostream& write_var_(std::ostream& os, model_info::var const* const v,
                           solver_loader::base const* const = nullptr) const {
    os << "v" << v->get_index();
    return os;
  }

  std::ostream& write_stdfcall(
      std::ostream& os, parser::stdfcallAST const* const et,
      solver_loader::base const* const sl = nullptr) const {
    if (et->get_stdf()->get_id() == parser::STDF_POW) {
      base::write_(os, et->get_args().front(), sl);
      os << "^";
      base::write_(os, et->get_args()[1], sl);
      return os;
    }
    return base::write_stdfcall(os, et, sl);
  }
};
}
}

#endif
