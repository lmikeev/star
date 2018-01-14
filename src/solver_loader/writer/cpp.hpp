/*
 *  cpp.hpp
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

#ifndef SOLVER_LOADER_WRITER_CPP_HPP_
#define SOLVER_LOADER_WRITER_CPP_HPP_

#include "base.hpp"
#include "../../experiment.hpp"

namespace solver_loader {
namespace writer {

class cpp : public base {
 protected:
  std::ostream& write_unary(
      std::ostream& os, parser::unaryAST const* const et,
      solver_loader::base const* const sl = nullptr) const {
    static std::string op_str[] = {"+", "-", "!"};
    os << op_str[et->get_op()];
    return write_(os, et->get_term(), sl);
  }

  std::ostream& write_binary(
      std::ostream& os, parser::binaryAST const* const et,
      solver_loader::base const* const sl = nullptr) const {
    static std::string op_str[] = {"*", "/",  "+",  "-",  "<",  "<=",
                                   ">", ">=", "==", "!=", "&&", "||"};
    write_(os, et->get_lterm(), sl) << op_str[et->get_op()];
    return write_(os, et->get_rterm(), sl);
  }

  std::ostream& write_cnst_(std::ostream& os, model_info::cnst const* const c,
                            solver_loader::base const* const = nullptr) const {
    if (c->get_index_p() < 0) {
      if (c->get_type()->is_boolean()) {
        os << (static_cast<bool>(c->get_value()->get()) ? "true" : "false");
      } else if (c->get_type()->is_integer()) {
        os << static_cast<int>(c->get_value()->get());
      } else {
        os << static_cast<double>(c->get_value()->get());
      }
    } else {
      os << "c" << c->get_index_p();
    }
    return os;
  }

  std::ostream& write_var_(std::ostream& os, model_info::var const* const v,
                           solver_loader::base const* const = nullptr) const {
    os << "v" << v->get_index();
    return os;
  }

  std::ostream& write_monom_el_(
      std::ostream& os, model_info::cnst const* const c, const int p,
      solver_loader::base const* const sl = nullptr) const {
    if (p <= 3) {
      write_cnst_(os, c, sl);
      for (int i = 1; i < p; i++) {
        os << "*";
        write_cnst_(os, c, sl);
      }
    } else {
      os << "pow(";
      write_cnst_(os, c, sl);
      os << ',' << p << ')';
    }
    return os;
  }

  std::ostream& write_monom_el_(
      std::ostream& os, model_info::var const* const v, const int p,
      solver_loader::base const* const sl = nullptr) const {
    if (p <= 3) {
      write_var_(os, v, sl);
      for (int i = 1; i < p; i++) {
        os << "*";
        write_var_(os, v, sl);
      }
    } else {
      os << "pow(";
      write_var_(os, v, sl);
      os << ',' << p << ')';
    }
    return os;
  }
};
}
}

#endif
