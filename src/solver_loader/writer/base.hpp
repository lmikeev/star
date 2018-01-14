/*
 *  base.hpp
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

#ifndef SOLVER_LOADER_WRITER_BASE_HPP_
#define SOLVER_LOADER_WRITER_BASE_HPP_

#include "../parser/ast.hpp"
#include "../base.hpp"

namespace solver_loader {
namespace writer {

class base {
 public:
  virtual ~base() {}

  virtual std::ostream& write_unary(
      std::ostream& os, parser::unaryAST const* const et,
      solver_loader::base const* const sl = nullptr) const {
    static std::string op_str[] = {"+", "-", "not"};
    os << op_str[et->get_op()] << " ";
    return write_(os, et->get_term(), sl);
  }

  virtual std::ostream& write_binary(
      std::ostream& os, parser::binaryAST const* const et,
      solver_loader::base const* const sl = nullptr) const {
    static std::string op_str[] = {"*", "/",  "+", "-",  "<",   "<=",
                                   ">", ">=", "=", "<>", "and", "or"};
    write_(os, et->get_lterm(), sl) << " " << op_str[et->get_op()] << " ";
    return write_(os, et->get_rterm(), sl);
  }

  virtual std::ostream& write_boolnum(
      std::ostream& os, parser::boolnumAST const* const et,
      solver_loader::base const* const = nullptr) const {
    os << (et->get_value() ? "true" : "false");
    return os;
  }

  virtual std::ostream& write_intnum(
      std::ostream& os, parser::intnumAST const* const et,
      solver_loader::base const* const = nullptr) const {
    os << et->get_value();
    return os;
  }

  virtual std::ostream& write_realnum(
      std::ostream& os, parser::realnumAST const* const et,
      solver_loader::base const* const = nullptr) const {
    os << et->get_value();
    return os;
  }

  std::ostream& write_cnst(
      std::ostream& os, parser::cnstAST const* const et,
      solver_loader::base const* const sl = nullptr) const {
    return write_cnst_(os, et->get_cnst(), sl);
  }

  std::ostream& write_var(std::ostream& os, parser::varAST const* const et,
                          solver_loader::base const* const sl = nullptr) const {
    return write_var_(os, et->get_var(), sl);
  }

  std::ostream& write_clock(
      std::ostream& os, parser::clockAST const* const et,
      solver_loader::base const* const sl = nullptr) const {
    return write_clock_(os, et->get_clock(), sl);
  }

  virtual std::ostream& write_polynom(
      std::ostream& os, parser::polynomAST const* const et,
      solver_loader::base const* const sl = nullptr) const {
    assert(!et->get_terms().empty());

    bool f = true;
    for (auto const& a : et->get_terms()) {
      assert(std::fabs(a.k) > std::numeric_limits<double>::epsilon());

      if (!f) {
        os << ((a.k > std::numeric_limits<double>::epsilon()) ? '+' : '-');
      }
      f = false;

      if (!a.c.empty() || !a.v.empty()) {
        bool g = true;
        if (std::fabs(std::fabs(a.k) - 1.0) >
            std::numeric_limits<double>::epsilon()) {
          os << std::fabs(a.k);
          g = false;
        }

        for (auto const& c : a.c) {
          if (!g) {
            os << '*';
          }
          g = false;

          assert(c.second > 0);

          write_monom_el_(os, c.first, c.second, sl);
        }

        for (auto const& v : a.v) {
          if (!g) {
            os << '*';
          }
          g = false;

          assert(v.second > 0);

          write_monom_el_(os, v.first, v.second, sl);
        }
      } else {
        os << std::fabs(a.k);
      }
    }
    return os;
  }

  virtual std::ostream& write_userfcall(
      std::ostream& os, parser::userfcallAST const* const et,
      solver_loader::base const* const sl = nullptr) const {
    os << et->get_f()->get_name() << "(";
    std::size_t ai = 0;
    for (auto const& a : et->get_args()) {
      if (ai) {
        os << ",";
      }
      a->write(os, this, sl);
      ai++;
    }
    os << ")";
    return os;
  }

  virtual std::ostream& write_userfarg(
      std::ostream& os, parser::userfargAST const* const et,
      solver_loader::base const* const = nullptr) const {
    os << (*et->get_fargs())[et->get_index()].name;
    return os;
  }

  virtual std::ostream& write_stdfcall(
      std::ostream& os, parser::stdfcallAST const* const et,
      solver_loader::base const* const sl = nullptr) const {
    if (et->get_stdf()->get_id() == parser::STDF_TIME) {
      os << "t";
    } else {
      os << et->get_stdf()->get_name() << "(";
      std::size_t ai = 0;
      for (auto const& a : et->get_args()) {
        if (ai) {
          os << ",";
        }
        a->write(os, this, sl);
        ai++;
      }
      os << ")";
    }
    return os;
  }

  virtual std::ostream& write_ifelse(
      std::ostream& os, parser::ifelseAST const* const et,
      solver_loader::base const* const sl = nullptr) const {
    et->get_cond()->write(os, this, sl) << " ? ";
    et->get_then_term()->write(os, this, sl) << " : ";
    return et->get_else_term()->write(os, this, sl);
  }

 protected:
  std::ostream& write_(std::ostream& os, parser::AST const* const et,
                       solver_loader::base const* const sl = nullptr) const {
    const bool leaf_ = et->is_leaf();
    if (!leaf_) {
      os << "(";
    }
    et->write(os, this, sl);
    if (!leaf_) {
      os << ")";
    }
    return os;
  }

  virtual std::ostream& write_cnst_(
      std::ostream& os, model_info::cnst const* const c,
      solver_loader::base const* const = nullptr) const {
    os << c->get_name();
    return os;
  }

  virtual std::ostream& write_var_(
      std::ostream& os, model_info::var const* const v,
      solver_loader::base const* const = nullptr) const {
    os << v->get_name();
    return os;
  }

  virtual std::ostream& write_clock_(
      std::ostream& os, model_info::clock const* const cl,
      solver_loader::base const* const = nullptr) const {
    os << cl->get_name();
    return os;
  }

  virtual std::ostream& write_monom_el_(
      std::ostream& os, model_info::cnst const* const c, const int p,
      solver_loader::base const* const sl = nullptr) const {
    write_cnst_(os, c, sl);
    if (p > 1) {
      os << "^" << p;
    }
    return os;
  }

  virtual std::ostream& write_monom_el_(
      std::ostream& os, model_info::var const* const v, const int p,
      solver_loader::base const* const sl = nullptr) const {
    write_var_(os, v, sl);
    if (p > 1) {
      os << "^" << p;
    }
    return os;
  }
};
}
}

#endif
