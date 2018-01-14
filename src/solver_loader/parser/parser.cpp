/*
 *  parser.cpp
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

#include <stack>
#include <cmath>
#include <boost/algorithm/string.hpp>
#include "parser.hpp"
#include "../base.hpp"
#include "../model_info/trsys.hpp"
#include "../writer/cpp.hpp"
#include "../writer/matlab.hpp"
#include "../../experiment.hpp"

namespace solver_loader {

namespace model_info {

type const* const trsys::typeBoolean(new type_boolean("boolean", nullptr));
type const* const trsys::typeInteger(new type_integer("integer", nullptr));
type const* const trsys::typeSpecies(new type_species("species", nullptr));
type const* const trsys::typeReal(new type_real("real", nullptr));
}

namespace parser {

const stdfnc sfncs[] = {
    {STDF_TIME, "time_", 0, model_info::T_REAL},
    {STDF_MASS_ACTION, "mass_action", 1, model_info::T_REAL},
    {STDF_ABS, "abs", 1, model_info::T_AUTO},
    {STDF_MIN, "min", -1, model_info::T_AUTO},
    {STDF_MAX, "max", -1, model_info::T_AUTO},
    {STDF_FLOOR, "floor", 1, model_info::T_INTEGER},
    {STDF_CEIL, "ceil", 1, model_info::T_INTEGER},
    {STDF_MOD, "mod", 1, model_info::T_AUTO},
    {STDF_EXP, "exp", 1, model_info::T_REAL},
    {STDF_SQRT, "sqrt", 1, model_info::T_REAL},
    {STDF_POW, "pow", 2, model_info::T_AUTO},
    {STDF_LN, "ln", 1, model_info::T_REAL},
    {STDF_LOG10, "log10", 1, model_info::T_REAL},
    {STDF_LOG, "log", 2, model_info::T_REAL},
    {STDF_SIN, "sin", 1, model_info::T_REAL},
    {STDF_COS, "cos", 1, model_info::T_REAL},
    {STDF_TAN, "tan", 1, model_info::T_REAL},
    {STDF_ASIN, "asin", 1, model_info::T_REAL},
    {STDF_ACOS, "acos", 1, model_info::T_REAL},
    {STDF_ATAN, "atan", 1, model_info::T_REAL},
    {STDF_ATAN2, "atan2", 2, model_info::T_REAL},
    {STDF_SINH, "sinh", 1, model_info::T_REAL},
    {STDF_COSH, "cosh", 1, model_info::T_REAL},
    {STDF_TANH, "tanh", 1, model_info::T_REAL},
    {STDF_ASINH, "asinh", 1, model_info::T_REAL},
    {STDF_ACOSH, "acosh", 1, model_info::T_REAL},
    {STDF_ATANH, "atanh", 1, model_info::T_REAL}};

AST const* const AST::et_0(new parser::polynomAST(0.0));
AST const* const AST::et_1(new parser::polynomAST(1.0));

writer::base const* const std_writer(new writer::base());
writer::base const* const cpp_writer(new writer::cpp());
writer::base const* const matlab_writer(new writer::matlab());

stdfnc const* AST::get_stdf(const stdfnc_id id) { return &sfncs[id]; }

std::ostream& unaryAST::write(std::ostream& os, writer::base const* const w,
                              solver_loader::base const* const sl) const {
  return w->write_unary(os, this, sl);
}

std::ostream& binaryAST::write(std::ostream& os, writer::base const* const w,
                               solver_loader::base const* const sl) const {
  return w->write_binary(os, this, sl);
}

std::ostream& boolnumAST::write(std::ostream& os, writer::base const* const w,
                                solver_loader::base const* const sl) const {
  return w->write_boolnum(os, this, sl);
}

std::ostream& intnumAST::write(std::ostream& os, writer::base const* const w,
                               solver_loader::base const* const sl) const {
  return w->write_intnum(os, this, sl);
}

std::ostream& realnumAST::write(std::ostream& os, writer::base const* const w,
                                solver_loader::base const* const sl) const {
  return w->write_realnum(os, this, sl);
}

std::ostream& cnstAST::write(std::ostream& os, writer::base const* const w,
                             solver_loader::base const* const sl) const {
  return w->write_cnst(os, this, sl);
}

std::ostream& varAST::write(std::ostream& os, writer::base const* const w,
                            solver_loader::base const* const sl) const {
  return w->write_var(os, this, sl);
}

std::ostream& clockAST::write(std::ostream& os, writer::base const* const w,
                              solver_loader::base const* const sl) const {
  return w->write_clock(os, this, sl);
}

std::ostream& polynomAST::write(std::ostream& os, writer::base const* const w,
                                solver_loader::base const* const sl) const {
  return w->write_polynom(os, this, sl);
}

std::ostream& userfcallAST::write(std::ostream& os, writer::base const* const w,
                                  solver_loader::base const* const sl) const {
  return w->write_userfcall(os, this, sl);
}

std::ostream& userfargAST::write(std::ostream& os, writer::base const* const w,
                                 solver_loader::base const* const sl) const {
  return w->write_userfarg(os, this, sl);
}

std::ostream& stdfcallAST::write(std::ostream& os, writer::base const* const w,
                                 solver_loader::base const* const sl) const {
  return w->write_stdfcall(os, this, sl);
}

std::ostream& ifelseAST::write(std::ostream& os, writer::base const* const w,
                               solver_loader::base const* const sl) const {
  return w->write_ifelse(os, this, sl);
}

bool set_error(state& p, char const* const err) {
  sprintf(p.err_buf, "%s", err);
  return false;
}

bool set_error_expected(state& p, const char c) {
  sprintf(p.err_buf, "'%c' expected", c);
  return false;
}

char* get_error(state& p) { return p.err_buf; }

int get_line_number(state& p) { return p.line_number; }

int get_col_number(state& p) { return p.col_number; }

void get_next(state& p) {
  p.last = *p.s++;
  p.line += p.last;
  p.col_number++;
}

void print_line(state& p) {
#ifndef STAR_WEB_INTERFACE

  printf("%d\t%s", p.line_number, p.line.c_str());
#endif
  p.line = "";
}

void chk_newline(state& p) {
  if (p.last == '\r' || p.last == '\n') {
    if (p.last == '\r' && *p.s == '\n') {
      get_next(p);
    }
    print_line(p);
    p.line_number++;
    p.col_number = 0;
  }
}

token_t get_token(
    state& p, model_info::trsys const* const context,
    std::vector<model_info::fnc_arg> const* const fargs = nullptr) {
  while (isspace(p.last)) {
    chk_newline(p);
    get_next(p);
  }

  if (p.last == '*') {
    get_next(p);
    return TOK_MULT;
  }

  if (p.last == '/') {
    get_next(p);
    if (p.last == '/') {
      get_next(p);
      while (p.last != '\r' && p.last != '\n') {
        get_next(p);
      }
      if (p.last == '\r' && *p.s == '\n') {
        get_next(p);
      }
      return get_token(p, context, fargs);
    }

    if (p.last == '*') {
      do {
        chk_newline(p);
        get_next(p);
      } while (p.last != '\0' && (p.last != '*' || *p.s != '/'));

      if (p.last == '*' && *p.s == '/') {
        get_next(p);
        get_next(p);
      } else {
        sprintf(p.err_buf, "unexpected end of file");
        return TOK_UNKNOWN;
      }

      return get_token(p, context, fargs);
    }
    return TOK_DIV;
  }

  if (p.last == '+') {
    get_next(p);
    return TOK_PLUS;
  }

  if (p.last == '-') {
    get_next(p);
    if (p.last == '>') {
      get_next(p);
      return TOK_RARROW;
    }
    return TOK_MINUS;
  }

  if (p.last == '<') {
    get_next(p);
    if (p.last == '=') {
      get_next(p);
      return TOK_LEQ;
    }
    if (p.last == '>') {
      get_next(p);
      return TOK_INEQ;
    }
    if (p.last == '-' && *p.s == '>') {
      get_next(p);
      get_next(p);
      return TOK_DARROW;
    }
    return TOK_LESS;
  }

  if (p.last == '>') {
    get_next(p);
    if (p.last == '=') {
      get_next(p);
      return TOK_GREQ;
    }
    return TOK_GR;
  }

  if (p.last == '=') {
    get_next(p);
    if (p.last == '<') {
      get_next(p);
      return TOK_LEQ;
    }
    if (p.last == '>') {
      get_next(p);
      return TOK_GREQ;
    }
    return TOK_EQ;
  }

  if (p.last == '!') {
    get_next(p);
    return TOK_EXCL;
  }

  if (p.last == '(') {
    get_next(p);
    return TOK_R_OP;
  }

  if (p.last == ')') {
    get_next(p);
    return TOK_R_CL;
  }

  if (p.last == '{') {
    get_next(p);
    return TOK_C_OP;
  }

  if (p.last == '}') {
    get_next(p);
    return TOK_C_CL;
  }

  if (p.last == '[') {
    get_next(p);
    return TOK_S_OP;
  }

  if (p.last == ']') {
    get_next(p);
    return TOK_S_CL;
  }

  if (p.last == '?') {
    get_next(p);
    return TOK_QUESTION;
  }

  if (p.last == ',') {
    get_next(p);
    return TOK_COMMA;
  }

  if (p.last == '.') {
    get_next(p);
    if (p.last == '.') {
      get_next(p);
      return TOK_DDOT;
    }
    return TOK_DOT;
  }

  if (p.last == ':') {
    get_next(p);
    return TOK_COLON;
  }

  if (p.last == ';') {
    get_next(p);
    return TOK_SEMICOLON;
  }

  if (p.last == '@') {
    get_next(p);
    return TOK_AT;
  }

  if (p.last == '~') {
    get_next(p);
    return TOK_TILDE;
  }

  if (isdigit(p.last)) {
    std::string num_str;
    do {
      num_str += p.last;
      get_next(p);
    } while (isdigit(p.last));

    if (p.last == '.' && *p.s != '.') {
      num_str += p.last;
      get_next(p);

      if (isdigit(p.last)) {
        do {
          num_str += p.last;
          get_next(p);
        } while (isdigit(p.last));

        if (p.last == 'e') {
          num_str += p.last;
          get_next(p);

          if (p.last == '-') {
            num_str += p.last;
            get_next(p);
          }

          if (isdigit(p.last)) {
            do {
              num_str += p.last;
              get_next(p);
            } while (isdigit(p.last));
          } else {
            set_error(p, "invalid number format");
            return TOK_UNKNOWN;
          }
        }
      } else {
        set_error(p, "invalid number format");
        return TOK_UNKNOWN;
      }

      p.realnum_val = strtod(num_str.c_str(), 0);
      return TOK_REALNUM;
    }

    if (p.last == 'e') {
      num_str += p.last;
      get_next(p);

      if (p.last == '-') {
        num_str += p.last;
        get_next(p);
      }

      if (isdigit(p.last)) {
        do {
          num_str += p.last;
          get_next(p);
        } while (isdigit(p.last));
      } else {
        set_error(p, "invalid number format");
        return TOK_UNKNOWN;
      }

      p.realnum_val = strtod(num_str.c_str(), 0);
      return TOK_REALNUM;
    }

    p.intnum_val = strtol(num_str.c_str(), 0, 10);
    return TOK_INTNUM;
  }

  if (isalpha(p.last)) {
    std::string in_str;
    do {
      in_str += p.last;
      get_next(p);
    } while (isalnum(p.last) || p.last == '_');

    std::string str = in_str;

    if (str == "true") {
      p.boolnum_val = true;
      return TOK_BOOLNUM;
    }

    if (str == "false") {
      p.boolnum_val = false;
      return TOK_BOOLNUM;
    }

    if (str == "boolean") {
      p.t = model_info::trsys::typeBoolean;
      return TOK_STDTYPE;
    }

    if (str == "integer") {
      p.t = model_info::trsys::typeInteger;
      return TOK_STDTYPE;
    }

    if (str == "species") {
      p.t = model_info::trsys::typeSpecies;
      return TOK_STDTYPE;
    }

    if (str == "real") {
      p.t = model_info::trsys::typeReal;
      return TOK_STDTYPE;
    }

    if (str == "and") {
      return TOK_AND;
    }
    if (str == "or") {
      return TOK_OR;
    }
    if (str == "not") {
      return TOK_NOT;
    }

    if (str == "type") {
      return TOK_TYPE;
    }
    if (str == "const") {
      return TOK_CONST;
    }
    if (str == "var") {
      return TOK_VAR;
    }
    if (str == "molecule") {
      return TOK_MOLECULE;
    }
    if (str == "clock") {
      return TOK_CLOCK;
    }
    if (str == "function") {
      return TOK_FUNCTION;
    }
    if (str == "guarded_commands") {
      return TOK_GUARDEDCMDS;
    }
    if (str == "chemical_reactions") {
      return TOK_CHEMREACTIONS;
    }
    if (str == "pta") {
      return TOK_PTA;
    }
    if (str == "init") {
      return TOK_INIT;
    }
    if (str == "end") {
      return TOK_END;
    }

    if (str == "time" || str == "mass_action" || str == "abs" || str == "min" ||
        str == "max" || str == "floor" || str == "ceil" || str == "mod" ||
        str == "exp" || str == "sqrt" || str == "pow" || str == "ln" ||
        str == "log10" || str == "log" || str == "sin" || str == "cos" ||
        str == "tan" || str == "asin" || str == "acos" || str == "atan" ||
        str == "atan2" || str == "sinh" || str == "cosh" || str == "tanh" ||
        str == "asinh" || str == "acosh" || str == "atanh") {
      auto const& sfi =
          std::find_if(sfncs, sfncs + sizeof(sfncs) / sizeof(stdfnc),
                       [&](const stdfnc& sf) { return sf.name == str; });
      p.stdf = sfi;
      return TOK_STDFNC;
    }

    p.id_str = in_str;

    if (p.last == '\'') {
      get_next(p);
      return TOK_ID_PRIME;
    }
    return TOK_ID;
  }

  if (p.last == '\0') {
    return TOK_EOF;
  }

  sprintf(p.err_buf, "unrecognized symbol '%c'", p.last);
  return TOK_UNKNOWN;
}

void get_next_token(
    state& p, model_info::trsys const* const context,
    std::vector<model_info::fnc_arg> const* const fargs = nullptr) {
  p.tok = get_token(p, context, fargs);
}

bool tok_is_cnst(state& p, model_info::trsys const* const context) {
  if (p.tok == TOK_ID) {
    p.c = context->findcnst(p.id_str);
    if (p.c != nullptr) {
      return true;
    }
  }
  return false;
}

bool tok_is_type(state& p, model_info::trsys const* const context) {
  if (p.tok == TOK_STDTYPE) {
    return true;
  }
  if (p.tok == TOK_ID) {
    p.t = context->findtype(p.id_str);
    if (p.t != nullptr) {
      return true;
    }
  }
  return false;
}

bool tok_is_fnc(state& p, model_info::trsys const* const context) {
  if (p.tok == TOK_ID) {
    p.f = context->findfnc(p.id_str);
    if (p.f != nullptr) {
      return true;
    }
  }
  return false;
}

bool tok_is_var(state& p, model_info::trsys const* const context) {
  if (p.tok == TOK_ID) {
    p.v = context->findvar(p.id_str);
    if (p.v != nullptr) {
      return true;
    }
  }
  return false;
}

bool tok_is_var_prime(state& p, model_info::trsys const* const context) {
  if (p.tok == TOK_ID_PRIME) {
    p.v = context->findvar(p.id_str);
    if (p.v != nullptr) {
      return true;
    }
  }
  return false;
}

bool tok_is_farg(
    state& p, model_info::trsys const* const,
    std::vector<model_info::fnc_arg> const* const fargs = nullptr) {
  if (p.tok == TOK_ID && fargs != nullptr) {
    auto const& ai = std::find_if(
        fargs->begin(), fargs->end(),
        [&](const model_info::fnc_arg& a) { return a.name == p.id_str; });

    if (ai != fargs->end()) {
      p.farg_i = ai - fargs->begin();
      return true;
    }
  }
  return false;
}

bool is_expr_const(AST const* const et) {
  if (et->is_number()) {
  } else if (et->is_cnst()) {
  } else if (et->is_unary()) {
    if (!is_expr_const(static_cast<unaryAST const*>(et)->get_term())) {
      return false;
    }
  } else if (et->is_binary()) {
    if (!is_expr_const(static_cast<binaryAST const*>(et)->get_lterm())) {
      return false;
    }
    if (!is_expr_const(static_cast<binaryAST const*>(et)->get_rterm())) {
      return false;
    }
  } else if (et->is_fcall()) {
    if (et->is_stdfcall() &&
        static_cast<stdfcallAST const*>(et)->get_stdf()->get_id() ==
            STDF_MASS_ACTION) {
      return false;
    }

    for (auto const& a : static_cast<callfAST const*>(et)->get_args()) {
      if (!is_expr_const(a)) {
        return false;
      }
    }
  } else if (et->is_ifelse()) {
    ifelseAST const* const iet = static_cast<ifelseAST const*>(et);
    if (!is_expr_const(iet->get_cond())) {
      return false;
    }
    if (!is_expr_const(iet->get_then_term())) {
      return false;
    }
    if (!is_expr_const(iet->get_else_term())) {
      return false;
    }
  } else {
    return false;
  }

  return true;
}

bool chk_expr_const(state& p, AST const* const et) {
  if (!is_expr_const(et)) {
    return set_error(p, "constant expression expected");
  }
  return true;
}

bool is_expr_boolean(AST const* const et) {
  if (et->is_boolnum()) {
  } else if (et->is_cnst() &&
             static_cast<cnstAST const*>(et)
                 ->get_cnst()
                 ->get_type()
                 ->is_boolean()) {
  } else if (et->is_var() &&
             static_cast<varAST const*>(et)
                 ->get_var()
                 ->get_type()
                 ->is_boolean()) {
  } else if (et->is_unary() &&
             static_cast<unaryAST const*>(et)->get_op_tok() == TOK_NOT) {
  } else if (et->is_binary()) {
    const token_t op_tok = static_cast<binaryAST const*>(et)->get_op_tok();

    if (op_tok != TOK_LESS && op_tok != TOK_LEQ && op_tok != TOK_GREQ &&
        op_tok != TOK_GR && op_tok != TOK_EQ && op_tok != TOK_INEQ &&
        op_tok != TOK_AND && op_tok != TOK_OR) {
      return false;
    }
  } else if (et->is_userfcall() &&
             static_cast<userfcallAST const*>(et)
                     ->get_f()
                     ->get_rtype()
                     ->get_type() == model_info::T_BOOLEAN) {
  } else if (et->is_stdfcall() &&
             static_cast<stdfcallAST const*>(et)->get_stdf()->get_rtype() ==
                 model_info::T_BOOLEAN) {
  } else if (et->is_ifelse()) {
    ifelseAST const* const iet = static_cast<ifelseAST const*>(et);
    if (!is_expr_boolean(iet->get_then_term())) {
      return false;
    }
    if (!is_expr_boolean(iet->get_else_term())) {
      return false;
    }
  } else {
    return false;
  }

  return true;
}

bool chk_expr_boolean(state&, AST const* const) { return true; }

bool is_expr_integer(AST const* const et) {
  if (et->is_intnum()) {
  } else if (et->is_cnst() &&
             static_cast<cnstAST const*>(et)
                 ->get_cnst()
                 ->get_type()
                 ->is_integer()) {
  } else if (et->is_var() &&
             static_cast<varAST const*>(et)
                 ->get_var()
                 ->get_type()
                 ->is_integer()) {
  } else if (et->is_unary() &&
             (static_cast<unaryAST const*>(et)->get_op_tok() == TOK_PLUS ||
              static_cast<unaryAST const*>(et)->get_op_tok() == TOK_MINUS)) {
  } else if (et->is_binary()) {
    const token_t op_tok = static_cast<binaryAST const*>(et)->get_op_tok();
    if (op_tok != TOK_MULT && op_tok != TOK_DIV && op_tok != TOK_PLUS &&
        op_tok != TOK_MINUS) {
      return false;
    }
  } else if (et->is_userfcall() &&
             static_cast<userfcallAST const*>(et)
                     ->get_f()
                     ->get_rtype()
                     ->get_type() == model_info::T_INTEGER) {
  } else if (et->is_stdfcall() &&
             static_cast<stdfcallAST const*>(et)->get_stdf()->get_rtype() ==
                 model_info::T_INTEGER) {
  } else if (et->is_ifelse()) {
    ifelseAST const* const iet = static_cast<ifelseAST const*>(et);
    if (!is_expr_integer(iet->get_then_term())) {
      return false;
    }
    if (!is_expr_integer(iet->get_else_term())) {
      return false;
    }
  } else {
    return false;
  }

  return true;
}

bool chk_expr_integer(state& p, AST const* const et) {
  if (!is_expr_integer(et)) {
    return set_error(p, "integer expression expected");
  }
  return true;
}

bool is_expr_real(AST const* const et) { return !is_expr_boolean(et); }

bool chk_expr_real(state& p, AST const* const et) {
  if (!is_expr_real(et)) {
    return set_error(p, "unexpected boolean expression");
  }
  return true;
}

bool chk_updates_(state& p, AST const* const et,
                  std::vector<model_info::transition_update_item>& updates) {
  if (et->is_binary()) {
    binaryAST const* const etb = static_cast<binaryAST const*>(et);
    if (etb->get_op_tok() == TOK_AND) {
      if (!chk_updates_(p, etb->get_lterm(), updates)) {
        return false;
      }
      if (!chk_updates_(p, etb->get_rterm(), updates)) {
        return false;
      }
    } else if (etb->get_op_tok() == TOK_EQ) {
      if (!etb->get_lterm()->is_var()) {
        return set_error(
            p, "variable expected on the left side of the update function");
      }

      parser::varAST const* etv =
          static_cast<parser::varAST const*>(etb->get_lterm());
      if (!etv->is_prime()) {
        return set_error(p, "' expected");
      }

      auto const& uit =
          std::find_if(updates.begin(), updates.end(),
                       [&](const model_info::transition_update_item& ui) {
            return ui.v == etv->get_var();
          });

      if (uit != updates.end()) {
        sprintf(p.err_buf, "duplicate update for '%s'",
                etv->get_var()->get_name().c_str());
        return false;
      }

      if (!chk_expr_integer(p, etb->get_rterm())) {
        return false;
      }

      model_info::transition_update_item u;
      u.v = etv->get_var();
      u.u = etb->get_rterm();
      updates.push_back(u);
    } else {
      return set_error(p, "invalid update function");
    }
  } else {
    return set_error(p, "invalid update function");
  }

  return true;
}

bool chk_type(state& p, model_info::type const* typ, AST const* const et) {
  if (typ->is_boolean()) {
    if (!chk_expr_boolean(p, et)) {
      return false;
    }
  } else if (typ->is_integer()) {
    if (!chk_expr_integer(p, et)) {
      return false;
    }
  } else {
    if (!chk_expr_real(p, et)) {
      return false;
    }
  }
  return true;
}

bool chk_auto_type(state& p, model_info::type const*& typ,
                   AST const* const et) {
  if (typ == nullptr) {
    if (is_expr_boolean(et)) {
      typ = model_info::trsys::typeBoolean;
    } else if (is_expr_integer(et)) {
      typ = model_info::trsys::typeInteger;
    } else {
      typ = model_info::trsys::typeReal;
    }
  } else {
    if (!chk_type(p, typ, et)) {
      return false;
    }
  }
  return true;
}

bool chk_id(state& p) {
  if (p.tok != TOK_ID) {
    return set_error(p, "identifier expected");
  }
  return true;
}

bool chk_var(state& p, model_info::trsys const* const context) {
  if (!tok_is_var(p, context)) {
    return set_error(p, "variable expected");
  }
  return true;
}

bool chk_updates(state& p, AST const* const et,
                 std::vector<model_info::transition_update_item>& updates) {
  if (et->is_boolnum() &&
      static_cast<boolnumAST const*>(et)->get_value() == true) {
    return true;
  }

  return chk_updates_(p, et, updates);
}

bool parse_expr(state& p, model_info::trsys const* const context,
                AST const*& et,
                std::vector<model_info::fnc_arg> const* const fargs);

bool parse_args(state& p, model_info::trsys const* const context,
                const int narg, std::vector<AST const*>& args,
                std::vector<model_info::fnc_arg> const* const fargs = nullptr) {
  if (narg || narg < 0) {
    if (p.tok != TOK_R_OP) {
      return set_error_expected(p, '(');
    }
    get_next_token(p, context, fargs);
  }

  for (int i = 0; i < narg || narg < 0; i++) {
    if (i) {
      if (p.tok != TOK_COMMA) {
        return set_error_expected(p, ',');
      }
      get_next_token(p, context, fargs);
    }
    AST const* arg;
    if (!parse_expr(p, context, arg, fargs)) {
      return false;
    }
    args.push_back(arg);

    if (narg < 0 && p.tok == TOK_R_CL) {
      break;
    }
  }

  if (narg || narg < 0) {
    if (p.tok != TOK_R_CL) {
      return set_error_expected(p, ')');
    }
    get_next_token(p, context, fargs);
  }

  return true;
}

bool parse_expr_term(
    state& p, model_info::trsys const* const context, AST const*& et,
    std::vector<model_info::fnc_arg> const* const fargs = nullptr) {
  switch (p.tok) {
    case TOK_BOOLNUM: {
      et = new boolnumAST(p.boolnum_val);
      get_next_token(p, context, fargs);
      break;
    }

    case TOK_INTNUM: {
      et = new intnumAST(p.intnum_val);
      get_next_token(p, context, fargs);
      break;
    }

    case TOK_REALNUM: {
      et = new realnumAST(p.realnum_val);
      get_next_token(p, context, fargs);
      break;
    }

    case TOK_ID: {
      if (tok_is_cnst(p, context)) {
        et = new cnstAST(p.c);
        get_next_token(p, context, fargs);
      } else if (tok_is_var(p, context)) {
        if (p.v->get_compartment() != nullptr) {
          et = new binaryAST(parser::TOK_DIV, new varAST(p.v),
                             new cnstAST(p.v->get_compartment()));
        } else {
          et = new varAST(p.v);
        }
        get_next_token(p, context, fargs);
      } else if (tok_is_fnc(p, context)) {
        model_info::fnc const* const f = p.f;
        get_next_token(p, context, fargs);
        std::vector<AST const*> args;
        if (!parse_args(p, context, f->get_narg(), args, fargs)) {
          return false;
        }
        et = new userfcallAST(f, args);
      } else if (tok_is_farg(p, context, fargs)) {
        et = new userfargAST(fargs, p.farg_i);
        get_next_token(p, context, fargs);
      } else {
        sprintf(p.err_buf, "unknown identifier '%s'", p.id_str.c_str());
        return false;
      }
      break;
    }

    case TOK_ID_PRIME: {
      if (tok_is_var_prime(p, context)) {
        et = new varAST(p.v, true);
        get_next_token(p, context, fargs);
      }
      break;
    }

    case TOK_STDFNC: {
      stdfnc const* const stdf = p.stdf;
      get_next_token(p, context, fargs);
      std::vector<AST const*> args;
      if (!parse_args(p, context, stdf->get_narg(), args, fargs)) {
        return false;
      }
      et = new stdfcallAST(stdf, args);
      break;
    }

    case TOK_R_OP: {
      get_next_token(p, context, fargs);
      if (!parse_expr(p, context, et, fargs)) {
        return false;
      }
      if (p.tok != TOK_R_CL) {
        return set_error_expected(p, ')');
      }
      get_next_token(p, context, fargs);
      break;
    }

    default: { return set_error(p, "unexpected token"); }
  }

  return true;
}

int get_unary_op(const token_t op_tok) {
  switch (op_tok) {
    case TOK_PLUS:
      return 0;

    case TOK_MINUS:
      return 1;

    case TOK_NOT:
      return 2;

    default:
      break;
  }

  return -1;
}

int get_binary_op_prec(const token_t op_tok) {
  switch (op_tok) {
    case TOK_MULT:
    case TOK_DIV:
      return 60;

    case TOK_PLUS:
    case TOK_MINUS:
      return 50;

    case TOK_LESS:
    case TOK_LEQ:
    case TOK_GR:
    case TOK_GREQ:
      return 40;

    case TOK_EQ:
    case TOK_INEQ:
      return 30;

    case TOK_AND:
      return 20;

    case TOK_OR:
      return 10;

    default:
      break;
  }

  return -1;
}

bool parse_expr_unary(
    state& p, model_info::trsys const* const context, AST const*& et,
    std::vector<model_info::fnc_arg> const* const fargs = nullptr) {
  const token_t op_tok = p.tok;
  const int op = get_unary_op(op_tok);
  if (op >= 0) {
    get_next_token(p, context, fargs);
  }

  if (!parse_expr_term(p, context, et, fargs)) {
    return false;
  }

  if (op >= 0) {
    if (op_tok == TOK_NOT) {
      if (!chk_expr_boolean(p, et)) {
        return false;
      }
    } else {
      if (!chk_expr_real(p, et)) {
        return false;
      }
    }

    et = new unaryAST(op_tok, et);
  }

  return true;
}

bool parse_expr_rterm(
    state& p, model_info::trsys const* const context, const int expr_prec,
    AST const* lterm, AST const*& et,
    std::vector<model_info::fnc_arg> const* const fargs = nullptr) {
  while (1) {
    const int op_prec = get_binary_op_prec(p.tok);

    if (op_prec < expr_prec || p.tok == TOK_EOF) {
      et = lterm;
      return true;
    }

    const token_t op_tok = p.tok;
    get_next_token(p, context, fargs);

    AST const* rterm;
    if (!parse_expr_unary(p, context, rterm, fargs)) {
      return false;
    }

    const int next_op_prec = get_binary_op_prec(p.tok);

    if (op_prec < next_op_prec) {
      if (!parse_expr_rterm(p, context, op_prec + 1, rterm, rterm, fargs)) {
        return false;
      }
    }

    switch (op_tok) {
      case TOK_AND:
      case TOK_OR: {
        if (!chk_expr_boolean(p, lterm) || !chk_expr_boolean(p, rterm)) {
          return false;
        }
        break;
      }

      case TOK_EQ:
      case TOK_INEQ: {
        if (is_expr_boolean(lterm)) {
          if (!chk_expr_boolean(p, rterm)) {
            return false;
          }
        } else if (is_expr_integer(lterm)) {
          if (!chk_expr_integer(p, rterm)) {
            return false;
          }
        } else {
          if (!chk_expr_real(p, rterm)) {
            return false;
          }
        }
        break;
      }

      default: {
        if (!chk_expr_real(p, lterm) || !chk_expr_real(p, rterm)) {
          return false;
        }
        break;
      }
    }

    lterm = new binaryAST(op_tok, lterm, rterm);
  }

  return true;
}

bool parse_expr_(
    state& p, model_info::trsys const* const context, AST const*& et,
    std::vector<model_info::fnc_arg> const* const fargs = nullptr) {
  AST const* lterm;
  if (!parse_expr_unary(p, context, lterm, fargs)) {
    return false;
  }

  if (!parse_expr_rterm(p, context, 0, lterm, et, fargs)) {
    return false;
  }

  return true;
}

bool parse_expr(state& p, model_info::trsys const* const context,
                AST const*& et,
                std::vector<model_info::fnc_arg> const* const fargs = nullptr) {
  if (!parse_expr_(p, context, et, fargs)) {
    return false;
  }

  if (p.tok == TOK_QUESTION) {
    if (!chk_expr_boolean(p, et)) {
      return false;
    }

    get_next_token(p, context, fargs);
    AST const* then_term;
    if (!parse_expr(p, context, then_term, fargs)) {
      return false;
    }
    if (p.tok != TOK_COLON) {
      return set_error_expected(p, ':');
    }
    get_next_token(p, context, fargs);
    AST const* else_term;
    if (!parse_expr(p, context, else_term, fargs)) {
      return false;
    }
    et = new ifelseAST(et, then_term, else_term);
  }

  return true;
}

bool parse_type_(state& p, model_info::trsys* ts, model_info::type const*& t,
                 const std::string& name = "") {
  if (tok_is_type(p, ts)) {
    t = p.t;
    get_next_token(p, ts);
  } else {
    AST const* minvalue;
    if (!parse_expr(p, ts, minvalue)) {
      return false;
    }

    if (!chk_expr_const(p, minvalue)) {
      return false;
    }

    if (!chk_expr_integer(p, minvalue)) {
      return false;
    }

    if (p.tok != TOK_DDOT) {
      return set_error(p, "'..' expected");
    }
    get_next_token(p, ts);

    AST const* maxvalue;
    if (!parse_expr(p, ts, maxvalue)) {
      return false;
    }

    if (!chk_expr_const(p, maxvalue)) {
      return false;
    }

    if (!chk_expr_integer(p, maxvalue)) {
      return false;
    }

    t = new model_info::type_subrange(name, ts, new model_info::val(minvalue),
                                      new model_info::val(maxvalue));
  }

  return true;
}

bool parse_cnst_definition(state& p, model_info::trsys* ts) {
  get_next_token(p, ts);

  bool f_ = true;
  while (f_ || p.tok == TOK_COMMA) {
    if (f_) {
      f_ = false;
    } else {
      get_next_token(p, ts);
    }

    if (!chk_id(p)) {
      return false;
    }

    const std::string cname = p.id_str;
    get_next_token(p, ts);

    model_info::type const* ctype = nullptr;
    if (p.tok == TOK_COLON) {
      get_next_token(p, ts);

      if (!parse_type_(p, ts, ctype)) {
        return false;
      }
    }

    if (p.tok != TOK_EQ) {
      return set_error_expected(p, '=');
    }
    get_next_token(p, ts);

    AST const* cval;
    if (!parse_expr(p, ts, cval)) {
      return false;
    }

    if (!chk_expr_const(p, cval)) {
      return false;
    }

    if (!chk_auto_type(p, ctype, cval)) {
      return false;
    }

    ts->add_cnst(
        new model_info::cnst(cname, ts, ctype, new model_info::val(cval)));
  }

  return true;
}

bool parse_var_definition(state& p, model_info::trsys* ts) {
  get_next_token(p, ts);

  std::vector<std::string> vnames;
  bool f_ = true;
  while (f_ || p.tok == TOK_COMMA) {
    if (f_) {
      f_ = false;
    } else {
      get_next_token(p, ts);
    }

    if (!chk_id(p)) {
      return false;
    }

    vnames.push_back(p.id_str);
    get_next_token(p, ts);

    if (p.tok == TOK_COLON) {
      get_next_token(p, ts);
      model_info::type const* vtype;
      if (!parse_type_(p, ts, vtype)) {
        return false;
      }

      for (auto const& vname : vnames) {
        ts->add_var(new model_info::var(vname, ts, vtype));
      }
      vnames.clear();
    }
  }

  for (auto const& vname : vnames) {
    ts->add_var(new model_info::var(vname, ts, model_info::trsys::typeInteger));
  }

  return true;
}

bool parse_clock_definition(state& p, model_info::trsys* ts) {
  get_next_token(p, ts);

  bool f_ = true;
  while (f_ || p.tok == TOK_COMMA) {
    if (f_) {
      f_ = false;
    } else {
      get_next_token(p, ts);
    }

    if (!chk_id(p)) {
      return false;
    }

    ts->add_clock(new model_info::clock(p.id_str, ts));
    get_next_token(p, ts);
  }

  return true;
}

bool parse_type_definition(state& p, model_info::trsys* ts) {
  get_next_token(p, ts);

  bool f_ = true;
  while (f_ || p.tok == TOK_COMMA) {
    if (f_) {
      f_ = false;
    } else {
      get_next_token(p, ts);
    }

    if (!chk_id(p)) {
      return false;
    }

    const std::string tname = p.id_str;
    get_next_token(p, ts);

    if (p.tok != TOK_EQ) {
      return set_error_expected(p, '=');
    }
    get_next_token(p, ts);

    const bool alias_ = tok_is_type(p, ts);
    model_info::type const* t_;

    if (!parse_type_(p, ts, t_, tname)) {
      return false;
    }

    if (alias_) {
      ts->add_type(new model_info::type_alias(tname, ts, t_));
    } else {
      ts->add_type(t_);
    }
  }

  return true;
}

bool parse_fnc_definition(state& p, model_info::trsys* ts) {
  get_next_token(p, ts);

  bool f_ = true;
  while (f_ || p.tok == TOK_COMMA) {
    if (f_) {
      f_ = false;
    } else {
      get_next_token(p, ts);
    }

    if (!chk_id(p)) {
      return false;
    }

    const std::string fname = p.id_str;
    get_next_token(p, ts);

    std::vector<model_info::fnc_arg>* const args =
        new std::vector<model_info::fnc_arg>;
    if (p.tok == TOK_R_OP) {
      get_next_token(p, ts);

      bool f_ = true;
      std::vector<std::string> argnames;
      while (f_ || p.tok == TOK_COMMA) {
        if (f_) {
          f_ = false;
        } else {
          get_next_token(p, ts);
        }

        if (!chk_id(p)) {
          return false;
        }

        auto const& ait = std::find_if(
            args->begin(), args->end(),
            [&](const model_info::fnc_arg& ai) { return ai.name == p.id_str; });

        if (ait != args->end()) {
          sprintf(p.err_buf, "duplicate argument '%s'", p.id_str.c_str());
          return false;
        }

        auto const& ait2 = std::find_if(
            argnames.begin(), argnames.end(),
            [&](const std::string& ani) { return ani == p.id_str; });

        if (ait2 != argnames.end()) {
          sprintf(p.err_buf, "duplicate argument '%s'", p.id_str.c_str());
          return false;
        }

        argnames.push_back(p.id_str);
        get_next_token(p, ts);

        if (p.tok == TOK_COLON) {
          get_next_token(p, ts);

          model_info::type const* argtype;
          if (!parse_type_(p, ts, argtype)) {
            return false;
          }

          for (auto const& argname : argnames) {
            model_info::fnc_arg a;
            a.name = argname;
            a.typ = argtype;
            args->push_back(a);
          }
          argnames.clear();
        }

        for (auto const& argname : argnames) {
          model_info::fnc_arg a;
          a.name = argname;
          a.typ = model_info::trsys::typeReal;
          args->push_back(a);
        }
      }

      if (p.tok != TOK_R_CL) {
        return set_error_expected(p, ')');
      }
      get_next_token(p, ts);
    }

    model_info::type const* ftype = model_info::trsys::typeReal;
    if (p.tok == TOK_COLON) {
      get_next_token(p, ts);
      if (!parse_type_(p, ts, ftype)) {
        return false;
      }
    }

    if (p.tok != TOK_EQ) {
      return set_error_expected(p, '=');
    }
    get_next_token(p, ts, args);

    AST const* body;
    if (!parse_expr(p, ts, body, args)) {
      return false;
    }

    ts->add_fnc(new model_info::fnc(fname, ts, ftype, args, body));
  }

  return true;
}

bool set_chemreaction_guard_updates(
    state& p, AST const* opt_guard,
    const std::vector<model_info::chemreaction_item>& reactants,
    const std::vector<model_info::chemreaction_item>& stoichiometry,
    AST const*& guard, std::vector<model_info::transition_update_item>& updates,
    base const* const sl) {
  std::stringstream last_err;
  if (!set_chemreaction_guard_updates(opt_guard, reactants, stoichiometry,
                                      guard, updates, last_err, sl)) {
    sprintf(p.err_buf, "%s", last_err.str().c_str());
    return false;
  }
  return true;
}

bool set_chemreaction_guard_updates(
    AST const* opt_guard,
    const std::vector<model_info::chemreaction_item>& reactants,
    const std::vector<model_info::chemreaction_item>& stoichiometry,
    AST const*& guard, std::vector<model_info::transition_update_item>& updates,
    std::stringstream& last_err, base const* const sl) {
  guard = nullptr;

  for (auto const& r : reactants) {
    if (sl == nullptr || sl->is_var_control(r.v)) {
      const int c = static_cast<int>(r.c + 0.50) - 1;
      AST const* const g =
          (r.v->get_type()->is_boolean() || !c)
              ? static_cast<AST const*>(new varAST(r.v))
              : static_cast<AST const*>(
                    new binaryAST(TOK_GR, new varAST(r.v), new intnumAST(c)));

      if (guard != nullptr) {
        guard = new binaryAST(TOK_AND, guard, g);
      } else {
        guard = g;
      }
    }
  }

  if (opt_guard != nullptr) {
    if (guard != nullptr) {
      guard = new binaryAST(TOK_AND, guard, opt_guard);
    } else {
      guard = opt_guard;
    }
  }

  if (guard == nullptr) {
    guard = new boolnumAST(true);
  }

  for (auto const& s : stoichiometry) {
    if (sl == nullptr || sl->is_var_control(s.v)) {
      model_info::transition_update_item u;
      u.v = s.v;
      const int c = static_cast<int>(s.c);
      if (s.v->get_type()->is_boolean()) {
        if (c < -1 || c > 1) {
          last_err << "invalid stoichiometry for '" << s.v->get_name().c_str()
                   << "'";
          return false;
        }
        if (c < 0) {
          u.u = new boolnumAST(false);
        } else {
          u.u = new boolnumAST(true);
        }
      } else {
        if (c < 0) {
          u.u = new binaryAST(TOK_MINUS, new varAST(s.v), new intnumAST(-c));
        } else {
          u.u = new binaryAST(TOK_PLUS, new varAST(s.v), new intnumAST(c));
        }
      }
      updates.push_back(u);
    }
  }

  return true;
}

bool parse_chemreaction(state& p, model_info::trsys* ts, base const* const sl) {
  std::string name = "";

  if (p.tok == TOK_S_OP) {
    get_next_token(p, ts);

    if (p.tok == TOK_ID) {
      name = p.id_str;
      get_next_token(p, ts);
    }

    if (p.tok != TOK_S_CL) {
      return set_error_expected(p, ']');
    }
    get_next_token(p, ts);
  }

  AST const* opt_guard = nullptr;

  if (p.tok == TOK_R_OP) {
    if (!parse_expr(p, ts, opt_guard)) {
      return false;
    }

    if (!chk_expr_boolean(p, opt_guard)) {
      return false;
    }
  }

  std::vector<model_info::chemreaction_item> reactants;

  bool f_ = true;
  while (f_ || p.tok == TOK_PLUS) {
    if (!f_) {
      get_next_token(p, ts);
    }

    double c = 1.0;
    if (p.tok == TOK_INTNUM) {
      c = static_cast<double>(p.intnum_val);
      get_next_token(p, ts);

      if (c < 0) {
        return set_error(p, "invalid reactants stoichiometry");
      } else if (c == 0) {
        if (f_) {
          break;
        }
        return set_error(p, "invalid reactants stoichiometry");
      }
    } else if (p.tok == TOK_REALNUM) {
      c = p.realnum_val;
      get_next_token(p, ts);

      if (c < 0) {
        return set_error(p, "invalid reactants stoichiometry");
      }
    }

    if (!chk_var(p, ts)) {
      return false;
    }

    auto const& crit = std::find_if(
        reactants.begin(), reactants.end(),
        [&](const model_info::chemreaction_item& ri) { return ri.v == p.v; });

    if (crit == reactants.end()) {
      const model_info::chemreaction_item cri = {p.v, c};
      reactants.push_back(cri);
    } else {
      crit->c += c;
    }

    get_next_token(p, ts);

    if (f_) {
      f_ = false;
    }
  }

  token_t rtype = p.tok;
  if (rtype != TOK_RARROW && rtype != TOK_DARROW) {
    return set_error(p, "'->' or '<->' expected");
  }
  get_next_token(p, ts);

  std::vector<model_info::chemreaction_item> stoichiometry;
  stoichiometry.assign(reactants.begin(), reactants.end());
  for (auto& ri : stoichiometry) {
    ri.c = -ri.c;
  }

  std::vector<model_info::chemreaction_item> products;
  f_ = true;
  while (f_ || p.tok == TOK_PLUS) {
    if (!f_) {
      get_next_token(p, ts);
    }

    double c = 1.0;
    if (p.tok == TOK_INTNUM) {
      c = static_cast<double>(p.intnum_val);
      get_next_token(p, ts);

      if (c < 0) {
        return set_error(p, "invalid products stoichiometry");
      } else if (c == 0) {
        if (f_) {
          break;
        }
        return set_error(p, "invalid products stoichiometry");
      }
    } else if (p.tok == TOK_REALNUM) {
      c = p.realnum_val;
      get_next_token(p, ts);

      if (c < 0) {
        return set_error(p, "invalid products stoichiometry");
      }
    }

    if (!chk_var(p, ts)) {
      return false;
    }

    auto const& crit = std::find_if(
        products.begin(), products.end(),
        [&](const model_info::chemreaction_item& ri) { return ri.v == p.v; });

    if (crit == products.end()) {
      const model_info::chemreaction_item cri = {p.v, c};
      products.push_back(cri);
    } else {
      crit->c += c;
    }

    auto const& crit2 = std::find_if(
        stoichiometry.begin(), stoichiometry.end(),
        [&](const model_info::chemreaction_item& ri) { return ri.v == p.v; });

    if (crit2 == stoichiometry.end()) {
      const model_info::chemreaction_item cri = {p.v, c};
      stoichiometry.push_back(cri);
    } else {
      crit2->c += c;
    }

    get_next_token(p, ts);

    if (f_) {
      f_ = false;
    }
  }

  stoichiometry.erase(
      std::remove_if(stoichiometry.begin(), stoichiometry.end(),
                     [](const model_info::chemreaction_item& ri) {
        return std::fabs(ri.c) < std::numeric_limits<double>::epsilon();
      }),
      stoichiometry.end());

  if (p.tok != TOK_AT) {
    return set_error_expected(p, '@');
  }
  get_next_token(p, ts);

  AST const* rate = nullptr;
  if (!parse_expr(p, ts, rate)) {
    return false;
  }

  if (!chk_expr_real(p, rate)) {
    return false;
  }

  AST const* guard;
  std::vector<model_info::transition_update_item> updates;

  if (!set_chemreaction_guard_updates(p, opt_guard, reactants, stoichiometry,
                                      guard, updates, sl)) {
    return false;
  }

  model_info::chemreaction const* const cr =
      new model_info::chemreaction(name, ts, guard, rate, updates, reactants,
                                   products, stoichiometry, opt_guard);

  ts->add_transition(cr);

  if (rtype == TOK_DARROW) {
    if (p.tok != TOK_COMMA) {
      return set_error_expected(p, ',');
    }
    get_next_token(p, ts);

    AST const* reverse_rate = nullptr;
    if (!parse_expr(p, ts, reverse_rate)) {
      return false;
    }

    if (!chk_expr_real(p, reverse_rate)) {
      return false;
    }

    for (auto& ri : stoichiometry) {
      ri.c = -ri.c;
    }

    AST const* reverse_guard;
    std::vector<model_info::transition_update_item> reverse_updates;

    if (!set_chemreaction_guard_updates(p, opt_guard, products, stoichiometry,
                                        reverse_guard, reverse_updates, sl)) {
      return false;
    }

    ts->add_transition(new model_info::chemreaction(
        name, ts, reverse_guard, reverse_rate, reverse_updates, products,
        reactants, stoichiometry, opt_guard, true, cr));
  }

  return true;
}

bool parse_guardedcmd(state& p, model_info::trsys* ts) {
  std::string name = "";

  if (p.tok == TOK_S_OP) {
    get_next_token(p, ts);

    if (!chk_id(p)) {
      return false;
    }

    name = p.id_str;
    get_next_token(p, ts);

    if (p.tok != TOK_S_CL) {
      return set_error_expected(p, ']');
    }
    get_next_token(p, ts);
  }

  AST const* guard;
  if (!parse_expr(p, ts, guard)) {
    return false;
  }

  if (!chk_expr_boolean(p, guard)) {
    return false;
  }

  if (p.tok != TOK_COLON) {
    return set_error_expected(p, ':');
  }
  get_next_token(p, ts);

  AST const* u_;
  if (!parse_expr(p, ts, u_)) {
    return false;
  }

  std::vector<model_info::transition_update_item> updates;
  if (!chk_updates(p, u_, updates)) {
    return false;
  }

  if (p.tok != TOK_AT) {
    return set_error_expected(p, '@');
  }
  get_next_token(p, ts);

  AST const* rate;
  if (!parse_expr(p, ts, rate)) {
    return false;
  }

  if (!chk_expr_real(p, rate)) {
    return false;
  }

  ts->add_transition(
      new model_info::transition(name, ts, guard, rate, updates));

  return true;
}

bool parse_init(state& p, model_info::trsys* ts) {
  if (!ts->get_ics().empty()) {
    return set_error(p, "duplicate initial conditions");
  }

  get_next_token(p, ts);

  std::vector<model_info::ic_s> ls;
  double psum = 0.0;
  int isdistr = -1;

  while (p.tok != TOK_END) {
    model_info::ic_s is;

    if (p.tok == TOK_INTNUM || p.tok == TOK_REALNUM) {
      if (isdistr < 0) {
        isdistr = 1;
      } else if (isdistr == 0) {
        return set_error(p, "'end' expected");
      }

      if (p.tok == TOK_INTNUM) {
        is.p = static_cast<double>(p.intnum_val);
      } else if (p.tok == TOK_REALNUM) {
        is.p = p.realnum_val;
      }

      if (is.p < 0.0 || is.p > 1.0) {
        return set_error(p, "invalid state probability");
      }
      get_next_token(p, ts);

      if (p.tok != TOK_COLON) {
        return set_error_expected(p, ':');
      }
      get_next_token(p, ts);
    } else {
      if (isdistr < 0) {
        isdistr = 0;
      } else if (isdistr == 1) {
        return set_error(p, "state probability expected");
      }

      is.p = 1.0;
    }

    bool f_ = true;
    while (f_ || p.tok == TOK_COMMA) {
      if (f_) {
        f_ = false;
      } else {
        get_next_token(p, ts);
      }

      if (!chk_var(p, ts)) {
        return false;
      }

      model_info::ic_s_i isi;
      isi.v = p.v;
      get_next_token(p, ts);

      if (p.tok != TOK_EQ) {
        return set_error_expected(p, '=');
      }
      get_next_token(p, ts);

      AST const* et;
      if (!parse_expr(p, ts, et)) {
        return false;
      }

      if (!chk_expr_const(p, et)) {
        return false;
      }

      isi.value = new model_info::val(et);
      is.li.push_back(isi);
    }

    if (p.tok == TOK_SEMICOLON) {
      get_next_token(p, ts);
    }

    ls.push_back(is);
    psum += is.p;
  }
  get_next_token(p, ts);

  if (std::fabs(1.0 - psum) > std::numeric_limits<double>::epsilon()) {
    return set_error(p, "state probabilities don't sum up to 1.0");
  }

  ts->add_ic(new model_info::ic("", ts, ls));

  return true;
}

bool parse_expr(char const* const str, model_info::trsys const* const context,
                AST const*& et, state& p) {
  static char err_buf[1024];

  p.init(str, err_buf);

  get_next_token(p, context);

  return parse_expr(p, context, et);
}

bool parse_trsys(state& p, std::stack<model_info::trsys*>& tss,
                 base const* const sl) {
  bool success = true;

  while (success) {
    if (p.tok == TOK_EOF) {
      if (tss.top()->get_trsys() == nullptr) {
        tss.pop();
      } else {
        success = set_error(p, "unexpected end of file");
      }
      break;
    } else if (p.tok == TOK_END) {
      if (tss.top()->get_trsys() != nullptr) {
        tss.pop();
        get_next_token(p, tss.top());
      } else {
        success = set_error(p, "unexpected 'end'");
      }
      break;
    } else {
      switch (p.tok) {
        case TOK_CONST: {
          success = parse_cnst_definition(p, tss.top());
          break;
        }

        case TOK_VAR: {
          success = parse_var_definition(p, tss.top());
          break;
        }

        case TOK_CLOCK: {
          success = parse_clock_definition(p, tss.top());
          break;
        }

        case TOK_TYPE: {
          success = parse_type_definition(p, tss.top());
          break;
        }

        case TOK_FUNCTION: {
          success = parse_fnc_definition(p, tss.top());
          break;
        }

        case TOK_CHEMREACTIONS: {
          get_next_token(p, tss.top());
          model_info::trsys* ts_ = new model_info::trsys(
              "", tss.top(), model_info::TS_CHEMREACTIONS);
          tss.top()->add_trsys(ts_);
          tss.push(ts_);
          success = parse_trsys(p, tss, sl);
          break;
        }

        case TOK_GUARDEDCMDS: {
          get_next_token(p, tss.top());
          model_info::trsys* ts_ =
              new model_info::trsys("", tss.top(), model_info::TS_GUARDEDCMDS);
          tss.top()->add_trsys(ts_);
          tss.push(ts_);
          success = parse_trsys(p, tss, sl);
          break;
        }

        case TOK_INIT: {
          success = parse_init(p, tss.top());
          break;
        }

        case TOK_SEMICOLON: {
          get_next_token(p, tss.top());
          break;
        }

        case TOK_UNKNOWN: {
          success = false;
          break;
        }

        default: {
          if (tss.top()->get_type() == model_info::TS_GUARDEDCMDS) {
            success = parse_guardedcmd(p, tss.top());
          } else {
            success = parse_chemreaction(p, tss.top(), sl);
          }
          break;
        }
      }
    }
  }

  return success;
}

bool parse_model(char const* const str, model_info::trsys*& ts, state& p,
                 base const* const sl) {
  static char err_buf[1024];

  p.init(str, err_buf);

  ts = new model_info::trsys("", nullptr);

  std::stack<model_info::trsys*> tss;
  tss.push(ts);

  get_next_token(p, tss.top());

  if (!parse_trsys(p, tss, sl)) {
    ts = nullptr;

    print_line(p);
#ifndef STAR_WEB_INTERFACE
    std::cout << std::endl;
#endif

    return false;
  }

  return true;
}

void write_model_description(std::ostream& os,
                             model_info::trsys const* const ts,
                             char const* tab) {
  if (!ts->get_cnsts().empty()) {
    os << std::endl;
    os << tab << "const" << std::endl;
    const std::size_t clen = ts->get_cnsts().size();
    for (std::size_t i = 0; i < clen; i++) {
      model_info::cnst const* const c = ts->get_cnsts()[i];
      os << tab << "  " << c->get_name() << " : " << c->get_type()->get_name()
         << " = ";
      c->get_value()->get_expr()->write(os, std_writer);
      if (i < clen - 1) {
        os << ",";
      }
      os << std::endl;
    }
  }

  if (!ts->get_vars().empty()) {
    os << std::endl;
    os << tab << "var" << std::endl;
    const std::size_t vlen = ts->get_vars().size();
    for (std::size_t i = 0; i < vlen; i++) {
      model_info::var const* const v = ts->get_vars()[i];
      os << tab << "  " << v->get_name() << " : " << v->get_type()->get_name();
      if (i < vlen - 1) {
        os << ",";
      }
      os << std::endl;
    }
  }

  if (!ts->get_fncs().empty()) {
    os << std::endl;
    os << tab << "function" << std::endl;
    const std::size_t flen = ts->get_fncs().size();
    for (std::size_t i = 0; i < flen; i++) {
      model_info::fnc const* const f = ts->get_fncs()[i];
      os << tab << "  " << f->get_name() << " : " << f->get_rtype()->get_name()
         << " = ";
      f->get_rexpr()->write(os, std_writer);
      if (i < flen - 1) {
        os << ",";
      }
      os << std::endl;
    }
  }

  if (!ts->get_transitions().empty()) {
    os << std::endl;
    for (auto t : ts->get_transitions()) {
      if (t->is_chemreaction()) {
        model_info::chemreaction const* const cr =
            static_cast<model_info::chemreaction const*>(t);
        if (!cr->is_reversible() || cr->get_reverse()) {
          os << tab;

          const std::vector<model_info::chemreaction_item>& reactants =
              cr->get_reactants();
          const std::vector<model_info::chemreaction_item>& products =
              cr->get_products();

          if (!reactants.empty()) {
            for (auto re = reactants.begin(); re != reactants.end(); re++) {
              if (re != reactants.begin()) {
                os << "+ ";
              }
              if (re->c != 1) {
                os << re->c << " ";
              }
              os << re->v->get_name() << " ";
            }
          } else {
            os << "0 ";
          }

          if (cr->get_reverse()) {
            os << "<->";
          } else {
            os << "->";
          }

          if (!products.empty()) {
            for (auto re = products.begin(); re != products.end(); re++) {
              if (re != products.begin()) {
                os << " +";
              }
              if (re->c != 1) {
                os << " " << re->c;
              }
              os << " " << re->v->get_name();
            }
          } else {
            os << " 0";
          }
          os << "  @  ";
          cr->get_rate()->write(os, std_writer);
          if (cr->get_reverse()) {
            os << ", ";
            cr->get_reverse()->get_rate()->write(os, std_writer);
          }
          os << std::endl;
        }
      }
    }
  }

  std::string tab_ = tab;
  tab_ += "  ";

  for (auto& tss : ts->get_trsyss()) {
    if (tss->get_type() == model_info::TS_CHEMREACTIONS) {
      os << tab << "chemical_reactions" << std::endl;
    }

    write_model_description(os, tss, tab_.c_str());

    if (tss->get_type() == model_info::TS_CHEMREACTIONS) {
      os << tab << "end" << std::endl;
    }
  }

  model_info::ic const* const init = ts->get_ics().front();

  os << std::endl;
  os << tab << "init" << std::endl;
  const std::vector<model_info::ic_s>& ls = init->get_states();
  const bool isdistr_ = ls.size() > 1;
  for (auto const& lsi : ls) {
    os << "  ";
    if (isdistr_) {
      os << lsi.p << " : ";
    }
    for (std::size_t i = 0; i < lsi.li.size(); i++) {
      const model_info::ic_s_i& si = lsi.li[i];

      if (i) {
        os << ", ";
        if (i % 5 == 0) {
          os << std::endl;
          if (isdistr_) {
            os << tab << "\t";
          } else {
            os << "  ";
          }
        }
      }
      os << si.v->get_name() << " = ";
      si.value->get_expr()->write(os, std_writer);
    }
    os << std::endl;
  }
  os << tab << "end" << std::endl;
}

int intAST_get_value(AST const* const et) {
  assert(et->is_integral());

  if (et->is_boolnum()) {
    return static_cast<int>(static_cast<boolnumAST const*>(et)->get_value());
  }
  return static_cast<intnumAST const*>(et)->get_value();
}

double numAST_get_value(AST const* const et) {
  assert(et->is_number());

  if (et->is_boolnum()) {
    return static_cast<double>(static_cast<boolnumAST const*>(et)->get_value());
  }
  if (et->is_intnum()) {
    return static_cast<double>(static_cast<intnumAST const*>(et)->get_value());
  }
  return static_cast<realnumAST const*>(et)->get_value();
}

bool same(AST const* const et1, AST const* const et2) {
  if (et1 == et2) {
    return true;
  }

  if (et1->is_cnst() && et2->is_cnst()) {
    return static_cast<cnstAST const*>(et1)->get_cnst() ==
           static_cast<cnstAST const*>(et2)->get_cnst();
  }

  if (et1->is_var() && et2->is_var()) {
    return static_cast<varAST const*>(et1)->get_var() ==
           static_cast<varAST const*>(et2)->get_var();
  }

  return false;
}

AST const* build_power(AST const* const et, AST const* const etp,
                       model_info::transition const* const tr,
                       userfcallAST const* const etf) {
  if (etp->is_number()) {
    if (etp->is_integral()) {
      const int p = intAST_get_value(etp);
      if (!p) {
        return AST::et_1;
      }
      if (p == 1) {
        return et;
      }
    }

    if (et->is_number()) {
      const double p = numAST_get_value(etp);
      return new realnumAST(std::pow(numAST_get_value(et), p));
    }
  }

  if (et->is_stdfcall()) {
    stdfcallAST const* etstdf = static_cast<parser::stdfcallAST const*>(et);
    if (etstdf->get_stdf()->get_id() == parser::STDF_POW) {
      return build_power(etstdf->get_args().front(),
                         build_product(etstdf->get_args()[1], etp, tr, etf), tr,
                         etf);
    }
  }

  const std::vector<parser::AST const*> args = {et, etp};
  return new stdfcallAST(AST::get_stdf(STDF_POW), args);
}

AST const* build_product(AST const* const et, model_info::var const* const v,
                         model_info::transition const* const tr,
                         userfcallAST const* const etf) {
  monom a;
  a.v[v] = 1.0;
  return build_product(et, new polynomAST({a}), tr, etf);
}

AST const* build_product(AST const* const etl, AST const* const etr,
                         model_info::transition const* const tr,
                         userfcallAST const* const etf) {
  if (etl != AST::et_1 && etr != AST::et_1) {
    if (etl == AST::et_0 || etr == AST::et_0) {
      return AST::et_0;
    }
    if (etl->is_number() && etr->is_number()) {
      if (etl->is_realnum() || etr->is_realnum()) {
        return new realnumAST(numAST_get_value(etl) * numAST_get_value(etr));
      }
      if (etl->is_intnum() && etr->is_intnum()) {
        return new intnumAST(static_cast<intnumAST const*>(etl)->get_value() *
                             static_cast<intnumAST const*>(etr)->get_value());
      }
    }

    if (etr->is_polynom()) {
      if (etl->is_polynom()) {
        std::vector<monom> terms;
        for (auto const& a :
             static_cast<parser::polynomAST const*>(etl)->get_terms()) {
          for (auto const& b :
               static_cast<parser::polynomAST const*>(etr)->get_terms()) {
            monom c(a, b);
            bool f = false;
            for (auto& d : terms) {
              if (d.cmp(c)) {
                d.k += c.k;
                f = true;
                break;
              }
            }
            if (!f) {
              terms.push_back(c);
            }
          }
        }
        polynomAST* const etp = new polynomAST();
        for (auto const& a : terms) {
          if (std::fabs(a.k) > std::numeric_limits<double>::epsilon()) {
            etp->get_terms().push_back(a);
          }
        }

        if (etp->get_terms().empty()) {
          return AST::et_0;
        }

        return etp;
      }

      if (etl->is_binary()) {
        binaryAST const* const etb = static_cast<binaryAST const*>(etl);
        if (etb->get_op_tok() == TOK_DIV && etb->get_lterm()->is_polynom()) {
          return build_fraction(build_product(etb->get_lterm(), etr, tr, etf),
                                etb->get_rterm(), tr, etf);
        } else if (etb->get_op_tok() == TOK_MULT &&
                   etb->get_lterm()->is_polynom()) {
          return build_product(build_product(etb->get_lterm(), etr, tr, etf),
                               etb->get_rterm(), tr, etf);
        }
      }
    }

    if (same(etl, etr)) {
      return build_power(etl, new intnumAST(2), tr, etf);
    }
    return new binaryAST(TOK_MULT, etl, etr);
  } else if (etl != AST::et_1) {
    return etl;
  }
  return etr;
}

AST const* build_fraction(AST const* const etn, AST const* const etd,
                          model_info::transition const* const,
                          userfcallAST const* const) {
  if (etn != AST::et_0 && etd != AST::et_1) {
    if (etn->is_polynom() && etd->is_polynom()) {
      polynomAST* const etpd =
          new polynomAST(static_cast<polynomAST const*>(etd));
      if (etpd->is_const()) {
        const double d = etpd->get_terms().front().k;
        polynomAST* const etp =
            new polynomAST(static_cast<polynomAST const*>(etn));
        for (auto& a : etp->get_terms()) {
          a.k /= d;
        }
        return etp;
      }
    }
    return new binaryAST(TOK_DIV, etn, etd);
  }
  return etn;
}

AST const* build_sum(AST const* const etl, AST const* const etr,
                     model_info::transition const* const,
                     userfcallAST const* const) {
  if (etl != AST::et_0 && etr != AST::et_0) {
    if (etl->is_polynom() && etr->is_polynom()) {
      std::vector<monom> terms(
          static_cast<polynomAST const*>(etl)->get_terms());
      for (auto const& b : static_cast<polynomAST const*>(etr)->get_terms()) {
        bool f = false;
        for (auto& a : terms) {
          if (b.cmp(a)) {
            a.k += b.k;
            f = true;
            break;
          }
        }
        if (!f) {
          terms.push_back(b);
        }
      }
      polynomAST* const etp = new polynomAST();
      for (auto const& a : terms) {
        if (std::fabs(a.k) > std::numeric_limits<double>::epsilon()) {
          etp->get_terms().push_back(a);
        }
      }
      return etp;
    }
    if (same(etl, etr)) {
      return new binaryAST(TOK_MULT, new intnumAST(2), etr);
    }
    return new binaryAST(TOK_PLUS, etl, etr);
  } else if (etl != AST::et_0) {
    return etl;
  }
  return etr;
}

AST const* build_difference(AST const* const etl, AST const* const etr,
                            model_info::transition const* const,
                            userfcallAST const* const) {
  if (etl != etr) {
    if (etl != AST::et_0 && etr != AST::et_0) {
      if (etl->is_polynom() && etr->is_polynom()) {
        std::vector<monom> terms(
            static_cast<polynomAST const*>(etl)->get_terms());
        for (auto const& b : static_cast<polynomAST const*>(etr)->get_terms()) {
          bool f = false;
          for (auto& a : terms) {
            if (b.cmp(a)) {
              a.k -= b.k;
              f = true;
              break;
            }
          }
          if (!f) {
            monom b_(b);
            b_.k = -b_.k;
            terms.push_back(b_);
          }
        }
        polynomAST* const etp = new polynomAST();
        for (auto const& a : terms) {
          if (std::fabs(a.k) > std::numeric_limits<double>::epsilon()) {
            etp->get_terms().push_back(a);
          }
        }
        return etp;
      }
      return new binaryAST(TOK_MINUS, etl, etr);
    } else if (etl != AST::et_0) {
      return etl;
    } else if (etr != AST::et_0) {
      return new unaryAST(TOK_MINUS, etr);
    }
  }
  return AST::et_0;
}

AST const* build_negative(AST const* const et,
                          model_info::transition const* const,
                          userfcallAST const* const) {
  if (et != AST::et_0) {
    if (et->is_polynom()) {
      polynomAST* const etp =
          new polynomAST(static_cast<polynomAST const*>(et));
      for (auto& a : etp->get_terms()) {
        a.k = -a.k;
      }
      return etp;
    }
    return new unaryAST(TOK_MINUS, et);
  }
  return AST::et_0;
}

AST const* build_expr(AST const* const et,
                      model_info::transition const* const tr,
                      userfcallAST const* const etf) {
  if (et->is_number()) {
    monom a(numAST_get_value(et));
    return new polynomAST({a});
  } else if (et->is_cnst()) {
    model_info::cnst const* const c =
        static_cast<cnstAST const*>(et)->get_cnst();
    if (c->get_index_p() >= 0) {
      monom a;
      a.c[c] = 1.0;
      return new polynomAST({a});
    }
    monom a(c->get_value()->get());
    return new polynomAST({a});
  } else if (et->is_var()) {
    monom a;
    a.v[static_cast<varAST const*>(et)->get_var()] = 1.0;
    return new polynomAST({a});
  } else if (et->is_unary()) {
    unaryAST const* const etu = static_cast<unaryAST const*>(et);
    AST const* const ett = build_expr(etu->get_term(), tr, etf);
    if (ett != et) {
      if (etu->get_op_tok() == TOK_MINUS) {
        return build_negative(etu->get_term(), tr, etf);
      }
      return new unaryAST(etu->get_op_tok(), ett);
    }
  } else if (et->is_binary()) {
    binaryAST const* const etb = static_cast<binaryAST const*>(et);
    AST const* const etl = build_expr(etb->get_lterm(), tr, etf);
    AST const* const etr = build_expr(etb->get_rterm(), tr, etf);
    if (etl != etb->get_lterm() || etr != etb->get_rterm()) {
      if (etb->get_op_tok() == TOK_MINUS) {
        return build_difference(etl, etr, tr, etf);
      } else if (etb->get_op_tok() == TOK_PLUS) {
        return build_sum(etl, etr, tr, etf);
      } else if (etb->get_op_tok() == TOK_MULT) {
        return build_product(etl, etr, tr, etf);
      } else if (etb->get_op_tok() == TOK_DIV) {
        return build_fraction(etl, etr, tr, etf);
      }
      return new binaryAST(etb->get_op_tok(), etl, etr);
    }
  } else if (et->is_userfarg()) {
    return build_expr(
        etf->get_args()[static_cast<userfargAST const*>(et)->get_index()], tr,
        etf);
  } else if (et->is_stdfcall()) {
    stdfcallAST const* const etu = static_cast<stdfcallAST const*>(et);
    if (etu->get_stdf()->get_id() == STDF_MASS_ACTION) {
      assert(tr != nullptr);

      AST const* etl = build_expr(etu->get_args().front(), tr, etf);
      for (auto const& r :
           static_cast<model_info::chemreaction const*>(tr)->get_reactants()) {
        AST const* etr = nullptr;
        const int c = static_cast<int>(r.c + 0.50);
        if (c == 1) {
          monom a;
          a.v[r.v] = 1.0;
          etr = new polynomAST({a});
        } else if (c == 2) {
          monom a1(0.50);
          a1.v[r.v] = 2.0;
          monom a2(-0.50);
          a2.v[r.v] = 1.0;
          etr = new polynomAST({a1, a2});
        } else {
          assert(false);
        }
        etl = build_product(etl, etr, tr, etf);
      }
      return etl;
    } else {
      assert(false);
    }
  } else if (et->is_userfcall()) {
    userfcallAST const* const etf_ = static_cast<userfcallAST const*>(et);
    return build_expr(etf_->get_f()->get_rexpr(), tr, etf_);
  }
  return et;
}

#if HAVE_LIBSBML

bool sbml_read_expr_stdf(ASTNode const* const anode,
                         const std::size_t stdf_index,
                         model_info::trsys const* const context,

                         const int, AST const*& et, std::stringstream& last_err,
                         std::vector<model_info::fnc_arg> const* const fargs,
                         std::vector<lparam_s> const* const lparams) {
  const stdfnc& sf = sfncs[stdf_index];
  const std::size_t narg = anode->getNumChildren();

  if (sf.narg >= 0 && (int)narg != sf.narg) {
    last_err << "wrong number of arguments";
    return false;
  }

  std::vector<AST const*> eta(narg);
  for (std::size_t i = 0; i < narg; i++) {
    if (!sbml_read_expr(anode->getChild(i), context, 0, eta[i], last_err, fargs,
                        lparams)) {
      return false;
    }
  }

  et = new stdfcallAST(&sf, eta);

  return true;
}

bool sbml_read_expr(ASTNode const* const anode,
                    model_info::trsys const* const context, const int lr,
                    AST const*& et, std::stringstream& last_err,
                    std::vector<model_info::fnc_arg> const* const fargs,
                    std::vector<lparam_s> const* const lparams) {
  if (anode->getType() == AST_PLUS) {
    AST const* etl;
    if (!sbml_read_expr(anode->getLeftChild(), context, lr, etl, last_err,
                        fargs, lparams)) {
      return false;
    }
    AST const* etr;
    if (!sbml_read_expr(anode->getRightChild(), context, lr, etr, last_err,
                        fargs, lparams)) {
      return false;
    }
    et = build_sum(etl, etr);
  } else if (anode->getType() == AST_MINUS) {
    if (anode->getRightChild() != nullptr) {
      if (!lr) {
        if (anode->getRightChild() != nullptr) {
          AST const* etl;
          if (!sbml_read_expr(anode->getLeftChild(), context, 0, etl, last_err,
                              fargs, lparams)) {
            return false;
          }
          AST const* etr;
          if (!sbml_read_expr(anode->getRightChild(), context, 0, etr, last_err,
                              fargs, lparams)) {
            return false;
          }
          et = build_difference(etl, etr);
        } else {
          AST const* etr;
          if (!sbml_read_expr(anode->getLeftChild(), context, 0, etr, last_err,
                              fargs, lparams)) {
            return false;
          }
          et = build_negative(etr);
        }
      } else if (lr < 0) {
        if (!sbml_read_expr(anode->getLeftChild(), context, 0, et, last_err,
                            fargs, lparams)) {
          return false;
        }
      } else {
        if (!sbml_read_expr(anode->getRightChild(), context, 0, et, last_err,
                            fargs, lparams)) {
          return false;
        }
      }
    } else {
      AST const* etr;
      if (!sbml_read_expr(anode->getLeftChild(), context, lr, etr, last_err,
                          fargs, lparams)) {
        return false;
      }
      et = build_negative(etr);
    }
  } else if (anode->getType() == AST_TIMES) {
    AST const* etl;
    if (!sbml_read_expr(anode->getLeftChild(), context, lr, etl, last_err,
                        fargs, lparams)) {
      return false;
    }
    AST const* etr;
    if (!sbml_read_expr(anode->getRightChild(), context, lr, etr, last_err,
                        fargs, lparams)) {
      return false;
    }
    et = build_product(etl, etr);
  } else if (anode->getType() == AST_DIVIDE) {
    AST const* etl;
    if (!sbml_read_expr(anode->getLeftChild(), context, lr, etl, last_err,
                        fargs, lparams)) {
      return false;
    }
    AST const* etr;
    if (!sbml_read_expr(anode->getRightChild(), context, lr, etr, last_err,
                        fargs, lparams)) {
      return false;
    }
    et = build_fraction(etl, etr);
  } else if (anode->getType() == AST_NAME) {
    et = nullptr;
    if (lparams != nullptr) {
      auto const& pi = std::find_if(
          lparams->begin(), lparams->end(),
          [&](lparam_s const& p) { return p.name == anode->getName(); });

      if (pi != lparams->end()) {
        et = new cnstAST(pi->c);
      }
    }

    if (et == nullptr) {
      if (fargs != nullptr) {
        auto const& ai = std::find_if(fargs->begin(), fargs->end(),
                                      [&](const model_info::fnc_arg& a) {
          return a.name == anode->getName();
        });

        if (ai != fargs->end()) {
          et = new userfargAST(fargs, ai - fargs->begin());
        }
      }

      if (et == nullptr) {
        model_info::var const* const v = context->findvar(anode->getName());
        if (v != nullptr) {
          if (v->get_compartment() != nullptr) {
            et = new binaryAST(TOK_DIV, new varAST(v),
                               new cnstAST(v->get_compartment()));
          } else {
            et = new varAST(v);
          }
        } else {
          model_info::cnst const* const c = context->findcnst(anode->getName());
          if (c != nullptr) {
            et = new cnstAST(c);
          } else {
            last_err << "unknown identifier '" << anode->getName() << "'";
            return false;
          }
        }
      }
    }
  } else if (anode->getType() == AST_NAME_TIME) {
    last_err << "time-dependent rate functions are not supported";
    return false;
  } else if (anode->getType() == AST_INTEGER) {
    et = new intnumAST(anode->getInteger());
  } else if (anode->getType() == AST_REAL || anode->getType() == AST_REAL_E) {
    et = new realnumAST(anode->getReal());
  } else if (anode->getType() == AST_CONSTANT_FALSE) {
    et = new boolnumAST(false);
  } else if (anode->getType() == AST_CONSTANT_TRUE) {
    et = new boolnumAST(true);
  } else if (anode->getType() == AST_CONSTANT_E) {
    et = new realnumAST(M_E);
  } else if (anode->getType() == AST_CONSTANT_PI) {
    et = new realnumAST(M_PI);
  } else if (anode->getType() == AST_FUNCTION) {
    model_info::fnc const* const f = context->findfnc(anode->getName());

    if (f == nullptr) {
      last_err << "unknown function '" << anode->getName() << "'";
      return false;
    }

    const std::size_t narg = anode->getNumChildren();
    std::vector<AST const*> eta(narg);
    for (std::size_t i = 0; i < narg; i++) {
      if (!sbml_read_expr(anode->getChild(i), context, 0, eta[i], last_err,
                          fargs, lparams)) {
        return false;
      }
    }

    et = new userfcallAST(f, eta);
  } else if (anode->getType() == AST_FUNCTION_ABS) {
    if (!sbml_read_expr_stdf(anode, STDF_ABS, context, lr, et, last_err, fargs,
                             lparams)) {
      return false;
    }
  } else if (anode->getType() == AST_FUNCTION_ARCCOS) {
    if (!sbml_read_expr_stdf(anode, STDF_ACOS, context, lr, et, last_err, fargs,
                             lparams)) {
      return false;
    }
  } else if (anode->getType() == AST_FUNCTION_ARCCOSH) {
    if (!sbml_read_expr_stdf(anode, STDF_ACOSH, context, lr, et, last_err,
                             fargs, lparams)) {
      return false;
    }
  } else if (anode->getType() == AST_FUNCTION_ARCCOT) {
  } else if (anode->getType() == AST_FUNCTION_ARCCOTH) {
  } else if (anode->getType() == AST_FUNCTION_ARCCSC) {
  } else if (anode->getType() == AST_FUNCTION_ARCCSCH) {
  } else if (anode->getType() == AST_FUNCTION_ARCSEC) {
  } else if (anode->getType() == AST_FUNCTION_ARCSECH) {
  } else if (anode->getType() == AST_FUNCTION_ARCSIN) {
    if (!sbml_read_expr_stdf(anode, STDF_ASIN, context, lr, et, last_err, fargs,
                             lparams)) {
      return false;
    }
  } else if (anode->getType() == AST_FUNCTION_ARCSINH) {
    if (!sbml_read_expr_stdf(anode, STDF_ASINH, context, lr, et, last_err,
                             fargs, lparams)) {
      return false;
    }
  } else if (anode->getType() == AST_FUNCTION_ARCTAN) {
    if (!sbml_read_expr_stdf(anode, STDF_ATAN, context, lr, et, last_err, fargs,
                             lparams)) {
      return false;
    }
  } else if (anode->getType() == AST_FUNCTION_ARCTANH) {
    if (!sbml_read_expr_stdf(anode, STDF_ATANH, context, lr, et, last_err,
                             fargs, lparams)) {
      return false;
    }
  } else if (anode->getType() == AST_FUNCTION_CEILING) {
    if (!sbml_read_expr_stdf(anode, STDF_CEIL, context, lr, et, last_err, fargs,
                             lparams)) {
      return false;
    }
  } else if (anode->getType() == AST_FUNCTION_COS) {
    if (!sbml_read_expr_stdf(anode, STDF_COS, context, lr, et, last_err, fargs,
                             lparams)) {
      return false;
    }
  } else if (anode->getType() == AST_FUNCTION_COSH) {
    if (!sbml_read_expr_stdf(anode, STDF_COSH, context, lr, et, last_err, fargs,
                             lparams)) {
      return false;
    }
  } else if (anode->getType() == AST_FUNCTION_COT) {
  } else if (anode->getType() == AST_FUNCTION_COTH) {
  } else if (anode->getType() == AST_FUNCTION_CSC) {
  } else if (anode->getType() == AST_FUNCTION_CSCH) {
  } else if (anode->getType() == AST_FUNCTION_DELAY) {
    last_err << "delays are not supported";
    return false;
  } else if (anode->getType() == AST_FUNCTION_EXP) {
    if (!sbml_read_expr_stdf(anode, STDF_EXP, context, lr, et, last_err, fargs,
                             lparams)) {
      return false;
    }
  } else if (anode->getType() == AST_FUNCTION_FLOOR) {
    if (!sbml_read_expr_stdf(anode, STDF_FLOOR, context, lr, et, last_err,
                             fargs, lparams)) {
      return false;
    }
  } else if (anode->getType() == AST_FUNCTION_LN) {
    if (!sbml_read_expr_stdf(anode, STDF_LN, context, lr, et, last_err, fargs,
                             lparams)) {
      return false;
    }
  } else if (anode->getType() == AST_FUNCTION_LOG) {
    if (!sbml_read_expr_stdf(anode, STDF_LOG, context, lr, et, last_err, fargs,
                             lparams)) {
      return false;
    }
  } else if (anode->getType() == AST_FUNCTION_PIECEWISE) {
    last_err << "piecewise functions are not supported";
    return false;
  } else if (anode->getType() == AST_FUNCTION_POWER) {
    if (!sbml_read_expr_stdf(anode, STDF_POW, context, lr, et, last_err, fargs,
                             lparams)) {
      return false;
    }
  } else if (anode->getType() == AST_FUNCTION_ROOT) {
    if (anode->isSqrt()) {
      if (!sbml_read_expr_stdf(anode, STDF_SQRT, context, lr, et, last_err,
                               fargs, lparams)) {
        return false;
      }
    } else {
      std::vector<AST const*> eta(2);
      if (!sbml_read_expr(anode, context, lr, eta[0], last_err, fargs,
                          lparams)) {
        return false;
      }
      if (!sbml_read_expr(anode, context, lr, eta[1], last_err, fargs,
                          lparams)) {
        return false;
      }

      eta[1] = build_fraction(AST::et_1, eta[1]);
      et = new stdfcallAST(&sfncs[STDF_POW], eta);
    }
  } else if (anode->getType() == AST_FUNCTION_SEC) {
  } else if (anode->getType() == AST_FUNCTION_SECH) {
  } else if (anode->getType() == AST_FUNCTION_SIN) {
    if (!sbml_read_expr_stdf(anode, STDF_SIN, context, lr, et, last_err, fargs,
                             lparams)) {
      return false;
    }
  } else if (anode->getType() == AST_FUNCTION_SINH) {
    if (!sbml_read_expr_stdf(anode, STDF_SINH, context, lr, et, last_err, fargs,
                             lparams)) {
      return false;
    }
  } else if (anode->getType() == AST_FUNCTION_TAN) {
    if (!sbml_read_expr_stdf(anode, STDF_TAN, context, lr, et, last_err, fargs,
                             lparams)) {
      return false;
    }
  } else if (anode->getType() == AST_FUNCTION_TANH) {
    if (!sbml_read_expr_stdf(anode, STDF_TANH, context, lr, et, last_err, fargs,
                             lparams)) {
      return false;
    }
  } else {
    last_err << "unhandled ASTNode type " << anode->getType() << std::endl;
    return false;
  }

  return true;
}

#endif
}

model_info::val_t compute_const_expr(
    parser::AST const* const et, const bool comp_init = false,
    const std::vector<model_info::val_t>& fargs =
        std::vector<model_info::val_t>()) {
  if (et->is_boolnum()) {
    return static_cast<model_info::val_t>(
        static_cast<parser::boolnumAST const*>(et)->get_value());
  } else if (et->is_intnum()) {
    return static_cast<model_info::val_t>(
        static_cast<parser::intnumAST const*>(et)->get_value());
  } else if (et->is_realnum()) {
    return static_cast<model_info::val_t>(
        static_cast<parser::realnumAST const*>(et)->get_value());
  } else if (et->is_cnst()) {
    return static_cast<parser::cnstAST const*>(et)
        ->get_cnst()
        ->get_value()
        ->get();
  } else if (et->is_var()) {
    if (comp_init) {
      parser::AST const* const eti =
          static_cast<parser::varAST const*>(et)->get_var()->get_init();
      if (eti != nullptr) {
        return compute_const_expr(eti, comp_init, fargs);
      }
      return 0.0;
    }
  } else if (et->is_unary()) {
    parser::unaryAST const* const etu =
        static_cast<parser::unaryAST const*>(et);
    const model_info::val_t v =
        compute_const_expr(etu->get_term(), comp_init, fargs);
    switch (etu->get_op_tok()) {
      case parser::TOK_NOT:
        return static_cast<model_info::val_t>(!static_cast<bool>(v));

      case parser::TOK_MINUS:
        return -v;

      case parser::TOK_PLUS:
        return v;

      default:
        break;
    }
  } else if (et->is_binary()) {
    parser::binaryAST const* const etb =
        static_cast<parser::binaryAST const*>(et);
    const model_info::val_t v1 =
        compute_const_expr(etb->get_lterm(), comp_init, fargs);
    const model_info::val_t v2 =
        compute_const_expr(etb->get_rterm(), comp_init, fargs);
    switch (etb->get_op_tok()) {
      case parser::TOK_MULT:
        return v1 * v2;

      case parser::TOK_DIV:
        return v1 / v2;

      case parser::TOK_PLUS:
        return v1 + v2;

      case parser::TOK_MINUS:
        return v1 - v2;

      case parser::TOK_LESS:
        return v1 < v2;

      case parser::TOK_LEQ:
        return v1 <= v2;

      case parser::TOK_GR:
        return v1 > v2;

      case parser::TOK_GREQ:
        return v1 >= v2;

      case parser::TOK_EQ:
        return v1 == v2;

      case parser::TOK_INEQ:
        return v1 != v2;

      case parser::TOK_AND:
        return static_cast<model_info::val_t>(static_cast<bool>(v1) &&
                                              static_cast<bool>(v2));

      case parser::TOK_OR:
        return static_cast<model_info::val_t>(static_cast<bool>(v1) ||
                                              static_cast<bool>(v2));

      default:
        break;
    }
  } else if (et->is_userfcall()) {
    parser::userfcallAST const* const etf =
        static_cast<parser::userfcallAST const*>(et);
    const std::size_t narg = etf->get_args().size();
    std::vector<model_info::val_t> x(narg);
    for (std::size_t i = 0; i < narg; i++) {
      x[i] = compute_const_expr(etf->get_args()[i], comp_init, fargs);
    }
    return compute_const_expr(etf->get_f()->get_rexpr(), comp_init, x);
  } else if (et->is_userfarg()) {
    return fargs[static_cast<parser::userfargAST const*>(et)->get_index()];
  } else if (et->is_stdfcall()) {
    parser::stdfcallAST const* const etf =
        static_cast<parser::stdfcallAST const*>(et);
    const std::size_t narg = etf->get_args().size();
    std::vector<model_info::val_t> x(narg);
    for (std::size_t i = 0; i < narg; i++) {
      x[i] = compute_const_expr(etf->get_args()[i], comp_init, fargs);
    }
    switch (etf->get_stdf()->get_id()) {
      case parser::STDF_ABS:
        return std::fabs(x[0]);

      case parser::STDF_MIN: {
        double xmin = x[0];
        for (std::size_t i = 1; i < narg; i++) {
          if (x[i] < xmin) {
            xmin = x[i];
          }
        }
        return xmin;
      }

      case parser::STDF_MAX: {
        double xmax = x[0];
        for (std::size_t i = 1; i < narg; i++) {
          if (x[i] > xmax) {
            xmax = x[i];
          }
        }
        return xmax;
      }

      case parser::STDF_FLOOR:
        return std::floor(x[0]);

      case parser::STDF_CEIL:
        return std::ceil(x[0]);

      case parser::STDF_MOD:
        return std::fmod(x[0], x[1]);

      case parser::STDF_EXP:
        return std::exp(x[0]);

      case parser::STDF_SQRT:
        return std::sqrt(x[0]);

      case parser::STDF_POW:
        return std::pow(x[0], x[1]);

      case parser::STDF_LN:
        return std::log(x[0]);

      case parser::STDF_LOG10:
        return std::log10(x[0]);

      case parser::STDF_LOG:
        return std::log(x[0]) / std::log(x[1]);

      case parser::STDF_SIN:
        return std::sin(x[0]);

      case parser::STDF_COS:
        return std::cos(x[0]);

      case parser::STDF_TAN:
        return std::tan(x[0]);

      case parser::STDF_ASIN:
        return std::asin(x[0]);

      case parser::STDF_ACOS:
        return std::acos(x[0]);

      case parser::STDF_ATAN:
        return std::atan(x[0]);

      case parser::STDF_ATAN2:
        return std::atan2(x[0], x[1]);

      case parser::STDF_SINH:
        return std::sinh(x[0]);

      case parser::STDF_COSH:
        return std::cosh(x[0]);

      case parser::STDF_TANH:
        return std::tanh(x[0]);

      case parser::STDF_ASINH:
        return std::asinh(x[0]);

      case parser::STDF_ACOSH:
        return std::acosh(x[0]);

      case parser::STDF_ATANH:
        return std::atanh(x[0]);

      default:
        break;
    }
  } else if (et->is_ifelse()) {
    parser::ifelseAST const* const eti =
        static_cast<parser::ifelseAST const*>(et);
    const bool b = static_cast<bool>(
        compute_const_expr(eti->get_cond(), comp_init, fargs));
    if (b) {
      return compute_const_expr(eti->get_then_term(), comp_init, fargs);
    } else {
      return compute_const_expr(eti->get_else_term(), comp_init, fargs);
    }
  }

  assert(false);
  return 0.0;
}

void compute_value(model_info::val* const value, const bool comp_init) {
  value->set(compute_const_expr(value->get_expr(), comp_init));
}
}
