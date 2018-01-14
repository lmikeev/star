/*
 *  parser.hpp
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

#ifndef SOLVER_LOADER_PARSER_PARSER_HPP_
#define SOLVER_LOADER_PARSER_PARSER_HPP_

#include "ast.hpp"

namespace solver_loader {

class base;

namespace model_info {

class trsys;
class cnst;
class var;
class fnc;
class type;

struct chemreaction_item;
struct transition_update_item;
}

namespace parser {

struct state {
  char const* s;
  int line_number;
  int col_number;

  char last;
  token_t tok;

  std::string line;

  model_info::cnst const* c;
  model_info::var const* v;
  model_info::fnc const* f;
  model_info::type const* t;

  stdfnc const* stdf;

  std::string id_str;

  bool boolnum_val;
  long int intnum_val;
  double realnum_val;

  int farg_i;

  char* err_buf;

  void init(char const* const str, char* err_buf_) {
    s = str;
    line_number = 1;
    last = ' ';
    err_buf = err_buf_;

    boolnum_val = false;
    intnum_val = 0;
    realnum_val = 0.0;
  }
};

bool parse_expr(char const* const str, model_info::trsys const* const context,
                AST const*& et, state& p);
bool parse_model(char const* const str, model_info::trsys*& ts, state& p,
                 base const* const sl = nullptr);

char* get_error(state& p);
int get_line_number(state& p);
int get_col_number(state& p);

bool set_chemreaction_guard_updates(
    AST const* opt_guard,
    const std::vector<model_info::chemreaction_item>& reactants,
    const std::vector<model_info::chemreaction_item>& stoichiometry,
    AST const*& guard, std::vector<model_info::transition_update_item>& updates,
    std::stringstream& last_err, base const* const sl = nullptr);

void write_model_description(std::ostream& os,
                             model_info::trsys const* const ts,
                             char const* tab = "");
}
}

#endif
