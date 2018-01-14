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

#ifndef SOLVER_LOADER_TASK_INFO_BASE_HPP_
#define SOLVER_LOADER_TASK_INFO_BASE_HPP_

#include <sstream>
#include "../parser/ast.hpp"

#define PROPS_DEF_PROP(type, name) \
  type name;                       \
  bool name##_set;

#define PROPS_SET_PROP(type, name)       \
  void set_##name(const type& name##_) { \
    name = name##_;                      \
    name##_set = true;                   \
  }

#define PROPS_GET_PROP(type, name)       \
  bool get_##name(type& name##_) const { \
    if (name##_set) {                    \
      name##_ = name;                    \
      return true;                       \
    }                                    \
    return false;                        \
  }

#define PROPS_SET_GET_PROP(type, name) \
  PROPS_SET_PROP(type, name)           \
  PROPS_GET_PROP(type, name)

namespace solver_loader {

class base;

namespace task_info {

enum type {
  TASK_BASE,

  TASK_DUMP_STATS,
  TASK_DUMP_MOMENTS,
  TASK_DUMP_DISTR,
  TASK_DUMP_SAMPLES,

  TASK_PLOT,
  TASK_PLOT_2D,
  TASK_PLOT_DISTR,
  TASK_PLOT_DISTR_2D,

  TASK_PLOT_OBJF,
  TASK_PLOT_OBJF_2D,
};

class base {
 protected:
  std::string name;

 public:
  base(const std::string& name) : name(name) {}
  base() {}
  virtual ~base() {}

  base(const base& ti) : base(ti.get_name()) {}

  char const* get_name() const { return name.c_str(); }

  virtual type get_type() const { return TASK_BASE; }

  virtual bool check(solver_loader::base const* const) { return true; }

  virtual bool post_check(solver_loader::base const* const) { return true; }

  virtual base* clone() const { return new base(*this); }

  bool check_var_names(solver_loader::base const* const sl,
                       const std::vector<std::string>& varnames,
                       std::vector<model_info::var const*>& vars);
};

class find_by_type : std::unary_function<type*, bool> {
 private:
  const type t;

 public:
  find_by_type(const type t) : t(t) {}

  bool operator()(task_info::base const* const ti) const {
    return t == ti->get_type();
  }
};

class cmp {
 private:
  const std::string cmp_path;

 public:
  cmp(const std::string& cmp_path = "") : cmp_path(cmp_path) {}

  bool do_cmp() const { return cmp_path != ""; }

  const std::string& get_cmp_path() const { return cmp_path; }
};

class cond {
 private:
  std::string cnd_str;

 protected:
  parser::AST const* cnd_et;
  int cnd_index;

 public:
  cond(const std::string& cnd_str = "")
      : cnd_str(cnd_str), cnd_et(nullptr), cnd_index(-1) {}

  ~cond() {
    if (cnd_et != nullptr) {
      delete cnd_et;
    }
  }

  void remove_cond() { cnd_str = ""; }

  const std::string& get_cond_str() const { return cnd_str; }

  parser::AST const* get_cond_et() const { return cnd_et; }

  int get_cond_index() const { return cnd_index; }

  bool check_cond(solver_loader::base const* const sl, std::string& name);
};

class plot_props {
 private:
  PROPS_DEF_PROP(std::string, title)
  PROPS_DEF_PROP(std::string, xlabel)
  PROPS_DEF_PROP(std::string, ylabel)

 public:
  plot_props() : title_set(false), xlabel_set(false), ylabel_set(false) {}

  PROPS_SET_GET_PROP(std::string, title)
  PROPS_SET_GET_PROP(std::string, xlabel)
  PROPS_SET_GET_PROP(std::string, ylabel)
};

class plot_line_props {
 private:
  PROPS_DEF_PROP(std::string, title)
  PROPS_DEF_PROP(unsigned int, color)
  PROPS_DEF_PROP(unsigned int, type)
  PROPS_DEF_PROP(unsigned int, width)

 public:
  plot_line_props()
      : title_set(false), color_set(false), type_set(false), width_set(false) {}

  PROPS_SET_GET_PROP(std::string, title)
  PROPS_SET_GET_PROP(unsigned int, color)
  PROPS_SET_GET_PROP(unsigned int, type)
  PROPS_SET_GET_PROP(unsigned int, width)
};

class plot_surface_props {
 private:
  PROPS_DEF_PROP(std::string, palette)

 public:
  plot_surface_props() : palette_set(false) {}

  PROPS_SET_GET_PROP(std::string, palette)
};

struct plot_expr {
  std::string expr;
  std::string dparam;
  std::string dparam2;
  plot_line_props line_props;

  parser::AST const* et;
  model_info::cnst const* dparami;
  model_info::cnst const* dparamj;

  int expr_index;

  bool chm;

  void free() {}

  plot_expr()
      : et(nullptr),
        dparami(nullptr),
        dparamj(nullptr),
        expr_index(-1),
        chm(false) {}

  ~plot_expr() { free(); }
};

class plot_base : public base, public cond {
 private:
  const plot_props props;
  const bool isdynamic_;

 public:
  plot_base(const plot_props& props, const std::string& cnd = "",
            const bool isdynamic_ = true)
      : cond(cnd), props(props), isdynamic_(isdynamic_) {}

  const plot_props& get_props() const { return props; }

  bool is_dynamic() const { return isdynamic_; }
};

class plot_2d_base : public plot_base {
 protected:
  const std::string expr1, expr2;
  const std::string dparam, dparam2;

  parser::AST const* et1;
  parser::AST const* et2;
  model_info::cnst const* dparami;
  model_info::cnst const* dparamj;

  int expr1_index, expr2_index;

 public:
  plot_2d_base(const plot_props& props, const std::string& expr1,
               const std::string& expr2, const std::string& cnd = "",
               const std::string& dparam = "", const std::string& dparam2 = "",
               const bool isdynamic_ = true)
      : plot_base(props, cnd, isdynamic_),
        expr1(expr1),
        expr2(expr2),
        dparam(dparam),
        dparam2(dparam2),
        et1(nullptr),
        et2(nullptr),
        dparami(nullptr),
        dparamj(nullptr),
        expr1_index(-1),
        expr2_index(-1) {}

  ~plot_2d_base() {}

  const std::string& get_expr1() const { return expr1; }

  const std::string& get_expr2() const { return expr2; }

  parser::AST const* get_et1() const { return et1; }

  parser::AST const* get_et2() const { return et2; }

  const std::string& get_dparam() const { return dparam; }

  const std::string& get_dparam2() const { return dparam2; }

  model_info::cnst const* get_dparami() const { return dparami; }

  model_info::cnst const* get_dparamj() const { return dparamj; }

  int get_expr1_index() const { return expr1_index; }

  int get_expr2_index() const { return expr2_index; }

  bool check_exprs(solver_loader::base const* const sl);
};
}
}

#endif
