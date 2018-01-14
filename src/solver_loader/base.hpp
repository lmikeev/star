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

#ifndef SOLVER_LOADER_BASE_HPP_
#define SOLVER_LOADER_BASE_HPP_

#define PAR2_LEN(N) ((N) * (N + 1) / 2)

#include <set>
#include <typeinfo>
#include <typeindex>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/function.hpp>
#include "model_info/trsys.hpp"
#include "task_info/all.hpp"
#include "parser/ast.hpp"
#include "options.hpp"
#include "../../config.hpp"

#ifdef WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

struct star_fnsi;
class experiment;

namespace solver {
class base;
class model;
}

namespace solver_loader {

namespace writer {

class base;
}

enum type {
  SL_BASE,

  SL_TRANSIENT,

  SL_SCAN,
  SL_OPTIMIZE,
  SL_FIT
};

class base : public model_info::trsys, public options {
 protected:
  experiment const* const E;
  const bool issubsloader;

  solver::base* S;
  solver::model* M;

  char dir[32];
  unsigned int id;

  std::ofstream fcout, fcerr;

  mutable boost::filesystem::path tmp_dir;
  mutable std::stringstream last_err;

  mutable std::vector<parser::AST const*> exprs;

  bool create_dir() const;
  bool load_solver();

  virtual bool preload_check_options() { return true; }
  virtual bool afterload_check_options() { return true; }

  bool check_model_options();

  virtual bool check_tasks() { return true; }
  virtual bool post_check_tasks() { return true; }

#ifdef STAR_CODEGEN

#ifdef WIN32
  HINSTANCE solverlibHandle;
#else
  void* solverlibHandle;
#endif

  mutable char model_class[16];

  bool open_log(std::ofstream& os) const;
  void close_log(std::ofstream& os) const;

  void write_set_vars_i(std::ostream& os,
                        const std::set<model_info::var const*>& ev,
                        char const* const tab = "    ") const;
  void write_set_vars_d(std::ostream& os,
                        const std::set<model_info::var const*>& ev,
                        char const* const tab = "    ") const;

  void write_set_params_(std::ostream& os,
                         const std::set<model_info::cnst const*>& ep,
                         char const* const tab = "    ") const;
  void write_set_params(std::ostream& os,
                        const std::set<model_info::cnst const*>& ep) const;

  void write_drate_h_dc(std::ostream& os, parser::AST const* const rate_h,
                        const std::size_t ci,
                        const std::size_t rates_h_len) const;

  void write_drate_h_dx(std::ostream& os, parser::AST const* const rate_h,
                        const std::size_t ci,
                        const std::size_t rates_h_len) const;

  void write_drate_h_dx(std::ostream& os, std::size_t* const I,
                        std::size_t* const I_, const std::size_t d,
                        const std::size_t ci0,
                        parser::AST const* const drate_h_dx,
                        const std::size_t rates_h_len) const;

  void write_rate_h(std::ostream& os, parser::AST const* const rate_h,
                    const std::size_t ci0, const std::size_t rates_h_len) const;

  void write_rate_h(std::ostream& os, parser::AST const* const rate_h,
                    std::size_t* const I, std::size_t* const I_,
                    const std::size_t d, const std::size_t ci0,
                    const std::size_t rates_h_len) const;

  bool write_usermodel() const;
  bool write_solverext() const;
  bool write_Makefile() const;

  virtual char const* solver_params_str() const { return "sl"; }

  void write_solver_def(std::ostream& os, char const* const solver_class,
                        const bool subsloader = false) const {
    os << "typedef " << solver_class << " "
       << (subsloader ? "Subsolver" : "Solver") << ";" << std::endl;
  }

  virtual bool write_solver_def(std::ostream&, const bool = false) const {
    return false;
  }

  virtual bool write_Makefile_defs(std::ostream&) const { return true; }

  virtual bool write_solver() {
    if (!write_usermodel()) {
      return false;
    }
    if (!write_solverext()) {
      return false;
    }
    if (!write_Makefile()) {
      return false;
    }
    return true;
  }

  bool load_solver_lib();
#else

#endif

  bool expr_contains_obj(
      parser::AST const* const et,
      std::function<bool(parser::AST const* const e)> found) const {
    if (et->is_unary()) {
      parser::unaryAST const* const etu =
          static_cast<parser::unaryAST const*>(et);
      return expr_contains_obj(etu->get_term(), found);
    } else if (et->is_binary()) {
      parser::binaryAST const* const etb =
          static_cast<parser::binaryAST const*>(et);
      return expr_contains_obj(etb->get_lterm(), found) ||
             expr_contains_obj(etb->get_rterm(), found);
    } else if (et->is_stdfcall()) {
      parser::stdfcallAST const* const etf =
          static_cast<parser::stdfcallAST const*>(et);
      for (auto a : etf->get_args()) {
        if (expr_contains_obj(a, found)) {
          return true;
        }
      }
    } else if (et->is_userfcall()) {
      parser::userfcallAST const* const etf =
          static_cast<parser::userfcallAST const*>(et);
      for (auto a : etf->get_args()) {
        if (expr_contains_obj(a, found)) {
          return true;
        }
      }
    } else {
    }

    return found(et);
  }

  bool expr_contains_clock(parser::AST const* const et) const {
    return expr_contains_obj(
        et, [](parser::AST const* const e) { return e->is_clock(); });
  }

  bool expr_contains_cnst(parser::AST const* const et) const {
    return expr_contains_obj(et, [](parser::AST const* const e) {
      if (e->is_polynom()) {
        parser::polynomAST const* const etp =
            static_cast<parser::polynomAST const*>(e);
        for (auto const& a : etp->get_terms()) {
          if (!a.c.empty()) {
            return true;
          }
        }
        return false;
      }
      return e->is_cnst();
    });
  }

  bool expr_contains_cnst(parser::AST const* const et,
                          model_info::cnst const* const c) const {
    return expr_contains_obj(et, [=](parser::AST const* const e) {
      if (e->is_polynom()) {
        parser::polynomAST const* const etp =
            static_cast<parser::polynomAST const*>(e);
        for (auto const& a : etp->get_terms()) {
          if (a.c.count(c)) {
            return true;
          }
        }
        return false;
      }
      return e->is_cnst() &&
             static_cast<parser::cnstAST const* const>(e)->get_cnst() == c;
    });
  }

  bool expr_contains_var(parser::AST const* const et) const {
    return expr_contains_obj(et, [this](parser::AST const* const e) {
      if (e->is_polynom()) {
        parser::polynomAST const* const etp =
            static_cast<parser::polynomAST const*>(e);
        for (auto const& a : etp->get_terms()) {
          if (!a.v.empty()) {
            return true;
          }
        }
        return false;
      }
      return e->is_var();
    });
  }

  bool expr_contains_pvar(parser::AST const* const et) const {
    return expr_contains_obj(et, [this](parser::AST const* const e) {
      if (e->is_polynom()) {
        parser::polynomAST const* const etp =
            static_cast<parser::polynomAST const*>(e);
        for (auto const& a : etp->get_terms()) {
          for (auto const& vi : a.v) {
            if (is_var_population(vi.first)) {
              return true;
            }
          }
        }
        return false;
      }
      return e->is_var() &&
             is_var_population(
                 static_cast<parser::varAST const* const>(e)->get_var());
    });
  }

  bool expr_contains_pvar(parser::AST const* const et,
                          model_info::var const* const v) const {
    return expr_contains_obj(et, [=](parser::AST const* const e) {
      if (e->is_polynom()) {
        parser::polynomAST const* const etp =
            static_cast<parser::polynomAST const*>(e);
        for (auto const& a : etp->get_terms()) {
          if (a.v.count(v)) {
            return true;
          }
        }
        return false;
      }
      return e->is_var() &&
             static_cast<parser::varAST const* const>(e)->get_var() == v;
    });
  }

  bool expr_contains_not_pvar(parser::AST const* const et) const {
    return expr_contains_obj(et, [this](parser::AST const* const e) {
      if (e->is_polynom()) {
        parser::polynomAST const* const etp =
            static_cast<parser::polynomAST const*>(e);
        for (auto const& a : etp->get_terms()) {
          if (!a.v.empty()) {
            for (auto const& vi : a.v) {
              if (!is_var_population(vi.first)) {
                return true;
              }
            }
          } else {
            return true;
          }
        }
        return false;
      }

      return !e->is_var() ||
             !is_var_population(
                 static_cast<parser::varAST const* const>(e)->get_var());
    });
  }

  template <class T>
  bool expr_eq_obj(parser::AST const* const, T const* const) const {
    return false;
  }

  template <class T>
  bool cnst_eq_obj(model_info::cnst const* const, T const* const) const {
    return false;
  }

  template <class T>
  bool var_eq_obj(model_info::var const* const, T const* const) const {
    return false;
  }

  template <class T>
  parser::AST const* build_expr_der(parser::AST const* const et,
                                    T const* const o) const {
    if (expr_eq_obj<T>(et, o)) {
      return parser::AST::et_1;
    } else if (et->is_unary()) {
      parser::unaryAST const* etu = static_cast<parser::unaryAST const*>(et);
      parser::AST const* const d = build_expr_der<T>(etu->get_term(), o);
      switch (etu->get_op_tok()) {
        case parser::TOK_MINUS:
          return build_negative(d);

        case parser::TOK_PLUS:
        default:
          return et;
      }
    } else if (et->is_binary()) {
      parser::binaryAST const* etb = static_cast<parser::binaryAST const*>(et);
      parser::AST const* const dl = build_expr_der<T>(etb->get_lterm(), o);
      parser::AST const* const dr = build_expr_der<T>(etb->get_rterm(), o);
      switch (etb->get_op_tok()) {
        case parser::TOK_MULT:
          return build_sum(build_product(dl, etb->get_rterm()),
                           build_product(dr, etb->get_lterm()));

        case parser::TOK_DIV: {
          parser::AST const* const d1 = build_product(dl, etb->get_rterm());
          parser::AST const* const d2 = build_product(dr, etb->get_lterm());
          if (d1 != parser::AST::et_0 && d2 != parser::AST::et_0) {
            return build_fraction(
                build_difference(d1, d2),
                build_product(etb->get_rterm(), etb->get_rterm()));
          } else if (d1 != parser::AST::et_0) {
            return build_fraction(dl, etb->get_rterm());
          } else if (d2 != parser::AST::et_0) {
            return build_negative(build_fraction(
                d2, build_product(etb->get_rterm(), etb->get_rterm())));
          }
          break;
        }

        case parser::TOK_PLUS:
          return build_sum(dl, dr);

        case parser::TOK_MINUS:
          return build_difference(dl, dr);

        default:
          assert(false);
          break;
      }
    } else if (et->is_polynom()) {
      std::vector<parser::monom> terms;
      for (auto const& a :
           static_cast<parser::polynomAST const*>(et)->get_terms()) {
        parser::monom b(a.k);
        bool f = false;
        for (auto const& ci : a.c) {
          if (strcmp(typeid(ci.first).name(),
                     typeid(model_info::cnst const*).name()) == 0 &&
              ci.first == dynamic_cast<model_info::cnst const*>(o)) {
            b.k *= ci.second;
            if (ci.second - 1.0 > std::numeric_limits<double>::epsilon()) {
              b.c[ci.first] = ci.second - 1.0;
            }
            f = true;
          } else {
            b.c[ci.first] = ci.second;
          }
        }
        for (auto const& vi : a.v) {
          if (var_eq_obj<T>(vi.first, o)) {
            b.k *= vi.second;
            if (vi.second - 1.0 > std::numeric_limits<double>::epsilon()) {
              b.v[vi.first] = vi.second - 1.0;
            }
            f = true;
          } else {
            b.v[vi.first] = vi.second;
          }
        }
        if (f) {
          terms.push_back(b);
        }
      }
      if (!terms.empty()) {
        return new parser::polynomAST(terms);
      }
    } else if (et->is_stdfcall()) {
      parser::stdfcallAST const* etf =
          static_cast<parser::stdfcallAST const*>(et);
      switch (etf->get_stdf()->get_id()) {
        case parser::STDF_EXP:
          return build_product(et,
                               build_expr_der<T>(etf->get_args().front(), o));

        case parser::STDF_POW: {
          return build_product(
              etf->get_args()[1],
              build_power(
                  etf->get_args().front(),
                  build_difference(etf->get_args()[1], parser::AST::et_1)));
        }

        default:
          assert(false);
          break;
      }
    } else {
    }
    return parser::AST::et_0;
  }

  void flatten(trsys const* const ts);

  bool load_model_src(char const* const src);
  bool load_sbml_model(const std::string& fname);
  bool load_custom_model(const int id);
  bool load_fnsi_model(star_fnsi const* const fnsi,
                       const std::string& model_name);

  void ic_dfs(const std::size_t u,
              const std::vector<std::vector<std::size_t>>& g,
              std::vector<bool>& dfs_num, std::vector<std::size_t>& topo) const;

 public:
  base(experiment* const e, const bool issubsloader = false);
  virtual ~base();

  bool load_model();

  bool compile_solver();
  bool run_solver();

  bool update_progress(const double progress) const;
  bool update_progress(char const* const progress) const;

  bool is_var_control(model_info::var const* const v) const;

  bool is_var_population(model_info::var const* const v) const {
    return !is_var_control(v);
  }

  experiment const* get_experiment() const { return E; }

  solver::base* get_solver() const { return S; }

  solver::model* get_model() const { return M; }

  void set_model(solver::model* const m) { M = m; }

  char const* get_dir() const { return dir; }

  boost::filesystem::path get_path() const;

  unsigned int get_id() const { return id; }

  std::string get_last_error() const { return last_err.str(); }

  std::stringstream& last_error() const { return last_err; }

  virtual enum type get_type() const { return SL_BASE; }

  virtual char* sprint_solver_description(char* const buf) const { return buf; }

  virtual char* sprint_description(char* const buf) const;

  std::vector<model_info::var const*> c_vars, p_vars;
  std::vector<model_info::transition const*> trs_stoch, trs_det;
  mutable std::vector<model_info::cnst const*> sa_params;
  mutable std::vector<model_info::var const*> sa_ivars;
  mutable std::vector<model_info::var const*> sa_evars;

  std::set<int> o_vars;
  std::vector<double> o_err;

  const std::vector<model_info::transition const*>& get_trs_stoch() const {
    return trs_stoch;
  }

  const std::vector<model_info::transition const*>& get_trs_det() const {
    return trs_det;
  }

  bool is_var_observable(model_info::var const* const v) const {
    return o_vars.empty() || o_vars.count(v->get_index());
  }

  int get_var_minvalue(model_info::var const* const v) const {
    if (v->get_type()->is_boolean() || v->get_type()->is_species()) {
      return 0;
    }
    if (v->get_type()->is_subrange()) {
      return static_cast<int>(
          static_cast<model_info::type_subrange const*>(v->get_type())
              ->get_minvalue()
              ->get());
    }
    return std::numeric_limits<int>::min();
  }

  int get_var_maxvalue(model_info::var const* const v) const {
    if (v->get_type()->is_boolean()) {
      return 1;
    }
    if (v->get_type()->is_subrange()) {
      return static_cast<int>(
          static_cast<model_info::type_subrange const*>(v->get_type())
              ->get_maxvalue()
              ->get());
    }
    return std::numeric_limits<int>::max();
  }

  void add_param(model_info::cnst const* const p) const;
  bool add_param(const std::string& dpar,
                 model_info::cnst const*& dparam) const;
  bool check_dparam(const std::string& dpar, const std::string& dpar2,
                    model_info::cnst const*& dparami,
                    model_info::cnst const*& dparamj) const;

  void add_ivar(model_info::var const* const iv) const;
  void add_evar(model_info::var const* const ev) const;

  std::size_t cov_index_(std::size_t const* const I, const std::size_t d,
                         const std::size_t len) const {
    uint64_t k = 1;
    for (std::size_t i = 0; i < len - d; i++) {
      k *= I[d] + i;
    }
    for (std::size_t i = 2; i <= len - d; i++) {
      k /= i;
    }
    if (d + 1 < len) {
      k += cov_index_(I, d + 1, len);
    }
    return static_cast<std::size_t>(k);
  }

  std::size_t cov_indexl(std::size_t const* const I,
                         const std::size_t len) const {
    return cov_index_(I, 0, len);
  }

  const std::vector<parser::AST const*>& get_exprs() const { return exprs; }

  void add_expr(parser::AST const* const et, int& et_index) const {
    et_index = exprs.size();
    exprs.push_back(et);
  }

  parser::AST const* build_rate_g(parser::AST const* const et) const {
    if (expr_contains_pvar(et)) {
      if (et->is_binary()) {
        parser::binaryAST const* etb =
            static_cast<parser::binaryAST const* const>(et);
        if (etb->get_op_tok() == parser::TOK_MULT) {
          parser::AST const* etl = parser::AST::et_1;
          if (expr_contains_not_pvar(etb->get_lterm())) {
            etl = build_rate_g(etb->get_lterm());
          }

          parser::AST const* etr = parser::AST::et_1;
          if (expr_contains_not_pvar(etb->get_rterm())) {
            etr = build_rate_g(etb->get_rterm());
          }

          return build_product(etl, etr);
        } else if (etb->get_op_tok() == parser::TOK_DIV) {
          parser::AST const* etl = parser::AST::et_1;
          if (expr_contains_not_pvar(etb->get_lterm())) {
            etl = build_rate_g(etb->get_lterm());
          }

          parser::AST const* etr = parser::AST::et_1;
          if (expr_contains_not_pvar(etb->get_rterm())) {
            etr = build_rate_g(etb->get_rterm());
          }

          return build_fraction(etl, etr);
        }
      }
      return parser::AST::et_1;
    }
    return et;
  }

  parser::AST const* build_rate_h(parser::AST const* const et) const {
    if (expr_contains_pvar(et)) {
      if (et->is_binary()) {
        parser::binaryAST const* etb =
            static_cast<parser::binaryAST const* const>(et);
        if (etb->get_op_tok() == parser::TOK_MULT) {
          parser::AST const* etl = parser::AST::et_1;
          if (expr_contains_pvar(etb->get_lterm())) {
            etl = build_rate_h(etb->get_lterm());
          }

          parser::AST const* etr = parser::AST::et_1;
          if (expr_contains_pvar(etb->get_rterm())) {
            etr = build_rate_h(etb->get_rterm());
          }

          return build_product(etl, etr);
        } else if (etb->get_op_tok() == parser::TOK_DIV) {
          parser::AST const* etl = parser::AST::et_1;
          if (expr_contains_pvar(etb->get_lterm())) {
            etl = build_rate_h(etb->get_lterm());
          }

          parser::AST const* etr = parser::AST::et_1;
          if (expr_contains_pvar(etb->get_rterm())) {
            etr = build_rate_h(etb->get_rterm());
          }

          return build_fraction(etl, etr);
        }
      }
      return et;
    }
    return parser::AST::et_1;
  }

  parser::AST const* build_expr_der(parser::AST const* const et,
                                    model_info::cnst const* const c) const {
    return build_expr_der<model_info::cnst>(et, c);
  }

  parser::AST const* build_expr_der(parser::AST const* const et,
                                    model_info::var const* const v) const {
    return build_expr_der<model_info::var>(et, v);
  }

  void expr_vars_params(parser::AST const* const et,
                        std::set<model_info::var const*>& ev,
                        std::set<model_info::cnst const*>& ep) const {
    if (et->is_var()) {
      ev.insert(static_cast<parser::varAST const*>(et)->get_var());
    } else if (et->is_cnst()) {
      model_info::cnst const* const c =
          static_cast<parser::cnstAST const*>(et)->get_cnst();
      if (c->get_index_p() >= 0) {
        ep.insert(c);
      }
    } else if (et->is_unary()) {
      expr_vars_params(static_cast<parser::unaryAST const*>(et)->get_term(), ev,
                       ep);
    } else if (et->is_binary()) {
      parser::binaryAST const* const etb =
          static_cast<parser::binaryAST const*>(et);
      expr_vars_params(etb->get_lterm(), ev, ep);
      expr_vars_params(etb->get_rterm(), ev, ep);
    } else if (et->is_polynom()) {
      parser::polynomAST const* const etp =
          static_cast<parser::polynomAST const*>(et);
      for (auto const& a : etp->get_terms()) {
        for (auto const& ci : a.c) {
          if (ci.first->get_index_p() >= 0) {
            ep.insert(ci.first);
          }
        }
        for (auto const& vi : a.v) {
          ev.insert(vi.first);
        }
      }
    } else if (et->is_stdfcall()) {
      parser::stdfcallAST const* const etf =
          static_cast<parser::stdfcallAST const*>(et);
      for (auto a : etf->get_args()) {
        expr_vars_params(a, ev, ep);
      }
    } else if (et->is_userfcall()) {
      parser::userfcallAST const* const etf =
          static_cast<parser::userfcallAST const*>(et);
      for (auto a : etf->get_args()) {
        expr_vars_params(a, ev, ep);
      }
    } else {
    }
  }

  bool expr_is_param(parser::AST const* const et) const {
    if (et->is_cnst()) {
      model_info::cnst const* const c =
          static_cast<parser::cnstAST const* const>(et)->get_cnst();
      auto const& pi = std::find_if(
          sa_params.begin(), sa_params.end(),
          [&](model_info::cnst const* const c_) { return c_ == c; });
      return pi != sa_params.end();
    }

    if (et->is_var()) {
      model_info::var const* const iv =
          static_cast<parser::varAST const* const>(et)->get_var();
      auto const& pi = std::find_if(
          sa_ivars.begin(), sa_ivars.end(),
          [&](model_info::var const* const v_) { return v_ == iv; });
      return pi != sa_ivars.end();
    }

    return false;
  }

  void update_err(const double p, const double q, double& d, double& d2,
                  double& d_max, double& r, double& r2, double& r_max,
                  double& d_r_max) const {
    const double d_ = std::fabs(p - q);
    d += d_;
    d2 += d_ * d_;
    d_max = std::max(d_max, d_);

    const double r_ = std::fabs(q);
    r += r_;
    r2 += r_ * r_;
    r_max = std::max(r_max, r_);

    d_r_max = std::max(d_r_max, d_ / r_);
  }

  std::ostream& write_nrm_header(std::ostream& os) const {
    const char csvSep = ',';
    os << "L1" << csvSep << "L2" << csvSep << "Linf" << csvSep << "L1(rel)"
       << csvSep << "L2(rel)" << csvSep << "Linf(rel)" << csvSep
       << "Linf(rel,c)";
    return os;
  }

  std::ostream& write_nrm(std::ostream& os, const double d, const double d2,
                          const double d_max, const double r, const double r2,
                          const double r_max, const double d_r_max) const {
    const char csvSep = ',';
    os << d << csvSep << std::sqrt(d2) << csvSep << d_max << csvSep << d / r
       << csvSep << std::sqrt(d2 / r2) << csvSep << d_max / r_max << csvSep
       << d_r_max;
    return os;
  }
};
}

#endif
