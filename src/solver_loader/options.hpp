/*
 *  options.hpp
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

#ifndef SOLVER_LOADER_OPTIONS_HPP_
#define SOLVER_LOADER_OPTIONS_HPP_

#include <cmath>
#include <map>
#include <random>
#include <set>

extern "C" {
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
}

#include "task_info/base.hpp"
#include "model_info/var.hpp"

enum fnc_id {
  FN_LINSPACE,
  FN_LOGSPACE,
  FN_UNIFORM,
  FN_BERNOULLI,
  FN_BINOMIAL,
  FN_NEG_BINOMIAL,
  FN_GEOMETRIC,
  FN_POISSON,
  FN_EXPONENTIAL,
  FN_GAMMA,
  FN_WEIBULL,
  FN_NORMAL,
  FN_LOGNORMAL,
};

class fnc {
 public:
  fnc() {}

  virtual ~fnc() {}

  virtual fnc_id get_id() const = 0;

  virtual bool is_distr() const { return false; }

  virtual bool get_values(std::vector<double>&) const { return false; }

  virtual double gen_sample() const { return 0.0; }
};

class fnc_linspace : public fnc {
 protected:
  const double a;
  const double b;
  const unsigned int n;

 public:
  fnc_linspace(const double a, const double b, const unsigned int n = 10)
      : a(a), b(b), n(n) {}

  double get_min() const { return a; }

  double get_max() const { return b; }

  unsigned int get_n() const { return n; }

  virtual fnc_id get_id() const { return FN_LINSPACE; }

  virtual bool get_values(std::vector<double>& vals) const {
    double x = a;
    const double dx = (b - a) / (double)n;
    for (std::size_t i = 0; i <= n; i++, x += dx) {
      vals.push_back(x);
    }
    return true;
  }
};

class fnc_logspace : public fnc_linspace {
 public:
  fnc_logspace(const double min, const double max, const unsigned int n = 10)
      : fnc_linspace(min, max, n) {}

  fnc_id get_id() const { return FN_LOGSPACE; }

  bool get_values(std::vector<double>& vals) const {
    double x = a;
    const double dx = (b - a) / (double)n;
    for (std::size_t i = 0; i <= n; i++, x += dx) {
      vals.push_back(std::pow(10.0, x));
    }
    return true;
  }
};

class fnc_distr_1 : public fnc {
 protected:
  mutable std::random_device rand_device;
  mutable std::mt19937 rand_engine;

  const double param1;

 public:
  fnc_distr_1(const double param1)
      : rand_engine(rand_device()), param1(param1) {}

  bool is_distr() const { return true; }

  double get_param1() const { return param1; }
};

class fnc_distr_2 : public fnc_distr_1 {
 protected:
  const double param2;

 public:
  fnc_distr_2(const double param1, const double param2)
      : fnc_distr_1(param1), param2(param2) {}

  double get_param2() const { return param2; }
};

class fnc_uniform : public fnc_distr_2 {
 private:
  mutable std::uniform_real_distribution<> distr;

 public:
  fnc_uniform(const double a, const double b)
      : fnc_distr_2(a, b), distr(a, b) {}

  fnc_id get_id() const { return FN_UNIFORM; }

  double gen_sample() const { return distr(rand_engine); }
};

class fnc_bernoulli : public fnc_distr_1 {
 private:
  mutable std::bernoulli_distribution distr;

 public:
  fnc_bernoulli(const double p) : fnc_distr_1(p), distr(p) {}

  fnc_id get_id() const { return FN_BERNOULLI; }

  double gen_sample() const { return distr(rand_engine); }
};

class fnc_binomial : public fnc_distr_2 {
 private:
  mutable std::binomial_distribution<> distr;

 public:
  fnc_binomial(const double k, const double p)
      : fnc_distr_2(k, p), distr(k, p) {}

  fnc_id get_id() const { return FN_BINOMIAL; }

  double gen_sample() const { return distr(rand_engine); }
};

class fnc_neg_binomial : public fnc_distr_2 {
 private:
  mutable std::negative_binomial_distribution<> distr;

 public:
  fnc_neg_binomial(const double k, const double p)
      : fnc_distr_2(k, p), distr(k, p) {}

  fnc_id get_id() const { return FN_NEG_BINOMIAL; }

  double gen_sample() const { return distr(rand_engine); }
};

class fnc_geometric : public fnc_distr_1 {
 private:
  mutable std::geometric_distribution<> distr;

 public:
  fnc_geometric(const double p) : fnc_distr_1(p), distr(p) {}

  fnc_id get_id() const { return FN_GEOMETRIC; }

  double gen_sample() const { return distr(rand_engine); }
};

class fnc_poisson : public fnc_distr_1 {
 private:
  mutable std::poisson_distribution<> distr;

 public:
  fnc_poisson(const double mu) : fnc_distr_1(mu), distr(mu) {}

  fnc_id get_id() const { return FN_POISSON; }

  double gen_sample() const { return distr(rand_engine); }
};

class fnc_exponential : public fnc_distr_1 {
 private:
  mutable std::exponential_distribution<> distr;

 public:
  fnc_exponential(const double lambda) : fnc_distr_1(lambda), distr(lambda) {}

  fnc_id get_id() const { return FN_EXPONENTIAL; }

  double gen_sample() const { return distr(rand_engine); }
};

class fnc_gamma : public fnc_distr_2 {
 private:
  mutable std::gamma_distribution<> distr;

 public:
  fnc_gamma(const double alpha, const double beta)
      : fnc_distr_2(alpha, beta), distr(alpha, beta) {}

  fnc_id get_id() const { return FN_GAMMA; }

  double gen_sample() const { return distr(rand_engine); }
};

class fnc_weibull : public fnc_distr_2 {
 private:
  mutable std::weibull_distribution<> distr;

 public:
  fnc_weibull(const double a, const double b)
      : fnc_distr_2(a, b), distr(a, b) {}

  fnc_id get_id() const { return FN_WEIBULL; }

  double gen_sample() const { return distr(rand_engine); }
};

class fnc_normal : public fnc_distr_2 {
 private:
  mutable std::normal_distribution<> distr;

 public:
  fnc_normal(const double mu, const double sigma)
      : fnc_distr_2(mu, sigma), distr(mu, sigma) {}

  virtual fnc_id get_id() const { return FN_NORMAL; }

  virtual double gen_sample() const { return distr(rand_engine); }
};

class fnc_lognormal : public fnc_distr_2 {
 private:
  mutable std::lognormal_distribution<> distr;

 public:
  fnc_lognormal(const double mu, const double sigma)
      : fnc_distr_2(mu, sigma), distr(mu, sigma) {}

  fnc_id get_id() const { return FN_LOGNORMAL; }

  virtual double gen_sample() const { return distr(rand_engine); }
};

enum option_t {
  OPT_BOOLNUM,
  OPT_INTNUM,
  OPT_REALNUM,
  OPT_STR,
  OPT_REALNUM_VEC,
  OPT_STR_VEC,
};

class option {
 public:
  virtual ~option() {}
  virtual option_t get_type() const = 0;

  virtual option* clone() const = 0;
};

class option_boolnum : public option {
 private:
  const bool value;

 public:
  option_boolnum(const bool value) : value(value) {}

  bool get_value() const { return value; }

  option_t get_type() const { return OPT_BOOLNUM; }

  option_boolnum* clone() const { return new option_boolnum(*this); }
};

class option_intnum : public option {
 private:
  const int value;

 public:
  option_intnum(const int value) : value(value) {}

  int get_value() const { return value; }

  option_t get_type() const { return OPT_INTNUM; }

  option_intnum* clone() const { return new option_intnum(*this); }
};

class option_realnum : public option {
 private:
  const double value;

 public:
  option_realnum(const double value) : value(value) {}

  double get_value() const { return value; }

  option_t get_type() const { return OPT_REALNUM; }

  option_realnum* clone() const { return new option_realnum(*this); }
};

class option_str : public option {
 private:
  const std::string value;

 public:
  option_str(const std::string& value) : value(value) {}

  option_str(const option_str& str) : value(str.get_value()) {}

  const std::string& get_value() const { return value; }

  option_t get_type() const { return OPT_STR; }

  option_str* clone() const { return new option_str(*this); }
};

class option_realnum_arr : public option {
 private:
  const std::vector<double> value;

 public:
  option_realnum_arr(const std::vector<double>& value) : value(value) {}

  option_realnum_arr(const option_realnum_arr& dvec)
      : option_realnum_arr(dvec.get_value()) {}

  const std::vector<double>& get_value() const { return value; }

  option_t get_type() const { return OPT_REALNUM_VEC; }

  option_realnum_arr* clone() const { return new option_realnum_arr(*this); }
};

class option_str_arr : public option {
 private:
  const std::vector<std::string> value;

 public:
  option_str_arr(const std::vector<std::string>& value) : value(value) {}

  option_str_arr(const option_str_arr& strvec)
      : option_str_arr(strvec.get_value()) {}

  const std::vector<std::string>& get_value() const { return value; }

  option_t get_type() const { return OPT_STR_VEC; }

  option_str_arr* clone() const { return new option_str_arr(*this); }
};

class set_param_info {
 private:
  const std::string param;

 public:
  set_param_info(const std::string& param) : param(param) {}
  virtual ~set_param_info() {}

  virtual option_t get_type() const = 0;

  const std::string& get_param() const { return param; }
};

class set_param_info_bool : public set_param_info, public option_boolnum {
 public:
  set_param_info_bool(const std::string& param, const bool value)
      : set_param_info(param), option_boolnum(value) {}

  option_t get_type() const { return option_boolnum::get_type(); }
};

class set_param_info_int : public set_param_info, public option_intnum {
 public:
  set_param_info_int(const std::string& param, const int value)
      : set_param_info(param), option_intnum(value) {}

  option_t get_type() const { return option_intnum::get_type(); }
};

class set_param_info_double : public set_param_info, public option_realnum {
 public:
  set_param_info_double(const std::string& param, const double value)
      : set_param_info(param), option_realnum(value) {}

  option_t get_type() const { return option_realnum::get_type(); }
};

class set_param_info_str : public set_param_info, public option_str {
 public:
  set_param_info_str(const std::string& param, const std::string& value)
      : set_param_info(param), option_str(value) {}

  option_t get_type() const { return option_str::get_type(); }
};

struct est_param_info {
  std::string param;
  double min, max, value0;

  union {
    solver_loader::model_info::cnst const* c;
    solver_loader::model_info::var const* v;
  };
};

struct sample_param_info {
  std::string param;
  fnc const* f;

  union {
    solver_loader::model_info::cnst const* c;
    solver_loader::model_info::var const* v;
  };
};

struct iterate_param_info {
  std::string param;
  std::vector<double> values;

  union {
    solver_loader::model_info::cnst const* c;
    solver_loader::model_info::var const* v;
  };
};

namespace solver_loader {

class options {
 private:
  std::map<std::string, option*> options_map;

 public:
  options()
      : output_dir(""),
        do_steady_state(false),
        do_rare_event(false),
        do_sa(false),
        sa_no_der(false),
        sa_no_2der(false),
        kinetics(0),
        custom_model_id(-1),
        nmoments(2),
        det_centered(false),
        hs_has_index(false),
        hs_has_flags(false),
        hs_has_exitrate(false),
        hsucc_has_tr(false),
        hsucc_has_rates_g(true),
        tol(1e-15),
        abs_tol(1e-15),
        rel_tol(1e-3),
        norm_control(false),
        use_transient_subsolver(false),
        estimate_init(false),
        estimate_param(false),
        estimate_obs_error(false) {
    tspan.push_back(0.0);
    tspan.push_back(1.0);
  }

  options(options const* const o) {
    subsloader_handler = o->subsloader_handler;

    output_dir = o->output_dir;

    do_steady_state = o->do_steady_state;
    do_rare_event = o->do_rare_event;
    do_sa = o->do_sa;

    sa_no_der = o->sa_no_der;
    sa_no_2der = o->sa_no_2der;

    kinetics = o->kinetics;

    custom_model_id = o->custom_model_id;
    model_name = o->model_name;
    sbml_model_name = o->sbml_model_name;
    init_src = o->init_src;

    nmoments = o->nmoments;

    det_centered = o->det_centered;

    hybrid_vars_force_det = o->hybrid_vars_force_det;
    hybrid_vars_force_stoch = o->hybrid_vars_force_stoch;

    hs_has_index = o->hs_has_index;
    hs_has_flags = o->hs_has_flags;
    hs_has_exitrate = o->hs_has_exitrate;

    hsucc_has_tr = o->hsucc_has_tr;
    hsucc_has_rates_g = o->hsucc_has_rates_g;

    tol = o->tol;
    abs_tol = o->abs_tol;
    rel_tol = o->rel_tol;
    norm_control = o->norm_control;

    use_transient_subsolver = o->use_transient_subsolver;

    method_name = o->method_name;
    opt_method_name = o->opt_method_name;

    objf = o->objf;

    tspan.assign(o->tspan.begin(), o->tspan.end());

    for (auto const& t : o->tasks) {
      tasks.push_back(t->clone());
    }

    set_params.assign(o->set_params.begin(), o->set_params.end());
    set_inits.assign(o->set_inits.begin(), o->set_inits.end());
    set_obs_errors.assign(o->set_obs_errors.begin(), o->set_obs_errors.end());

    time_series_data_src = o->time_series_data_src;
    steady_state_data_src = o->steady_state_data_src;

    observables.assign(o->observables.begin(), o->observables.end());

    estimate_init = o->estimate_init;
    estimate_param = o->estimate_param;
    estimate_obs_error = o->estimate_obs_error;

    estimate_inits.assign(o->estimate_inits.begin(), o->estimate_inits.end());
    estimate_params.assign(o->estimate_params.begin(),
                           o->estimate_params.end());
    estimate_obs_errors.assign(o->estimate_obs_errors.begin(),
                               o->estimate_obs_errors.end());

    iterate_inits.assign(o->iterate_inits.begin(), o->iterate_inits.end());
    iterate_params.assign(o->iterate_params.begin(), o->iterate_params.end());
    iterate_obs_errors.assign(o->iterate_obs_errors.begin(),
                              o->iterate_obs_errors.end());

    sample_inits.assign(o->sample_inits.begin(), o->sample_inits.end());
    sample_params.assign(o->sample_params.begin(), o->sample_params.end());
    sample_obs_errors.assign(o->sample_obs_errors.begin(),
                             o->sample_obs_errors.end());

    for (auto const& o_ : o->get_options()) {
      options_map[o_.first] = o_.second->clone();
    }
  }

  virtual ~options() {
    unset();
    clear_tasks();
  }

  const std::map<std::string, option*>& get_options() const {
    return options_map;
  }

  void set_stoch() { kinetics = 0; }

  void set_det() { kinetics = 1; }

  void set_hybrid() { kinetics = 2; }

  bool is_stoch() const { return kinetics == 0; }

  bool is_det() const { return kinetics == 1; }

  bool is_hybrid() const { return kinetics == 2; }

  void set(const std::string& name) { set(name, true); }

  void set(const std::string& name, const bool value) {
    options_map[name] = new option_intnum(value);
  }

  void set(const std::string& name, const int value) {
    options_map[name] = new option_intnum(value);
  }

  void set(const std::string& name, const unsigned int value) {
    options_map[name] = new option_intnum(value);
  }

  void set(const std::string& name, const long int value) {
    options_map[name] = new option_intnum(value);
  }

  void set(const std::string& name, const double value) {
    options_map[name] = new option_realnum(value);
  }

  void set(const std::string& name, const char* value) {
    options_map[name] = new option_str(value);
  }

  void set(const std::string& name, const std::vector<double>& value) {
    options_map[name] = new option_realnum_arr(value);
  }

  void set(const std::string& name, const std::vector<std::string>& value) {
    options_map[name] = new option_str_arr(value);
  }

  bool get(const std::string& name, bool& value) const {
    auto const& vi = options_map.find(name);
    if (vi != options_map.end()) {
      if (vi->second->get_type() == OPT_INTNUM) {
        value = static_cast<bool>(
            static_cast<option_intnum*>(vi->second)->get_value());
      } else {
        return false;
      }
      return true;
    }
    return false;
  }

  bool get(const std::string& name, int& value) const {
    auto const& vi = options_map.find(name);
    if (vi != options_map.end()) {
      if (vi->second->get_type() == OPT_INTNUM) {
        value = static_cast<option_intnum*>(vi->second)->get_value();
      } else if (vi->second->get_type() == OPT_REALNUM) {
        value = static_cast<int>(
            static_cast<option_realnum*>(vi->second)->get_value());
      } else {
        return false;
      }
      return true;
    }
    return false;
  }

  bool get(const std::string& name, unsigned int& value) const {
    auto const& vi = options_map.find(name);
    if (vi != options_map.end()) {
      if (vi->second->get_type() == OPT_INTNUM) {
        value = static_cast<option_intnum*>(vi->second)->get_value();
      } else if (vi->second->get_type() == OPT_REALNUM) {
        value = static_cast<unsigned int>(
            static_cast<option_realnum*>(vi->second)->get_value());
      } else {
        return false;
      }
      return true;
    }
    return false;
  }

  bool get(const std::string& name, long int& value) const {
    auto const& vi = options_map.find(name);
    if (vi != options_map.end()) {
      if (vi->second->get_type() == OPT_INTNUM) {
        value = static_cast<option_intnum*>(vi->second)->get_value();
      } else if (vi->second->get_type() == OPT_REALNUM) {
        value = static_cast<long int>(
            static_cast<option_realnum*>(vi->second)->get_value());
      } else {
        return false;
      }
      return true;
    }
    return false;
  }

  bool get(const std::string& name, double& value) const {
    auto const& vi = options_map.find(name);
    if (vi != options_map.end()) {
      if (vi->second->get_type() == OPT_REALNUM) {
        value = static_cast<option_realnum*>(vi->second)->get_value();
      } else if (vi->second->get_type() == OPT_INTNUM) {
        value = static_cast<double>(
            static_cast<option_intnum*>(vi->second)->get_value());
      } else {
        return false;
      }
      return true;
    }
    return false;
  }

  bool get(const std::string& name, std::string& value) const {
    auto const& vi = options_map.find(name);
    if (vi != options_map.end() && vi->second->get_type() == OPT_STR) {
      value = static_cast<option_str*>(vi->second)->get_value();
      return true;
    }
    return false;
  }

  bool get(const std::string& name, std::vector<double>& value) const {
    auto const& vi = options_map.find(name);
    if (vi != options_map.end() && vi->second->get_type() == OPT_REALNUM_VEC) {
      value = static_cast<option_realnum_arr*>(vi->second)->get_value();
      return true;
    }
    return false;
  }

  bool get(const std::string& name, std::vector<std::string>& value) const {
    auto const& vi = options_map.find(name);
    if (vi != options_map.end() && vi->second->get_type() == OPT_STR_VEC) {
      value = static_cast<option_str_arr*>(vi->second)->get_value();
      return true;
    }
    return false;
  }

  bool isset(const std::string& name) const {
    return options_map.find(name) != options_map.end();
  }

  void unset(const std::string& name) { options_map.erase(name); }

  void unset() {
    for (auto& s : options_map) {
      s.second = nullptr;
    }
    options_map.clear();
  }

  void clear_tasks() {
    for (auto& t : tasks) {
    }
    tasks.clear();
  }

  int sprint_description_param(char* const s, char const* const pname,
                               bool const pval,
                               const bool highlight = false) const {
    return sprint_description_param(s, pname, pval ? "true" : "false",
                                    highlight);
  }

  int sprint_description_param(char* const s, char const* const pname,
                               int const pval,
                               const bool highlight = false) const {
    return sprintf(s, highlight ? "\033[1;32m%s:\033[0m \033[0;34m%d;\033[0m "
                                : "%s: %d; ",
                   pname, pval);
  }

  int sprint_description_param(char* const s, char const* const pname,
                               double const pval,
                               const bool highlight = false) const {
    return sprintf(s, highlight ? "\033[1;32m%s:\033[0m \033[0;34m%g;\033[0m "
                                : "%s: %g; ",
                   pname, pval);
  }

  int sprint_description_param(char* const s, char const* const pname,
                               char const* const pval,
                               const bool highlight = false) const {
    return sprintf(s, highlight ? "\033[1;32m%s:\033[0m \033[0;34m%s;\033[0m "
                                : "%s: %s; ",
                   pname, pval);
  }

  int sprint_description_int_setting(char* const s, char const* const sname,
                                     char const* const stitle,
                                     const bool highlight = false) const {
    int value;
    if (get(sname, value)) {
      return sprint_description_param(s, stitle, value, highlight);
    }
    return 0;
  }

  int sprint_description_double_setting(char* const s, char const* const sname,
                                        char const* const stitle,
                                        const bool highlight = false) const {
    double value;
    if (get(sname, value)) {
      return sprint_description_param(s, stitle, value, highlight);
    }
    return 0;
  }

  int sprint_description_string_setting(char* const s, char const* const sname,
                                        char const* const stitle,
                                        const bool highlight = false) const {
    std::string value;
    if (get(sname, value)) {
      return sprint_description_param(s, stitle, value.c_str(), highlight);
    }
    return 0;
  }

  int sprint_description_solver_name(char* const s, char const* const name,
                                     const bool highlight = false) const {
    return sprint_description_param(s, "solver", name, highlight);
  }

  int sprint_description_method_name(char* s,
                                     const bool highlight = false) const {
    if (method_name != "") {
      return sprint_description_param(s, "method", method_name.c_str(),
                                      highlight);
    }
    return 0;
  }

  int sprint_description_model_name(char* s,
                                    const bool highlight = false) const {
    if (model_name != "") {
      return sprint_description_param(s, "model", model_name.c_str(),
                                      highlight);
    }
    if (sbml_model_name != "") {
      return sprint_description_param(s, "sbml_model", sbml_model_name.c_str(),
                                      highlight);
    }
    return 0;
  }

  int sprint_description_vars(char* s,
                              const std::vector<model_info::var const*>& vars,
                              char const* const stitle) const {
    if (vars.empty()) {
      return 0;
    }

    char const* const t = s;
    char* c = s;

    c += sprintf(c, "%s: { ", stitle);

    bool f = true;
    for (auto const& v : vars) {
      if (!f) {
        c += sprintf(c, ", ");
      } else {
        f = false;
      }
      c += sprintf(c, "%s", v->get_name().c_str());
    }

    c += sprintf(c, " }; ");

    return c - t;
  }

  int sprint_description_set_params(
      char* s, const std::vector<set_param_info const*>& params,
      char const* const stitle, const bool highlight = false) const {
    if (params.empty()) {
      return 0;
    }

    char* c = s;
    int r = 0;

    int r_ = sprintf(c, "%s: { ", stitle);
    c += r_;
    r += r_;

    for (auto const& p : params) {
      switch (p->get_type()) {
        case OPT_BOOLNUM:
          r_ = sprint_description_param(
              c, p->get_param().c_str(),
              static_cast<set_param_info_bool const*>(p)->get_value(),
              highlight);
          break;

        case OPT_INTNUM:
          r_ = sprint_description_param(
              c, p->get_param().c_str(),
              static_cast<set_param_info_int const*>(p)->get_value(),
              highlight);
          break;

        case OPT_REALNUM:
          r_ = sprint_description_param(
              c, p->get_param().c_str(),
              static_cast<set_param_info_double const*>(p)->get_value(),
              highlight);
          break;

        case OPT_STR:
          r_ = sprint_description_param(
              c, p->get_param().c_str(),
              static_cast<set_param_info_str const*>(p)->get_value().c_str(),
              highlight);
          break;

        default:
          r_ = 0;
          break;
      }
      c += r_;
      r += r_;
    }

    r_ = sprintf(c, "}; ");
    c += r_;
    r += r_;

    return r;
  }

  lua_CFunction subsloader_handler;

  std::string output_dir;

  bool do_steady_state;
  bool do_rare_event;
  bool do_sa;

  bool sa_no_der;
  bool sa_no_2der;

  int kinetics;

  int custom_model_id;
  std::string model_name, sbml_model_name;
  std::string init_src;

  std::string method_name;
  std::string opt_method_name;

  std::string objf;
  std::vector<parser::AST const*> objf_et;
  std::vector<int> objf_index;

  std::vector<double> tspan;

  std::vector<solver_loader::task_info::base*> tasks;

  unsigned int nmoments;

  bool det_centered;

  std::set<std::string> hybrid_vars_force_det;
  std::set<std::string> hybrid_vars_force_stoch;

  bool hs_has_index;
  bool hs_has_flags;
  bool hs_has_exitrate;

  bool hsucc_has_tr;
  bool hsucc_has_rates_g;

  double tol, abs_tol, rel_tol;
  bool norm_control;

  bool use_transient_subsolver;

  std::vector<set_param_info const*> set_params;
  std::vector<set_param_info const*> set_inits;
  std::vector<set_param_info const*> set_obs_errors;

  std::string time_series_data_src;
  std::string steady_state_data_src;

  std::vector<std::string> observables;

  bool estimate_init;
  bool estimate_param;
  bool estimate_obs_error;

  std::vector<est_param_info> estimate_inits;
  std::vector<est_param_info> estimate_params;
  std::vector<est_param_info> estimate_obs_errors;

  std::vector<sample_param_info> sample_inits;
  std::vector<sample_param_info> sample_params;
  std::vector<sample_param_info> sample_obs_errors;

  std::vector<iterate_param_info> iterate_inits;
  std::vector<iterate_param_info> iterate_params;
  std::vector<iterate_param_info> iterate_obs_errors;
};
}

#endif
