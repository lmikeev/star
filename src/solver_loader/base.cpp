/*
 *  base.cpp
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

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "base.hpp"
#include "model_info/trsys.hpp"
#include "task_info/all.hpp"
#include "parser/parser.hpp"
#include "parser/sbml.hpp"
#include "writer/cpp.hpp"
#include "writer/matlab.hpp"
#include "../experiment.hpp"
#include "../solver/model.hpp"
#include "../solver/base.hpp"

namespace solver_loader {

base::base(experiment* const e, const bool issubsloader)
    : model_info::trsys("", nullptr),
      options(e->get_options()),
      E(e),
      issubsloader(issubsloader),
      S(nullptr),
      id(0)
#ifdef STAR_CODEGEN
      ,
      solverlibHandle(nullptr)
#endif
{
}

base::~base() {}

bool base::is_var_control(model_info::var const* const v) const {
  if (is_det()) {
    return false;
  }

  if (is_stoch()) {
    return true;
  }

  if (hybrid_vars_force_stoch.count(v->get_name())) {
    return true;
  }

  if (hybrid_vars_force_det.count(v->get_name())) {
    return false;
  }

  if (v->get_type()->is_boolean()) {
    return true;
  }

  if (v->get_type()->is_subrange()) {
    auto const t = static_cast<model_info::type_subrange const*>(v->get_type());
    if (static_cast<int>(t->get_maxvalue()->get()) -
            static_cast<int>(t->get_minvalue()->get()) <
        20) {
      return true;
    }
  }

  return false;
}

char* base::sprint_description(char* const buf) const {
  char* c = buf;
  c += sprint_description_model_name(c);
  c += sprint_description_set_params(c, set_params, "params");
  c += sprint_description_set_params(c, set_inits, "init");
  return sprint_solver_description(c);
}

bool base::create_dir() const {
  namespace fs = boost::filesystem;
  try {
    fs::path p = E->get_dir() / dir;
    fs::create_directory(p);

    tmp_dir = p / "temp";
    fs::create_directory(tmp_dir.c_str());
  } catch (const fs::filesystem_error& e) {
    last_err << e.what();
    return false;
  }
  return true;
}

#ifdef STAR_CODEGEN
bool base::load_solver_lib() {
  char cmd[1024];
#ifndef WIN32
  sprintf(cmd, "make -C %s all", tmp_dir.c_str());
#else

  sprintf(cmd, "mingw32-make -C %s all", tmp_dir.c_str());
#endif
  if (system(cmd)) {
    last_err << "could not compile solver";
    return false;
  }

  boost::filesystem::path p = tmp_dir / "libsolver.so";

  solver::base* (*fnc_gs)(base* const);
#ifdef WIN32
  solverlibHandle = LoadLibrary(E->convws(p.wstring()));
  if (solverlibHandle == nullptr) {
    last_err << "could not load solver library";
    return false;
  }

  *(void**)(&fnc_gs) = (void*)GetProcAddress(solverlibHandle, "get_solver");
  if (fnc_gs == nullptr) {
    last_err << "could not get function get_solver";
    return false;
  }
#else
  solverlibHandle = dlopen(p.c_str(), RTLD_NOW);
  if (solverlibHandle == nullptr) {
    last_err << dlerror();
    return false;
  }

  *(void**)(&fnc_gs) = dlsym(solverlibHandle, "get_solver");
  if (fnc_gs == nullptr) {
    last_err << dlerror();
    return false;
  }
#endif

  S = (*fnc_gs)(this);

  return true;
}
#endif

template <>
bool base::expr_eq_obj<model_info::cnst>(
    parser::AST const* const et, model_info::cnst const* const c) const {
  return et->is_cnst() &&
         static_cast<parser::cnstAST const* const>(et)->get_cnst() == c;
}

template <>
bool base::expr_eq_obj<model_info::var>(parser::AST const* const et,
                                        model_info::var const* const v) const {
  return et->is_var() &&
         static_cast<parser::varAST const* const>(et)->get_var() == v;
}

template <>
bool base::cnst_eq_obj<model_info::cnst>(
    model_info::cnst const* const c, model_info::cnst const* const o) const {
  return c == o;
}

template <>
bool base::var_eq_obj<model_info::var>(model_info::var const* const v,
                                       model_info::var const* const o) const {
  return v == o;
}

void base::ic_dfs(const std::size_t u,
                  const std::vector<std::vector<std::size_t>>& g,
                  std::vector<bool>& dfs_vis,
                  std::vector<std::size_t>& topo) const {
  dfs_vis[u] = true;
  for (auto const v : g[u]) {
    if (!dfs_vis[v]) {
      ic_dfs(v, g, dfs_vis, topo);
    }
  }
  topo.push_back(u);
}

bool base::check_model_options() {
  for (auto const& t : get_types()) {
    if (t->is_subrange()) {
      auto const& t_ = static_cast<model_info::type_subrange const*>(t);
      compute_value(t_->get_minvalue());
      compute_value(t_->get_maxvalue());
    }
  }

  for (auto const& sp : set_params) {
    model_info::cnst const* const c = findcnst(sp->get_param());
    if (c == nullptr) {
      last_err << "unknown parameter '" << sp->get_param() << "'";
      return false;
    }

    if (c->get_type()->is_boolean()) {
      if (sp->get_type() != OPT_BOOLNUM) {
        last_err << "couldn't set '" << c->get_name() << "', boolean expected";
        return false;
      }
      c->get_value()->set(static_cast<model_info::val_t>(
          static_cast<set_param_info_bool const* const>(sp)->get_value()));
    } else if (c->get_type()->is_integer()) {
      if (sp->get_type() != OPT_INTNUM) {
        last_err << "couldn't set '" << c->get_name() << "', integer expected";
        return false;
      }

      const model_info::val_t v_ = static_cast<model_info::val_t>(
          static_cast<set_param_info_int const* const>(sp)->get_value());

      if (c->get_type()->is_subrange()) {
        auto const& t_ =
            static_cast<model_info::type_subrange const*>(c->get_type());
        const int minv_ = static_cast<int>(t_->get_minvalue()->get());
        const int maxv_ = static_cast<int>(t_->get_maxvalue()->get());
        const int vi_ = static_cast<int>(v_);
        if (vi_ < minv_ || vi_ > maxv_) {
          last_err << "couldn't set '" << c->get_name()
                   << "', the value is out of range";
          return false;
        }
      }

      c->get_value()->set(v_);
    } else {
      if (sp->get_type() != OPT_REALNUM) {
        last_err << "couldn't set '" << c->get_name() << "', real expected";
        return false;
      }
      c->get_value()->set(static_cast<model_info::val_t>(
          static_cast<set_param_info_double const* const>(sp)->get_value()));
    }
  }

  for (std::size_t i = 0; i < get_cnsts().size(); i++) {
    auto const& c = get_cnsts()[i];

    if (c->get_type()->is_subrange()) {
      auto const& t_ =
          static_cast<model_info::type_subrange const*>(c->get_type());
      compute_value(t_->get_minvalue());
      compute_value(t_->get_maxvalue());
    }

    if (!c->get_value()->is_set()) {
      compute_value(c->get_value());
    }

    c->set_index(i);
  }

  for (auto const& v : get_vars()) {
    if (v->get_type()->is_subrange()) {
      auto const& t_ =
          static_cast<model_info::type_subrange const*>(v->get_type());
      compute_value(t_->get_minvalue());
      compute_value(t_->get_maxvalue());
    }

    if (is_var_control(v)) {
      c_vars.push_back(v);
    } else {
      p_vars.push_back(v);
    }
  }

  for (std::size_t i = 0; i < p_vars.size(); i++) {
    p_vars[i]->set_index(i);
  }
  for (std::size_t i = 0; i < c_vars.size(); i++) {
    c_vars[i]->set_index(p_vars.size() + i);
  }

  for (auto const& tr : get_transitions()) {
    if (!tr->get_updates().empty()) {
      trs_stoch.push_back(tr);
    } else {
      trs_det.push_back(tr);
    }
  }

  for (auto const& o : observables) {
    solver_loader::model_info::var const* const v = findvar(o);
    if (v == nullptr) {
      last_err << "unknown variable '" << o << "'";
      return false;
    }
    o_vars.insert(v->get_index());
  }

  o_err.assign(get_vars().size(), 0.0);
  for (auto& oe : set_obs_errors) {
    solver_loader::model_info::var const* const v = findvar(oe->get_param());
    if (v == nullptr) {
      last_err << "unknown variable '" << oe->get_param() << "'";
      return false;
    }

    switch (oe->get_type()) {
      case OPT_BOOLNUM:
        if (static_cast<set_param_info_bool const*>(oe)->get_value()) {
          o_err[v->get_index()] = 1.0;
        }
        break;

      case OPT_INTNUM:
        o_err[v->get_index()] =
            static_cast<set_param_info_int const*>(oe)->get_value();
        break;

      case OPT_REALNUM:
        o_err[v->get_index()] =
            static_cast<set_param_info_double const*>(oe)->get_value();
        break;

      default:
        last_err << "invalid value of an observation error";
        return false;
    }
  }

  for (auto const& ic : get_ics()) {
    for (auto const& ics : ic->get_states()) {
      const int V = (int)ics.li.size();

      std::vector<std::vector<std::size_t>> g(V);
      std::vector<bool> dfs_vis(V, false);
      std::vector<std::size_t> topo;

      for (auto const& icsv : ics.li) {
        std::set<model_info::var const*> ev;
        std::set<model_info::cnst const*> ep;

        expr_vars_params(icsv.value->get_expr(), ev, ep);

        for (auto const& dv : ev) {
          g[dv->get_index()].push_back(icsv.v->get_index());
        }
      }

      for (int i = 0; i < V; i++) {
        if (!dfs_vis[i]) {
          ic_dfs(i, g, dfs_vis, topo);
        }
      }

      for (int i = V - 1; i >= 0; i--) {
        compute_value(ics.li[topo[i]].value, true);
      }
    }
  }

  sort(tasks.begin(), tasks.end(),
       [=](solver_loader::task_info::base const* const ti1,
           solver_loader::task_info::base const* const ti2) {
    return ti1->get_type() < ti2->get_type() ||
           (ti1->get_type() == ti2->get_type() &&
            strcmp(ti1->get_name(), ti2->get_name()) <= 0);
  });

  if (!check_tasks()) {
    return false;
  }

  return true;
}

bool base::load_solver() {
#ifdef STAR_CODEGEN
  if (!write_solver()) {
    return false;
  }

  if (!load_solver_lib()) {
    return false;
  }
#else
  last_err << "not yet supported";
  return false;
#endif

  return true;
}

bool base::open_log(std::ofstream& os) const {
  boost::filesystem::path p = E->get_dir() / get_dir() / "log.txt";
  os.open(E->convws(p.wstring()), std::ios::out | std::ios::app);
  if (!os.is_open()) {
    return false;
  }
  time_t rawtime;
  time(&rawtime);
  os << "=== " << ctime(&rawtime);
  return true;
}

void base::close_log(std::ofstream& os) const {
  os << std::endl;
  os.close();
}

bool base::compile_solver() {
  bool success = preload_check_options() && load_model() &&
                 check_model_options() && afterload_check_options() &&
                 post_check_tasks();

  static char tmp[1024];
  sprint_description(tmp);

#ifdef STAR_WEB_INTERFACE
  id = E->get_dbconnector()->create_result(E->get_id(), tmp,
                                           "Compiling solver...");
  if (!id) {
    last_err << "Error creating result";
    return false;
  }
#else
  static unsigned int id_ = 1;
  id = id_++;

  E->print_status("Compiling solver...");
  printf("%s\n", tmp);
#endif

  sprintf(dir, "%d", id);

  std::ofstream log;
  if (!create_dir() || !load_solver() || !open_log(log)) {
    success = false;
  }

  if (log.is_open()) {
    log << "PID : " << getpid() << std::endl << std::endl;
    log << tmp << std::endl;
    close_log(log);
  }

  {
    boost::filesystem::path p = E->get_dir() / get_dir() / "cout.txt";
    fcout.open(E->convws(p.wstring()));
    if (!fcout.is_open()) {
      last_err << "Couldn't redirect cout";
      return false;
    }
    std::cout.rdbuf(fcout.rdbuf());
  }

  {
    boost::filesystem::path p = E->get_dir() / get_dir() / "cerr.txt";
    fcerr.open(E->convws(p.wstring()));
    if (!fcerr.is_open()) {
      last_err << "Couldn't redirect cerr";
      return false;
    }
    std::cerr.rdbuf(fcerr.rdbuf());
  }

#ifdef STAR_WEB_INTERFACE
  if (success) {
    if (!E->get_dbconnector()->set_result_progress(id, "Queued")) {
      last_err << "Error setting result progress";
      return false;
    }
  } else {
    if (!E->get_dbconnector()->end_result(id, false, last_err.str().c_str())) {
      last_err << "Error finishing result";
      return false;
    }
  }
#endif

  return success;
}

bool base::run_solver() {
#ifdef STAR_WEB_INTERFACE
  if (!E->get_dbconnector()->start_result(id, "Starting...")) {
    last_err << "Error starting solver";
    return false;
  }
#endif

  bool success = false;

  if (S->run()) {
    success = true;
  }

#ifdef STAR_WEB_INTERFACE
  if (!E->get_dbconnector()->end_result(
          id, success, success ? "Completed" : last_err.str().c_str())) {
    last_err << "Error finishing result";
    return false;
  }
#endif

#ifdef STAR_CODEGEN
#else
#endif
  return success;
}

bool base::update_progress(const double progress) const {
  static char tmp[128];
  sprintf(tmp, "%.1f%c completed", progress, '%');
#ifdef STAR_WEB_INTERFACE
  return get_experiment()->get_dbconnector()->set_result_progress(id, tmp);
#else
  std::ofstream log;
  if (open_log(log)) {
    log << tmp << std::endl;
    close_log(log);
  }
  std::cout << tmp << std::endl;
  return true;
#endif
}

bool base::update_progress(char const* const progress) const {
#ifdef STAR_WEB_INTERFACE
  return get_experiment()->get_dbconnector()->set_result_progress(id, progress);
#else
  std::cout << progress << std::endl;
  return true;
#endif
}

#ifdef STAR_CODEGEN
void base::write_set_params_(std::ostream& os,
                             const std::set<model_info::cnst const*>& ep,
                             char const* const tab) const {
  for (auto const& p : ep) {
    os << tab << "const ";
    if (p->get_type()->is_ordinal()) {
      os << "int ";
    } else {
      os << "double ";
    }
    os << "c" << p->get_index_p() << " = c[" << p->get_index_p() << "];"
       << std::endl;
  }
}

void base::write_set_params(std::ostream& os,
                            const std::set<model_info::cnst const*>& ep) const {
  if (!ep.empty()) {
    os << "    double const* const c = o->param_vals.data();" << std::endl;
    write_set_params_(os, ep);
  }
}

void base::write_set_vars_i(std::ostream& os,
                            const std::set<model_info::var const*>& ev,
                            char const* const tab) const {
  for (auto const& v : ev) {
    if (is_var_control(v)) {
      os << tab << "const int v" << v->get_index() << " = s->v"
         << v->get_index() << ";" << std::endl;
    }
  }
}

void base::write_set_vars_d(std::ostream& os,
                            const std::set<model_info::var const*>& ev,
                            char const* const tab) const {
  for (auto const& v : ev) {
    os << tab << "const double v" << v->get_index() << " = ";
    if (is_var_control(v)) {
      os << "s->v" << v->get_index();
    } else {
      os << "x[" << v->get_index() << "]";
    }
    os << ";" << std::endl;
  }
}

void base::write_drate_h_dc(std::ostream& os, parser::AST const* const rate_h,
                            const std::size_t ci,
                            const std::size_t rates_h_len) const {
  if (do_sa && !sa_no_der) {
    os << std::endl;

    for (std::size_t i = 0, ij = 0; i < sa_params.size(); i++) {
      parser::AST const* const drate_h_dc =
          build_expr_der(rate_h, sa_params[i]);
      if (drate_h_dc != parser::AST::et_0) {
        os << "      h[" << rates_h_len*(1 + i) + ci << "] = ";
        drate_h_dc->write(os, parser::cpp_writer, this) << ";" << std::endl;

        if (!sa_no_2der) {
          for (std::size_t j = 0; j <= i; j++, ij++) {
            parser::AST const* const d2rate_h_dc2 =
                build_expr_der(drate_h_dc, sa_params[j]);
            if (d2rate_h_dc2 != parser::AST::et_0) {
              os << "      h[" << rates_h_len*(1 + sa_params.size() + ij) + ci
                 << "] = ";
              d2rate_h_dc2->write(os, parser::cpp_writer, this) << ";"
                                                                << std::endl;
            }
          }
        }
      }
    }
  }
}

void base::write_drate_h_dx(std::ostream& os, parser::AST const* const rate_h,
                            const std::size_t ci0,
                            const std::size_t rates_h_len) const {
  if (nmoments > 1 && !p_vars.empty()) {
    os << std::endl;

    for (std::size_t i = 0; i < p_vars.size(); i++) {
      std::vector<std::size_t> I(nmoments, 0);
      std::vector<std::size_t> I_(p_vars.size(), 0);

      I[0] = i;
      I_[i]++;

      assert(I_[i] == 1);

      parser::AST const* const drate_h_dx = build_expr_der(rate_h, p_vars[i]);

      if (drate_h_dx != parser::AST::et_0) {
        const std::size_t ci = ci0 + 1 + i;
        os << "      h[" << ci << "] = ";
        drate_h_dx->write(os, parser::cpp_writer, this) << ";" << std::endl;

        write_drate_h_dc(os, drate_h_dx, ci, rates_h_len);

        write_drate_h_dx(os, I.data(), I_.data(), 1, ci0 + p_vars.size(),
                         drate_h_dx, rates_h_len);
      }

      I_[i]--;

      assert(I_[i] == 0);
    }
  }
}

void base::write_drate_h_dx(std::ostream& os, std::size_t* const I,
                            std::size_t* const I_, const std::size_t d,
                            const std::size_t ci0,
                            parser::AST const* const drate_h_dx,
                            const std::size_t rates_h_len) const {
  const std::size_t cov_len_ = solver::covlen_c(p_vars.size(), d + 1);

  for (std::size_t i = 0; i <= I[d - 1]; i++) {
    I[d] = i;
    I_[i]++;

    parser::AST const* const d2rate_h_dx2 =
        build_expr_der(drate_h_dx, p_vars[i]);

    if (d2rate_h_dx2 != parser::AST::et_0) {
      const std::size_t ci = ci0 + 1 + cov_indexl(I, d + 1);

      os << "      h[" << ci << "] = ";
      d2rate_h_dx2->write(os, parser::cpp_writer, this) << ";" << std::endl;

      write_drate_h_dc(os, d2rate_h_dx2, ci, rates_h_len);

      if (d + 1 < nmoments) {
        write_drate_h_dx(os, I, I_, d + 1, ci0 + cov_len_, d2rate_h_dx2,
                         rates_h_len);
      }
    }

    I_[i]--;
  }
}

void base::write_rate_h(std::ostream& os, parser::AST const* const rate_h,
                        const std::size_t ci0,
                        const std::size_t rates_h_len) const {
  os << std::endl;
  os << "      h[" << ci0 << "] = ";
  rate_h->write(os, parser::cpp_writer, this) << ";" << std::endl;

  if ((!p_vars.empty() && nmoments > 1) || (do_sa && !sa_no_der)) {
    write_drate_h_dc(os, rate_h, ci0, rates_h_len);
    write_drate_h_dx(os, rate_h, ci0, rates_h_len);
  }
}

void base::write_rate_h(std::ostream& os, parser::AST const* const rate_h,
                        std::size_t* const I, std::size_t* const I_,
                        const std::size_t d, const std::size_t ci0,
                        const std::size_t rates_h_len) const {
  const std::size_t x_len_ = solver::xlen_c(p_vars.size(), nmoments);
  const std::size_t rate_h_len = 1 + ((nmoments > 1) ? x_len_ : 0);

  const std::size_t cov_len_ = solver::covlen_c(p_vars.size(), d + 1);

  for (std::size_t i = 0; i <= I[d - 1]; i++) {
    I[d] = i;
    I_[i]++;

    parser::AST const* const rate_h_x = build_product(rate_h, p_vars[i]);
    if (rate_h_x != parser::AST::et_0) {
      const std::size_t ci = ci0 + cov_indexl(I, d + 1);
      write_rate_h(os, rate_h_x, ci * rate_h_len, rates_h_len);

      if (d + 1 < nmoments) {
        write_rate_h(os, rate_h_x, I, I_, d + 1, ci0 + cov_len_, rates_h_len);
      }
    }

    I_[i]--;
  }
}

bool base::write_usermodel() const {
  sprintf(model_class, "model_%d", id);

  std::ofstream os;
  boost::filesystem::path p = tmp_dir / "usermodel.hpp";
  os.open(E->convws(p.wstring()));
  if (!os.is_open()) {
    last_err << "Error creating " << p.c_str();
    return false;
  }
  os << "#ifndef USERMODEL_HPP_" << std::endl;
  os << "#define USERMODEL_HPP_" << std::endl;
  os << std::endl;
  os << "#include \"../../../../src/solver/options.hpp\"" << std::endl;

  os << std::endl;
  os << "struct cstate {" << std::endl;
  for (auto const& v : c_vars) {
    if (v->get_type()->is_boolean()) {
      os << "  unsigned int v" << v->get_index() << " : 1;" << std::endl;
    } else {
      if (v->get_type()->is_subrange()) {
        auto const t =
            static_cast<model_info::type_subrange const*>(v->get_type());
        unsigned int dv = static_cast<int>(t->get_maxvalue()->get()) -
                          static_cast<int>(t->get_minvalue()->get()) + 1;

        int sizebits = 0;
        do {
          sizebits++;
        } while (dv >>= 1);
        os << "  unsigned int v" << v->get_index() << " : " << sizebits << ";"
           << std::endl;
      } else {
        os << "  int v" << v->get_index() << ";" << std::endl;
      }
    }
  }
  os << "};" << std::endl;
  os << std::endl;
  os << "namespace solver {" << std::endl;
  os << std::endl;
  std::size_t tri = 0;

  for (auto const& tr : get_transitions()) {
    std::set<model_info::var const*> ev;
    std::set<model_info::cnst const*> ep;

    parser::AST const* const rate = build_expr(tr->get_rate(), tr);
    parser::AST const* const rate_g = build_rate_g(rate);
    parser::AST const* const rate_h = build_rate_h(rate);

    os << std::endl;
    os << "class tr" << tri << " : public transition {" << std::endl;
    os << "public:" << std::endl;
    os << "  tr" << tri << "(options const* const o, const std::size_t index) "
                           ": transition(o, index) {}" << std::endl;
    if (!is_det()) {
      os << std::endl;
      os << "  bool is_enabled(cstate const* const s, const double /*t*/ = "
            "0.0)const {" << std::endl;
      expr_vars_params(tr->get_guard(), ev, ep);
      write_set_params(os, ep);
      write_set_vars_i(os, ev);
      ev.clear();
      ep.clear();
      os << "    return ";
      tr->get_guard()->write(os, parser::cpp_writer, this) << ";" << std::endl;
      os << "  }" << std::endl;
      os << std::endl;
      os << "  void update(cstate* const s, const double /*t*/ = 0.0) const {"
         << std::endl;
      for (auto const& u : tr->get_updates()) {
        expr_vars_params(u.u, ev, ep);
      }
      write_set_params(os, ep);
      write_set_vars_i(os, ev);
      ev.clear();
      ep.clear();
      for (auto const& u : tr->get_updates()) {
        if (is_var_control(u.v)) {
          os << "    s->v" << u.v->get_index() << " = ";
          u.u->write(os, parser::cpp_writer, this) << ";" << std::endl;
        }
      }
      os << "  }" << std::endl;
    }
    if (tr->is_chemreaction()) {
      os << std::endl;
      os << "  double const* change() const {" << std::endl;
      os << "    static double ch[" << get_vars().size() << "];" << std::endl;
      for (auto const& re : static_cast<model_info::chemreaction const*>(tr)
                                ->get_stoichiometry()) {
        os << "    ch[" << re.v->get_index() << "] = " << re.c << ";"
           << std::endl;
      }
      os << "    return ch;" << std::endl;
      os << "  }" << std::endl;
    }
    os << std::endl;
    os << "  void rates_g(cstate const* const s, double* const g, const "
          "double /*t*/ = 0.0)"
          "const {" << std::endl;
    expr_vars_params(rate_g, ev, ep);
    write_set_params(os, ep);
    write_set_vars_d(os, ev);
    ev.clear();
    ep.clear();

    os << "    g[0] = std::max(static_cast<double>(";
    rate_g->write(os, parser::cpp_writer, this) << "), 0.0);" << std::endl;

    if (do_sa && !sa_no_der) {
      os << "    if (g[0] > 0.0) {" << std::endl;
      for (std::size_t i = 0, ij = 0; i < sa_params.size(); i++) {
        parser::AST const* const drate_g_dc =
            build_expr_der(rate_g, sa_params[i]);
        if (drate_g_dc != parser::AST::et_0) {
          os << "      g[" << 1 + i << "] = ";
          drate_g_dc->write(os, parser::cpp_writer, this) << ";" << std::endl;

          if (!sa_no_2der) {
            for (std::size_t j = 0; j <= i; j++, ij++) {
              parser::AST const* const d2rate_g_dc2 =
                  build_expr_der(drate_g_dc, sa_params[j]);
              if (d2rate_g_dc2 != parser::AST::et_0) {
                os << "      g[" << 1 + sa_params.size() + ij << "] = ";
                d2rate_g_dc2->write(os, parser::cpp_writer, this) << ";"
                                                                  << std::endl;
              }
            }
          }
        }
      }
      os << "    }" << std::endl;
    }
    os << "  }" << std::endl;

    if (!is_stoch()) {
      const std::size_t x_len_ = solver::xlen_c(p_vars.size(), nmoments);
      const std::size_t rate_h_len = 1 + ((nmoments > 1) ? x_len_ : 0);
      const std::size_t rate_h_cnt = 1 + ((!det_centered) ? x_len_ : 0);
      const std::size_t rates_h_len = rate_h_len * rate_h_cnt;

      os << std::endl;
      os << "  void rates_h(cstate const* const s, double const* const x, "
            "double* const h, const double /*t*/ = 0.0) const {" << std::endl;
      expr_vars_params(rate_h, ev, ep);
      if (!det_centered) {
        for (auto const v : get_vars()) {
          ev.insert(v);
        }
      }
      write_set_params(os, ep);
      write_set_vars_d(os, ev);
      ev.clear();
      ep.clear();

      os << "    h[0] = std::max(static_cast<double>(";
      rate_h->write(os, parser::cpp_writer, this) << "), 0.0);" << std::endl;

      os << "    if (h[0] > 0.0) {" << std::endl;

      write_drate_h_dc(os, rate_h, 0, rates_h_len);

      write_drate_h_dx(os, rate_h, 0, rates_h_len);

      if (!det_centered) {
        std::vector<std::size_t> I(nmoments, 0);
        std::vector<std::size_t> I_(p_vars.size(), 0);

        for (std::size_t i = 0; i < p_vars.size(); i++) {
          I[0] = i;
          I_[i]++;

          assert(I_[i] == 1);

          parser::AST const* const rate_h_x = build_product(rate_h, p_vars[i]);
          if (rate_h_x != parser::AST::et_0) {
            write_rate_h(os, rate_h_x, (1 + i) * rate_h_len, rates_h_len);

            if (nmoments > 1) {
              write_rate_h(os, rate_h_x, I.data(), I_.data(), 1,
                           1 + p_vars.size(), rates_h_len);
            }
          }

          I_[i]--;

          assert(I_[i] == 0);
        }
      }

      os << "    }" << std::endl;

      os << "  }" << std::endl;
    }
    os << std::endl;
    os << "};" << std::endl;
    tri++;
  }
  os << std::endl;
  os << "class " << model_class << " : public model {" << std::endl;
  os << "public:" << std::endl;
  os << "  " << model_class << "(options const* const o)" << std::endl;
  os << "    : model(o) {" << std::endl;
  for (std::size_t i = 0; i < get_transitions().size(); i++) {
    os << "    add_transition(new tr" << i << "(o, " << i << "));" << std::endl;
  }
  os << "  }" << std::endl;
  os << std::endl;
  os << "  void cs_set(cstate* const cs, int const* const x) const {"
     << std::endl;
  for (auto const& v : c_vars) {
    os << "    cs->v" << v->get_index() << " = x[" << v->get_index() << "];"
       << std::endl;
  }
  os << "  }" << std::endl;
  os << std::endl;
  os << "  void cs_get(int* const r, cstate const* const cs) const {"
     << std::endl;
  for (auto const& v : c_vars) {
    os << "    r[" << v->get_index() << "] = cs->v" << v->get_index() << ";"
       << std::endl;
  }
  os << "  }" << std::endl;
  os << "};" << std::endl;
  os << std::endl;
  os << "} //solver" << std::endl;
  os << std::endl;
  os << "#endif /* USERMODEL_HPP_ */" << std::endl;
  os.close();

  return true;
}

bool base::write_solverext() const {
  std::ofstream os;
  boost::filesystem::path p =
      tmp_dir / (!issubsloader ? "solverext.cpp" : "subsolverext.cpp");
  os.open(E->convws(p.wstring()));
  if (!os.is_open()) {
    last_err << "Error creating " << p.c_str();
    return false;
  }
  os << "#include \"solverext.hpp\"" << std::endl;
  std::set<model_info::var const*> ev;
  std::set<model_info::cnst const*> ep;
  for (std::size_t i = 0; i < exprs.size(); i++) {
    os << std::endl;
    os << "double expr_" << i << "(cstate const* const s, double const* const "
                                 "x, double const* const c) {" << std::endl;
    expr_vars_params(exprs[i], ev, ep);
    write_set_params_(os, ep, "  ");
    if (is_stoch()) {
      write_set_vars_i(os, ev, "  ");
    } else {
      write_set_vars_d(os, ev, "  ");
    }
    ev.clear();
    ep.clear();
    os << "  return ";
    exprs[i]->write(os, parser::cpp_writer, this) << ";" << std::endl;
    os << "}" << std::endl;
  }
  os << std::endl;
  os << "solver::Solver* g_S;" << std::endl;
  os << std::endl;
  os << "extern \"C\" solver::base* get_solver(solver_loader::base* const sl) {"
     << std::endl;
  os << "  g_S = new solver::Solver(" << solver_params_str() << ");"
     << std::endl;
  os << "  sl->set_model(new solver::" << model_class << "(g_S));" << std::endl;
  for (std::size_t i = 0; i < exprs.size(); i++) {
    os << "  g_S->get_exprs()[" << i << "] = expr_" << i << ";" << std::endl;
  }
  os << "  return g_S;" << std::endl;
  os << "}" << std::endl;
  os << std::endl;
  os << "extern \"C\" void free_solver() {" << std::endl;
  os << "  delete g_S;" << std::endl;
  os << "}" << std::endl;
  os.close();

  p = tmp_dir / "solverext.hpp";
  os.open(E->convws(p.wstring()));
  if (!os.is_open()) {
    last_err << "Error creating " << p.c_str();
    return false;
  }
  os << "#ifndef SOLVEREXT_HPP_" << std::endl;
  os << "#define SOLVEREXT_HPP_" << std::endl;
  os << std::endl;
  os << "#include \"usermodel.hpp\"" << std::endl;
  if (!write_solver_def(os)) {
    return false;
  }
  os << std::endl;
  os << "#endif /* SOLVEREXT_HPP_ */" << std::endl;
  os.close();

  return true;
}

bool base::write_Makefile() const {
  std::ofstream os;
  boost::filesystem::path p =
      tmp_dir / (!issubsloader ? "Makefile" : "Makefile_subsolver");
  os.open(E->convws(p.wstring()));
  if (!os.is_open()) {
    last_err << "Error creating " << p.c_str();
    return false;
  }

  if (!issubsloader) {
    os << "include ../../../../src/solver/Makefile" << std::endl;
  } else {
  }

#ifdef NDEBUG
  os << "SOLVERLIB_CXX_FLAGS  += -D'NDEBUG'" << std::endl;
#endif

  os << "SOLVERLIB_CXX_FLAGS  += -D'STAR_CODEGEN'" << std::endl;

  if (is_det()) {
    os << "SOLVERLIB_CXX_FLAGS  += -D'STAR_DET'" << std::endl;
  } else if (is_hybrid()) {
    os << "SOLVERLIB_CXX_FLAGS  += -D'STAR_HYBRID'" << std::endl;
  }

  if (!is_stoch() && det_centered) {
    os << "SOLVERLIB_CXX_FLAGS  += -D'STAR_DET_CENTERED'" << std::endl;
  }

  if (hsucc_has_tr) {
    os << "SOLVERLIB_CXX_FLAGS  += -D'STAR_HSUCC_HAS_TR'" << std::endl;
  }

  if (hsucc_has_rates_g) {
    os << "SOLVERLIB_CXX_FLAGS  += -D'STAR_HSUCC_HAS_RATES_G'" << std::endl;
  }

  if (do_steady_state) {
    os << "SOLVERLIB_CXX_FLAGS  += -D'STAR_DO_STEADY_STATE'" << std::endl;
  }

  if (do_rare_event) {
    os << "SOLVERLIB_CXX_FLAGS  += -D'STAR_DO_RARE_EVENT'" << std::endl;
  }

  if (do_sa) {
    os << "SOLVERLIB_CXX_FLAGS  += -D'STAR_DO_SA'" << std::endl;

    if (sa_no_der) {
      os << "SOLVERLIB_CXX_FLAGS  += -D'STAR_SA_NO_DER'" << std::endl;
    } else {
      if (sa_no_2der) {
        os << "SOLVERLIB_CXX_FLAGS  += -D'STAR_SA_NO_2DER'" << std::endl;
      }
    }
  }

  if (hs_has_index) {
    os << "SOLVERLIB_CXX_FLAGS  += -D'STAR_HS_HAS_INDEX'" << std::endl;
  }

  if (hs_has_flags) {
    os << "SOLVERLIB_CXX_FLAGS  += -D'STAR_HS_HAS_FLAGS'" << std::endl;
  }

  if (hs_has_exitrate) {
    os << "SOLVERLIB_CXX_FLAGS  += -D'STAR_HS_HAS_EXITRATE'" << std::endl;
  }

  if (!write_Makefile_defs(os)) {
    return false;
  }

  os.close();

  return true;
}
#endif

bool base::load_model() {
  if (model_name != "") {
    if (E->get_fnsi() != nullptr) {
      if (!load_fnsi_model(E->get_fnsi(), model_name)) {
        return false;
      }
    } else {
      char* buf = nullptr;
#ifdef STAR_WEB_INTERFACE
      if (!E->get_dbconnector()->load_model_script(buf, E->get_user_id(),
                                                   model_name.c_str())) {
        last_err << "couldn't open model '" << model_name << "'";
        return false;
      }
#else
      std::ifstream in(model_name.c_str(), std::ios::in | std::ios::binary);
      if (!in.is_open()) {
        last_err << "couldn't open model '" << model_name << "'";
        return false;
      }
      std::string src;
      in.seekg(0, std::ios::end);
      src.resize(in.tellg());
      in.seekg(0, std::ios::beg);
      in.read(&src[0], src.size());
      in.close();
      buf = (char*)src.c_str();
#endif

      if (!load_model_src(buf)) {
        return false;
      }
    }
  } else {
    if (sbml_model_name != "") {
#if HAVE_LIBSBML
      if (!load_sbml_model(sbml_model_name)) {
        return false;
      }
#else
      last_err << "SBML models not supported";
      return false;
#endif
    } else {
      if (custom_model_id < 0) {
        last_err << "model is not specified";
        return false;
      }

      return load_custom_model(custom_model_id);
    }
  }

  return true;
}

bool base::load_model_src(char const* const src) {
  parser::state ps;
  trsys* ts;
  if (!parser::parse_model(src, ts, ps, this)) {
    last_err << parser::get_line_number(ps) << ":" << parser::get_col_number(ps)
             << ": " << parser::get_error(ps);
    return false;
  }

  flatten(ts);
  ts = nullptr;

  return true;
}

void base::flatten(trsys const* const ts) {
  for (auto& tss : ts->get_trsyss()) {
    flatten(tss);
  }

  for (auto& c : ts->get_cnsts()) {
    add_cnst(c);
  }

  for (auto& v : ts->get_vars()) {
    add_var(v);
  }

  for (auto& f : ts->get_fncs()) {
    add_fnc(f);
  }

  for (auto& tr : ts->get_transitions()) {
    add_transition(tr);
  }

  for (auto& t : ts->get_types()) {
    add_type(t);
  }

  for (auto& ic : ts->get_ics()) {
    add_ic(ic);
  }
}

bool base::load_sbml_model(const std::string& fname) {
  trsys* ts;
  if (!parser::sbml_read(fname.c_str(), ts, last_err, this)) {
    return false;
  }

  flatten(ts);
  ts = nullptr;

  parser::write_model_description(std::cout, this);

  return true;
}

bool base::load_custom_model(const int id) {
  if (id == 0) {
  }

  last_err << "[custom_model=" << id << "]: "
           << "unknown model id";
  return false;
}

bool base::load_fnsi_model(star_fnsi const* const fnsi,
                           const std::string& model_name) {
  if (fnsi->load_model(model_name.c_str())) {
    const int Nvars = fnsi->get_variables_count();
    for (int i = 0; i < Nvars; i++) {
      static char buf[128];
      fnsi->get_variable_name(i, buf, sizeof(buf));
    }
  } else {
    static char buf[1024];
    fnsi->get_last_error(buf, sizeof(buf));
    last_err << model_name << ": " << buf;
    return false;
  }

  return true;
}

void base::add_param(model_info::cnst const* const p) const {
  auto const& pi = std::find_if(sa_params.begin(), sa_params.end(),
                                [&](model_info::cnst const* const c) {
    return c->get_name() == p->get_name();
  });

  if (pi == sa_params.end()) {
    p->set_index_p(sa_params.size());
    sa_params.push_back(p);
  }
}

bool base::add_param(const std::string& dpar,
                     model_info::cnst const*& dparam) const {
  if (dpar != "") {
    auto const& pi = std::find_if(
        sa_params.begin(), sa_params.end(),
        [&](model_info::cnst const* const c) { return c->get_name() == dpar; });

    if (pi == sa_params.end()) {
      dparam = findcnst(dpar);
      if (dparam == nullptr) {
        last_err << "unknown parameter '" << dpar << "' [2]";
        return false;
      }
      dparam->set_index_p(sa_params.size());
      sa_params.push_back(dparam);
    } else {
      dparam = *pi;
    }
  } else {
    dparam = nullptr;
  }
  return true;
}

bool base::check_dparam(const std::string& dpar, const std::string& dpar2,
                        model_info::cnst const*& dparami,
                        model_info::cnst const*& dparamj) const {
  if (!add_param(dpar, dparami)) {
    return false;
  }

  if (!add_param(dpar2, dparamj)) {
    return false;
  }

  if (dparami == nullptr) {
    dparami = dparamj;
  }

  if (dparamj != nullptr) {
    if (dparami->get_index_p() < dparamj->get_index_p()) {
      std::swap(dparami, dparamj);
    }
  }

  return true;
}

void base::add_ivar(model_info::var const* const iv) const {
  auto const& pi = std::find_if(sa_ivars.begin(), sa_ivars.end(),
                                [&](model_info::var const* const v) {
    return v->get_name() == iv->get_name();
  });

  if (pi == sa_ivars.end()) {
    iv->set_index_i(sa_ivars.size());
    sa_ivars.push_back(iv);
  }
}

void base::add_evar(model_info::var const* const ev) const {
  auto const& pi = std::find_if(sa_evars.begin(), sa_evars.end(),
                                [&](model_info::var const* const v) {
    return v->get_name() == ev->get_name();
  });

  if (pi == sa_evars.end()) {
    ev->set_index_e(sa_evars.size());
    sa_evars.push_back(ev);
  }
}

boost::filesystem::path base::get_path() const {
  return get_experiment()->get_dir() / get_dir();
}

namespace task_info {

bool base::check_var_names(solver_loader::base const* const sl,
                           const std::vector<std::string>& var_names,
                           std::vector<model_info::var const*>& vars) {
  if (!var_names.empty()) {
    vars.resize(var_names.size());
    for (std::size_t i = 0; i < vars.size(); i++) {
      model_info::var const* const v = sl->findvar(var_names[i]);
      if (v != nullptr) {
        vars[i] = v;
        name += " " + v->get_name();
      } else {
        sl->last_error() << "unknown identifier '" << var_names[i] << "'";
        return false;
      }
    }
  } else {
    vars = sl->get_vars();
  }
  return true;
}

bool cond::check_cond(solver_loader::base const* const sl, std::string& name) {
  if (cnd_str != "") {
    parser::state ps;
    if (!parser::parse_expr(cnd_str.c_str(), sl, cnd_et, ps)) {
      sl->last_error() << parser::get_error(ps);
      return false;
    }
    sl->add_expr(cnd_et, cnd_index);

    name += ", " + cnd_str;
  }
  return true;
}

bool plot::check_exprs(solver_loader::base const* const sl) {
  if (!exprs.empty()) {
    name += " {";

    parser::state ps;
    for (auto& e : exprs) {
      if (!sl->check_dparam(e.dparam, e.dparam2, e.dparami, e.dparamj)) {
        return false;
      }

      if (!parser::parse_expr(e.expr.c_str(), sl, e.et, ps)) {
        sl->last_error() << parser::get_error(ps);
        return false;
      }
      name += " " + e.expr;

      if (e.dparam != "") {
        if (!e.et->is_var()) {
          sl->last_error() << "Derivatives of expressions are not supported, "
                              "please specify a variable";
          return false;
        }

        name += " [d{" + e.dparam + "}";
        if (e.dparam2 != "") {
          name += "d{" + e.dparam2 + "}";
        }
        name += "]";
      }

      sl->add_expr(e.et, e.expr_index);
    }

    name += " }";
  }

  if (!check_cond(sl, name)) {
    return false;
  }

  return true;
}

bool plot_2d_base::check_exprs(solver_loader::base const* const sl) {
  if (!sl->check_dparam(dparam, dparam2, dparami, dparamj)) {
    return false;
  }

  parser::state ps;
  if (!parser::parse_expr(expr1.c_str(), sl, et1, ps)) {
    sl->last_error() << parser::get_error(ps);
    return false;
  }

  if (!parser::parse_expr(expr2.c_str(), sl, et2, ps)) {
    sl->last_error() << parser::get_error(ps);
    return false;
  }

  sl->add_expr(et1, expr1_index);
  sl->add_expr(et2, expr2_index);

  name += " { " + expr1 + " " + expr2;

  if (dparam != "") {
    name += " [d{" + dparam + "}";
    if (dparam2 != "") {
      name += "d{" + dparam2 + "}";
    }
    name += "]";
  }

  name += " }";

  return true;
}

bool plot_distr_2d::check(solver_loader::base const* const sl) {
  name = "plot_distr_2d";
  if (!check_exprs(sl)) {
    return false;
  }

  if (!check_cond(sl, name)) {
    return false;
  }

  if (!et1->is_var()) {
    sl->last_error() << "variable expected: '" << expr1 << "'";
    return false;
  }
  if (!et2->is_var()) {
    sl->last_error() << "variable expected: '" << expr2 << "'";
    return false;
  }

  return true;
}

bool plot_objf::check(solver_loader::base const* const sl) {
  name = "plot_objf";
  if (!check_exprs(sl)) {
    return false;
  }

  return true;
}

bool plot_objf::post_check(solver_loader::base const* const sl) {
  for (auto& e : exprs) {
    if (!sl->expr_is_param(e.et)) {
      sl->last_error() << "parameter expected: '" << e.expr << "'";
      return false;
    }
  }

  return true;
}

bool plot_objf_2d::check(solver_loader::base const* const sl) {
  name = "plot_objf_2d";
  if (!check_exprs(sl)) {
    return false;
  }

  return true;
}

bool plot_objf_2d::post_check(solver_loader::base const* const sl) {
  if (!sl->expr_is_param(et1)) {
    sl->last_error() << "parameter expected: '" << expr1 << "'";
    return false;
  }
  if (!sl->expr_is_param(et2)) {
    sl->last_error() << "parameter expected: '" << expr2 << "'";
    return false;
  }

  return true;
}
}
}
