/*
 *  trsys.hpp
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

#ifndef SOLVER_LOADER_MODEL_INFO_TRSYS_HPP_
#define SOLVER_LOADER_MODEL_INFO_TRSYS_HPP_

#include <algorithm>
#include "cnst.hpp"
#include "var.hpp"
#include "clock.hpp"
#include "fnc.hpp"
#include "transition.hpp"
#include "type.hpp"
#include "ic.hpp"

namespace solver_loader {
namespace model_info {

enum trsys_t { TS_NONE, TS_GUARDEDCMDS, TS_CHEMREACTIONS, TS_MIXED };

class trsys : public obj {
 private:
  mutable trsys_t tr_type;

  std::vector<trsys const*> trsyss;
  std::vector<cnst const*> cnsts;
  std::vector<var const*> vars;
  std::vector<clock const*> clocks;
  std::vector<fnc const*> fncs;
  std::vector<transition const*> transitions;
  std::vector<type const*> types;
  std::vector<ic const*> ics;

 protected:
  template <typename T>
  T const* findobj(const std::vector<T const*>& objs,
                   const std::string& obj_name) const {
    auto const& oi = std::find_if(
        objs.begin(), objs.end(),
        [&](T const* const o) { return o->get_name() == obj_name; });

    if (oi == objs.end()) {
      return nullptr;
    }
    return *oi;
  }

 public:
  trsys(const std::string& name, trsys const* const ts,
        trsys_t tr_type = TS_NONE)
      : obj(name, ts), tr_type(tr_type) {}

  ~trsys() {}

  void update_type(const trsys_t new_type) const {
    if (tr_type == TS_NONE) {
      tr_type = new_type;
    } else if (tr_type != TS_MIXED && tr_type != new_type) {
      tr_type = TS_MIXED;
    }

    if (ts != nullptr) {
      ts->update_type(tr_type);
    }
  }

  trsys_t get_type() const { return tr_type; }

  const std::vector<trsys const*>& get_trsyss() const { return trsyss; }

  const std::vector<cnst const*>& get_cnsts() const { return cnsts; }

  const std::vector<var const*>& get_vars() const { return vars; }

  const std::vector<clock const*>& get_clocks() const { return clocks; }

  const std::vector<fnc const*>& get_fncs() const { return fncs; }

  const std::vector<transition const*>& get_transitions() const {
    return transitions;
  }

  const std::vector<type const*>& get_types() const { return types; }

  const std::vector<ic const*>& get_ics() const { return ics; }

  void add_trsys(trsys const* const ts) { trsyss.push_back(ts); }

  void add_cnst(cnst const* const c) { cnsts.push_back(c); }

  void add_var(var const* const v) { vars.push_back(v); }

  void add_clock(clock const* const cl) { clocks.push_back(cl); }

  void add_fnc(fnc const* const f) { fncs.push_back(f); }

  void add_transition(transition const* const tr) {
    transitions.push_back(tr);

    if (tr->is_chemreaction()) {
      update_type(TS_CHEMREACTIONS);
    } else {
      update_type(TS_GUARDEDCMDS);
    }
  }

  void add_type(type const* const t) { types.push_back(t); }

  void add_ic(ic const* const i) { ics.push_back(i); }

  var const* findvar(const std::string& v) const {
    if (ts != nullptr) {
      var const* const v_ = ts->findvar(v);
      if (v_ != nullptr) {
        return v_;
      }
    }
    return findobj<var>(vars, v);
  }

  clock const* findclock(const std::string& cl) const {
    if (ts != nullptr) {
      clock const* const cl_ = ts->findclock(cl);
      if (cl_ != nullptr) {
        return cl_;
      }
    }
    return findobj<clock>(clocks, cl);
  }

  cnst const* findcnst(const std::string& c) const {
    if (ts != nullptr) {
      cnst const* const c_ = ts->findcnst(c);
      if (c_ != nullptr) {
        return c_;
      }
    }
    return findobj<cnst>(cnsts, c);
  }

  type const* findtype(const std::string& t) const {
    if (ts != nullptr) {
      type const* const t_ = ts->findtype(t);
      if (t_ != nullptr) {
        return t_;
      }
    }
    return findobj<type>(types, t);
  }

  fnc const* findfnc(const std::string& f) const {
    if (ts != nullptr) {
      fnc const* const f_ = ts->findfnc(f);
      if (f_ != nullptr) {
        return f_;
      }
    }
    return findobj<fnc>(fncs, f);
  }

  static type const* const typeBoolean;
  static type const* const typeInteger;
  static type const* const typeSpecies;
  static type const* const typeReal;
};
}
}

#endif
