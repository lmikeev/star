/*
 *  experiment.cpp
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

#include "experiment.hpp"
#include "solver_loader/task_info/all.hpp"
#include "solver_loader/all.hpp"
#include "solver/util/html_helper.hpp"
#include "solver/base.hpp"

#ifdef STAR_HAVE_MATLAB
#include "solver/matlab/libmatlabsd.hpp"
#endif

#ifdef STAR_WEB_INTERFACE

bool db_connect(mysqlconnector& dbconnector_,
                const boost::property_tree::ptree& propTree) {
  const std::string opt_host_name = propTree.get("MYSQL.opt_host_name", "");
  const std::string opt_user_name = propTree.get("MYSQL.opt_user_name", "");
  const std::string opt_password = propTree.get("MYSQL.opt_password", "");
  const unsigned int opt_port_num = propTree.get("MYSQL.opt_port_num", 0);
  const std::string opt_socket_name = propTree.get("MYSQL.opt_socket_name", "");
  const std::string opt_db_name = propTree.get("MYSQL.opt_db_name", "");
  const unsigned int opt_flags = propTree.get("MYSQL.opt_flags", 0);

  if (!dbconnector_.connect("star", opt_host_name, opt_user_name, opt_password,
                            opt_port_num, opt_socket_name, opt_db_name,
                            opt_flags)) {
    return false;
  }
  return true;
}

#endif

bool load_config(boost::property_tree::ptree& propTree,
                 char const* const config_ini) {
  namespace pt = boost::property_tree;
  try {
    read_ini(config_ini, propTree);
  } catch (pt::ini_parser_error& e) {
#ifdef STAR_WEB_INTERFACE
    BOOST_PROPERTY_TREE_THROW(
        pt::ini_parser_error(e.message(), config_ini, e.line()));
    return false;
#endif
  }
  return true;
}

experiment* get_experiment(lua_State* L) {
  lua_getglobal(L, "experiment");
  experiment* E = (experiment*)lua_topointer(L, -1);
  lua_pop(L, 1);
  return E;
}

bool chk_zero_arg(lua_State* L, const int narg) {
  if (narg != 0) {
    luaL_error(L, "No arguments expected");
    return false;
  }
  return true;
}

bool chk_one_arg(lua_State* L, const int narg) {
  if (narg != 1) {
    luaL_error(L, "One argument expected");
    return false;
  }
  return true;
}

bool chk_two_args(lua_State* L, const int narg) {
  if (narg != 2) {
    luaL_error(L, "Two arguments expected");
    return false;
  }
  return true;
}

bool chk_bool(lua_State* L, int index, const char* setting = nullptr) {
  if (!lua_isboolean(L, index)) {
    if (setting != nullptr) {
      luaL_error(L, "Boolean expected [%s]", setting);
    } else {
      luaL_error(L, "Boolean expected [%d]", index);
    }
    return false;
  }
  return true;
}

bool chk_number(lua_State* L, int index, const char* setting = nullptr) {
  if (!lua_isnumber(L, index)) {
    if (setting != nullptr) {
      luaL_error(L, "Number expected [%s]", setting);
    } else {
      luaL_error(L, "Number expected [%d]", index);
    }
    return false;
  }
  return true;
}

bool chk_str(lua_State* L, int index, const char* setting = nullptr) {
  if (!lua_isstring(L, index) || lua_isnumber(L, index)) {
    if (setting != nullptr) {
      luaL_error(L, "String expected [%s]", setting);
    } else {
      luaL_error(L, "String expected [%d]", index);
    }
    return false;
  }
  return true;
}

bool chk_table(lua_State* L, int index, const char* setting = nullptr) {
  if (!lua_istable(L, index)) {
    if (setting != nullptr) {
      luaL_error(L, "Table expected [%s]", setting);
    } else {
      luaL_error(L, "Table expected [%d]", index);
    }
    return false;
  }
  return true;
}

bool chk_dvec(lua_State* L, const int index, const char* setting = nullptr) {
  for (lua_pushnil(L); lua_next(L, -2); lua_pop(L, 1)) {
    if (!lua_isnumber(L, -2)) {
      if (setting != nullptr) {
        luaL_error(L, "Array expected [%s]", setting);
      } else {
        luaL_error(L, "Array expected [%d]", index);
      }
      return false;
    }

    if (!chk_number(L, -1)) {
      return false;
    }
  }

  return true;
}

bool chk_strvec(lua_State* L, const int index, const char* setting = nullptr) {
  for (lua_pushnil(L); lua_next(L, -2); lua_pop(L, 1)) {
    if (!lua_isnumber(L, -2)) {
      if (setting != nullptr) {
        luaL_error(L, "Array expected [%s]", setting);
      } else {
        luaL_error(L, "Array expected [%d]", index);
      }
      return false;
    }

    if (!chk_str(L, -1)) {
      return false;
    }
  }

  return true;
}

int linspace(lua_State* L) {
  const int narg = lua_gettop(L);
  if (narg < 2 || narg > 3) {
    luaL_error(L, "linspace: invalid number of arguments");
    return 0;
  }

  if (!chk_number(L, 1)) {
    return 0;
  }
  if (!chk_number(L, 2)) {
    return 0;
  }
  if (narg > 2 && !chk_number(L, 3)) {
    return 0;
  }

  const double a = lua_tonumber(L, 1);
  const double b = lua_tonumber(L, 2);
  const int n = (narg < 3) ? 10 : lua_tonumber(L, 3);

  static char tmp[64];

  const double a_min = -1000000.0;
  const double b_max = 1000000.0;
  const int n_min = 1;
  const int n_max = 1000;

  if (a >= b) {
    sprintf(tmp, "linspace: invalid arguments (%g >= %g)", a, b);
    luaL_error(L, tmp);
    return 0;
  }

  if (a < a_min) {
    sprintf(tmp, "linspace: invalid arguments (%g < %g)", a, a_min);
    luaL_error(L, tmp);
    return 0;
  }

  if (b > b_max) {
    sprintf(tmp, "linspace: invalid arguments (%g > %g)", b, b_max);
    luaL_error(L, tmp);
    return 0;
  }

  if (n < n_min || n > n_max) {
    sprintf(tmp, "logspace: invalid arguments (%d < %d || %d > %d)", n, n_min,
            n, n_max);
    luaL_error(L, tmp);
    return 0;
  }

  lua_pushlightuserdata(L, new fnc_linspace(a, b, n));

  return 1;
}

int logspace(lua_State* L) {
  const int narg = lua_gettop(L);
  if (narg < 2 || narg > 3) {
    luaL_error(L, "logspace: invalid number of arguments");
    return 0;
  }

  if (!chk_number(L, 1)) {
    return 0;
  }
  if (!chk_number(L, 2)) {
    return 0;
  }
  if (narg > 2 && !chk_number(L, 3)) {
    return 0;
  }

  const double a = lua_tonumber(L, 1);
  const double b = lua_tonumber(L, 2);
  const int n = (narg < 3) ? 10 : lua_tonumber(L, 3);

  static char tmp[64];

  const double a_min = -30.0;
  const double b_max = 6.0;
  const int n_min = 1;
  const int n_max = 1000;

  if (a >= b) {
    sprintf(tmp, "logspace: invalid arguments (%g >= %g)", a, b);
    luaL_error(L, tmp);
    return 0;
  }

  if (a < a_min) {
    sprintf(tmp, "logspace: invalid arguments (%g < %g)", a, a_min);
    luaL_error(L, tmp);
    return 0;
  }

  if (b > b_max) {
    sprintf(tmp, "logspace: invalid arguments (%g > %g)", b, b_max);
    luaL_error(L, tmp);
    return 0;
  }

  if (n < n_min || n > n_max) {
    sprintf(tmp, "logspace: invalid arguments (%d < %d || %d > %d)", n, n_min,
            n, n_max);
    luaL_error(L, tmp);
    return 0;
  }

  lua_pushlightuserdata(L, new fnc_logspace(a, b, n));

  return 1;
}

template <class T>
int distr_1_(lua_State* L, char const* const name) {
  const int narg = lua_gettop(L);

  if (!narg) {
    lua_pushlightuserdata(L, new T(0.50));
  } else if (narg == 1) {
    if (!chk_number(L, 1)) {
      return 0;
    }
    lua_pushlightuserdata(L, new T(lua_tonumber(L, 1)));
  } else {
    static char tmp[64];
    sprintf(tmp, "%s: invalid number of arguments", name);
    luaL_error(L, tmp);
    return 0;
  }

  return 1;
}

template <class T>
int distr_2_(lua_State* L, char const* const name) {
  const int narg = lua_gettop(L);

  if (!narg) {
    lua_pushlightuserdata(L, new T(0.0, 1.0));
  } else if (narg == 2) {
    if (!chk_number(L, 1)) {
      return 0;
    }
    if (!chk_number(L, 2)) {
      return 0;
    }
    lua_pushlightuserdata(L, new T(lua_tonumber(L, 1), lua_tonumber(L, 2)));
  } else {
    static char tmp[64];
    sprintf(tmp, "%s: invalid number of arguments", name);
    luaL_error(L, tmp);
    return 0;
  }

  return 1;
}

int uniform(lua_State* L) { return distr_2_<fnc_uniform>(L, "uniform"); }

int bernoulli(lua_State* L) { return distr_1_<fnc_bernoulli>(L, "bernoulli"); }

int binomial(lua_State* L) { return distr_2_<fnc_binomial>(L, "binomial"); }

int neg_binomial(lua_State* L) {
  return distr_2_<fnc_neg_binomial>(L, "neg_binomial");
}

int geometric(lua_State* L) { return distr_2_<fnc_binomial>(L, "geometric"); }

int poisson(lua_State* L) { return distr_1_<fnc_poisson>(L, "poisson"); }

int exponential(lua_State* L) {
  return distr_1_<fnc_exponential>(L, "exponential");
}

int gamma(lua_State* L) { return distr_2_<fnc_gamma>(L, "gamma"); }

int weibull(lua_State* L) { return distr_2_<fnc_weibull>(L, "weibull"); }

int normal(lua_State* L) { return distr_2_<fnc_normal>(L, "normal"); }

int lognormal(lua_State* L) { return distr_2_<fnc_lognormal>(L, "lognormal"); }

bool load_number(lua_State* L, const int index, int& i) {
  i = lua_tonumber(L, index);
  return true;
}

bool load_number(lua_State* L, const int index, unsigned int& ui) {
  ui = lua_tonumber(L, index);
  return true;
}

bool load_number(lua_State* L, const int index, double& d) {
  d = lua_tonumber(L, index);
  return true;
}

bool load_number(lua_State* L, const int index, const char* key) {
  experiment_options* const eo = get_experiment(L)->get_options();
  eo->set(key, lua_tonumber(L, index));
  return true;
}

bool load_bool(lua_State* L, const int index, bool& b) {
  b = lua_toboolean(L, index);
  return true;
}

bool load_bool(lua_State* L, const int index, const char* key) {
  experiment_options* const eo = get_experiment(L)->get_options();
  eo->set(key, lua_toboolean(L, index));
  return true;
}

bool load_str(lua_State* L, const int index, std::string& str) {
  str = lua_tostring(L, index);
  return true;
}

bool load_str(lua_State* L, const int index, const char* key) {
  experiment_options* const eo = get_experiment(L)->get_options();
  eo->set(key, lua_tostring(L, index));
  return true;
}

int set_str(lua_State* L, std::string& str) {
  const int narg = lua_gettop(L);
  if (!chk_one_arg(L, narg)) {
    return 0;
  }
  if (!chk_str(L, 1)) {
    return 0;
  }
  if (!load_str(L, 1, str)) {
    return 0;
  }
  return 0;
}

bool load_dvec(lua_State* L, const int, std::vector<double>& dvec) {
  dvec.clear();
  for (lua_pushnil(L); lua_next(L, -2); lua_pop(L, 1)) {
    if (!lua_isnumber(L, -2)) {
      luaL_error(L, "Array expected");
      return false;
    }

    if (!chk_number(L, -1)) {
      return false;
    }

    dvec.push_back(lua_tonumber(L, -1));
  }

  return true;
}

int set_dvec(lua_State* L, std::vector<double>& dvec) {
  const int narg = lua_gettop(L);
  if (!chk_one_arg(L, narg)) {
    return 0;
  }
  if (!chk_table(L, 1)) {
    return 0;
  }
  if (!load_dvec(L, 1, dvec)) {
    return 0;
  }
  return 0;
}

void stack_dump(lua_State* L) {
  int i = lua_gettop(L);
  printf(" ----------------  Stack Dump ----------------\n");
  while (i) {
    int t = lua_type(L, i);
    switch (t) {
      case LUA_TSTRING:
        printf("%d:`%s'\n", i, lua_tostring(L, i));
        break;

      case LUA_TBOOLEAN:
        printf("%d: %s\n", i, lua_toboolean(L, i) ? "true" : "false");
        break;

      case LUA_TNUMBER:
        printf("%d: %g\n", i, lua_tonumber(L, i));
        break;

      default:
        printf("%d: %s\n", i, lua_typename(L, t));
        break;
    }
    i--;
  }
  printf("--------------- Stack Dump Finished ---------------\n");
}

bool load_strvec(lua_State* L, const int index,
                 std::vector<std::string>& strvec) {
  strvec.clear();
  for (lua_pushnil(L); lua_next(L, index - 1); lua_pop(L, 1)) {
    if (!lua_isnumber(L, -2)) {
      luaL_error(L, "Array expected");
      return false;
    }

    if (!chk_str(L, -1)) {
      return false;
    }

    strvec.push_back(lua_tostring(L, -1));
  }

  return true;
}

bool load_strset(lua_State* L, const int index, std::set<std::string>& strset) {
  strset.clear();
  for (lua_pushnil(L); lua_next(L, index - 1); lua_pop(L, 1)) {
    if (!lua_isnumber(L, -2)) {
      luaL_error(L, "Array expected");
      return false;
    }

    if (!chk_str(L, -1)) {
      return false;
    }

    strset.insert(lua_tostring(L, -1));
  }
  return true;
}

int set_strvec(lua_State* L, std::vector<std::string>& strvec) {
  const int narg = lua_gettop(L);
  if (!chk_one_arg(L, narg)) {
    return 0;
  }
  if (!chk_table(L, 1)) {
    return 0;
  }
  if (!load_strvec(L, 1, strvec)) {
    return 0;
  }
  return 0;
}

bool load_se_params(lua_State* L, const int index,
                    std::vector<set_param_info const*>& setparams,
                    bool& estparam, std::vector<est_param_info>& estparams,
                    std::vector<sample_param_info>& sampleparams,
                    std::vector<iterate_param_info>& iterateparams,

                    const bool = false, const bool = false) {
  for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
       lua_pop(L, 1)) {
    std::vector<std::string> params;

    if (lua_istable(L, -2)) {
      if (!load_strvec(L, -2, params)) {
        return false;
      }
    } else {
      if (!chk_str(L, -2)) {
        return false;
      }

      params.push_back(lua_tostring(L, -2));
    }

    if (lua_isboolean(L, -1)) {
      const bool value = lua_toboolean(L, -1);
      for (auto const& p : params) {
        setparams.push_back(new set_param_info_bool(p, value));
      }
    }

    else if (lua_isnumber(L, -1)) {
      const double value = lua_tonumber(L, -1);
      for (auto const& p : params) {
        setparams.push_back(new set_param_info_double(p, value));
      }
    } else if (lua_isstring(L, -1)) {
      const char* value = lua_tostring(L, -1);
      for (auto const& p : params) {
        setparams.push_back(new set_param_info_str(p, value));
      }
    } else if (lua_istable(L, -1)) {
      bool isvec = true;

      for (lua_pushnil(L); lua_next(L, -2); lua_pop(L, 1)) {
        if (!lua_isnumber(L, -2) || !lua_isnumber(L, -1)) {
          isvec = false;
        }
      }

      if (isvec) {
        iterate_param_info ip;
        for (lua_pushnil(L); lua_next(L, -2); lua_pop(L, 1)) {
          ip.values.push_back(lua_tonumber(L, -1));
        }

        for (auto const& p : params) {
          ip.param = p;
          iterateparams.push_back(ip);
        }
      } else {
        double smin = 1e-5, smax = 1000, sinit = 0.50 * (smin + smax);
        bool fmin = false, fmax = false, finit = false;

        for (lua_pushnil(L); lua_next(L, -2); lua_pop(L, 1)) {
          if (!chk_str(L, -2)) {
            return false;
          }

          const char* sname = lua_tostring(L, -2);

          if (!strcmp(sname, "min")) {
            if (!chk_number(L, -1)) {
              return false;
            }
            smin = lua_tonumber(L, -1);
            fmin = true;
          } else if (!strcmp(sname, "max")) {
            if (!chk_number(L, -1)) {
              return false;
            }
            smax = lua_tonumber(L, -1);
            fmax = true;
          } else if (!strcmp(sname, "init")) {
            if (!chk_number(L, -1)) {
              return false;
            }
            sinit = lua_tonumber(L, -1);
            finit = true;
          } else {
            luaL_error(L, "Invalid setting");
            return false;
          }
        }

        if (!fmin) {
          smin = 1e-5;
        }
        if (!fmax) {
          smax = 1000.0;
        }
        if (!finit) {
          sinit = 0.50 * (smin + smax);
        }

        for (auto const& p : params) {
          est_param_info ep;
          ep.param = p;
          ep.min = smin;
          ep.max = smax;
          ep.value0 = sinit;
          estparams.push_back(ep);
        }
        estparam = true;
      }
    } else if (lua_islightuserdata(L, index)) {
      iterate_param_info ip;

      fnc* f = static_cast<fnc*>(lua_touserdata(L, index));
      if (f->get_values(ip.values)) {
        for (auto const& p : params) {
          ip.param = p;
          iterateparams.push_back(ip);
        }
      } else if (f->is_distr()) {
        for (auto const& p : params) {
          sample_param_info sp;
          sp.param = p;
          sp.f = f;
          sampleparams.push_back(sp);
        }
      } else {
        luaL_error(L, "Invalid function");
        return false;
      }
    } else {
      luaL_error(L, "Number or table expected");
      return false;
    }
  }

  return true;
}

bool set_se_params(lua_State* L, std::vector<set_param_info const*>& setparams,
                   bool& estparam, std::vector<est_param_info>& estparams,
                   std::vector<sample_param_info>& sampleparams,
                   std::vector<iterate_param_info>& iterateparams,
                   const bool clear = true, const bool ignoreZeros = false) {
  const int narg = lua_gettop(L);
  if (!chk_one_arg(L, narg)) {
    return 0;
  }
  if (!chk_table(L, 1)) {
    return false;
  }
  if (!load_se_params(L, 1, setparams, estparam, estparams, sampleparams,
                      iterateparams, clear, ignoreZeros)) {
    return 0;
  }
  return 0;
}

#ifndef WIN32
#define RGB(r, g, b)                                  \
  ((uint32_t)(((uint8_t)(r) | ((uint16_t)(g) << 8)) | \
              (((uint32_t)(uint8_t)(b)) << 16)))
#endif

bool load_color(lua_State* L, const int, unsigned int& color) {
  if (lua_isboolean(L, -1)) {
    const bool value = lua_toboolean(L, -1);
    if (value) {
      color = RGB(255, 255, 255);
    } else {
      color = RGB(0, 0, 0);
    }
  } else if (lua_isnumber(L, -1)) {
    const double value = lua_tonumber(L, -1);
    if (value < 0.0) {
      luaL_error(L, "Invalid color");
      return false;
    }
    double i_;
    if (modf(value, &i_) < std::numeric_limits<double>::epsilon()) {
      color = static_cast<unsigned int>(value);
    } else {
      if (value > 1.0) {
        luaL_error(L, "Invalid color");
        return false;
      }
      color = static_cast<unsigned int>(
          static_cast<double>(RGB(255, 255, 255)) * value);
    }
  } else if (lua_isstring(L, -1)) {
    const char* value = lua_tostring(L, -1);
    if (!strcmp(value, "black")) {
      color = RGB(0, 0, 0);
    } else if (!strcmp(value, "blue")) {
      color = RGB(0, 0, 255);
    } else if (!strcmp(value, "green")) {
      color = RGB(0, 128, 0);
    } else if (!strcmp(value, "red")) {
      color = RGB(0, 128, 0);
    } else if (!strcmp(value, "white")) {
      color = RGB(255, 255, 255);
    } else {
    }
  } else if (lua_istable(L, -1)) {
  }
  return true;
}

bool load_init(lua_State* L, const int index) {
  experiment_options* const eo = get_experiment(L)->get_options();

  if (lua_istable(L, index)) {
    if (!load_se_params(L, index, eo->set_inits, eo->estimate_init,
                        eo->estimate_inits, eo->sample_inits,
                        eo->iterate_inits)) {
      return false;
    }
    eo->init_src = "";
  } else {
    if (!chk_str(L, index, "init")) {
      return false;
    }
    if (!load_str(L, index, eo->init_src)) {
      return false;
    }
  }
  return true;
}

bool load_tspan(lua_State* L, const int index) {
  experiment_options* const eo = get_experiment(L)->get_options();

  if (lua_istable(L, index)) {
    if (!load_dvec(L, index, eo->tspan)) {
      return false;
    }
  } else if (lua_islightuserdata(L, index)) {
    eo->tspan.clear();

    fnc* fn_ = static_cast<fnc*>(lua_touserdata(L, index));
    if (!fn_->get_values(eo->tspan)) {
      luaL_error(L, "Invalid function");
      return false;
    }
  } else {
    luaL_error(L, "Invalid format");
    return false;
  }
  return true;
}

bool load_dump_moments_(lua_State* L, const int index) {
  experiment_options* const eo = get_experiment(L)->get_options();

  std::vector<std::string> var_names;
  int nmoments_ = -1;
  bool raw = false;
  std::string cond = "";
  std::string cmp_path = "";

  bool ists = true;
  bool ists_ = false;
  bool isb = false;
  bool isbv = true;

  for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
       lua_pop(L, 1)) {
    if (lua_isnumber(L, -2)) {
      if (!ists || isb) {
        luaL_error(L, "dump_moments: invalid format");
        return false;
      }

      if (lua_isboolean(L, -1)) {
        isb = true;
        isbv = lua_toboolean(L, -1);
      } else {
        if (!chk_str(L, -1)) {
          return false;
        }

        var_names.push_back(lua_tostring(L, -1));

        ists_ = true;
      }
    } else if (lua_isstring(L, -2)) {
      if (ists_ || isb) {
        luaL_error(L, "dump_moments: invalid format");
        return false;
      }

      ists = false;

      const char* key = lua_tostring(L, -2);

      if (strcmp(key, "vars") == 0) {
        if (lua_istable(L, -1)) {
          for (lua_pushnil(L); lua_next(L, -2); lua_pop(L, 1)) {
            if (lua_isnumber(L, -2)) {
              if (!chk_str(L, -1)) {
                return false;
              }

              var_names.push_back(lua_tostring(L, -1));
            } else if (lua_isstring(L, -2)) {
              var_names.push_back(lua_tostring(L, -2));
            } else {
              luaL_error(L, "dump_moments: variable name expected");
              return false;
            }
          }
        } else {
          if (!chk_str(L, -1)) {
            return false;
          }

          var_names.push_back(lua_tostring(L, -1));
        }
      } else if (strcmp(key, "cond") == 0) {
        if (!chk_str(L, -1)) {
          return false;
        }
        if (!load_str(L, -1, cond)) {
          return false;
        }
      } else if (strcmp(key, "raw") == 0) {
        if (!chk_bool(L, -1)) {
          return false;
        }
        if (!load_bool(L, -1, raw)) {
          return false;
        }
      } else if (strcmp(key, "nmoments") == 0) {
        if (!chk_number(L, -1)) {
          return false;
        }
        if (!load_number(L, -1, nmoments_)) {
          return false;
        }
      } else if (strcmp(key, "cmp_path") == 0) {
        if (!chk_str(L, -1)) {
          return false;
        }
        if (!load_str(L, -1, cmp_path)) {
          return false;
        }
      }
    } else {
      luaL_error(L, "dump_moments: invalid format");
      return false;
    }
  }

  if (isbv) {
    eo->tasks.push_back(new solver_loader::task_info::dump_moments(
        var_names, nmoments_, raw, cond, cmp_path));
  }

  return true;
}

bool load_dump_moments(lua_State* L, const int index) {
  experiment_options* const eo = get_experiment(L)->get_options();

  if (lua_istable(L, index)) {
    bool isempty = true;
    bool istt = true;

    for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
         lua_pop(L, 1)) {
      isempty = false;
      if (!lua_istable(L, -1)) {
        istt = false;
      }
    }

    if (isempty) {
      eo->tasks.push_back(new solver_loader::task_info::dump_moments());
    } else if (istt) {
      for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
           lua_pop(L, 1)) {
        if (!load_dump_moments_(L, -1)) {
          return false;
        }
      }
    } else {
      if (!load_dump_moments_(L, index)) {
        return false;
      }
    }
  } else if (lua_isboolean(L, index)) {
    if (lua_toboolean(L, index)) {
      eo->tasks.push_back(new solver_loader::task_info::dump_moments());
    }
  } else if (lua_isstring(L, index)) {
    std::vector<std::string> var_names;
    var_names.push_back(lua_tostring(L, index));
    eo->tasks.push_back(new solver_loader::task_info::dump_moments(var_names));
  } else {
    luaL_error(L, "dump_moments: invalid format");
    return false;
  }

  return true;
}

bool load_dump_distr_(lua_State* L, const int index) {
  experiment_options* const eo = get_experiment(L)->get_options();

  std::vector<std::string> var_names;
  std::string cond = "";
  std::string cmp_path = "";

  bool ists = true;
  bool ists_ = false;
  bool isb = false;
  bool isbv = true;

  for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
       lua_pop(L, 1)) {
    if (lua_isnumber(L, -2)) {
      if (!ists || isb) {
        luaL_error(L, "dump_distr: invalid format");
        return false;
      }

      if (lua_isboolean(L, -1)) {
        isb = true;
        isbv = lua_toboolean(L, -1);
      } else {
        if (!chk_str(L, -1)) {
          return false;
        }

        var_names.push_back(lua_tostring(L, -1));

        ists_ = true;
      }
    } else if (lua_isstring(L, -2)) {
      if (ists_ || isb) {
        luaL_error(L, "dump_distr: invalid format");
        return false;
      }

      ists = false;

      const char* key = lua_tostring(L, -2);

      if (strcmp(key, "vars") == 0) {
        if (lua_istable(L, -1)) {
          for (lua_pushnil(L); lua_next(L, -2); lua_pop(L, 1)) {
            if (lua_isnumber(L, -2)) {
              if (!chk_str(L, -1)) {
                return false;
              }

              var_names.push_back(lua_tostring(L, -1));
            } else if (lua_isstring(L, -2)) {
              var_names.push_back(lua_tostring(L, -2));
            } else {
              luaL_error(L, "dump_distr: variable name expected");
              return false;
            }
          }
        } else {
          if (!chk_str(L, -1)) {
            return false;
          }

          var_names.push_back(lua_tostring(L, -1));
        }
      } else if (strcmp(key, "cond") == 0) {
        if (!chk_str(L, -1)) {
          return false;
        }
        if (!load_str(L, -1, cond)) {
          return false;
        }
      } else if (strcmp(key, "cmp_path") == 0) {
        if (!chk_str(L, -1)) {
          return false;
        }
        if (!load_str(L, -1, cmp_path)) {
          return false;
        }
      }
    } else {
      luaL_error(L, "dump_distr: invalid format");
      return false;
    }
  }

  if (isbv) {
    eo->tasks.push_back(
        new solver_loader::task_info::dump_distr(var_names, cond, cmp_path));
  }

  return true;
}

bool load_dump_distr(lua_State* L, const int index) {
  experiment_options* const eo = get_experiment(L)->get_options();

  if (lua_istable(L, index)) {
    bool isempty = true;
    bool istt = true;

    for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
         lua_pop(L, 1)) {
      isempty = false;
      if (!lua_istable(L, -1)) {
        istt = false;
      }
    }

    if (isempty) {
      eo->tasks.push_back(new solver_loader::task_info::dump_distr());
    } else if (istt) {
      for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
           lua_pop(L, 1)) {
        if (!load_dump_distr_(L, -1)) {
          return false;
        }
      }
    } else {
      if (!load_dump_distr_(L, index)) {
        return false;
      }
    }
  } else if (lua_isboolean(L, index)) {
    if (lua_toboolean(L, index)) {
      eo->tasks.push_back(new solver_loader::task_info::dump_distr());
    }
  } else if (lua_isstring(L, index)) {
    std::vector<std::string> var_names;
    var_names.push_back(lua_tostring(L, index));
    eo->tasks.push_back(new solver_loader::task_info::dump_distr(var_names));
  } else {
    luaL_error(L, "dump_distr: invalid format");
    return false;
  }

  return true;
}

bool load_plot_(lua_State* L, const int index, const char* = "") {
  experiment_options* const eo = get_experiment(L)->get_options();

  solver_loader::task_info::plot_props props;
  std::vector<solver_loader::task_info::plot_expr> exprs;
  std::string dparam = "";
  std::string dparam2 = "";
  bool chm = false;
  std::string cond = "";
  bool plot_stddev = false;
  std::string copy_data_path = "";

  bool ists = true;
  bool ists_ = false;
  bool isb = false;
  bool isbv = true;

  for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
       lua_pop(L, 1)) {
    if (lua_isnumber(L, -2)) {
      if (!ists || isb) {
        luaL_error(L, "plot: invalid format");
        return false;
      }

      if (lua_isboolean(L, -1)) {
        isb = true;
        isbv = lua_toboolean(L, -1);
      } else {
        if (!chk_str(L, -1)) {
          return false;
        }

        solver_loader::task_info::plot_expr e;
        e.expr = lua_tostring(L, -1);
        e.dparam = "";
        e.dparam2 = "";
        e.chm = false;
        exprs.push_back(e);

        ists_ = true;
      }
    } else if (lua_isstring(L, -2)) {
      if (ists_ || isb) {
        luaL_error(L, "plot: invalid format");
        return false;
      }

      ists = false;

      const char* key = lua_tostring(L, -2);

      if (strcmp(key, "expr") == 0) {
        if (lua_istable(L, -1)) {
          for (lua_pushnil(L); lua_next(L, -2); lua_pop(L, 1)) {
            if (lua_isnumber(L, -2)) {
              if (!chk_str(L, -1)) {
                return false;
              }

              solver_loader::task_info::plot_expr e;
              e.expr = lua_tostring(L, -1);
              e.dparam = "";
              e.dparam2 = "";
              e.chm = false;
              exprs.push_back(e);
            } else if (lua_isstring(L, -2)) {
              solver_loader::task_info::plot_expr e;
              e.expr = lua_tostring(L, -2);
              e.dparam = "";
              e.dparam2 = "";
              e.chm = false;
              exprs.push_back(e);
            } else {
              luaL_error(L, "plot: expression expected");
              return false;
            }
          }
        } else {
          if (!chk_str(L, -1)) {
            return false;
          }

          solver_loader::task_info::plot_expr e;
          e.expr = lua_tostring(L, -1);
          e.dparam = "";
          e.dparam2 = "";
          e.chm = false;
          exprs.push_back(e);
        }
      } else if (strcmp(key, "title") == 0) {
        if (!chk_str(L, -1)) {
          return false;
        }
        std::string title;
        if (!load_str(L, -1, title)) {
          return false;
        }
        props.set_title(title);
      } else if (strcmp(key, "dparam") == 0) {
        if (!chk_str(L, -1)) {
          return false;
        }
        if (!load_str(L, -1, dparam)) {
          return false;
        }
      } else if (strcmp(key, "dparam2") == 0) {
        if (!chk_str(L, -1)) {
          return false;
        }
        if (!load_str(L, -1, dparam2)) {
          return false;
        }
      } else if (strcmp(key, "chm") == 0) {
        if (!chk_bool(L, -1)) {
          return false;
        }
        if (!load_bool(L, -1, chm)) {
          return false;
        }
      } else if (strcmp(key, "plot_stddev") == 0) {
        if (!chk_bool(L, -1)) {
          return false;
        }
        if (!load_bool(L, -1, plot_stddev)) {
          return false;
        }
      } else if (strcmp(key, "copy_data_path") == 0) {
        if (!chk_str(L, -1)) {
          return false;
        }
        if (!load_str(L, -1, copy_data_path)) {
          return false;
        }
      }
    } else {
      luaL_error(L, "plot: invalid format");
      return false;
    }
  }

  if (!ists) {
    for (auto& e : exprs) {
      e.dparam = dparam;
      e.dparam2 = dparam2;
      e.chm = chm;
    }
  }

  if (isbv) {
    eo->tasks.push_back(new solver_loader::task_info::plot(
        props, exprs, cond, copy_data_path, plot_stddev));
  }

  return true;
}

bool load_plot(lua_State* L, const int index) {
  experiment_options* const eo = get_experiment(L)->get_options();

  if (lua_istable(L, index)) {
    bool isempty = true;
    bool istt = true;

    for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
         lua_pop(L, 1)) {
      isempty = false;
      if (!lua_istable(L, -1)) {
        istt = false;
      }
    }

    if (isempty) {
      solver_loader::task_info::plot_props props;
      std::vector<solver_loader::task_info::plot_expr> exprs;
      eo->tasks.push_back(new solver_loader::task_info::plot(props, exprs));
    } else if (istt) {
      for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
           lua_pop(L, 1)) {
        const char* title_ = (lua_isstring(L, -2) && !lua_isnumber(L, -2))
                                 ? lua_tostring(L, -2)
                                 : "";
        if (!load_plot_(L, -1, title_)) {
          return false;
        }
      }
    } else {
      if (!load_plot_(L, index)) {
        return false;
      }
    }
  } else if (lua_isboolean(L, index)) {
    if (lua_toboolean(L, index)) {
      solver_loader::task_info::plot_props props;
      std::vector<solver_loader::task_info::plot_expr> exprs;
      eo->tasks.push_back(new solver_loader::task_info::plot(props, exprs));
    }
  } else if (lua_isstring(L, index)) {
    solver_loader::task_info::plot_props props;
    std::vector<solver_loader::task_info::plot_expr> exprs;
    solver_loader::task_info::plot_expr e;
    e.expr = lua_tostring(L, index);
    exprs.push_back(e);
    eo->tasks.push_back(new solver_loader::task_info::plot(props, exprs));
  } else {
    luaL_error(L, "plot: invalid format");
    return false;
  }

  return true;
}

bool load_plot_2d_(lua_State* L, const int index, const char* = "") {
  experiment_options* const eo = get_experiment(L)->get_options();

  solver_loader::task_info::plot_props props;
  solver_loader::task_info::plot_line_props line_props;
  std::string expr1;
  std::string expr2;
  std::string dparam = "";
  std::string dparam2 = "";

  for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
       lua_pop(L, 1)) {
    if (!chk_str(L, -2)) {
      return false;
    }

    const char* key = lua_tostring(L, -2);

    if (strcmp(key, "expr1") == 0) {
      if (!chk_str(L, -1)) {
        return false;
      }
      if (!load_str(L, -1, expr1)) {
        return false;
      }
    } else if (strcmp(key, "expr2") == 0) {
      if (!chk_str(L, -1)) {
        return false;
      }
      if (!load_str(L, -1, expr2)) {
        return false;
      }
    } else if (strcmp(key, "dparam") == 0) {
      if (!chk_str(L, -1)) {
        return false;
      }
      if (!load_str(L, -1, dparam)) {
        return false;
      }
    } else if (strcmp(key, "dparam2") == 0) {
      if (!chk_str(L, -1)) {
        return false;
      }
      if (!load_str(L, -1, dparam2)) {
        return false;
      }
    } else if (strcmp(key, "title") == 0) {
      if (!chk_str(L, -1)) {
        return false;
      }
      std::string title;
      if (!load_str(L, -1, title)) {
        return false;
      }
      props.set_title(title);
    } else {
    }
  }

  if (expr1 == "") {
    luaL_error(L, "plot_2d: 'expr1' expected");
    return false;
  }

  if (expr2 == "") {
    luaL_error(L, "plot_2d: 'expr2' expected");
    return false;
  }

  eo->tasks.push_back(new solver_loader::task_info::plot_2d(
      props, expr1, expr2, line_props, dparam, dparam2));

  return true;
}

bool load_plot_2d(lua_State* L, const int index) {
  if (lua_istable(L, index)) {
    bool istt = true;
    for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
         lua_pop(L, 1)) {
      if (!lua_istable(L, -1)) {
        istt = false;
      }
    }

    if (istt) {
      for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
           lua_pop(L, 1)) {
        const char* title_ = (lua_isstring(L, -2) && !lua_isnumber(L, -2))
                                 ? lua_tostring(L, -2)
                                 : "";
        if (!load_plot_2d_(L, -1, title_)) {
          return false;
        }
      }
    } else {
      if (!load_plot_2d_(L, index)) {
        return false;
      }
    }
  } else {
    luaL_error(L, "plot_2d: invalid format");
    return false;
  }

  return true;
}

bool load_plot_distr_(lua_State* L, const int index, const char* = "") {
  experiment_options* const eo = get_experiment(L)->get_options();

  solver_loader::task_info::plot_props props;
  std::vector<solver_loader::task_info::plot_expr> exprs;
  std::string cnd = "";
  std::string cmp_path = "";
  std::string dparam = "";
  std::string dparam2 = "";
  bool chm = false;

  bool ists = true;
  bool ists_ = false;
  bool isb = false;
  bool isbv = true;

  for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
       lua_pop(L, 1)) {
    if (lua_isnumber(L, -2)) {
      if (!ists || isb) {
        luaL_error(L, "plot_distr: invalid format");
        return false;
      }

      if (lua_isboolean(L, -1)) {
        isb = true;
        isbv = lua_toboolean(L, -1);
      } else {
        if (!chk_str(L, -1)) {
          return false;
        }

        solver_loader::task_info::plot_expr e;
        e.expr = lua_tostring(L, -1);
        e.dparam = "";
        e.dparam2 = "";
        e.chm = false;
        exprs.push_back(e);

        ists_ = true;
      }
    } else if (lua_isstring(L, -2)) {
      if (ists_ || isb) {
        luaL_error(L, "plot_distr: invalid format");
        return false;
      }

      ists = false;

      const char* key = lua_tostring(L, -2);

      if (strcmp(key, "vars") == 0) {
        if (lua_istable(L, -1)) {
          for (lua_pushnil(L); lua_next(L, -2); lua_pop(L, 1)) {
            if (lua_isnumber(L, -2)) {
              if (!chk_str(L, -1)) {
                return false;
              }

              solver_loader::task_info::plot_expr e;
              e.expr = lua_tostring(L, -1);
              e.dparam = "";
              e.dparam2 = "";
              e.chm = false;
              exprs.push_back(e);
            } else if (lua_isstring(L, -2)) {
              solver_loader::task_info::plot_expr e;
              e.expr = lua_tostring(L, -2);
              e.dparam = "";
              e.dparam2 = "";
              e.chm = false;
              exprs.push_back(e);
            } else {
              luaL_error(L, "plot_distr: variable name expected");
              return false;
            }
          }
        } else {
          if (!chk_str(L, -1)) {
            return false;
          }

          solver_loader::task_info::plot_expr e;
          e.expr = lua_tostring(L, -1);
          e.dparam = "";
          e.dparam2 = "";
          e.chm = false;
          exprs.push_back(e);
        }
      } else if (strcmp(key, "title") == 0) {
        if (!chk_str(L, -1)) {
          return false;
        }
        std::string title;
        if (!load_str(L, -1, title)) {
          return false;
        }
        props.set_title(title);
      } else if (strcmp(key, "cond") == 0) {
        if (!chk_str(L, -1)) {
          return false;
        }
        if (!load_str(L, -1, cnd)) {
          return false;
        }
      } else if (strcmp(key, "cmp_path") == 0) {
        if (!chk_str(L, -1)) {
          return false;
        }
        if (!load_str(L, -1, cmp_path)) {
          return false;
        }
      } else if (strcmp(key, "dparam") == 0) {
        if (!chk_str(L, -1)) {
          return false;
        }
        if (!load_str(L, -1, dparam)) {
          return false;
        }
      } else if (strcmp(key, "dparam2") == 0) {
        if (!chk_str(L, -1)) {
          return false;
        }
        if (!load_str(L, -1, dparam2)) {
          return false;
        }
      } else if (strcmp(key, "chm") == 0) {
        if (!chk_bool(L, -1)) {
          return false;
        }
        if (!load_bool(L, -1, chm)) {
          return false;
        }
      }
    } else {
      luaL_error(L, "plot_distr: invalid format");
      return false;
    }
  }

  if (!ists) {
    for (auto& e : exprs) {
      e.dparam = dparam;
      e.dparam2 = dparam2;
      e.chm = chm;
    }
  }

  if (isbv) {
    eo->tasks.push_back(
        new solver_loader::task_info::plot_distr(props, exprs, cnd, cmp_path));
  }

  return true;
}

bool load_plot_distr(lua_State* L, const int index) {
  experiment_options* const eo = get_experiment(L)->get_options();

  if (lua_istable(L, index)) {
    bool isempty = true;
    bool istt = true;

    for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
         lua_pop(L, 1)) {
      isempty = false;
      if (!lua_istable(L, -1)) {
        istt = false;
      }
    }

    if (isempty) {
      solver_loader::task_info::plot_props props;
      std::vector<solver_loader::task_info::plot_expr> exprs;
      eo->tasks.push_back(
          new solver_loader::task_info::plot_distr(props, exprs));
    } else if (istt) {
      for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
           lua_pop(L, 1)) {
        const char* title_ = (lua_isstring(L, -2) && !lua_isnumber(L, -2))
                                 ? lua_tostring(L, -2)
                                 : "";
        if (!load_plot_distr_(L, -1, title_)) {
          return false;
        }
      }
    } else {
      if (!load_plot_distr_(L, index)) {
        return false;
      }
    }
  } else if (lua_isboolean(L, index)) {
    if (lua_toboolean(L, index)) {
      solver_loader::task_info::plot_props props;
      std::vector<solver_loader::task_info::plot_expr> exprs;
      eo->tasks.push_back(
          new solver_loader::task_info::plot_distr(props, exprs));
    }
  } else if (lua_isstring(L, index)) {
    solver_loader::task_info::plot_props props;
    std::vector<solver_loader::task_info::plot_expr> exprs;
    solver_loader::task_info::plot_expr e;
    e.expr = lua_tostring(L, index);
    exprs.push_back(e);
    eo->tasks.push_back(new solver_loader::task_info::plot_distr(props, exprs));
  } else {
    luaL_error(L, "plot_distr: invalid format");
    return false;
  }

  return true;
}

bool load_plot_distr_2d_(lua_State* L, const int index, const char* = "") {
  experiment_options* const eo = get_experiment(L)->get_options();

  solver_loader::task_info::plot_props props;
  solver_loader::task_info::plot_surface_props surface_props;
  std::string expr1;
  std::string expr2;
  std::string cnd = "";
  std::string dparam = "";
  std::string dparam2 = "";

  for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
       lua_pop(L, 1)) {
    if (!chk_str(L, -2)) {
      return false;
    }

    const char* key = lua_tostring(L, -2);

    if (strcmp(key, "var1") == 0) {
      if (!chk_str(L, -1)) {
        return false;
      }
      if (!load_str(L, -1, expr1)) {
        return false;
      }
    } else if (strcmp(key, "var2") == 0) {
      if (!chk_str(L, -1)) {
        return false;
      }
      if (!load_str(L, -1, expr2)) {
        return false;
      }
    } else if (strcmp(key, "cond") == 0) {
      if (!chk_str(L, -1)) {
        return false;
      }
      if (!load_str(L, -1, cnd)) {
        return false;
      }
    } else if (strcmp(key, "dparam") == 0) {
      if (!chk_str(L, -1)) {
        return false;
      }
      if (!load_str(L, -1, dparam)) {
        return false;
      }
    } else if (strcmp(key, "dparam2") == 0) {
      if (!chk_str(L, -1)) {
        return false;
      }
      if (!load_str(L, -1, dparam2)) {
        return false;
      }
    } else if (strcmp(key, "title") == 0) {
      if (!chk_str(L, -1)) {
        return false;
      }
      std::string title;
      if (!load_str(L, -1, title)) {
        return false;
      }
      props.set_title(title);
    }
  }

  if (expr1 == "") {
    luaL_error(L, "plot_distr_2d: 'var1' expected");
    return false;
  }

  if (expr2 == "") {
    luaL_error(L, "plot_distr_2d: 'var2' expected");
    return false;
  }

  eo->tasks.push_back(new solver_loader::task_info::plot_distr_2d(
      props, expr1, expr2, surface_props, cnd, dparam, dparam2));

  return true;
}

bool load_plot_distr_2d(lua_State* L, const int index) {
  if (lua_istable(L, index)) {
    bool istt = true;
    for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
         lua_pop(L, 1)) {
      if (!lua_istable(L, -1)) {
        istt = false;
      }
    }

    if (istt) {
      for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
           lua_pop(L, 1)) {
        const char* title_ = (lua_isstring(L, -2) && !lua_isnumber(L, -2))
                                 ? lua_tostring(L, -2)
                                 : "";
        if (!load_plot_distr_2d_(L, -1, title_)) {
          return false;
        }
      }
    } else {
      if (!load_plot_distr_2d_(L, index)) {
        return false;
      }
    }
  } else {
    luaL_error(L, "plot_distr_2d: invalid format");
    return false;
  }

  return true;
}

bool load_plot_objf_(lua_State* L, const int index, const char* = "") {
  experiment_options* const eo = get_experiment(L)->get_options();

  solver_loader::task_info::plot_props props;
  std::vector<solver_loader::task_info::plot_expr> exprs;

  bool ists = true;
  bool ists_ = false;
  bool isb = false;
  bool isbv = true;

  for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
       lua_pop(L, 1)) {
    if (lua_isnumber(L, -2)) {
      if (!ists || isb) {
        luaL_error(L, "plot_objf: invalid format");
        return false;
      }

      if (lua_isboolean(L, -1)) {
        isb = true;
        isbv = lua_toboolean(L, -1);
      } else {
        if (!chk_str(L, -1)) {
          return false;
        }

        solver_loader::task_info::plot_expr e;
        e.expr = lua_tostring(L, -1);
        e.dparam = "";
        e.dparam2 = "";
        e.chm = false;
        exprs.push_back(e);

        ists_ = true;
      }
    } else if (lua_isstring(L, -2)) {
      if (ists_ || isb) {
        luaL_error(L, "plot_objf: invalid format");
        return false;
      }

      ists = false;

      const char* key = lua_tostring(L, -2);

      if (strcmp(key, "param") == 0) {
        if (lua_istable(L, -1)) {
          for (lua_pushnil(L); lua_next(L, -2); lua_pop(L, 1)) {
            if (lua_isnumber(L, -2)) {
              if (!chk_str(L, -1)) {
                return false;
              }

              solver_loader::task_info::plot_expr e;
              e.expr = lua_tostring(L, -1);
              e.dparam = "";
              e.dparam2 = "";
              e.chm = false;
              exprs.push_back(e);
            } else if (lua_isstring(L, -2)) {
              solver_loader::task_info::plot_expr e;
              e.expr = lua_tostring(L, -2);
              e.dparam = "";
              e.dparam2 = "";
              e.chm = false;
              exprs.push_back(e);
            } else {
              luaL_error(L, "plot_objf: expression expected");
              return false;
            }
          }
        } else {
          if (!chk_str(L, -1)) {
            return false;
          }

          solver_loader::task_info::plot_expr e;
          e.expr = lua_tostring(L, -1);
          e.dparam = "";
          e.dparam2 = "";
          e.chm = false;
          exprs.push_back(e);
        }
      } else if (strcmp(key, "title") == 0) {
        if (!chk_str(L, -1)) {
          return false;
        }
        std::string title;
        if (!load_str(L, -1, title)) {
          return false;
        }
        props.set_title(title);
      } else {
      }
    } else {
      luaL_error(L, "plot_objf: invalid format");
      return false;
    }
  }

  if (isbv) {
    eo->tasks.push_back(new solver_loader::task_info::plot_objf(props, exprs));
  }

  return true;
}

bool load_plot_objf(lua_State* L, const int index) {
  experiment_options* const eo = get_experiment(L)->get_options();

  if (lua_istable(L, index)) {
    bool isempty = true;
    bool istt = true;

    for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
         lua_pop(L, 1)) {
      isempty = false;
      if (!lua_istable(L, -1)) {
        istt = false;
      }
    }

    if (isempty) {
      solver_loader::task_info::plot_props props;
      std::vector<solver_loader::task_info::plot_expr> exprs;
      eo->tasks.push_back(
          new solver_loader::task_info::plot_objf(props, exprs));
    } else if (istt) {
      for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
           lua_pop(L, 1)) {
        const char* title_ = (lua_isstring(L, -2) && !lua_isnumber(L, -2))
                                 ? lua_tostring(L, -2)
                                 : "";
        if (!load_plot_objf_(L, -1, title_)) {
          return false;
        }
      }
    } else {
      if (!load_plot_objf_(L, index)) {
        return false;
      }
    }
  } else if (lua_isboolean(L, index)) {
    if (lua_toboolean(L, index)) {
      solver_loader::task_info::plot_props props;
      std::vector<solver_loader::task_info::plot_expr> exprs;
      eo->tasks.push_back(
          new solver_loader::task_info::plot_objf(props, exprs));
    }
  } else if (lua_isstring(L, index)) {
    solver_loader::task_info::plot_props props;
    std::vector<solver_loader::task_info::plot_expr> exprs;
    solver_loader::task_info::plot_expr e;
    e.expr = lua_tostring(L, index);
    exprs.push_back(e);
    eo->tasks.push_back(new solver_loader::task_info::plot_objf(props, exprs));
  } else {
    luaL_error(L, "plot_objf: invalid format");
    return false;
  }

  return true;
}

bool load_plot_objf_2d_(lua_State* L, const int index, const char* = "") {
  experiment_options* const eo = get_experiment(L)->get_options();

  solver_loader::task_info::plot_props props;
  solver_loader::task_info::plot_surface_props surface_props;
  std::string expr1;
  std::string expr2;

  for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
       lua_pop(L, 1)) {
    if (!chk_str(L, -2)) {
      return false;
    }

    const char* key = lua_tostring(L, -2);

    if (strcmp(key, "param1") == 0) {
      if (!chk_str(L, -1)) {
        return false;
      }
      if (!load_str(L, -1, expr1)) {
        return false;
      }
    } else if (strcmp(key, "param2") == 0) {
      if (!chk_str(L, -1)) {
        return false;
      }
      if (!load_str(L, -1, expr2)) {
        return false;
      }
    } else if (strcmp(key, "title") == 0) {
      if (!chk_str(L, -1)) {
        return false;
      }
      std::string title;
      if (!load_str(L, -1, title)) {
        return false;
      }
      props.set_title(title);
    } else {
    }
  }

  if (expr1 == "") {
    luaL_error(L, "plot_objf_2d: 'param1' expected");
    return false;
  }

  if (expr2 == "") {
    luaL_error(L, "plot_objf_2d: 'param2' expected");
    return false;
  }

  eo->tasks.push_back(new solver_loader::task_info::plot_objf_2d(
      props, expr1, expr2, surface_props));

  return true;
}

bool load_plot_objf_2d(lua_State* L, const int index) {
  if (lua_istable(L, index)) {
    bool istt = true;
    for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
         lua_pop(L, 1)) {
      if (!lua_istable(L, -1)) {
        istt = false;
      }
    }

    if (istt) {
      for (lua_pushnil(L); lua_next(L, (index < 0) ? index - 1 : index);
           lua_pop(L, 1)) {
        const char* title_ = (lua_isstring(L, -2) && !lua_isnumber(L, -2))
                                 ? lua_tostring(L, -2)
                                 : "";
        if (!load_plot_objf_2d_(L, -1, title_)) {
          return false;
        }
      }
    } else {
      if (!load_plot_objf_2d_(L, index)) {
        return false;
      }
    }
  } else {
    luaL_error(L, "plot_objf_2d: invalid format");
    return false;
  }

  return true;
}

template <class T>
bool load_bool_task(lua_State* L, const int vali, const char* key,
                    const solver_loader::task_info::type tit) {
  experiment_options* const eo = get_experiment(L)->get_options();

  if (!chk_bool(L, vali, key)) {
    return false;
  }

  bool f_ = false;
  if (!load_bool(L, vali, f_)) {
    return false;
  }

  auto const& ti = std::find_if(eo->tasks.begin(), eo->tasks.end(),
                                solver_loader::task_info::find_by_type(tit));
  if (f_) {
    if (ti == eo->tasks.end()) {
      eo->tasks.push_back(new T());
    }
  } else {
    if (ti != eo->tasks.end()) {
      eo->tasks.erase(ti);
    }
  }
  return true;
}

bool set(lua_State* L, const char* key, const int vali) {
  experiment_options* const eo = get_experiment(L)->get_options();

  if (!strcmp(key, "output_dir")) {
    if (!chk_str(L, vali, key)) {
      return false;
    }
    if (!load_str(L, vali, eo->output_dir)) {
      return false;
    }
  } else if (!strcmp(key, "method")) {
    if (!chk_str(L, vali, key)) {
      return false;
    }
    if (!load_str(L, vali, eo->method_name)) {
      return false;
    }
  } else if (!strcmp(key, "opt_method")) {
    if (!chk_str(L, vali, key)) {
      return false;
    }
    if (!load_str(L, vali, eo->opt_method_name)) {
      return false;
    }
  }

  else if (!strcmp(key, "model")) {
    if (!chk_str(L, vali, key)) {
      return false;
    }
    if (!load_str(L, vali, eo->model_name)) {
      return false;
    }
    eo->sbml_model_name = "";
  } else if (!strcmp(key, "model_sbml")) {
    if (!chk_str(L, vali, key)) {
      return false;
    }
    if (!load_str(L, vali, eo->sbml_model_name)) {
      return false;
    }
    eo->model_name = "";
  } else if (!strcmp(key, "kinetics")) {
    if (!chk_str(L, vali, key)) {
      return false;
    }
    std::string kinetics;
    if (!load_str(L, vali, kinetics)) {
      return false;
    }
    if (kinetics == "stochastic" || kinetics == "stoch") {
      eo->set_stoch();
    } else if (kinetics == "deterministic" || kinetics == "det") {
      eo->set_det();
    } else if (kinetics == "hybrid") {
      eo->set_hybrid();
    } else {
      luaL_error(
          L, "stochastic(stoch)/deterministic(det)/hybrid kinetics expected");
      return false;
    }
  } else if (!strcmp(key, "tspan")) {
    if (!load_tspan(L, vali)) {
      return false;
    }
  } else if (!strcmp(key, "init")) {
    if (!load_init(L, vali)) {
      return false;
    }
  } else if (!strcmp(key, "param")) {
    if (!chk_table(L, vali, key)) {
      return false;
    }
    if (!load_se_params(L, vali, eo->set_params, eo->estimate_param,
                        eo->estimate_params, eo->sample_params,
                        eo->iterate_params)) {
      return false;
    }
  } else if (!strcmp(key, "hybrid_vars_force_det")) {
    if (!chk_strvec(L, vali, key)) {
      return false;
    }
    if (!load_strset(L, vali, eo->hybrid_vars_force_det)) {
      return false;
    }
  } else if (!strcmp(key, "hybrid_vars_force_stoch")) {
    if (!chk_strvec(L, vali, key)) {
      return false;
    }
    if (!load_strset(L, vali, eo->hybrid_vars_force_stoch)) {
      return false;
    }
  } else if (!strcmp(key, "obs_err")) {
    if (!chk_table(L, vali, key)) {
      return false;
    }
    if (!load_se_params(L, vali, eo->set_obs_errors, eo->estimate_obs_error,
                        eo->estimate_obs_errors, eo->sample_obs_errors,
                        eo->iterate_obs_errors, true, true)) {
      return false;
    }
  } else if (!strcmp(key, "objf")) {
    if (!chk_str(L, vali, key)) {
      return false;
    }
    if (!load_str(L, vali, eo->objf)) {
      return false;
    }
  } else if (!strcmp(key, "nruns")) {
    if (!chk_number(L, vali, key)) {
      return false;
    }
    if (!load_number(L, vali, key)) {
      return false;
    }
  } else if (!strcmp(key, "ntrajectories")) {
    if (!chk_number(L, vali, key)) {
      return false;
    }
    if (!load_number(L, vali, key)) {
      return false;
    }
  } else if (!strcmp(key, "nrepeat")) {
    if (!chk_number(L, vali, key)) {
      return false;
    }
    if (!load_number(L, vali, key)) {
      return false;
    }
  } else if (!strcmp(key, "nmoments")) {
    if (!chk_number(L, vali, key)) {
      return false;
    }
    if (!load_number(L, vali, eo->nmoments)) {
      return false;
    }
  } else if (!strcmp(key, "det_centered")) {
    if (!chk_bool(L, vali, key)) {
      return false;
    }
    if (!load_bool(L, vali, eo->det_centered)) {
      return false;
    }
  } else if (!strcmp(key, "time_series_data")) {
    if (!chk_str(L, vali, key)) {
      return false;
    }
    if (!load_str(L, vali, eo->time_series_data_src)) {
      return false;
    }
  } else if (!strcmp(key, "steady_state_data")) {
    if (!chk_str(L, vali, key)) {
      return false;
    }
    if (!load_str(L, vali, eo->steady_state_data_src)) {
      return false;
    }
  } else if (!strcmp(key, "observables")) {
    if (!chk_strvec(L, vali, key)) {
      return false;
    }
    if (!load_strvec(L, vali, eo->observables)) {
      return false;
    }
  } else if (!strcmp(key, "stats")) {
    if (!load_bool_task<solver_loader::task_info::dump_stats>(
            L, vali, key, solver_loader::task_info::TASK_DUMP_STATS)) {
      return false;
    }
  } else if (!strcmp(key, "dump_moments")) {
    if (!load_dump_moments(L, vali)) {
      return false;
    }
  } else if (!strcmp(key, "dump_distr")) {
    if (!load_dump_distr(L, vali)) {
      return false;
    }
  } else if (!strcmp(key, "dump_samples")) {
    if (!load_bool_task<solver_loader::task_info::dump_samples>(
            L, vali, key, solver_loader::task_info::TASK_DUMP_SAMPLES)) {
      return false;
    }
  } else if (!strcmp(key, "plot")) {
    if (!load_plot(L, vali)) {
      return false;
    }
  } else if (!strcmp(key, "plot_2d")) {
    if (!load_plot_2d(L, vali)) {
      return false;
    }
  } else if (!strcmp(key, "plot_distr")) {
    if (!load_plot_distr(L, vali)) {
      return false;
    }
  } else if (!strcmp(key, "plot_distr_2d")) {
    if (!load_plot_distr_2d(L, vali)) {
      return false;
    }
  } else if (!strcmp(key, "plot_objf")) {
    if (!load_plot_objf(L, vali)) {
      return false;
    }
  } else if (!strcmp(key, "plot_objf_2d")) {
    if (!load_plot_objf_2d(L, vali)) {
      return false;
    }
  } else if (!strcmp(key, "rare_event")) {
    if (!chk_str(L, vali, key)) {
      return false;
    }
    if (!load_str(L, vali, key)) {
      return false;
    }
  } else if (!strcmp(key, "uniformization_rate")) {
    if (!chk_number(L, vali, key)) {
      return false;
    }
    if (!load_number(L, vali, key)) {
      return false;
    }
  } else if (!strcmp(key, "uniformization_method")) {
    if (!chk_number(L, vali, key)) {
      return false;
    }
    if (!load_number(L, vali, key)) {
      return false;
    }
  } else if (!strcmp(key, "bpcumulativeprobeps")) {
    if (!chk_number(L, vali, key)) {
      return false;
    }
    if (!load_number(L, vali, key)) {
      return false;
    }
  } else if (!strcmp(key, "steady_state")) {
    if (!chk_bool(L, vali, key)) {
      return false;
    }
    if (!load_bool(L, vali, key)) {
      return false;
    }
  } else if (!strcmp(key, "tol")) {
    if (!chk_number(L, vali, key)) {
      return false;
    }
    if (!load_number(L, vali, eo->tol)) {
      return false;
    }
  } else if (!strcmp(key, "abs_tol")) {
    if (!chk_number(L, vali, key)) {
      return false;
    }
    if (!load_number(L, vali, eo->abs_tol)) {
      return false;
    }
  } else if (!strcmp(key, "rel_tol")) {
    if (!chk_number(L, vali, key)) {
      return false;
    }
    if (!load_number(L, vali, eo->rel_tol)) {
      return false;
    }
  } else {
    if (lua_isnumber(L, vali)) {
      eo->set(key, lua_tonumber(L, vali));
    } else if (lua_isboolean(L, vali)) {
      eo->set(key, (bool)lua_toboolean(L, vali));
    } else if (lua_isstring(L, vali)) {
      const char* val = lua_tostring(L, vali);
      eo->set(key, val);
    } else if (lua_istable(L, vali)) {
      if (chk_dvec(L, vali, key)) {
        std::vector<double> dvec;

        if (!load_dvec(L, vali, dvec)) {
          return false;
        }

        eo->set(key, dvec);
      } else if (chk_strvec(L, vali, key)) {
        std::vector<std::string> strvec;

        if (!load_strvec(L, vali, strvec)) {
          return false;
        }

        eo->set(key, strvec);
      }
    } else {
      luaL_error(L, "Invalid format");
      return false;
    }
  }

  return true;
}

bool load_options(lua_State* L, int index) {
  if (!lua_istable(L, index)) {
    luaL_error(L, "Invalid format");
    return false;
  }

  for (lua_pushnil(L); lua_next(L, -2); lua_pop(L, 1)) {
    if (!chk_str(L, -2)) {
      return false;
    }

    if (!::set(L, lua_tostring(L, -2), -1)) {
      return false;
    }
  }

  return true;
}

int set(lua_State* L) {
  const int narg = lua_gettop(L);
  if (!chk_two_args(L, narg)) {
    return 0;
  }
  if (!chk_str(L, 1)) {
    return 0;
  }
  if (!::set(L, lua_tostring(L, 1), 2)) {
    return 0;
  }
  return 0;
}

int unset(lua_State* L) {
  experiment_options* const eo = get_experiment(L)->get_options();
  const int narg = lua_gettop(L);
  if (narg) {
    for (int i = 1; i <= narg; i++) {
      if (!chk_str(L, i)) {
        return 0;
      }
      eo->unset(lua_tostring(L, i));
    }
  } else {
    eo->unset();
  }
  return 0;
}

int set_model(lua_State* L) {
  experiment_options* const eo = get_experiment(L)->get_options();
  eo->sbml_model_name = "";
  return set_str(L, eo->model_name);
}

int get_model_list(lua_State*) { return 1; }

int set_sbml_model(lua_State* L) {
  experiment_options* const eo = get_experiment(L)->get_options();
  eo->model_name = "";
  return set_str(L, eo->sbml_model_name);
}

int set_custom_model(lua_State* L) {
  experiment_options* const eo = get_experiment(L)->get_options();
  eo->model_name = "";
  eo->sbml_model_name = "";
  eo->custom_model_id = luaL_checknumber(L, 1);
  return 0;
}

int set_init(lua_State* L) {
  const int narg = lua_gettop(L);
  if (!chk_one_arg(L, narg)) {
    return 0;
  }
  if (!load_init(L, 1)) {
    return 0;
  }
  return 0;
}

int get_init(lua_State*) { return 1; }

int set_stoch(lua_State* L) {
  const int narg = lua_gettop(L);
  if (!chk_zero_arg(L, narg)) {
    return 0;
  }
  experiment_options* const eo = get_experiment(L)->get_options();
  eo->set_stoch();
  return 0;
}

int set_det(lua_State* L) {
  const int narg = lua_gettop(L);
  if (!chk_zero_arg(L, narg)) {
    return 0;
  }
  experiment_options* const eo = get_experiment(L)->get_options();
  eo->set_det();
  return 0;
}

int set_hybrid(lua_State* L) {
  const int narg = lua_gettop(L);
  if (!chk_zero_arg(L, narg)) {
    return 0;
  }
  experiment_options* const eo = get_experiment(L)->get_options();
  eo->set_hybrid();
  return 0;
}

int set_param(lua_State* L) {
  experiment_options* const eo = get_experiment(L)->get_options();
  return set_se_params(L, eo->set_params, eo->estimate_param,
                       eo->estimate_params, eo->sample_params,
                       eo->iterate_params);
}

int get_param(lua_State*) { return 1; }

int get_param_list(lua_State*) { return 1; }

int reset_params(lua_State* L) {
  experiment_options* const eo = get_experiment(L)->get_options();
  eo->set_params.clear();
  return 0;
}

int set_objf(lua_State* L) {
  experiment_options* const eo = get_experiment(L)->get_options();
  return set_str(L, eo->objf);
}

int set_tspan(lua_State* L) {
  const int narg = lua_gettop(L);
  if (!chk_one_arg(L, narg)) {
    return 0;
  }
  if (!load_tspan(L, 1)) {
    return 0;
  }
  return 0;
}

int get_tspan(lua_State*) { return 1; }

int plot(lua_State* L) {
  const int narg = lua_gettop(L);
  if (!chk_one_arg(L, narg)) {
    return 0;
  }
  if (!load_plot(L, 1)) {
    return 0;
  }
  return 0;
}

int plot_2d(lua_State* L) {
  const int narg = lua_gettop(L);
  if (!chk_one_arg(L, narg)) {
    return 0;
  }
  if (!load_plot_2d(L, 1)) {
    return 0;
  }
  return 0;
}

int plot_distr(lua_State* L) {
  const int narg = lua_gettop(L);
  if (!chk_one_arg(L, narg)) {
    return 0;
  }
  if (!load_plot_distr(L, 1)) {
    return 0;
  }
  return 0;
}

int plot_distr_2d(lua_State* L) {
  const int narg = lua_gettop(L);
  if (!chk_one_arg(L, narg)) {
    return 0;
  }
  if (!load_plot_distr_2d(L, 1)) {
    return 0;
  }
  return 0;
}

int plot_objf(lua_State* L) {
  const int narg = lua_gettop(L);
  if (!chk_one_arg(L, narg)) {
    return 0;
  }
  if (!load_plot_objf(L, 1)) {
    return 0;
  }
  return 0;
}

int plot_objf_2d(lua_State* L) {
  const int narg = lua_gettop(L);
  if (!chk_one_arg(L, narg)) {
    return 0;
  }
  if (!load_plot_objf_2d(L, 1)) {
    return 0;
  }
  return 0;
}

int clear_tasks(lua_State* L) {
  experiment_options* const eo = get_experiment(L)->get_options();
  eo->clear_tasks();
  return 0;
}

int set_observables(lua_State* L) {
  experiment_options* const eo = get_experiment(L)->get_options();
  return set_strvec(L, eo->observables);
}

int set_observation_error(lua_State* L) {
  experiment_options* const eo = get_experiment(L)->get_options();
  return set_se_params(L, eo->set_obs_errors, eo->estimate_obs_error,
                       eo->estimate_obs_errors, eo->sample_obs_errors,
                       eo->iterate_obs_errors, true, true);
}

int compose_transitions(lua_State*) { return 0; }

int set_params(lua_State* L) {
  experiment_options* const eo = get_experiment(L)->get_options();
  const int narg = lua_gettop(L);
  if (narg) {
    for (int i = 0; i < narg; i++) {
      est_param_info ep;
      ep.param = luaL_checkstring(L, 1);
      eo->estimate_params.push_back(ep);
    }
    eo->estimate_param = true;
  }
  return 0;
}

int add_sl(lua_State* L, solver_loader::base* const sl) {
  experiment* const E = get_experiment(L);

  if (sl->compile_solver()) {
    E->add_sl(sl);
  } else {
    E->print_error(sl->get_last_error());
  }
  return 0;
}

int chk_model(lua_State* L) {
  experiment* const E = get_experiment(L);
  set_str(L, E->get_options()->model_name);
  solver_loader::base sl(E);
  std::stringstream err;
  err << "Loading '" << E->get_options()->model_name << "'";
  E->print_status(err.str());
  err.str("");
  if (!sl.load_model()) {
    E->print_error(err.str());
  }
  return 0;
}

int chk_model_err(lua_State* L) {
  experiment* const E = get_experiment(L);
  set_str(L, E->get_options()->model_name);
  solver_loader::base sl(E);
  std::stringstream err;
  err << "Loading '" << E->get_options()->model_name << "'";
  E->print_status(err.str());
  err.str("");
  if (sl.load_model()) {
    E->print_error("did not fail");
  }
  return 0;
}

template <class T>
int add_sl(lua_State* L) {
  experiment* const E = get_experiment(L);

  E->push_options();

  const int narg = lua_gettop(L);
  if (narg) {
    if (!chk_one_arg(L, narg)) {
      return 0;
    }
    if (!load_options(L, 1)) {
      return 0;
    }
  }

  const int r = add_sl(L, new T(E));

  E->pop_options();

  return r;
}

int transient(lua_State* L) { return add_sl<solver_loader::transient>(L); }

int scan(lua_State* L) { return add_sl<solver_loader::scan>(L); }

int optimize(lua_State* L) { return add_sl<solver_loader::optimize>(L); }

int fit(lua_State* L) { return add_sl<solver_loader::fit>(L); }

char* lua_get_error(lua_State* L, char* const buf) {
  char const* s = lua_tostring(L, -1);

  if (*s == '[') {
    ++s;
    while (*s != ']') {
      ++s;
    }
    ++s;
    ++s;
  }

  char* p = buf;
  for (const char* c = s; *c != '\0'; ++c, ++p) {
    if (*c == '"') {
      *p = '\\';
      p++;
    }
    *p = *c;
  }
  *p = '\0';

  lua_pop(L, 1);

  return buf;
}

int panicf(lua_State* L) {
  char err[1024];
  lua_get_error(L, err);

  experiment* const E = get_experiment(L);
  E->print_error(err);

#ifdef STAR_WEB_INTERFACE
  if (!E->get_dbconnector()->end_experiment(E->get_id(),
                                            EXPERIMENT_STATUS_FAILED, err)) {
    sprintf(err, "Error stopping experiment #%d", E->get_id());
    E->print_error_(err);
  }
#endif  // STAR_WEB_INTERFACE

  return EXIT_FAILURE;
}

experiment::~experiment() {
  lua_close(L);

  for (auto& sl : SL) {
    // TODO(memory)
    /*delete sl;*/
  }
}

bool experiment::create_dir(char* const err) const {
  namespace fs = boost::filesystem;
  try {
    fs::remove_all(dir);
    fs::create_directory(dir);
  } catch (const fs::filesystem_error& e) {
    sprintf(err, "%s", e.what());
    return false;
  }

  boost::filesystem::path p = output_dir / "last.txt";
  std::ofstream os(convws(p.wstring()));
  if (!os.is_open()) {
    sprintf(err, "Couldn't create 'last.txt'");
    return false;
  }
  os << convws(dir.wstring());
  os.close();

  return true;
}

bool experiment::run_solvers(char* const err) const {
  bool success = true;

#ifdef STAR_HAVE_MCR
/*if (useMATLAB) {
  const char* args[] = {"-nodesktop", "-nodisplay", "-nosplash"};
  if (!mclInitializeApplication(args, 3)) {
    std::cerr << "mclInitializeApplication = false";
    return false;
  }
  if (!libmatlabsdInitialize()) {
    std::cerr << "libmatlabsdInitialize = false";
    return false;
  }
}*/
#endif  // STAR_HAVE_MCR

  for (auto& sl : SL) {
    if (!sl->run_solver()) {
      strcpy(err, sl->get_last_error().c_str());
      success = false;
      break;
    }
  }

#ifdef STAR_HAVE_MCR
/*if (useMATLAB) {
  // libmatlabsdTerminate();
  mclTerminateApplication();
}*/
#endif  // STAR_HAVE_MCR

  return success;
}

bool experiment::process(char const* const src) const {
  char err[1024];

#ifdef STAR_WEB_INTERFACE
  if (!dbconnector->start_experiment(id)) {
    sprintf(err, "Error starting experiment #%d", id);
    print_error_(err);
    return false;
  }
#else   // STAR_WEB_INTERFACE
  print_status("Starting experiment...");
  std::cout << src << std::endl;
#endif  // STAR_WEB_INTERFACE

  bool success = false;

  if (create_dir(err)) {
    if (!luaL_dostring(L, src)) {
      // TODO
      // get("norun", norun);

      if (run_solvers(err)) {
        success = true;
      }
    } else {
      lua_get_error(L, err);
    }
  }

  // TODO ?
  /*if (success) {
    bool notify;
    if (get("notify", notify) && notify) {
      static char email[32];
#ifdef STAR_WEB_INTERFACE
// TODO
#else   // STAR_WEB_INTERFACE
      sprintf(email, "linar.mikeev@gmail.com");
#endif  // STAR_WEB_INTERFACE

      static char tmp[128];
      sprintf(tmp, "echo \"Your job has completed.\" | mailx
                       - s
  \"SHAVE notification\" %s",
              email);
      if (system(tmp)) {
        sprintf(err, "Failed to send a notification");
        success = false;
      }
    }
  }*/

  if (!success) {
    print_error(err);
  }

#ifdef STAR_WEB_INTERFACE
  if (!dbconnector->end_experiment(
          id, success ? EXPERIMENT_STATUS_COMPLETED : EXPERIMENT_STATUS_FAILED,
          err)) {
    sprintf(err, "Error stopping experiment #%d", id);
    print_error_(err);
    return false;
  }
#else   // STAR_WEB_INTERFACE
  if (success) {
    write_report();
  }
#endif  // STAR_WEB_INTERFACE

  return success;
}

void experiment::lua_setup() {
  lua_atpanic(L, panicf);

  luaopen_base(L);
  luaopen_table(L);
  luaopen_string(L);
  luaopen_math(L);

  lua_pushlightuserdata(L, this);
  lua_setglobal(L, "experiment");

  lua_register(L, "linspace", linspace);
  lua_register(L, "logspace", logspace);

  lua_register(L, "uniform", uniform);
  lua_register(L, "bernoulli", bernoulli);
  lua_register(L, "binomial", binomial);
  lua_register(L, "neg_binomial", neg_binomial);
  lua_register(L, "geometric", geometric);
  lua_register(L, "poisson", poisson);
  lua_register(L, "exponential", exponential);
  lua_register(L, "gamma", gamma);
  lua_register(L, "weibull", weibull);
  lua_register(L, "normal", normal);
  lua_register(L, "lognormal", lognormal);

  lua_register(L, "set", ::set);
  lua_register(L, "unset", ::unset);

  lua_register(L, "set_model", set_model);
// lua_register(L, "get_model_list", get_model_list);

#ifndef STAR_WEB_INTERFACE
  lua_register(L, "set_sbml_model", set_sbml_model);
#endif  // STAR_WEB_INTERFACE

  lua_register(L, "set_custom_model", set_custom_model);

  lua_register(L, "set_stoch", set_stoch);
  lua_register(L, "set_det", set_det);
  lua_register(L, "set_hybrid", set_hybrid);

  lua_register(L, "set_param", set_param);
  // lua_register(L, "get_param", get_param);
  // lua_register(L, "get_param_list", get_param_list);
  lua_register(L, "reset_params", reset_params);

  lua_register(L, "set_init", set_init);
  lua_register(L, "get_init", get_init);

  lua_register(L, "set_tspan", set_tspan);
  lua_register(L, "get_tspan", get_tspan);

  lua_register(L, "set_objf", set_objf);

  lua_register(L, "plot", plot);
  lua_register(L, "plot_2d", plot_2d);

  lua_register(L, "plot_distr", plot_distr);
  lua_register(L, "plot_distr_2d", plot_distr_2d);

  lua_register(L, "plot_objf", plot_objf);
  lua_register(L, "plot_objf_2d", plot_objf_2d);

  lua_register(L, "clear_tasks", ::clear_tasks);

  lua_register(L, "set_observables", set_observables);

  lua_register(L, "set_observation_error", set_observation_error);

  lua_register(L, "transient", transient);

  lua_register(L, "scan", scan);
  lua_register(L, "optimize", optimize);
  lua_register(L, "fit", fit);

  /*lua_register(L, "steady_state", steady_state);
  lua_register(L, "verify", verify);

  lua_register(L, "profile", profile);

#ifndef STAR_WEB_INTERFACE
  lua_register(L, "compose_transitions", compose_transitions);
  lua_register(L, "set_pta", set_pta);
  lua_register(L, "set_clock_regions", set_clock_regions);

  lua_register(L, "set_params", set_params);
  lua_register(L, "set_levels", set_levels);

  lua_register(L, "multilevel_ce", multilevel_ce);

  lua_register(L, "dmpcmp", dmpcmp);
#endif  // STAR_WEB_INTERFACE*/

  // TODO
  luaopen_package(L);
  lua_register(L, "chk_model", chk_model);
  lua_register(L, "chk_model_err", chk_model_err);
}

void experiment::write_task_report(solver::task const* const tsk) const {
  solver_loader::base const* const sl = tsk->get_solver()->get_solver_loader();

  boost::filesystem::path p =
      dir / sl->get_dir() / tsk->get_dir() / "index.html";
  std::ofstream out(convws(p.wstring()));

  static char tmp[128];
  sprintf(tmp, "Task - %s", tsk->get_info()->get_name());
  // TODO: path
  solver::util::htmlHelper::write_beginning(out, tmp, "../../");

  std::size_t npt = 0, npnt = 0;
  std::size_t ndt = 0, ndnt = 0;

  for (auto const& d : tsk->get_plots()) {
    if (d->get_timepoint_index() < 0) {
      npnt++;
    } else {
      npt++;
    }
  }
  for (auto const& d : tsk->get_dumps()) {
    if (d->get_timepoint_index() < 0) {
      ndnt++;
    } else {
      ndt++;
    }
  }

  if (npnt) {
    solver::util::htmlHelper::write_begin_table(out);
    out << "<tr><th>Plot</th></tr>" << std::endl;
    for (auto const& d : tsk->get_plots()) {
      if (d->get_timepoint_index() < 0) {
        out << "<tr><td>";
        solver::util::htmlHelper::write_image(out, d->get_local_fname(),
                                              d->get_name());
        out << "</td></tr>" << std::endl;
      }
    }
    solver::util::htmlHelper::write_end_table(out);
  }
  if (ndnt) {
    solver::util::htmlHelper::write_begin_table(out);
    out << "<tr><th>Data</th></tr>" << std::endl;
    for (auto const& d : tsk->get_dumps()) {
      if (d->get_timepoint_index() < 0) {
        out << "<tr><td>";
        solver::util::htmlHelper::write_link(out, d->get_local_fname(),
                                             d->get_name());
        out << "</td></tr>";
      }
    }
    solver::util::htmlHelper::write_end_table(out);
  }
  if (npt) {
    solver::util::htmlHelper::write_begin_table(out);
    out << "<tr><th>Timepoint</th><th>Plot</th></tr>" << std::endl;
    for (auto const& d : tsk->get_plots()) {
      if (d->get_timepoint_index() >= 0) {
        double t;
        sl->get_solver()->get_timepoint(d->get_timepoint_index(), t);

        out << "<tr><td class='row-caption'>" << t << "</td><td>";
        solver::util::htmlHelper::write_image(out, d->get_local_fname(),
                                              d->get_name());
        out << "</td></tr>" << std::endl;
      }
    }
    solver::util::htmlHelper::write_end_table(out);
  }
  if (ndt) {
    solver::util::htmlHelper::write_begin_table(out);
    out << "<tr><th>Timepoint</th><th>Data</th></tr>" << std::endl;
    for (auto const& d : tsk->get_dumps()) {
      if (d->get_timepoint_index() >= 0) {
        double t;
        sl->get_solver()->get_timepoint(d->get_timepoint_index(), t);

        out << "<tr><td class='row-caption'>" << t << "</td><td>";
        solver::util::htmlHelper::write_link(out, d->get_local_fname(),
                                             d->get_name());
        out << "</td></tr>";
      }
    }
    solver::util::htmlHelper::write_end_table(out);
  }
  solver::util::htmlHelper::write_ending(out);
  out.close();
}

void experiment::write_solution_report(
    solver_loader::base const* const sl) const {
  boost::filesystem::path p = dir / sl->get_dir() / "index.html";
  std::ofstream out(convws(p.wstring()));

  static char desc[1024], tmp[1024];
  sl->sprint_description(desc);
  sprintf(tmp, "Solution - %s", desc);
  // TODO: path
  solver::util::htmlHelper::write_beginning(out, tmp, "../");
  solver::util::htmlHelper::write_begin_table(out);
  out << "<tr><th>Task</th></tr>" << std::endl;
  for (auto const& t : sl->get_solver()->get_tasks()) {
    out << "<tr><td>";
    p = boost::filesystem::path(t->get_dir()) / "index.html";
    solver::util::htmlHelper::write_link(out, p.c_str(),
                                         t->get_info()->get_name());
    out << "</td></tr>" << std::endl;

    write_task_report(t);
  }
  solver::util::htmlHelper::write_end_table(out);
  solver::util::htmlHelper::write_ending(out);
  out.close();
}

void experiment::write_report() const {
  boost::filesystem::path p = dir / "report.html";
  std::ofstream out(convws(p.wstring()));

  std::string cmd = "cp report.css ";
  cmd += convws(dir.wstring());
  system(cmd.c_str());

  static char tmp[128];
  sprintf(tmp, "Report - %s", dir.c_str());
  solver::util::htmlHelper::write_beginning(out, tmp, ".");
  solver::util::htmlHelper::write_begin_table(out);
  out << "<tr><th>Results</th></tr>" << std::endl;
  for (auto const& sl : SL) {
    static char desc[1024];
    sl->sprint_description(desc);

    out << "<tr><td>";
    p = boost::filesystem::path(sl->get_dir()) / "index.html";
    solver::util::htmlHelper::write_link(out, p.c_str(), desc);
    out << "</td></tr>" << std::endl;

    write_solution_report(sl);
  }
  solver::util::htmlHelper::write_end_table(out);
  solver::util::htmlHelper::write_ending(out);
  out.close();
}
