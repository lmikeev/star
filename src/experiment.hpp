/*
 *  experiment.hpp
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

#ifndef EXPERIMENT_HPP_
#define EXPERIMENT_HPP_

#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "experiment_options.hpp"
#include "solver_loader/task_info/base.hpp"
#include "../star.hpp"

namespace solver_loader {
class base;
}

namespace solver {
class task;
}

class experiment {
 private:
  experiment_options* options_;
  experiment_options* options;

  star_fnsi* const fnsi;
  lua_State* const L;

  const boost::filesystem::path output_dir;

  boost::filesystem::path dir;

  std::vector<solver_loader::base*> SL;

#ifdef STAR_WEB_INTERFACE
  mysqlconnector* dbconnector;
  unsigned int id;
  unsigned int user_id;
#endif

  void set_local_dir() {
    static char tmp[32];
    time_t rawtime;
    time(&rawtime);
    struct tm* timeinfo = localtime(&rawtime);
    strftime(tmp, sizeof(tmp), "%y-%m-%d_%H-%M-%S/", timeinfo);

    dir = output_dir / tmp;
  }

  bool create_dir(char* const err) const;
  bool run_solvers(char* const err) const;

 protected:
  void lua_setup();
  bool process(char const* const src) const;

  void write_task_report(solver::task const* const tsk) const;
  void write_solution_report(solver_loader::base const* const sl) const;
  void write_report() const;

 public:
  experiment(star_fnsi* const fnsi = nullptr)
      : options_(nullptr),
        options(new experiment_options()),
        fnsi(fnsi),
        L(luaL_newstate()),
        output_dir("output") {
    lua_setup();
  }
  ~experiment();

#ifdef STAR_WEB_INTERFACE

  bool run(mysqlconnector* const dbconnector_, const unsigned int user_id_,
           const unsigned int id_) {
    dbconnector = dbconnector_;
    id = id_;
    user_id = user_id_;

    static char tmp[16];
    sprintf(tmp, "%d", id);
    dir = boost::filesystem::path("output") / tmp;

    char* src = nullptr;
    if (!dbconnector->load_experiment_script(src, id)) {
      print_error_("Error loading experiment script");
      return false;
    }

    const bool success = process(src);

    return success;
  }

#else

  bool run_file(const std::string& fname) {
    set_local_dir();

    std::ifstream in(fname.c_str(), std::ios::in | std::ios::binary);
    if (!in.is_open()) {
      static char tmp[1024];
      sprintf(tmp, "Couldn't open '%s'", fname.c_str());
      print_error(tmp);
      return false;
    }

    std::string src;
    in.seekg(0, std::ios::end);
    src.resize(in.tellg());
    in.seekg(0, std::ios::beg);
    in.read(&src[0], src.size());
    in.close();

    return process(src.c_str());
  }

  bool run_src(const std::string& src) {
    set_local_dir();

    return process(src.c_str());
  }
#endif

  void push_options() {
    options_ = options;
    options = new experiment_options(options);
  }

  void pop_options() {
    assert(options_ != nullptr);

    options = options_;
    options_ = nullptr;
  }

  experiment_options* get_options() const { return options; }

  const boost::filesystem::path& get_dir() const { return dir; }

  star_fnsi* get_fnsi() const { return fnsi; }

#ifdef STAR_WEB_INTERFACE
  mysqlconnector* get_dbconnector() const { return dbconnector; }

  unsigned int get_id() const { return id; }

  unsigned int get_user_id() const { return user_id; }
#endif

  const std::vector<solver_loader::base*>& get_SL() const { return SL; }

  std::vector<solver_loader::base*>& get_SL() { return SL; }

  void add_sl(solver_loader::base* const sl) { SL.push_back(sl); }

  void print_status(const std::string& status) const {
    printf("%s\n", status.c_str());
  }

  void print_error_(const std::string& err) const {
    printf("ERROR: %s\n", err.c_str());
  }

  void print_error(const std::string& err) const {
#ifndef STAR_WEB_INTERFACE
    print_error_(err);
#endif
  }

  char const* convws(const std::wstring& ws) const {
    static char tmp[1024];
    const std::size_t l = ws.length();
    std::wcstombs(tmp, ws.c_str(), l);
    tmp[l] = '\0';
    return tmp;
  }
};

bool load_config(boost::property_tree::ptree& propTree,
                 char const* const config_ini = "star.ini");

#ifdef STAR_WEB_INTERFACE
bool db_connect(mysqlconnector& dbconnector_,
                const boost::property_tree::ptree& propTree);
#endif

#endif
