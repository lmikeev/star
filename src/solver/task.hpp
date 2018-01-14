/*
 *  task.hpp
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

#ifndef SOLVER_TASK_HPP_
#define SOLVER_TASK_HPP_

#include <vector>
#include <random>
#include <boost/filesystem.hpp>
#include "dump.hpp"
#include "../solver_loader/task_info/base.hpp"

namespace solver {

class base;

class task {
 private:
  base const* const S;
  unsigned int const id;
  solver_loader::task_info::base const* const info;

  const boost::filesystem::path path;

  char dir[32];

  std::vector<dump*> dumps;
  std::vector<dump*> plots;

  unsigned int dyn_plot_id;

  void set_dump_paths(std::wstring& local_path, std::string& local_fname,
                      std::string& web_fname, char const* const dump_name,
                      char const* const fname, char const* const fext,
                      const int index = -1) {
    static char tmp[128];
    if (index < 0) {
      sprintf(tmp, "%s.%s", fname, fext);
    } else {
      sprintf(tmp, "%s_%d.%s", fname, index, fext);
    }

    local_fname = tmp;

    const boost::filesystem::path p = path / dir / local_fname;
    local_path = p.wstring();

#ifdef STAR_WEB_INTERFACE
    static char buf[16];
    if (index < 0) {
      sprintf(tmp, "%d_%s_%s.%s", id, fname, rand_string(buf, 7), fext);
    } else {
      sprintf(tmp, "%d_%s_%d_%s.%s", id, fname, index, rand_string(buf, 7),
              fext);
    }
    web_fname = tmp;
#else

    web_fname = web_fname;
    srand(*dump_name);
#endif
  }

  char* rand_string(char* const s, const std::size_t len) {
    static std::random_device rand_device;
    static std::mt19937 rand_engine(rand_device());

    static const char alphanum[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";

    static std::uniform_int_distribution<uint32_t> rand_index(
        0, sizeof(alphanum) - 2);

    for (std::size_t i = 0; i < len; i++) {
      s[i] = alphanum[rand_index(rand_engine)];
    }
    s[len] = 0;

    return s;
  }

 public:
  task(base const* const s, const unsigned int id,
       solver_loader::task_info::base const* const info,
       const boost::filesystem::path& path)
      : S(s), id(id), info(info), path(path), dyn_plot_id(0) {
    sprintf(dir, "%d", id);
  }

  ~task() {}

  base const* get_solver() const { return S; }

  unsigned int get_id() const { return id; }

  solver_loader::task_info::base const* get_info() const { return info; }

  char const* get_dir() const { return dir; }

  char const* get_local_dir() const {
    const boost::filesystem::path p = path / dir;
    const std::wstring local_dir = p.wstring();

    static char tmp[1024];
    const std::size_t l = local_dir.length();
    std::wcstombs(tmp, local_dir.c_str(), l);
    tmp[l] = '\0';
    return tmp;
  }

  dump const* get_dump(const unsigned int i) const { return dumps[i]; }

  const std::vector<dump*>& get_dumps() const { return dumps; }

  dump const* get_plot(const unsigned int i) const { return plots[i]; }

  const std::vector<dump*>& get_plots() const { return plots; }

  bool add_dump(char const* const dump_name, char const* const fname,
                char const* const fext = "csv", const int index = -1) {
    std::wstring local_path;
    std::string local_fname, web_fname;
    set_dump_paths(local_path, local_fname, web_fname, dump_name, fname, fext,
                   index);

    dumps.push_back(
        new dump(this, dump_name, local_path, local_fname, web_fname));

    return true;
  }

  bool add_dump(const std::size_t timepoint_index, char const* const dump_name,
                char const* const fname, char const* const fext = "csv") {
    std::wstring local_path;
    std::string local_fname, web_fname;
    set_dump_paths(local_path, local_fname, web_fname, dump_name, fname, fext,
                   timepoint_index);

    dumps.push_back(new dump(this, dump_name, local_path, local_fname,
                             web_fname, timepoint_index));

    return true;
  }

  bool add_plot(char const* const dump_name, char const* const fname,
                char const* const fext = "png", const int index = -1) {
    std::wstring local_path;
    std::string local_fname, web_fname;
    set_dump_paths(local_path, local_fname, web_fname, dump_name, fname, fext,
                   index);

    plots.push_back(
        new dump(this, dump_name, local_path, local_fname, web_fname));

    return true;
  }

  bool add_plot(const std::size_t timepoint_index, char const* const dump_name,
                char const* const fname, char const* const fext = "png") {
    std::wstring local_path;
    std::string local_fname, web_fname;
    set_dump_paths(local_path, local_fname, web_fname, dump_name, fname, fext,
                   timepoint_index);

    plots.push_back(new dump(this, dump_name, local_path, local_fname,
                             web_fname, timepoint_index));

    return true;
  }

  void set_dyn_plot_id(const unsigned int dyn_plot_id_) {
    dyn_plot_id = dyn_plot_id_;
  }

  unsigned int get_dyn_plot_id() const { return dyn_plot_id; }

  bool update_web_copies() const {
#ifdef STAR_WEB_INTERFACE
    for (auto const& d : dumps) {
      if (!d->update_web_copy()) {
        return false;
      }
    }
    for (auto const& p : plots) {
      if (!p->update_web_copy()) {
        return false;
      }
    }
#endif

    return true;
  }
};
}

#endif
