/*
 *  dump.hpp
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

#ifndef SOLVER_DUMP_HPP_
#define SOLVER_DUMP_HPP_

#include <clocale>

namespace solver {

class task;

class dump {
 private:
  task const* const tsk;
  const std::string name;
  const int timepoint_index;

  std::wstring local_path;
  std::string local_fname;
  std::string web_fname;

 public:
  dump(task const* const tsk, char const* const name,
       const std::wstring& local_path, const std::string& local_fname,
       const std::string& web_fname, const int timepoint_index = -1)
      : tsk(tsk),
        name(name),
        timepoint_index(timepoint_index),
        local_path(local_path),
        local_fname(local_fname),
        web_fname(web_fname) {}

  virtual ~dump() {}

  task const* get_task() const { return tsk; }

  char const* get_name() const { return name.c_str(); }

  int get_timepoint_index() const { return timepoint_index; }

  char const* get_local_path() const {
    static char tmp[1024];
    const std::size_t l = local_path.length();
    std::wcstombs(tmp, local_path.c_str(), l);
    tmp[l] = '\0';
    return tmp;
  }

  char const* get_local_fname() const { return local_fname.c_str(); }

  char const* get_web_fname() const { return web_fname.c_str(); }

  bool update_web_copy() const {
#ifdef STAR_WEB_INTERFACE

    static const char* web_output_dir = "/var/www/star/app/webroot/r/";

    static char tmp[1024];
    sprintf(tmp, "cp %s %s%s", get_local_path(), web_output_dir,
            web_fname.c_str());
    if (system(tmp)) {
      return false;
    }
#endif

    return true;
  }
};
}

#endif
