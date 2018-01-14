/*
 *  star.hpp
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

#ifndef STAR_HPP_
#define STAR_HPP_

#include <inttypes.h>

#ifdef STAR_WEB_INTERFACE
#include "src/mysqlconnector/mysqlconnector.hpp"
#endif

typedef double fnsi_rate_t;
typedef uint32_t fnsi_state_t;

struct fnsi_transition {
  fnsi_rate_t rate;
  fnsi_state_t dst;
};

struct star_fnsi {
  int32_t (*load_model)(const char* fname);
  int32_t (*get_variables_count)();
  int32_t (*get_variable_name)(const uint32_t index, char* buf,
                               const uint32_t buf_len);
  int32_t (*get_max_transitions_count)();
  int32_t (*get_transitions)(const fnsi_state_t s, fnsi_transition* buf,
                             const uint32_t buf_len);
  int32_t (*get_variable_value)(const fnsi_state_t s, const uint32_t index);
  int32_t (*get_last_error)(char* buf, const uint32_t buf_len);
};

extern "C" int star_run_file(const char* experiment_fname, star_fnsi* fnsi);
extern "C" int star_run_src(const char* experiment_src, star_fnsi* fnsi);

#endif
