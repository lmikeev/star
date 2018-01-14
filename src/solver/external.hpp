/*
 *  external.hpp
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

#ifndef SOLVER_EXTERNAL_HPP_
#define SOLVER_EXTERNAL_HPP_

#include "base.hpp"

extern void external_set_solver(solver::base* S);

extern void lsoda_f(double t, double* y, double* ydot, void* data);

extern double nlopt_f(const std::vector<double>& x, std::vector<double>& grad,
                      void* data);

#if HAVE_DLIB
extern double dlib_opt_f(const column_vector& x);
extern column_vector dlib_opt_g(const column_vector& x);
extern general_matrix dlib_opt_h(const column_vector& x);
#endif

#endif
