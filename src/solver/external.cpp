/*
 *  external.cpp
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

#include "external.hpp"

solver::base* S;

void external_set_solver(solver::base* s) { S = s; }

void lsoda_f(double t, double* y, double* ydot, void* data) {
  S->lsoda_f(t, y, ydot, data);
}

double nlopt_f(const std::vector<double>& x, std::vector<double>& grad,
               void* data) {
  return S->nlopt_f(x, grad, data);
}

#if HAVE_DLIB

double dlib_opt_f(const column_vector& x) { return S->dlib_opt_f(x); }

column_vector dlib_opt_g(const column_vector& x) { return S->dlib_opt_g(x); }

general_matrix dlib_opt_h(const column_vector& x) { return S->dlib_opt_h(x); }

#endif

#if HAVE_MCR

#include "matlab/libmatlabsd.hpp"
#include "matlab/matlab_ft_external.h"
#include "matlab/matlab_f_external.h"
#include "matlab/matlab_fg_external.h"
#include "matlab/matlab_fgh_external.h"

bool MW_CALL_CONV
    matlab_ft(int nlhs, mxArray* plhs[], int nrhs, mxArray* prhs[]) {
  return S->matlab_ft(nlhs, plhs, nrhs, prhs);
}

bool MW_CALL_CONV
    matlab_f(int nlhs, mxArray* plhs[], int nrhs, mxArray* prhs[]) {
  return S->matlab_fgh(nlhs, plhs, nrhs, prhs);
}

bool MW_CALL_CONV
    matlab_fg(int nlhs, mxArray* plhs[], int nrhs, mxArray* prhs[]) {
  return S->matlab_fgh(nlhs, plhs, nrhs, prhs);
}

bool MW_CALL_CONV
    matlab_fgh(int nlhs, mxArray* plhs[], int nrhs, mxArray* prhs[]) {
  return S->matlab_fgh(nlhs, plhs, nrhs, prhs);
}

#endif
