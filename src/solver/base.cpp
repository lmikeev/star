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

#include "base.hpp"
#include "lib/gnuplot_i/gnuplot_i.hpp"

namespace solver {

bool base::plot(solver_loader::task_info::plot const* const ti,
                dump const* const dmp, dump const* const plt) const {
  try {
    Gnuplot gp("lp");

    gp << "set terminal png size 800,600 transparent truecolor";
    gp << "set key outside";
    gp << "set key right top";

    sprintf(tmp, "set datafile separator '%c'", csvSep);
    gp.cmd(tmp);

    sprintf(tmp, "set output '%s'", plt->get_local_path());
    gp.cmd(tmp);

    std::string tmp_;
    if (ti->get_props().get_xlabel(tmp_)) {
      gp.set_xlabel(tmp_);
    } else {
      gp.set_xlabel("t");
    }
    if (ti->get_props().get_ylabel(tmp_)) {
      gp.set_ylabel(tmp_);
    } else if (ti->get_exprs().size() == 1) {
      gp.set_ylabel(ti->get_exprs().front().expr);
    }

    char* c = tmp;
    c += sprintf(c, "plot ");
    if (ti->get_exprs().empty()) {
      for (std::size_t i = 0; i < nvars; i++) {
        if (i) {
          c += sprintf(c, "%c ", csvSep);
        }

        c += sprintf(c, "'%s' u 1:%lu every ::1 w lp t '%s'",
                     dmp->get_local_path(), i + 2,
                     sl->get_vars()[i]->get_name().c_str());
      }

      if (ti->plot_stddev()) {
        for (std::size_t i = 0; i < nvars; i++) {
          c += sprintf(
              c, ", '%s' u 1:%lu:%lu every ::1 notitle w yerrorbars ls %lu",
              dmp->get_local_path(), i + 2, nvars + i + 2, i + 1);
        }
      }
    } else {
      const std::size_t nexprs = ti->get_exprs().size();
      for (std::size_t i = 0; i < nexprs; i++) {
        auto const& e = ti->get_exprs()[i];

        if (i) {
          c += sprintf(c, "%c ", csvSep);
        }

        c += sprintf(c, "'%s' u 1:%lu every ::1 w lp", dmp->get_local_path(),
                     i + 2);

        if (e.line_props.get_title(tmp_)) {
          c += sprintf(c, " t '%s'", tmp_.c_str());
        } else {
          c += sprintf(c, " t '%s'", e.expr.c_str());
        }
      }

      if (ti->plot_stddev()) {
        for (std::size_t i = 0; i < nexprs; i++) {
          c += sprintf(
              c, ", '%s' u 1:%lu:%lu every ::1 notitle w yerrorbars ls %lu",
              dmp->get_local_path(), i + 2, nexprs + i + 2, i + 1);
        }
      }
    }
    gp.cmd(tmp);
  } catch (GnuplotException ge) {
    sl->last_error() << ge.what() << std::endl;
    return false;
  }

  return true;
}

bool base::plot_2d(solver_loader::task_info::plot_2d const* const ti,
                   dump const* const dmp, dump const* const plt) const {
  try {
    Gnuplot gp("lp");

    gp << "set terminal png size 600,600 transparent truecolor";
    gp << "set size square";
    gp << "unset key";

    sprintf(tmp, "set datafile separator '%c'", csvSep);
    gp.cmd(tmp);

    sprintf(tmp, "set output '%s'", plt->get_local_path());
    gp.cmd(tmp);

    std::string tmp_;
    if (ti->get_props().get_xlabel(tmp_)) {
      gp.set_xlabel(tmp_);
    } else {
      gp.set_xlabel(ti->get_expr1());
    }
    if (ti->get_props().get_ylabel(tmp_)) {
      gp.set_ylabel(tmp_);
    } else {
      gp.set_ylabel(ti->get_expr2());
    }

    sprintf(tmp, "plot '%s' u 1:2:3 every ::1 w lp", dmp->get_local_path());
    gp.cmd(tmp);
  } catch (GnuplotException ge) {
    sl->last_error() << ge.what() << std::endl;
    return false;
  }

  return true;
}

bool base::plot_distr(solver_loader::task_info::plot_distr const* const ti,
                      dump const* const dmp, dump const* const plt) const {
  try {
    Gnuplot gp("lp");

    gp << "set terminal png size 800,600 transparent truecolor";
    gp << "set style fill transparent solid 0.5";
    gp << "set key outside";
    gp << "set key right top";

    sprintf(tmp, "set datafile separator '%c'", csvSep);
    gp.cmd(tmp);

    sprintf(tmp, "set output '%s'", plt->get_local_path());
    gp.cmd(tmp);

    std::string tmp_;
    if (ti->get_props().get_xlabel(tmp_)) {
      gp.set_xlabel(tmp_);
    } else if (ti->get_exprs().size() == 1) {
      gp.set_xlabel(ti->get_exprs().front().expr);
    }
    if (ti->get_props().get_ylabel(tmp_)) {
      gp.set_ylabel(tmp_);
    }

    char* c = tmp;
    c += sprintf(c, "plot ");
    if (ti->get_exprs().empty()) {
      for (std::size_t i = 0; i < nvars; i++) {
        if (i) {
          c += sprintf(c, "%c ", csvSep);
        }

        c += sprintf(c, "'%s' u 1:%lu every ::1 w boxes t '%s'",
                     dmp->get_local_path(), i + 2,
                     sl->get_vars()[i]->get_name().c_str());
      }
    } else {
      for (std::size_t i = 0; i < ti->get_exprs().size(); i++) {
        auto const& e = ti->get_exprs()[i];

        if (i) {
          c += sprintf(c, "%c ", csvSep);
        }

        c += sprintf(c, "'%s' u 1:%lu every ::1 w boxes", dmp->get_local_path(),
                     i + 2);

        if (e.line_props.get_title(tmp_)) {
          c += sprintf(c, " t '%s'", tmp_.c_str());
        } else {
          c += sprintf(c, " t '%s'", e.expr.c_str());
        }
      }
    }
    gp.cmd(tmp);
  } catch (GnuplotException ge) {
    sl->last_error() << ge.what() << std::endl;
    return false;
  }

  return true;
}

bool base::plot_distr_2d(
    solver_loader::task_info::plot_distr_2d const* const ti,
    dump const* const dmp, dump const* const plt) const {
  try {
    Gnuplot gp("lp");

    gp << "set terminal png size 600,600 transparent truecolor";
    gp << "set size square";
    gp << "set auto fix";

    gp << "unset key";

    sprintf(tmp, "set datafile separator '%c'", csvSep);
    gp.cmd(tmp);

    sprintf(tmp, "set output '%s'", plt->get_local_path());
    gp.cmd(tmp);

    std::string tmp_;
    if (ti->get_props().get_xlabel(tmp_)) {
      gp.set_xlabel(tmp_);
    } else {
      gp.set_xlabel(ti->get_expr1());
    }
    if (ti->get_props().get_ylabel(tmp_)) {
      gp.set_ylabel(tmp_);
    } else {
      gp.set_ylabel(ti->get_expr2());
    }

    sprintf(tmp, "plot '%s' matrix w image", dmp->get_local_path());
    gp.cmd(tmp);
  } catch (GnuplotException ge) {
    sl->last_error() << ge.what() << std::endl;
    return false;
  }

  return true;
}

bool base::plot_objf(task const* const tsk) const {
  solver_loader::task_info::plot_objf const* const ti =
      static_cast<solver_loader::task_info::plot_objf const*>(tsk->get_info());

  dump const* const plt = tsk->get_plots().back();

  try {
    Gnuplot gp("lp");

    gp << "set terminal png size 600,600 transparent truecolor";
    gp << "set style fill transparent solid 0.5";

    sprintf(tmp, "set datafile separator '%c'", csvSep);
    gp.cmd(tmp);

    sprintf(tmp, "set output '%s'", plt->get_local_path());
    gp.cmd(tmp);

    std::string tmp_;
    if (ti->get_props().get_xlabel(tmp_)) {
      gp.set_xlabel(tmp_);
    } else if (ti->get_exprs().size() == 1) {
      gp.set_xlabel(ti->get_exprs().front().expr);
    }
    if (ti->get_props().get_ylabel(tmp_)) {
      gp.set_ylabel(tmp_);
    } else {
      gp.set_ylabel(sl->objf);
    }

    char* c = tmp;
    c += sprintf(c, "plot ");
    if (ti->get_exprs().empty()) {
      bool f = false;

      for (std::size_t i = 0; i < nparams; i++) {
        if (f) {
          c += sprintf(c, "%c ", csvSep);
        } else {
          f = true;
        }

        c += sprintf(c, "'%s' u 1:2 every ::1 w lp t '%s'",
                     tsk->get_dumps()[i]->get_local_path(),
                     sl->sa_params[i]->get_name().c_str());
      }

      for (std::size_t i = 0; i < nivars; i++) {
        if (f) {
          c += sprintf(c, "%c ", csvSep);
        } else {
          f = true;
        }

        c += sprintf(c, "'%s' u 1:2 every ::1 w lp t '%s'",
                     tsk->get_dumps()[nparams + i]->get_local_path(),
                     sl->sa_ivars[i]->get_name().c_str());
      }
    } else {
      for (std::size_t i = 0; i < ti->get_exprs().size(); i++) {
        auto const& e = ti->get_exprs()[i];

        if (i) {
          c += sprintf(c, "%c ", csvSep);
        }

        c += sprintf(c, "'%s' u 1:2 every ::1 w lp",
                     tsk->get_dumps()[i]->get_local_path());

        if (e.line_props.get_title(tmp_)) {
          c += sprintf(c, " t '%s'", tmp_.c_str());
        } else {
          c += sprintf(c, " t '%s'", e.expr.c_str());
        }
      }
    }
    gp.cmd(tmp);
  } catch (GnuplotException ge) {
    sl->last_error() << ge.what() << std::endl;
    return false;
  }

  return true;
}

bool base::plot_objf_2d(task const* const tsk) const {
  solver_loader::task_info::plot_objf_2d const* const ti =
      static_cast<solver_loader::task_info::plot_objf_2d const*>(
          tsk->get_info());

  dump const* const dmp = tsk->get_dumps().back();
  dump const* const plt = tsk->get_plots().back();

  try {
    Gnuplot gp("lp");

    gp << "set terminal png size 600,600 transparent truecolor";
    gp << "set size square";
    gp << "set auto fix";
    gp << "unset key";

    sprintf(tmp, "set datafile separator '%c'", csvSep);
    gp.cmd(tmp);

    sprintf(tmp, "set output '%s'", plt->get_local_path());
    gp.cmd(tmp);

    std::string tmp_;
    if (ti->get_props().get_xlabel(tmp_)) {
      gp.set_xlabel(tmp_);
    } else {
      gp.set_xlabel(ti->get_expr1());
    }
    if (ti->get_props().get_ylabel(tmp_)) {
      gp.set_ylabel(tmp_);
    } else {
      gp.set_ylabel(ti->get_expr2());
    }

    if (!do_sa() || sa_no_der()) {
      sprintf(tmp, "plot '%s' u 1:2:3 every ::1 w p palette pt 13",
              dmp->get_local_path());
    } else {
      sprintf(tmp, "plot '%s' u 1:2:4:5:3 every ::1 w vec filled lc palette",
              dmp->get_local_path());
    }
    gp.cmd(tmp);
  } catch (GnuplotException ge) {
    sl->last_error() << ge.what() << std::endl;
    return false;
  }

  return true;
}
}
