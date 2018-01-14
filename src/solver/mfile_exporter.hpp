/*
 *  mfile_exporter.hpp
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

#ifndef SOLVER_MFILE_EXPORTER_HPP_
#define SOLVER_MFILE_EXPORTER_HPP_

#include "../solver_loader/parser/parser.hpp"
#include "../solver_loader/writer/matlab.hpp"

namespace solver {

template <class lstate_data_accessor_t = lstate_data_accessor,
          class hstate_data_accessor_t = hstate_data_accessor,
          class hstate_succ_data_accessor_t = hstate_succ_data_accessor,
          typename lstate_t = lstate_h, typename hstate_t = hstate_h,
          typename hstate_succ_t = hstate_succ_h>
class mfileExporterT {
 public:
  typedef dataAccessorT<lstate_data_accessor_t, hstate_data_accessor_t,
                        hstate_succ_data_accessor_t, lstate_t, hstate_t,
                        hstate_succ_t> dataAccessor;

  typedef container::abstractStateHashListT<lstate_t, hstate_t, hstate_succ_t>
      abstractStateHashList;

  typedef lstate_t state;
  typedef hstate_t hstate;
  typedef hstate_succ_t hstate_succ;

  mfileExporterT(options const* const o, dataAccessor const* const da,
                 abstractStateHashList const* const states)
      : o(o),
        da(da),
        states(states),
        sI(o->nmoments),
        sIk(o->nmoments),
        sJ(o->nmoments),
        sJk(o->nmoments),
        sK(o->nmoments),
        sI_(o->nvars),
        sIk_(o->nvars),
        sJ_(o->nvars),
        sK_(o->nvars) {}

  bool write_mfile(char const* const file_name) const {
    if (!o->hs_has_index()) {
      std::cerr << "no indexing" << std::endl;
      return false;
    }

    std::ofstream os(file_name);

    const std::size_t nmodes = states->get_nactive();

    os << std::endl;
    os << " % " << std::endl;
    os << " %  The state of the system is stored as" << std::endl;
    os << " % [ p_1 x_1^1 .. x_1^N C_1^11 .. C_1^NN C_1^111 .. C_1^NNN .. "
       << std::endl;
    os << " %   p_M x_M^1 .. x_M^N C_M^11 .. C_M^NN C_M^111 .. C_M^NNN .. ]"
       << std::endl;
    os << " % " << std::endl;
    os << std::endl;
    os << " % " << std::endl;
    os << " % T - time points" << std::endl;
    os << " % Y - solutions" << std::endl;
    os << " % E - expectations" << std::endl;
    os << " % " << std::endl;
    os << std::endl;
    os << "function [T,Y,E] = output" << std::endl;
    os << std::endl;
    os << "  nreactions = " << o->ntransition_classes << ";" << std::endl;
    os << "  nvars = " << o->nvars << ";" << std::endl;
    if (!o->is_stoch()) {
      os << "  npvars = " << o->npvars << ";" << std::endl;
      os << "  nmoments = " << o->nmoments << ";" << std::endl;
    }
    os << "  nmodes = " << nmodes << ";" << std::endl;
    os << "  nmeqs = " << o->s_y_len << ";" << std::endl;
    os << "  neqs = nmodes * nmeqs;" << std::endl;
    os << "  rates_g_len = " << o->tr_rates_g_len << ";" << std::endl;
    os << "  rates_h_len = " << o->tr_rates_h_len << ";" << std::endl;
    os << std::endl;
    os << "  modes = [";
    for (state* s = states->get_first(); s != nullptr; s = s->lnext) {
      if (s != states->get_first()) {
        os << "    ";
      }

      cstate const* const cs = da->hs_cs(s->hs);
      int const* const v = da->cs_get(cs);
      for (std::size_t i = 0; i < o->nvars; i++) {
        if (i < o->npvars) {
          os << " 0";
        } else {
          os << " " << v[i];
        }
      }

      if (s->lnext != nullptr) {
        os << " ;" << std::endl;
      }
    }
    os << " ];" << std::endl;

    if (!o->is_stoch()) {
      os << std::endl;
      os << "  dvar = transpose([" << std::endl;
      for (std::size_t i = 0; i < o->ntransition_classes; i++) {
        os << "    ";

        double const* const ch =
            o->sl->get_model()->get_transitions()[i]->change();
        for (std::size_t j = 0; j < o->npvars; j++) {
          os << " " << ch[j];
        }

        if (i < o->ntransition_classes - 1) {
          os << " ;" << std::endl;
        }
      }
      os << " ]);" << std::endl;
    }

    os << std::endl;

    solver_loader::writer::base const* const matlab_writer =
        new solver_loader::writer::matlab();
    std::vector<std::vector<solver_loader::parser::AST const*>> rates_g(
        o->ntransition_classes,
        std::vector<solver_loader::parser::AST const*>(
            o->tr_rates_g_len, solver_loader::parser::AST::et_0));
    std::vector<std::vector<solver_loader::parser::AST const*>> rates_h(
        o->ntransition_classes,
        std::vector<solver_loader::parser::AST const*>(
            o->tr_rates_h_len, solver_loader::parser::AST::et_0));
    std::set<solver_loader::model_info::var const*> ev;
    std::set<solver_loader::model_info::cnst const*> ep;

    for (std::size_t i = 0; i < o->ntransition_classes; i++) {
      solver_loader::model_info::transition const* const tr =
          o->sl->get_transitions()[i];
      solver_loader::parser::AST const* const rate =
          solver_loader::parser::build_expr(tr->get_rate(), tr);

      rates_g[i][0] = o->sl->build_rate_g(rate);
      o->sl->expr_vars_params(rates_g[i][0], ev, ep);

      rates_h[i][0] = o->sl->build_rate_h(rate);
      o->sl->expr_vars_params(rates_h[i][0], ev, ep);

      if (o->nmoments > 1) {
        for (std::size_t j = 0; j < o->npvars; j++) {
          solver_loader::parser::AST const* const drate_h_dx =
              o->sl->build_expr_der(rates_h[i][0], o->sl->p_vars[j]);

          sI[0] = j;
          sI_[j]++;

          assert(sI_[j] == 1);

          if (drate_h_dx != solver_loader::parser::AST::et_0) {
            rates_h[i][1 + j] = drate_h_dx;

            fill_drate_h_dx(rates_h[i], sI.data(), sI_.data(), 1, 1 + o->npvars,
                            drate_h_dx);
          }

          sI_[j]--;

          assert(sI_[j] == 0);
        }
      }
    }

    os << "  function [g,h] = get_rates(m,x)" << std::endl;
    if (!o->sl->get_cnsts().empty()) {
      os << std::endl;
      for (auto const& c : o->sl->get_cnsts()) {
        os << "    c" << c->get_index() << " = ";
        if (c->get_type()->is_boolean()) {
          os << (static_cast<bool>(c->get_value()->get()) ? "true" : "false");
        } else if (c->get_type()->is_integer()) {
          os << static_cast<int>(c->get_value()->get());
        } else {
          os << static_cast<double>(c->get_value()->get());
        }
        os << ";   % " << c->get_name() << std::endl;
      }
    }
    if (!ev.empty()) {
      os << std::endl;
      for (auto const& v : ev) {
        const std::size_t vi = v->get_index();
        os << "    v" << vi << " = ";
        if (vi < o->npvars) {
          os << "x(" << vi + 1 << ")";
        } else {
          os << "m(" << vi + 1 << ")";
        }
        os << ";   % " << v->get_name() << std::endl;
      }
    }
    os << std::endl;
    os << "    g = zeros(nreactions,rates_g_len);" << std::endl;
    os << std::endl;
    for (std::size_t i = 0; i < o->ntransition_classes; i++) {
      for (std::size_t j = 0; j < o->tr_rates_g_len; j++) {
        if (rates_g[i][j] != solver_loader::parser::AST::et_0) {
          os << "    g(" << i + 1 << "," << j + 1 << ") = ";
          rates_g[i][j]->write(os, matlab_writer) << ";" << std::endl;
        }
      }
    }
    os << std::endl;
    os << "    h = zeros(nreactions,rates_h_len);" << std::endl;
    for (std::size_t i = 0; i < o->ntransition_classes; i++) {
      if (o->nmoments > 1) {
        os << std::endl;
      }
      for (std::size_t j = 0; j < o->tr_rates_h_len; j++) {
        if (rates_h[i][j] != solver_loader::parser::AST::et_0) {
          os << "    h(" << i + 1 << "," << j + 1 << ") = ";
          rates_h[i][j]->write(os, matlab_writer) << ";" << std::endl;
        }
      }
    }
    os << std::endl;
    os << "  end" << std::endl;

    os << std::endl;
    os << "  tspan = [";
    for (auto ti : o->sl->tspan) {
      os << " " << ti;
    }
    os << " ];" << std::endl;
    os << std::endl;
    os << "  y0 = zeros(neqs,1);" << std::endl;
    for (state* s = states->get_first(); s != nullptr; s = s->lnext) {
      hstate const* const hs = s->hs;
      void const* const hs_d = da->get_hstate_data(hs);
      const std::size_t s_i = da->hs_d_index(hs_d);

      const double p = da->s_p(s);
      double const* const x = da->s_x(s, p);

      const std::size_t mpi = s_i * o->s_y_len + 1;
      if (p > 0.0) {
        os << "  y0(" << mpi << ") = " << p << ";" << std::endl;
      }

      for (std::size_t i = 0; i < o->npvars; i++) {
        if (x[i] > std::numeric_limits<double>::epsilon()) {
          const std::size_t mvi = s_i * o->s_y_len + 2 + i;
          os << "  y0(" << mvi << ") = " << x[i] << ";" << std::endl;
        }
      }
    }

    if (o->is_stoch()) {
      os << std::endl;
      os << "  Qi = [];" << std::endl;
      os << "  Qj = [];" << std::endl;
      os << "  Qs = [];" << std::endl;
      os << "  x = [];" << std::endl;
      for (state* s = this->states->get_first(); s != nullptr; s = s->lnext) {
        hstate const* const hs = s->hs;
        if (!da->hs_is_unexplored(hs) && hs->succ != nullptr) {
          void const* const hs_d = da->get_hstate_data(hs);
          const std::size_t s_i = da->hs_d_index(hs_d);
          os << std::endl;
          os << "  [rates_g,~] = get_rates(modes(" << s_i + 1 << ",:),x);"
             << std::endl;
          hstate_succ const* hsucc = hs->succ;
          while (hsucc != nullptr) {
            const std::size_t tri = da->hsucc_tri(hsucc);
            const std::size_t succ_i = da->hs_index(hsucc->hs);
            os << "  Qi = [Qi," << succ_i + 1 << "];" << std::endl;
            os << "  Qj = [Qj," << s_i + 1 << "];" << std::endl;
            os << "  Qs = [Qs,rates_g(" << tri + 1 << ")];" << std::endl;
            hsucc = hsucc->lnext;
          }
          os << "  Qi = [Qi," << s_i + 1 << "];" << std::endl;
          os << "  Qj = [Qj," << s_i + 1 << "];" << std::endl;
          os << "  Qs = [Qs,-sum(rates_g)];" << std::endl;
        }
      }
      os << std::endl;
      os << "  Qt = sparse(Qi,Qj,Qs,length(y0),length(y0));" << std::endl;
    }

    os << std::endl;
    os << "  options = odeset('Stats','on', 'Refine',1, 'AbsTol',"
       << o->sl->abs_tol << ", 'RelTol'," << o->sl->rel_tol << ");"
       << std::endl;

    const bool jacobian = false;
    if (jacobian) {
      os << "  options = odeset(options, 'Jacobian',@FJac);" << std::endl;
    }

    os << std::endl;
    os << "  [T,Y] = ode23s(@yprime, tspan, y0, options);" << std::endl;

    os << std::endl;
    os << "   % compute expectations" << std::endl;
    os << "  E = zeros(length(T),nvars);" << std::endl;

    os << "  for ti=1:length(T)" << std::endl;
    os << "    for i=1:nmodes" << std::endl;
    os << "      mp = Y(ti,(i-1)*nmeqs + 1);" << std::endl;
    os << "      mpx = zeros(1,nvars);" << std::endl;
    for (std::size_t j = 0; j < o->nvars; j++) {
      os << "      mpx(" << j + 1 << ") = ";
      if (j < o->npvars) {
        os << "Y(ti, (i-1)*nmeqs + 2 + " << j << ")";
      } else {
        os << "mp * modes(i," << 1 + j << ")";
      }
      os << ";" << std::endl;
    }
    os << "      E(ti,:) = E(ti,:) + mpx;" << std::endl;
    os << "    end" << std::endl;
    os << "  end" << std::endl;
    os << std::endl;
    os << "  function dy = yprime(~,y)" << std::endl;
    os << std::endl;
    if (!o->is_stoch()) {
      os << "    dy = zeros(size(y));" << std::endl;

      for (state* s = states->get_first(); s != nullptr; s = s->lnext) {
        hstate const* const hs = s->hs;
        void const* const hs_d = da->get_hstate_data(hs);
        cstate const* const cs = da->hs_d_cs(hs_d);
        const std::size_t s_i = da->hs_d_index(hs_d);
        const std::size_t mi = s_i + 1;
        const std::size_t s_pi = s_i * o->s_y_len + 1;

        os << std::endl;
        os << "     % MODE " << mi << std::endl;
        os << "    mp = y(" << s_pi << "); % mode prob." << std::endl;
        os << "    if (mp > 1e-16)" << std::endl;
        os << "      m = modes(" << mi << ",:);" << std::endl;
        os << "      x = y(" << s_pi + 1 << ":" << s_pi + o->s_x_len
           << ") / mp;" << std::endl;
        os << std::endl;
        os << "      [rates_g,rates_h] = get_rates(m,x);" << std::endl;

        if (!da->hs_is_unexplored(hs)) {
          hstate_succ const* hsucc = hs->succ;
          while (hsucc != nullptr) {
            const std::size_t tri = da->hsucc_tri(hsucc);
            const std::size_t succ_i = da->hs_index(hsucc->hs);
            const std::size_t succ_pi = succ_i * o->s_y_len + 1;

            os << std::endl;
            os << "        % -> mode " << succ_i + 1 << std::endl;

            write_comp_succ(os, s_pi, succ_pi, tri, rates_h[tri]);

            hsucc = hsucc->lnext;
          }
        }

        if (!o->is_stoch()) {
          if (o->is_det_centered()) {
            write_comp_c(os, cs, s_pi, rates_h);
          } else {
            write_comp_c_raw(os, cs, s_pi, rates_h);
          }
        }

        os << "    end" << std::endl;
      }

      if (!o->is_stoch() && o->is_det_centered() && o->nmoments > 1) {
        for (state* s = states->get_first(); s != nullptr; s = s->lnext) {
          hstate const* const hs = s->hs;
          void const* const hs_d = da->get_hstate_data(hs);
          cstate const* const cs = da->hs_d_cs(hs_d);
          const std::size_t s_i = da->hs_d_index(hs_d);
          const std::size_t mi = s_i + 1;
          const std::size_t s_pi = s_i * o->s_y_len + 1;

          os << std::endl;
          os << "     % MODE " << mi << " [C]" << std::endl;
          os << "    mp = y(" << s_pi << "); % mode prob." << std::endl;
          os << "    if (mp > 1e-16)" << std::endl;
          os << "      m = modes(" << mi << ",:);" << std::endl;
          os << "      x = y(" << s_pi + 1 << ":" << s_pi + o->s_x_len
             << ") / mp;" << std::endl;
          os << std::endl;
          os << "      [rates_g,rates_h] = get_rates(m,x);" << std::endl;

          write_comp_c_cov(os, cs, s_pi, rates_h);

          if (o->is_hybrid() && !da->hs_is_unexplored(hs)) {
            hstate_succ const* hsucc = hs->succ;
            while (hsucc != nullptr) {
              const std::size_t tri = da->hsucc_tri(hsucc);
              const std::size_t succ_i = da->hs_index(hsucc->hs);
              const std::size_t succ_pi = succ_i * o->s_y_len + 1;

              os << std::endl;
              os << "       % -> mode " << succ_i + 1 << std::endl;
              os << "      smp = y(" << succ_pi << ");" << std::endl;
              os << "      if (smp > 1e-16)" << std::endl;
              os << "        dmu = x(1:" << o->npvars << ") - y(" << succ_pi + 1
                 << ":" << succ_pi + o->npvars << ") / smp;" << std::endl;
              os << "      else" << std::endl;
              os << "        dmu = x(1:" << o->npvars << ");" << std::endl;
              os << "      end" << std::endl;

              write_comp_succ_cov(os, s_pi, succ_pi, tri, rates_h[tri]);

              hsucc = hsucc->lnext;
            }
          }

          os << "    end" << std::endl;
        }
      }
    } else {
      os << "    dy = Qt * y;" << std::endl;
    }
    os << std::endl;
    os << "  end" << std::endl;
    os << std::endl;
    os << "end" << std::endl;
    os.close();

    return true;
  }

 private:
  options const* const o;
  dataAccessor const* const da;
  abstractStateHashList const* const states;

  mutable std::vector<std::size_t> sI, sIk, sJ, sJk, sK;
  mutable std::vector<std::size_t> sI_, sIk_, sJ_, sK_;

 protected:
  void fill_drate_h_dx(
      std::vector<solver_loader::parser::AST const*>& rates_h,
      std::size_t* const I, std::size_t* const I_, const std::size_t d,
      const std::size_t ci0,
      solver_loader::parser::AST const* const drate_h_dx) const {
    const std::size_t cov_len_ = o->covlen(o->npvars, d + 1);

    for (std::size_t i = 0; i <= I[d - 1]; i++) {
      I[d] = i;
      I_[i]++;

      solver_loader::parser::AST const* const d2rate_h_dx2 =
          o->sl->build_expr_der(drate_h_dx, o->sl->p_vars[i]);
      if (d2rate_h_dx2 != solver_loader::parser::AST::et_0) {
        rates_h[ci0 + o->cov_indexl(I, d + 1)] = d2rate_h_dx2;

        if (d + 1 < o->nmoments) {
          fill_drate_h_dx(rates_h, I, I_, d + 1, ci0 + cov_len_, d2rate_h_dx2);
        }
      }

      I_[i]--;
    }
  }

  std::ostream& write_ch(std::ostream& os, const double ch) const {
    if (ch > 0.0) {
      os << "+ " << ch;
    } else {
      os << "- " << -ch;
    }
    return os;
  }

  void write_ez_(std::ostream& os, const std::size_t tri,
                 const std::vector<solver_loader::parser::AST const*>& rates_h,
                 std::size_t* const K_, std::size_t* const K,
                 std::size_t* const Jk_, const std::size_t d,
                 const std::size_t kci0, const std::size_t ci0, bool& p_,
                 bool* nz = nullptr) const {
    std::size_t Jk_len;
    o->cov_I_I(o->npvars, sJk.data(), Jk_len, Jk_);

    if (Jk_len != 1) {
      const std::size_t kci = kci0 + o->cov_indexl(K, d);

      if (rates_h[kci] != solver_loader::parser::AST::et_0) {
        double c = 1.0;
        for (std::size_t i = 0; i < o->npvars; i++) {
          c *= o->_1_fac[K_[i]];
        }

        if (!Jk_len) {
          if (nz == nullptr) {
            if (p_) {
              os << " ..." << std::endl << "        + ";
            }

            os << c << " * rates_h(" << tri + 1 << "," << kci + 1 << ")";

            p_ = true;
          } else {
            *nz = true;
          }
        } else if (Jk_len > 1) {
          if (nz == nullptr) {
            if (p_) {
              os << " ..." << std::endl << "        + ";
            }

            const std::size_t ci = ci0 + o->cov_indexl(sJk.data(), Jk_len);

            os << c << " * rates_h(" << tri + 1 << "," << kci + 1 << ") * x("
               << ci + 1 << ")";

            p_ = true;
          } else {
            *nz = true;
          }
        }
      }
    }

    if (Jk_len + 1 <= o->nmoments) {
      const std::size_t kcov_len_ = o->covlen(o->npvars, d);
      const std::size_t cov_len_ = Jk_len ? o->covlen(o->npvars, Jk_len) : 0;

      for (std::size_t i = 0; i <= K[d - 1]; i++) {
        K[d] = i;
        K_[i]++;

        Jk_[i]++;

        write_ez_(os, tri, rates_h, K_, K, Jk_, d + 1, kci0 + kcov_len_,
                  ci0 + cov_len_, p_, nz);

        Jk_[i]--;

        K_[i]--;
      }
    }
  }

  void write_ez_(std::ostream& os, const std::size_t tri,
                 const std::vector<solver_loader::parser::AST const*>& rates_h,
                 std::size_t* const J_, std::size_t* const J,
                 const std::size_t J_len, bool* nz = nullptr) const {
    std::size_t ci0 = J_len ? o->cov_index0(o->npvars, J_len) : 0;
    bool p_ = false;
    if (!J_len) {
      if (nz == nullptr) {
        if (p_) {
          os << " ..." << std::endl << "        + ";
        }

        os << "rates_h(" << tri + 1 << ",1)";

        p_ = true;
      } else {
        *nz = true;
      }
    } else if (J_len > 1) {
      if (nz == nullptr) {
        if (p_) {
          os << " ..." << std::endl << "        + ";
        }

        const std::size_t ci = ci0 + o->cov_indexl(J, J_len);

        os << "rates_h(" << tri + 1 << ",1) * x(" << ci + 1 << ")";

        p_ = true;
      } else {
        *nz = true;
      }
    }

    if (J_len + 1 <= o->nmoments) {
      ci0 += (J_len ? o->covlen(o->npvars, J_len) : 0);

      for (std::size_t i = 0; i < o->npvars; i++) {
        sK[0] = i;
        sK_[i]++;

        assert(sK_[i] == 1);

        J_[i]++;

        write_ez_(os, tri, rates_h, sK_.data(), sK.data(), J_, 1, 1, ci0, p_,
                  nz);

        J_[i]--;

        sK_[i]--;

        assert(sK_[i] == 0);
      }
    }

    if (nz == nullptr) {
      os << ";" << std::endl;
    }
  }

  void write_ez(std::ostream& os, const std::size_t tri,
                const std::vector<solver_loader::parser::AST const*>& rates_h,
                std::size_t* const J_, std::size_t* const J,
                const std::size_t J_len, bool* nz = nullptr) const {
    if (nz == nullptr) {
      os << "      ez = ";
    }
    write_ez_(os, tri, rates_h, J_, J, J_len, nz);
  }

  void write_ez_(std::ostream& os, const std::size_t tri,
                 const std::vector<solver_loader::parser::AST const*>& rates_h,
                 std::size_t* const J_, bool* nz = nullptr) const {
    std::size_t J_len;
    o->cov_I_I(o->npvars, sJ.data(), J_len, J_);
    write_ez_(os, tri, rates_h, J_, sJ.data(), J_len, nz);
  }

  void write_ez(std::ostream& os, const std::size_t tri,
                const std::vector<solver_loader::parser::AST const*>& rates_h,
                std::size_t* const J_, bool* nz = nullptr) const {
    if (nz == nullptr) {
      os << "      ez = ";
    }
    write_ez_(os, tri, rates_h, J_, nz);
  }

  void write_eh(std::ostream& os, const std::size_t tri,
                const std::vector<solver_loader::parser::AST const*>& rates_h,
                bool* nz = nullptr) const {
    if (nz == nullptr) {
      os << "      eh = ";
    }
    write_ez_(os, tri, rates_h, sJ_.data(), nz);
  }

  void write_cov_z(
      std::ostream& os, const std::size_t tri,
      const std::vector<solver_loader::parser::AST const*>& rates_h,
      double const* const ch, std::size_t const* const I_,
      std::size_t* const J_, const std::size_t d, const bool cnt,
      const bool fmx, bool* nz = nullptr) const {
    assert(J_[d] == 0);

    const std::size_t I_mx = I_[d] + ((!fmx || d + 1 < o->npvars) ? 1 : 0);

    for (std::size_t i = 0; i < I_mx; i++) {
      J_[d] = i;

      if (d + 1 < o->npvars) {
        write_cov_z(os, tri, rates_h, ch, I_, J_, d + 1, cnt, fmx && i == I_[d],
                    nz);
      } else {
        bool nz_ = false;
        write_ez(os, tri, rates_h, J_, &nz_);

        if (nz_) {
          double k = 1.0;
          for (std::size_t j = 0; j < o->npvars; j++) {
            k *= o->cmb[I_[j]][J_[j]];
          }

          if (cnt) {
            for (std::size_t j = 0; j < o->npvars; j++) {
              const double ch_ = ch[j];
              for (std::size_t l = 0; l < I_[j] - J_[j]; l++) {
                k *= ch_;
              }
            }
          }

          if (std::fabs(k) > 0.0) {
            if (nz == nullptr) {
              os << std::endl;
              os << "         % k = (";
              for (std::size_t j = 0; j < o->npvars; j++) {
                os << " " << J_[j];
              }
              os << " )" << std::endl;

              write_ez(os, tri, rates_h, J_);

              os << "      z = z ";
              write_ch(os, k) << " * ez";

              if (!cnt) {
                if (o->npvars > 2) {
                  os << " * prod((dmu + dvar(:," << tri + 1 << ")) .^ [";
                  for (std::size_t j = 0; j < o->npvars; j++) {
                    if (j) {
                      os << ";";
                    }
                    os << " " << I_[j] - J_[j];
                  }
                  os << " ])";
                } else {
                  for (std::size_t j = 0; j < o->npvars; j++) {
                    const std::size_t d_ = I_[j] - J_[j];
                    if (d_) {
                      const double ch_ = ch[j];

                      os << " * ";
                      if (std::fabs(ch_) > 0.0) {
                        os << "(";
                      }
                      os << "dmu(" << j + 1 << ") ";
                      if (std::fabs(ch_) > 0.0) {
                        write_ch(os, ch_) << ")";
                      }
                      if (d_ > 1) {
                        os << "^" << d_;
                      }
                    }
                  }
                }
              }

              os << ";" << std::endl;
            } else {
              *nz = true;
            }
          }
        }
      }
    }
    J_[d] = 0;
  }

  void write_dcov(std::ostream& os, const std::size_t s_pi,
                  const std::size_t succ_pi, const std::size_t tri,
                  const std::vector<solver_loader::parser::AST const*>& rates_h,
                  double const* const ch, std::size_t* const I,
                  std::size_t* const I_, const std::size_t d,
                  const std::size_t ci0) const {
    const std::size_t cov_len_ = o->covlen(o->npvars, d + 1);

    for (std::size_t i = 0; i <= I[d - 1]; i++) {
      I[d] = i;
      I_[i]++;

      bool nz1 = false;
      write_ez(os, tri, rates_h, I_, I, d + 1, &nz1);

      bool nz2 = false;
      write_cov_z(os, tri, rates_h, ch, I_, sJ_.data(), 0, false, false, &nz2);

      if (nz1 || nz2) {
        os << std::endl;
        os << "        % C ^ (";
        for (std::size_t j = 0; j < o->npvars; j++) {
          os << " " << I_[j];
        }
        os << " ) [s]" << std::endl;

        const std::size_t ci = ci0 + o->cov_indexl(I, d + 1);

        if (nz1) {
          write_ez(os, tri, rates_h, I_, I, d + 1);
          os << "      dy(" << s_pi + 1 + ci << ") = dy(" << s_pi + 1 + ci
             << ") - mp * rates_g(" << tri + 1 << ",1) * ez;" << std::endl;
        }

        if (nz2) {
          os << "      z = 0;" << std::endl;
          write_cov_z(os, tri, rates_h, ch, I_, sJ_.data(), 0, false, false);

          os << "      dy(" << succ_pi + 1 + ci << ") = dy(" << succ_pi + 1 + ci
             << ") + mp * rates_g(" << tri + 1 << ",1) * z;" << std::endl;
        }
      }

      if (d + 1 < o->nmoments) {
        write_dcov(os, s_pi, succ_pi, tri, rates_h, ch, I, I_, d + 1,
                   ci0 + cov_len_);
      }

      I_[i]--;
    }
  }

  void write_dcov_c(
      std::ostream& os, const std::size_t s_pi, const std::size_t tri,
      const std::vector<solver_loader::parser::AST const*>& rates_h,
      double const* const ch, std::size_t* const I, std::size_t* const I_,
      const std::size_t d, const std::size_t ci0) const {
    const std::size_t cov_len_ = o->covlen(o->npvars, d + 1);

    for (std::size_t i = 0; i <= I[d - 1]; i++) {
      I[d] = i;
      I_[i]++;

      bool nz = false;
      write_cov_z(os, tri, rates_h, ch, I_, sJ_.data(), 0, true, true, &nz);

      if (nz) {
        os << std::endl;
        os << "        % C ^ (";
        for (std::size_t j = 0; j < o->npvars; j++) {
          os << " " << I_[j];
        }
        os << " ) [c]" << std::endl;
        os << "      z = 0;" << std::endl;
        write_cov_z(os, tri, rates_h, ch, I_, sJ_.data(), 0, true, true);

        const std::size_t ci = ci0 + o->cov_indexl(I, d + 1);

        os << "      dy(" << s_pi + 1 + ci << ") = dy(" << s_pi + 1 + ci
           << ") + mp * rates_g(" << tri + 1 << ",1) * z;" << std::endl;
      }

      if (d + 1 < o->nmoments) {
        write_dcov_c(os, s_pi, tri, rates_h, ch, I, I_, d + 1, ci0 + cov_len_);
      }

      I_[i]--;
    }
  }

  void write_dcov_2(std::ostream& os, const std::size_t s_pi,
                    const std::size_t vi, std::size_t* const I,
                    std::size_t* const I_, const std::size_t d,
                    const std::size_t ci0) const {
    const std::size_t cov_len_ = o->covlen(o->npvars, d + 1);

    for (std::size_t i = 0; i <= I[d - 1]; i++) {
      I[d] = i;
      I_[i]++;

      if (I_[vi]) {
        I_[vi]--;

        std::size_t Ik_len;
        o->cov_I_I(o->npvars, sIk.data(), Ik_len, I_);

        I_[vi]++;

        if (Ik_len != 1) {
          const std::size_t ci = ci0 + o->cov_indexl(I, d + 1);

          os << std::endl;
          os << "        % C ^ (";
          for (std::size_t j = 0; j < o->npvars; j++) {
            os << " " << I_[j];
          }
          os << " ) [2]" << std::endl;
          os << "      dy(" << s_pi + 1 + ci << ") = dy(" << s_pi + 1 + ci
             << ") - " << I_[vi] << " * (dy(" << s_pi + 1 + vi << ") - x("
             << vi + 1 << ") * dy(" << s_pi << "))";

          if (Ik_len > 1) {
            os << " * x(" << o->cov_indexg(o->npvars, sIk.data(), Ik_len) + 1
               << ")";
          }

          os << ";" << std::endl;
        }
      }

      if (d + 1 < o->nmoments) {
        write_dcov_2(os, s_pi, vi, I, I_, d + 1, ci0 + cov_len_);
      }

      I_[i]--;
    }
  }

  void write_comp_succ(
      std::ostream& os, const std::size_t s_pi, const std::size_t succ_pi,
      const std::size_t tri,
      const std::vector<solver_loader::parser::AST const*>& rates_h) const {
    write_eh(os, tri, rates_h);
    os << "        % p [s]" << std::endl;
    os << "      dp = mp * rates_g(" << tri + 1 << ",1) * eh;" << std::endl;
    os << "      dy(" << s_pi << ") = dy(" << s_pi << ") - dp;" << std::endl;
    os << "      dy(" << succ_pi << ") = dy(" << succ_pi << ") + dp;"
       << std::endl;

    if (o->is_hybrid()) {
      double const* const ch =
          o->sl->get_model()->get_transitions()[tri]->change();

      for (std::size_t i = 0; i < o->npvars; i++) {
        sI[0] = i;
        sI_[i]++;

        assert(sI_[i] == 1);

        bool nz = false;
        write_ez(os, tri, rates_h, sI_.data(), sI.data(), 1, &nz);

        os << std::endl;
        os << "        % x" << 1 + i << " [s]" << std::endl;
        if (nz) {
          write_ez(os, tri, rates_h, sI_.data(), sI.data(), 1);
          os << "      dpx = mp * rates_g(" << tri + 1 << ",1) * (ez + x("
             << 1 + i << ") * eh)";
        } else {
          os << "      dpx = mp * rates_g(" << tri + 1 << ",1) * x(" << 1 + i
             << ") * eh";
        }
        os << ";" << std::endl;

        os << "      dy(" << s_pi + 1 + i << ") = dy(" << s_pi + 1 + i
           << ") - dpx;" << std::endl;
        os << "      dy(" << succ_pi + 1 + i << ") = dy(" << succ_pi + 1 + i
           << ") + dpx";
        const double ch_ = ch[i];
        if (std::fabs(ch_) > 0.0) {
          os << " ";
          write_ch(os, ch_) << " * mp * rates_g(" << tri + 1 << ",1) * eh";
        }
        os << ";" << std::endl;

        sI_[i]--;

        assert(sI_[i] == 0);
      }
    }
  }

  void write_comp_c(
      std::ostream& os, cstate const* const cs, const std::size_t s_pi,
      const std::vector<std::vector<solver_loader::parser::AST const*>>&
          rates_h) const {
    for (auto const& tr : o->trs_det) {
      if (tr->is_enabled(cs)) {
        const std::size_t tri = tr->get_index();
        double const* const ch = tr->change();

        bool nz = false;
        write_eh(os, tri, rates_h[tri], &nz);

        if (nz) {
          os << std::endl;
          write_eh(os, tri, rates_h[tri]);
          os << "      dpx = mp * rates_g(" << tri + 1 << ",1) * eh;"
             << std::endl;

          for (std::size_t i = 0; i < o->npvars; i++) {
            const double ch_ = ch[i];
            if (std::fabs(ch_) > 0.0) {
              os << "        % x" << 1 + i << " [c]" << std::endl;
              os << "      dy(" << s_pi + 1 + i << ") = dy(" << s_pi + 1 + i
                 << ") ";
              write_ch(os, ch_) << " * dpx;" << std::endl;
            }
          }
        }
      }
    }
  }

  void write_comp_c_raw(
      std::ostream&, cstate const* const, const std::size_t,
      const std::vector<std::vector<solver_loader::parser::AST const*>>&)
      const {}

  void write_comp_succ_cov(
      std::ostream& os, const std::size_t s_pi, const std::size_t succ_pi,
      const std::size_t tri,
      const std::vector<solver_loader::parser::AST const*>& rates_h) const {
    double const* const ch =
        o->sl->get_model()->get_transitions()[tri]->change();

    for (std::size_t i = 0; i < o->npvars; i++) {
      sI[0] = i;
      sI_[i]++;

      assert(sI_[i] == 1);

      write_dcov(os, s_pi, succ_pi, tri, rates_h, ch, sI.data(), sI_.data(), 1,
                 o->npvars);

      sI_[i]--;

      assert(sI_[i] == 0);
    }
  }

  void write_comp_c_cov(
      std::ostream& os, cstate const* const cs, const std::size_t s_pi,
      const std::vector<std::vector<solver_loader::parser::AST const*>>&
          rates_h) const {
    for (std::size_t i = 0; i < o->npvars; i++) {
      sI[0] = i;
      sI_[i]++;

      assert(sI_[i] == 1);

      for (std::size_t j = 0; j < o->npvars; j++) {
        write_dcov_2(os, s_pi, j, sI.data(), sI_.data(), 1, o->npvars);
      }

      for (auto const& tr : o->trs_det) {
        if (tr->is_enabled(cs)) {
          const std::size_t tri = tr->get_index();
          double const* const ch = tr->change();

          write_dcov_c(os, s_pi, tri, rates_h[tri], ch, sI.data(), sI_.data(),
                       1, o->npvars);
        }
      }

      sI_[i]--;

      assert(sI_[i] == 0);
    }
  }
};
}

#endif
