/*
 *  state_hash_list.hpp
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

#ifndef SOLVER_STATE_HASH_LIST_HPP_
#define SOLVER_STATE_HASH_LIST_HPP_

#include "container/abstract_state_hash_list.hpp"
#include "util/bjenkins_lookup3_wrapper.hpp"
#include "task.hpp"
#include "options.hpp"

namespace solver {

class lstate_data_accessor : public container::data_accessor {
 protected:
  enum lstate_data_offset { do_y = 0 };

  const std::size_t len;
  const std::size_t size_;

 public:
  lstate_data_accessor(options const* const o, const std::size_t n = 1)
      : container::data_accessor(o, n * o->get_s_y_len() * sizeof(double)),
        len(o->get_s_y_len()),
        size_(len * sizeof(double)) {}

  double* y(void* const d) const { return static_cast<double*>(d); }

  double const* y(void const* const d) const {
    return static_cast<double const*>(d);
  }

  virtual void zero(void* const d) const { memset(y(d), 0, size_); }

  virtual double* f0(void* const) const { return nullptr; }

  virtual double const* f0(void const* const) const { return nullptr; }
};

struct hstate_data_accessor_helper {
  cstate cs;
  int foo;
};

class hstate_data_accessor : public container::hs_data_accessor {
 private:
  const std::size_t cs_offset;
  const std::size_t exitrate_offset;
  const std::size_t index_offset;
  const std::size_t flags_offset;

  const std::size_t cs_size;

 public:
  hstate_data_accessor(options const* const o)
      : container::hs_data_accessor(
            o, (o->hs_has_exitrate() ? sizeof(double) : 0) +
#ifdef STAR_CODEGEN
                   offsetof(hstate_data_accessor_helper, foo)
#else
                   TODO
#endif
                   + (o->hs_has_index() ? sizeof(std::size_t) : 0) +
                   (o->hs_has_flags() ? sizeof(std::size_t) : 0)),

        cs_offset(0),

        exitrate_offset(cs_offset +
#ifdef STAR_CODEGEN
                        offsetof(hstate_data_accessor_helper, foo)
#else
                        TODO
#endif
                        ),

        index_offset(exitrate_offset +
                     (o->hs_has_exitrate() ? sizeof(double) : 0)),

        flags_offset(index_offset +
                     (o->hs_has_index() ? sizeof(std::size_t) : 0)),

        cs_size(exitrate_offset - cs_offset) {
#ifndef NDEBUG
    std::cout << "cs_offset = " << cs_offset << std::endl;
    std::cout << "exitrate_offset = " << exitrate_offset << std::endl;
    std::cout << "index_offset = " << index_offset << std::endl;
#endif
  }

  std::size_t get_cs_size() const { return cs_size; }

  bool cs_equal(void const* const d1, void const* const d2) const {
    return memcmp(d1, d2, cs_size) == 0;
  }

  cstate* cs(void* const d) const {
    return static_cast<cstate*>(
        static_cast<void*>(static_cast<char*>(d) + cs_offset));
  }

  cstate const* cs(void const* const d) const {
    return static_cast<cstate const*>(
        static_cast<void const*>(static_cast<char const*>(d) + cs_offset));
  }

  double& exitrate(void* const d) const {
    assert(o->hs_has_exitrate());
    return *static_cast<double*>(
               static_cast<void*>(static_cast<char*>(d) + exitrate_offset));
  }

  double exitrate(void const* const d) const {
    assert(o->hs_has_exitrate());
    return *static_cast<double const*>(static_cast<void const*>(
        static_cast<char const*>(d) + exitrate_offset));
  }

  std::size_t& index(void* const d) const {
    assert(o->hs_has_index());
    return *static_cast<std::size_t*>(
               static_cast<void*>(static_cast<char*>(d) + index_offset));
  }

  std::size_t index(void const* const d) const {
    assert(o->hs_has_index());
    return *static_cast<std::size_t const*>(static_cast<void const*>(
        static_cast<char const*>(d) + index_offset));
  }

  std::size_t& flags(void* const d) const {
    assert(o->hs_has_flags());
    return *static_cast<std::size_t*>(
               static_cast<void*>(static_cast<char*>(d) + flags_offset));
  }

  std::size_t flags(void const* const d) const {
    assert(o->hs_has_flags());
    return *static_cast<std::size_t const*>(static_cast<void const*>(
        static_cast<char const*>(d) + flags_offset));
  }
};

class hstate_succ_data_accessor : public container::data_accessor {
 private:
  const std::size_t tri_offset;
  const std::size_t rates_g_offset;

 public:
  hstate_succ_data_accessor(options const* const o)
      : container::data_accessor(
            o, (o->hsucc_has_tr() ? sizeof(std::size_t) : 0) +
                   (o->hsucc_has_rates_g()
                        ? o->get_tr_rates_g_len() * sizeof(double)
                        : 0)),
        tri_offset(0),
        rates_g_offset(tri_offset +
                       (o->hsucc_has_tr() ? sizeof(std::size_t) : 0)) {}

  std::size_t& tri(void* const d) const {
    return *static_cast<std::size_t*>(
               static_cast<void*>(static_cast<char*>(d) + tri_offset));
  }

  std::size_t tri(void const* const d) const {
    return *static_cast<std::size_t const*>(static_cast<void const*>(
        static_cast<char const*>(d) + tri_offset));
  }

  double* rates_g(void* const d) const {
    return static_cast<double*>(
        static_cast<void*>(static_cast<char*>(d) + rates_g_offset));
  }

  double const* rates_g(void const* const d) const {
    return static_cast<double const*>(
        static_cast<void const*>(static_cast<char const*>(d) + rates_g_offset));
  }
};

struct hstate_h;

struct lstate_h {
  hstate_h* hs;

  lstate_h* lprev;
  lstate_h* lnext;
};

struct hstate_succ_h;

struct hstate_h {
  uint64_t hash;

  lstate_h* ls;

  hstate_h* hnext;
  hstate_h* lnext;

  hstate_succ_h* succ;
};

struct hstate_succ_h {
  hstate_h* hs;
  hstate_succ_h* lnext;
};

template <class lstate_data_accessor_t = lstate_data_accessor,
          class hstate_data_accessor_t = hstate_data_accessor,
          class hstate_succ_data_accessor_t = hstate_succ_data_accessor,
          typename lstate_t = lstate_h, typename hstate_t = hstate_h,
          typename hstate_succ_t = hstate_succ_h>
class dataAccessorT {
 public:
  typedef lstate_t state;
  typedef hstate_t hstate;
  typedef hstate_succ_t hstate_succ;

  dataAccessorT(
      options const* const o, lstate_data_accessor_t const* const lstate_da,
      hstate_data_accessor_t const* const hstate_da,
      hstate_succ_data_accessor_t const* const hstate_succ_da,
      container::allocator<lstate_t> const* const lstate_allocator,
      container::allocator<hstate_t> const* const hstate_allocator,
      container::allocator<hstate_succ_t> const* const hstate_succ_allocator)
      : o(o),
        lstate_da(lstate_da),
        hstate_da(hstate_da),
        hstate_succ_da(hstate_succ_da),
        lstate_allocator(lstate_allocator),
        hstate_allocator(hstate_allocator),
        hstate_succ_allocator(hstate_succ_allocator),
        rates_g_(o->get_tr_rates_g_len()),
        rates_h_(o->get_tr_rates_h_len()),
        sx_(o->s_x_len * o->mlt_u_len),
        cov_(o->xlen(o->nvars, o->nmoments)),
        iv_(o->nvars),
        phsucc_unexplored(&hsucc_unexplored) {}

  lstate_data_accessor_t const* get_lstate_da() const { return lstate_da; }

  hstate_data_accessor_t const* get_hstate_da() const { return hstate_da; }

  hstate_succ_data_accessor_t const* get_hstate_succ_da() const {
    return hstate_succ_da;
  }

  void* get_lstate_data(lstate_t* const ls) const {
    return lstate_allocator->get_data(ls);
  }

  void const* get_lstate_data(lstate_t const* const ls) const {
    return lstate_allocator->get_data(ls);
  }

  void* get_hstate_data(hstate_t* const hs) const {
    return hstate_allocator->get_data(hs);
  }

  void const* get_hstate_data(hstate_t const* const hs) const {
    return hstate_allocator->get_data(hs);
  }

  void* get_hstate_succ_data(hstate_succ_t* const hsucc) const {
    return hstate_succ_allocator->get_data(hsucc);
  }

  void const* get_hstate_succ_data(hstate_succ_t const* const hsucc) const {
    return hstate_succ_allocator->get_data(hsucc);
  }

  cstate* hs_d_cs(void* const hs_d) const { return hstate_da->cs(hs_d); }

  cstate const* hs_d_cs(void const* const hs_d) const {
    return hstate_da->cs(hs_d);
  }

  std::size_t& hs_d_index(void* const hs_d) const {
    return hstate_da->index(hs_d);
  }

  std::size_t hs_d_index(void const* const hs_d) const {
    return hstate_da->index(hs_d);
  }

  std::size_t& hs_index(hstate* const hs) const {
    return hs_d_index(get_hstate_data(hs));
  }

  std::size_t hs_index(hstate const* const hs) const {
    return hs_d_index(get_hstate_data(hs));
  }

  bool hs_is_erased(hstate const* const hs) const { return hs->ls == nullptr; }

  bool hs_is_temp(hstate const* const hs) const {
    return hs_flags(hs) & (1 << HS_FLAGS_TEMP);
  }

  void hs_set_temp(hstate* const hs) const {
    assert(!hs_is_temp(hs));
    hs_flags(hs) |= (1 << HS_FLAGS_TEMP);
  }

  void hs_unset_temp(hstate* const hs) const {
    assert(hs_is_temp(hs));
    hs_flags(hs) &= ~(1 << HS_FLAGS_TEMP);
  }

  bool hs_is_unexplored(hstate const* const hs) const {
    return hs->succ == phsucc_unexplored;
  }

  cstate* hs_cs(hstate* const hs) const { return hs_d_cs(get_hstate_data(hs)); }

  cstate const* hs_cs(hstate const* const hs) const {
    return hs_d_cs(get_hstate_data(hs));
  }

  double& hs_exitrate(hstate* const hs) const {
    return hstate_da->exitrate(get_hstate_data(hs));
  }

  double hs_exitrate(hstate const* const hs) const {
    return hstate_da->exitrate(get_hstate_data(hs));
  }

  std::size_t& hs_flags(hstate* const hs) const {
    return hstate_da->flags(get_hstate_data(hs));
  }

  std::size_t hs_flags(hstate const* const hs) const {
    return hstate_da->flags(get_hstate_data(hs));
  }

  void hs_init(hstate* const hs) const { hs->succ = phsucc_unexplored; }

  void hs_update_hash(hstate* const hs) const {
    hs->hash = util::BJenkins_lookup3_wrapper::datahash(
        hs_cs(hs), hstate_da->get_cs_size());
  }

  hstate_succ* hs_add_succ(hstate* const hs_succ, cstate const* const cs,
                           transition const* const tr) const {
    hstate_succ* const hsucc = hstate_succ_allocator->malloc();

    hsucc->hs = hs_succ;

    void* const hsucc_d = this->get_hstate_succ_data(hsucc);

    if (o->hsucc_has_tr()) {
      hstate_succ_da->tri(hsucc_d) = tr->get_index();
    }

    if (o->hsucc_has_rates_g()) {
      tr->rates_g(cs, hstate_succ_da->rates_g(hsucc_d));
    }

    return hsucc;
  }

  std::size_t hsucc_tri(void const* const hsucc_d) const {
    return hstate_succ_da->tri(hsucc_d);
  }

  std::size_t hsucc_tri(hstate_succ const* const hsucc) const {
    return hsucc_tri(this->get_hstate_succ_data(hsucc));
  }

  double const* tr_change(void const* const hsucc_d) const {
    return o->sl->get_model()->get_transitions()[hsucc_tri(hsucc_d)]->change();
  }

  void tr_rates(transition const* const tr, cstate const* const cs,
                double const* const x, double const*& rates_g,
                double const*& rates_h) const {
    tr->rates_g(cs, rates_g_.data());
    rates_g = rates_g_.data();

    if (!o->is_stoch()) {
      rates_h = tr_rates_h(tr, cs, x);
    }
  }

  void tr_rates(void const* const hsucc_d, cstate const* const cs,
                double const* const x, double const*& rates_g,
                double const*& rates_h) const {
    transition const* const tr =
        (!o->hsucc_has_rates_g() || o->is_hybrid())
            ? o->sl->get_model()
                  ->get_transitions()[hstate_succ_da->tri(hsucc_d)]
            : nullptr;

    if (o->hsucc_has_rates_g()) {
      rates_g = hstate_succ_da->rates_g(hsucc_d);
    } else {
      tr->rates_g(cs, rates_g_.data());
      rates_g = rates_g_.data();
    }

    if (!o->is_stoch()) {
      rates_h = tr_rates_h(tr, cs, x);
    }
  }

  double const* tr_rates_h(transition const* const tr, cstate const* const cs,
                           double const* const x) const {
    std::fill(rates_h_.begin(), rates_h_.end(), 0.0);
    tr->rates_h(cs, x, rates_h_.data());
    return rates_h_.data();
  }

  int const* cs_get(cstate const* const cs) const {
    o->sl->get_model()->cs_get(iv_.data(), cs);
    return iv_.data();
  }

  double* s_y(state* const s) const {
    return lstate_da->y(this->get_lstate_data(s));
  }

  double const* s_y(state const* const s) const {
    return lstate_da->y(this->get_lstate_data(s));
  }

  double& s_p(state* const s) const { return s_p(get_lstate_data(s)); }

  double s_p(state const* const s) const { return s_p(get_lstate_data(s)); }

  double& s_p(void* const ls_d) const { return o->v_p(lstate_da->y(ls_d)); }

  double s_p(void const* const ls_d) const {
    return o->v_p(lstate_da->y(ls_d));
  }

  double const* s_x(state const* const s, const double p) const {
    double const* const px = o->v_px(lstate_da->y(get_lstate_data(s)));
    if (p > 0.0) {
      for (std::size_t i = 0; i < o->npvars; i++) {
        sx_[i] = px[i] / p;
      }
      return sx_.data();
    }
    return px;
  }

  double const* s_x_cov(double const* const y) const {
    const double p = o->v_p(y);
    if (p > 0.0) {
      for (std::size_t i = 0; i < o->s_x_len; i++) {
        sx_[i] = o->v_px(y)[i] / p;

        if (o->do_sa() && !o->sa_no_der()) {
          for (std::size_t j = 0, jk = 0; j < o->nparams; j++) {
            sx_[(o->mlt_du_dc + j) * o->s_x_len + i] =
                (o->v_dpx_dc(y, j)[i] - o->v_dp_dc(y, j) * sx_[i]) / p;

            if (!o->sa_no_2der()) {
              for (std::size_t k = 0; k <= j; k++, jk++) {
                sx_[(o->mlt_d2u_dc2 + jk) * o->s_x_len + i] =
                    (o->v_d2px_dc2(y, jk)[i] -
                     o->v_dp_dc(y, j) *
                         sx_[(o->mlt_du_dc + k) * o->s_x_len + i] -
                     o->v_dp_dc(y, k) *
                         sx_[(o->mlt_du_dc + j) * o->s_x_len + i]) /
                    p;
              }
            }
          }
        }
      }
    } else {
      std::fill(sx_.begin(), sx_.end(), 0.0);
    }

    return sx_.data();
  }

  double* s_px(state* const s) const {
    return o->v_px(lstate_da->y(this->get_lstate_data(s)));
  }

  double const* s_px(state const* const s) const {
    return o->v_px(lstate_da->y(this->get_lstate_data(s)));
  }

 private:
  enum {
    HS_FLAGS_TEMP,
  };

  options const* const o;

  lstate_data_accessor_t const* const lstate_da;
  hstate_data_accessor_t const* const hstate_da;
  hstate_succ_data_accessor_t const* const hstate_succ_da;

  container::allocator<lstate_t> const* const lstate_allocator;
  container::allocator<hstate_t> const* const hstate_allocator;
  container::allocator<hstate_succ_t> const* const hstate_succ_allocator;

  mutable std::vector<double> rates_g_;
  mutable std::vector<double> rates_h_;

  mutable std::vector<double> sx_;
  mutable std::vector<double> cov_;

  mutable std::vector<int> iv_;

  hstate_succ hsucc_unexplored;
  hstate_succ* const phsucc_unexplored;
};

template <class lstate_data_accessor_t = lstate_data_accessor,
          class hstate_data_accessor_t = hstate_data_accessor,
          class hstate_succ_data_accessor_t = hstate_succ_data_accessor,
          typename lstate_t = lstate_h, typename hstate_t = hstate_h,
          typename hstate_succ_t = hstate_succ_h>
class stateHashListT : public container::abstractStateHashListT<
                           lstate_t, hstate_t, hstate_succ_t> {
 public:
  typedef container::abstractStateHashListT<lstate_t, hstate_t, hstate_succ_t>
      abstractStateHashList;

  typedef dataAccessorT<lstate_data_accessor_t, hstate_data_accessor_t,
                        hstate_succ_data_accessor_t, lstate_t, hstate_t,
                        hstate_succ_t> dataAccessor;

  typedef lstate_t state;
  typedef hstate_t hstate;
  typedef hstate_succ_t hstate_succ;

  stateHashListT(
      options const* const o, lstate_data_accessor_t const* const lstate_da,
      hstate_data_accessor_t const* const hstate_da,
      hstate_succ_data_accessor_t const* const hstate_succ_da,
      container::allocator<lstate_t>* const lstate_allocator,
      container::allocator<hstate_t>* const hstate_allocator,
      container::allocator<hstate_succ_t>* const hstate_succ_allocator,
      dataAccessor const* const da)
      : abstractStateHashList(lstate_allocator, hstate_allocator,
                              hstate_succ_allocator, hstate_da),
        o(o),
        lstate_da(lstate_da),
        hstate_da(hstate_da),
        hstate_succ_da(hstate_succ_da),
        lstate_allocator(lstate_allocator),
        hstate_allocator(hstate_allocator),
        hstate_succ_allocator(hstate_succ_allocator),
        da(da),
        cov_(o->xlen(o->nvars, o->nmoments)),
        iv_(o->nvars),
        sI(o->nmoments),
        sK(o->nmoments),
        sIk(o->nmoments),
        sI_(o->nvars),
        sJ_(o->nvars),
        sK_(o->nvars),
        sIk_(o->nvars) {}

  void update_hsucc_data() {
    if (o->hsucc_has_rates_g()) {
      assert(o->hsucc_has_tr());

      for (hstate* hs = this->get_hs_first(); hs != nullptr; hs = hs->lnext) {
        cstate const* const cs = da->hs_cs(hs);

        if (!da->hs_is_unexplored(hs)) {
          hstate_succ* hsucc = hs->succ;
          while (hsucc != nullptr) {
            void* const hsucc_d = da->get_hstate_succ_data(hsucc);

            transition const* const tr =
                o->sl->get_model()->get_transitions()[da->hsucc_tri(hsucc_d)];

            tr->rates_g(cs, hstate_succ_da->rates_g(hsucc_d));

            hsucc = hsucc->lnext;
          }
        }
      }
    }
  }

  void hs_restore_temp(hstate* const hs) {
    da->hs_unset_temp(hs);
    this->restore_temp(hs->ls);
  }

  void hs_restore(hstate* const hs) {
    assert(da->hs_is_erased(hs));
    this->add(hs);
  }

  void hs_restore_erased(hstate* const hs) {
    if (da->hs_is_erased(hs)) {
      this->add(hs);
    }
  }

  void s_make_temp(state* const s) {
    da->hs_set_temp(s->hs);
    this->make_temp(s);
  }

  void s_erase(state* const s) {
    lstate_da->zero(da->get_lstate_data(s));
    this->erase(s);
  }

  void erase_states() {
    for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
      s->hs->ls = nullptr;
      s->hs = nullptr;
      lstate_da->zero(da->get_lstate_data(s));
    }
    this->erase_all();
  }

  void s_init(state* const s) const { lstate_da->init(da->get_lstate_data(s)); }

  void s_zero_tmp(void* const ls_d) const { lstate_da->zero_tmp(ls_d); }

  void s_zero_tmp(state* const s) const { s_zero_tmp(da->get_lstate_data(s)); }

  void s_reset(state* const s) const {
    lstate_da->reset(da->get_lstate_data(s));
  }

  void init_states() {
    for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
      s_init(s);
    }
  }

  void reset_states() {
    for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
      s_reset(s);
    }
  }

  state* s_succ(cstate const* const cs, transition const* const tr, bool& fnew,
                const double t = 0.0) {
    hstate* hsucc = this->get_hstate_allocator()->get_last();
    cstate* csucc = da->hs_cs(hsucc);

    memcpy(csucc, cs, hstate_da->get_cs_size());

    tr->update(csucc, t);

    da->hs_update_hash(hsucc);

    state* const s = this->add(hsucc, fnew);

    if (fnew) {
      da->hs_init(s->hs);
      this->get_hstate_allocator()->malloc();
    }

    return s;
  }

  void s_explore_unexplored(state* const s) {
    hstate* const hs = s->hs;
    if (da->hs_is_unexplored(hs)) {
      cstate const* const cs = da->hs_cs(hs);

      hs->succ = nullptr;
      for (auto const& tr : o->trs_stoch) {
        if (tr->is_enabled(cs)) {
          bool fnew;
          state* const succ = s_succ(cs, tr, fnew);

          if (succ != s) {
            if (fnew) {
              s_init(succ);
            }

            hstate_succ* const hsucc = da->hs_add_succ(succ->hs, cs, tr);

            hsucc->lnext = hs->succ;
            hs->succ = hsucc;
          }
        }
      }
    }
  }

  void cs_print(cstate const* const cs) const {
    int const* const v = da->cs_get(cs);
    for (std::size_t i = 0; i < o->nvars; i++) {
      printf("%d ", v[i]);
    }
    printf("\n");
  }

  void s_print(state const* const s) const { cs_print(hs_cs(s->hs)); }

  void add_init_state(const solver_loader::model_info::ic_s& s0,
                      int const* const iv) {
    hstate* hs = this->get_hstate_allocator()->get_last();
    cstate* const cs = da->hs_cs(hs);
    o->sl->get_model()->cs_set(cs, iv);

    da->hs_update_hash(hs);

    bool fnew;
    state* const s = this->add(hs, fnew);
    double* const y = da->s_y(s);

    if (fnew) {
      da->hs_init(s->hs);
      this->get_hstate_allocator()->malloc();

      o->v_p(y) = s0.p;

      if (o->do_sa() && !o->sa_no_der()) {
        for (std::size_t i = 0, ij = 0; i < o->nivars; i++) {
          o->v_dp_di(y, i) = s0.dp_di[i];

          if (!o->sa_no_2der()) {
            for (std::size_t j = 0; j <= i; j++, ij++) {
              const std::size_t k = o->nparams + o->nevars + i;
              const std::size_t l = o->nparams + o->nevars + j;
              const std::size_t kl = o->CI(k, l);

              o->v_d2p_dc2(y, kl) = s0.d2p_di2[ij];
            }
          }
        }
      }

      if (o->is_hybrid()) {
        for (std::size_t i = 0; i < o->npvars; i++) {
          o->v_px(y)[i] = s0.p * iv[i];
        }
      }
    } else {
      o->v_p(y) += s0.p;

      if (o->do_sa() && !o->sa_no_der()) {
        for (std::size_t i = 0, ij = 0; i < o->nivars; i++) {
          o->v_dp_di(y, i) += s0.dp_di[i];

          if (!o->sa_no_2der()) {
            for (std::size_t j = 0; j <= i; j++, ij++) {
              const std::size_t k = o->nparams + o->nevars + i;
              const std::size_t l = o->nparams + o->nevars + j;
              const std::size_t kl = o->CI(k, l);

              o->v_d2p_dc2(y, kl) += s0.d2p_di2[ij];
            }
          }
        }
      }

      if (o->is_hybrid()) {
        for (std::size_t i = 0; i < o->npvars; i++) {
          o->v_px(y)[i] += s0.p * iv[i];
        }
      }
    }
  }

  void load_gaussian_init_helper(int* const iv, const std::size_t ivi,
                                 solver_loader::model_info::ic_s& s0,
                                 double& psum, const double esum = 0.0) {
    const double d = 1.0;
    const double C = -0.50 / d / d;

    if (ivi < o->nivars) {
      const double x = o->init_vals[ivi];
      const double r = 3.0 * d + 0.50;

      for (int v = std::max(0, (int)(x - r)); v <= (int)(x + r); v++) {
        iv[o->sl->sa_ivars[ivi]->get_index()] = v;

        load_gaussian_init_helper(iv, ivi + 1, s0, psum,
                                  esum + (x - v) * (x - v));
      }
    } else {
      s0.p = std::exp(C * esum);

      if (s0.p > 0.0) {
        if (!o->sa_no_der()) {
          for (std::size_t i = 0, ij = 0; i < o->nivars; i++) {
            const double f = C * 2.0 * (o->init_vals[i] -
                                        iv[o->sl->sa_ivars[i]->get_index()]);

            s0.dp_di[i] = s0.p * f;

            if (!o->sa_no_2der()) {
              for (std::size_t j = 0; j <= i; j++, ij++) {
                s0.d2p_di2[ij] = s0.dp_di[j] * f;
              }
            }
          }
        }

        psum += s0.p;

        add_init_state(s0, iv);
      }
    }
  }

  void load_gaussian_init(solver_loader::model_info::ic_s& s0) {
    erase_states();

    std::vector<int> iv(o->nvars, 0);
    for (auto const& s0i : s0.li) {
      iv[s0i.v->get_index()] = s0i.value->get();
    }

    double psum = 0.0;
    load_gaussian_init_helper(iv.data(), 0, s0, psum);

    const double _1_p_sum = 1.0 / psum;

    for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
      double* const y = da->s_y(s);
      for (std::size_t i = 0; i < o->s_y_len; i++) {
        y[i] *= _1_p_sum;
      }
      s_init(s);
    }
  }

  void load_init(solver_loader::model_info::ic const* const init) {
    erase_states();

    if (o->is_det()) {
      hstate* const hs = this->get_hstate_allocator()->get_last();
      o->sl->get_model()->cs_set(da->hs_cs(hs), iv_.data());
      da->hs_update_hash(hs);
      bool fnew;
      state* const s = this->add(hs, fnew);
      da->hs_init(s->hs);
      this->get_hstate_allocator()->malloc();

      da->s_p(s) = 1.0;
      for (auto const& s0 : init->get_states()) {
        for (auto const& s0i : s0.li) {
          da->s_px(s)[s0i.v->get_index()] += s0.p * s0i.value->get();
        }
      }
    } else {
      for (auto const& s0 : init->get_states()) {
        std::fill(iv_.begin(), iv_.end(), 0);
        for (auto const& s0i : s0.li) {
          iv_[s0i.v->get_index()] = s0i.value->get();
        }
        add_init_state(s0, iv_.data());
      }
    }

    if (!o->is_stoch() && !o->is_det_centered()) {
      for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
        void const* const ls_d = da->get_lstate_data(s);
        double const* const y = this->lstate_da->y(ls_d);
        double const* const x = da->s_x_cov(y);

        for (std::size_t i = 0; i < o->npvars; i++) {
          cov_[i] = x[i];

          if (o->nmoments > 1) {
            sI[0] = i;
            sI_[i]++;

            assert(sI_[i] == 1);

            raw_moments_(x, cov_.data(), o->nmoments, o->sl->p_vars, sI.data(),
                         sI_.data(), 1, o->npvars, o->npvars);

            sI_[i]--;

            assert(sI_[i] == 0);
          }
        }

        const double p = o->v_p(y);
        for (std::size_t i = 0; i < o->s_x_len; i++) {
          da->s_px(s)[i] = p * cov_[i];
        }
      }
    }
  }

  void raw_moments_(
      double const* const COV, double* const X, const std::size_t nmoments_,
      const std::vector<solver_loader::model_info::var const*>& vars_,
      std::size_t* const I, std::size_t* const I_, const std::size_t d,
      const std::size_t ci0, const std::size_t xci0) const {
    const std::size_t cov_len_ = o->covlen(o->nvars, d + 1);

    const std::size_t nvars_ = vars_.size();
    const std::size_t xcov_len_ = o->covlen(nvars_, d + 1);

    for (std::size_t i = 0; i <= I[d - 1]; i++) {
      I[d] = i;
      I_[i]++;

      X[xci0 + o->cov_indexl(I, d + 1)] =
          raw_moment_(I_, COV, vars_, 1.0, sJ_.data(), 0);

      if (d + 1 < nmoments_) {
        raw_moments_(COV, X, nmoments_, vars_, I, I_, d + 1, ci0 + cov_len_,
                     xci0 + xcov_len_);
      }

      I_[i]--;
    }
  }

  double raw_moment_(
      std::size_t const* const I_, double const* const x,
      const std::vector<solver_loader::model_info::var const*>& vars_,
      const double r_, std::size_t* const J_, const std::size_t d) const {
    const std::size_t nvars_ = vars_.size();

    if (d < nvars_) {
      double s = 0.0;
      double r = 1.0;

      assert(J_[d] == 0);

      for (std::size_t j = 0; j <= I_[d]; j++) {
        J_[d] = j;

        if (j) {
          r *= x[vars_[d]->get_index()];
        }

        s += raw_moment_(I_, x, vars_, r_ * r * o->cmb[I_[d]][j], J_, d + 1);
      }
      J_[d] = 0;
      return s;
    }

    for (std::size_t i = 0; i < nvars_; i++) {
      assert(sK_[i] == 0);

      sK_[i] = I_[i] - J_[i];
    }

    std::size_t K_len;
    o->cov_I_I(nvars_, sK.data(), K_len, sK_.data());

    for (std::size_t i = 0; i < nvars_; i++) {
      sK_[i] = 0;
    }

    if (K_len > 1) {
      return r_ * x[o->cov_indexg(o->npvars, sK.data(), K_len, vars_)];
    }

    if (K_len == 1) {
      return 0.0;
    }

    return r_;
  }

  void copy_moments_(
      double const* const COV, double* const X, const std::size_t nmoments_,
      const std::vector<solver_loader::model_info::var const*>& vars_,
      std::size_t* const I, std::size_t* const I_, const std::size_t d,
      const std::size_t ci0, const std::size_t xci0) const {
    const std::size_t cov_len_ = o->covlen(o->nvars, d + 1);

    const std::size_t nvars_ = vars_.size();
    const std::size_t xcov_len_ = o->covlen(nvars_, d + 1);

    for (std::size_t i = 0; i <= I[d - 1]; i++) {
      I[d] = i;
      I_[i]++;

      X[xci0 + o->cov_indexl(I, d + 1)] =
          COV[ci0 + o->cov_indexl(I, d + 1, vars_)];

      if (d + 1 < nmoments_) {
        copy_moments_(COV, X, nmoments_, vars_, I, I_, d + 1, ci0 + cov_len_,
                      xci0 + xcov_len_);
      }

      I_[i]--;
    }
  }

  void center_moments_(double const* const COV, double* const X,
                       const std::size_t nmoments_, const std::size_t nvars_,
                       std::size_t* const I, std::size_t* const I_,
                       const std::size_t d, const std::size_t ci0) const {
    const std::size_t cov_len_ = o->covlen(nvars_, d + 1);

    for (std::size_t i = 0; i <= I[d - 1]; i++) {
      I[d] = i;
      I_[i]++;

      X[ci0 + o->cov_indexl(I, d + 1)] =
          central_moment_(I_, X, COV, nvars_, 1.0, sJ_.data(), 0);

      if (d + 1 < nmoments_) {
        center_moments_(COV, X, nmoments_, nvars_, I, I_, d + 1,
                        ci0 + cov_len_);
      }

      I_[i]--;
    }
  }

  double central_moment_(std::size_t const* const I_, double const* const X,
                         double const* const COV, const std::size_t nvars_,
                         const double r_, std::size_t* const J_,
                         const std::size_t d) const {
    if (d < nvars_) {
      assert(J_[d] == 0);

      double s = 0.0;
      double r = -1.0;
      for (std::size_t j = 0; j <= I_[d]; j++) {
        J_[d] = j;

        if (j) {
          r *= X[d];
        }

        r = -r;

        s += central_moment_(I_, X, COV, nvars_, r_ * r * o->cmb[I_[d]][j], J_,
                             d + 1);
      }
      J_[d] = 0;
      return s;
    }

    for (std::size_t i = 0; i < nvars_; i++) {
      assert(sK_[i] == 0);

      sK_[i] = I_[i] - J_[i];
    }

    std::size_t K_len;
    o->cov_I_I(nvars_, sK.data(), K_len, sK_.data());

    for (std::size_t i = 0; i < nvars_; i++) {
      sK_[i] = 0;
    }

    if (K_len > 1) {
      return r_ * COV[o->cov_indexg(nvars_, sK.data(), K_len)];
    }

    if (K_len > 0) {
      return r_ * X[sK[0]];
    }

    return r_;
  }

  void s_raw_moments_stoch_(
      int const* const v, const double px_, double* const rX,
      const std::size_t nmoments_,
      const std::vector<solver_loader::model_info::var const*>& vars_,
      std::size_t* const I, std::size_t* const I_, const std::size_t d,
      const std::size_t ci0) const {
    const std::size_t nvars_ = vars_.size();
    const std::size_t cov_len_ = o->covlen(nvars_, d + 1);

    for (std::size_t i = 0; i <= I[d - 1]; i++) {
      I[d] = i;
      I_[i]++;

      const double px__ = px_ * v[vars_[i]->get_index()];

      rX[ci0 + o->cov_indexl(I, d + 1)] += px__;

      if (d + 1 < nmoments_) {
        s_raw_moments_stoch_(v, px__, rX, nmoments_, vars_, I, I_, d + 1,
                             ci0 + cov_len_);
      }

      I_[i]--;
    }
  }

  void s_raw_moments_hybrid_(
      int const* const v, double const* const cX, double* const rX,
      const std::size_t nmoments_,
      const std::vector<solver_loader::model_info::var const*>& vars_,
      std::size_t* const I, std::size_t* const I_, const std::size_t d,
      const double vv,
      std::vector<solver_loader::model_info::var const*>& dvars,
      const std::size_t ci0) const {
    const std::size_t nvars_ = vars_.size();
    const std::size_t cov_len_ = o->covlen(nvars_, d + 1);

    for (std::size_t i = 0; i <= I[d - 1]; i++) {
      I[d] = i;
      I_[i]++;

      double vv_ = vv;
      const std::size_t vi = vars_[i]->get_index();

      if (vi < o->npvars) {
        dvars.push_back(vars_[i]);
      } else {
        vv_ *= v[vi];
      }

      double c_ = vv_;

      if (!dvars.empty()) {
        std::vector<solver_loader::model_info::var const*> dvars_;
        std::size_t j = dvars.size() - 1;
        dvars_.push_back(dvars[j]);

        assert(sIk_[dvars_.size() - 1] == 0);

        sIk_[dvars_.size() - 1]++;

        while (j) {
          j--;
          if (dvars[j] != dvars[j + 1]) {
            dvars_.push_back(dvars[j]);

            assert(sIk_[dvars_.size() - 1] == 0);
          }
          sIk_[dvars_.size() - 1]++;
        }

        c_ *= raw_moment_(sIk_.data(), cX, dvars_, 1.0, sJ_.data(), 0);

        for (j = 0; j < dvars_.size(); j++) {
          sIk_[j] = 0;
        }
      }

      rX[ci0 + o->cov_indexl(I, d + 1)] = c_;

      if (d + 1 < nmoments_ && std::fabs(vv) > 0.0) {
        s_raw_moments_hybrid_(v, cX, rX, nmoments_, vars_, I, I_, d + 1, vv_,
                              dvars, ci0 + cov_len_);
      }

      if (vi < o->npvars) {
        dvars.pop_back();
      }
      I_[i]--;
    }
  }

  void s_copy_moments_hybrid_(
      int const* const v, double const* const srcX, double* const dstX,
      const std::size_t nmoments_,
      const std::vector<solver_loader::model_info::var const*>& vars_,
      std::size_t* const I, std::size_t* const I_, const std::size_t d,
      const double vv,
      std::vector<solver_loader::model_info::var const*>& dvars,
      const std::size_t ci0) const {
    const std::size_t nvars_ = vars_.size();
    const std::size_t cov_len_ = o->covlen(nvars_, d + 1);

    for (std::size_t i = 0; i <= I[d - 1]; i++) {
      I[d] = i;
      I_[i]++;

      double vv_ = vv;
      const std::size_t vi = vars_[i]->get_index();

      if (vi < o->npvars) {
        dvars.push_back(vars_[i]);
      } else {
        vv_ *= v[vi];
      }

      double c_ = vv_;

      if (!dvars.empty()) {
        std::vector<solver_loader::model_info::var const*> dvars_;
        std::size_t j = dvars.size() - 1;
        dvars_.push_back(dvars[j]);

        assert(sIk_[dvars_.size() - 1] == 0);

        sIk_[dvars_.size() - 1]++;

        while (j) {
          j--;
          if (dvars[j] != dvars[j + 1]) {
            dvars_.push_back(dvars[j]);

            assert(sIk_[dvars_.size() - 1] == 0);
          }
          sIk_[dvars_.size() - 1]++;
        }

        std::size_t sIk_len;
        o->cov_I_I(dvars_.size(), sIk.data(), sIk_len, sIk_.data());

        c_ *= srcX[o->cov_indexg(o->npvars, sIk.data(), sIk_len, dvars_)];

        std::fill(sIk.begin(), sIk.begin() + sIk_len, 0);
        std::fill(sIk_.begin(), sIk_.begin() + dvars_.size(), 0);
      }

      dstX[ci0 + o->cov_indexl(I, d + 1)] = c_;

      if (d + 1 < nmoments_ && std::fabs(vv) > 0.0) {
        s_copy_moments_hybrid_(v, srcX, dstX, nmoments_, vars_, I, I_, d + 1,
                               vv_, dvars, ci0 + cov_len_);
      }

      if (vi < o->npvars) {
        dvars.pop_back();
      }
      I_[i]--;
    }
  }

  void bfs(const std::size_t max_states = 1024) {
    std::queue<state*> Q;

    std::size_t si = 0;
    for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
      da->hs_index(s->hs) = si++;
      Q.push(s);
    }

    while (!Q.empty() && si < max_states) {
      state* const s = Q.front();
      Q.pop();

      hstate* const hs = s->hs;

      if (da->hs_is_unexplored(hs)) {
        cstate const* const cs = da->hs_cs(hs);

        hs->succ = nullptr;
        for (auto const& tr : o->trs_stoch) {
          if (tr->is_enabled(cs)) {
            bool fnew;
            state* const succ = s_succ(cs, tr, fnew);

            if (succ != s) {
              if (fnew) {
                da->hs_index(succ->hs) = si++;
                Q.push(succ);
              }

              hstate_succ* const hsucc = da->hs_add_succ(succ->hs, cs, tr);

              hsucc->lnext = hs->succ;
              hs->succ = hsucc;
            }
          }
        }
      }
    }
  }

  void comp_objf(double* const f, const double p_scale = 1.0) const {
    if (o->is_det()) {
      state const* const s = this->get_first();
      cstate const* const cs = da->hs_cs(s->hs);
      double const* const x = da->s_px(s);

      o->v_f(f) =
          o->exprs[o->sl->objf_index[o->mlt_u]](cs, x, o->param_vals.data());

    } else {
      std::fill(f, f + o->f_len, 0.0);

      if (o->is_hybrid()) {
        for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
          const double p = da->s_p(s) * p_scale;
          cstate const* const cs = da->hs_cs(s->hs);
          double const* const x = da->s_x(s, p);

          o->v_f(f) += p *
                       o->exprs[o->sl->objf_index[o->mlt_u]](
                           cs, x, o->param_vals.data());
        }
      } else {
        for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
          cstate const* const cs = da->hs_cs(s->hs);
          double const* const x = nullptr;
          double const* const y = da->s_y(s);

          const double p = o->v_p(y);
          const double f_ = o->exprs[o->sl->objf_index[o->mlt_u]](
              cs, x, o->param_vals.data());

          o->v_f(f) += p * f_;
          if (!o->sa_no_der()) {
            for (std::size_t i = 0; i < o->nparams; i++) {
              const double df_ = o->exprs[o->sl->objf_index[o->mlt_du_dc + i]](
                  cs, x, o->param_vals.data());
              o->v_df_dc(f)[i] += (o->v_dp_dc(y, i) * f_ + p * df_) * p_scale;
            }
          }
        }
      }
    }
  }

  void comp_exprs(double* const X, double* const F, double* const Xstd,
                  double* const Fstd, const double p_scale = 1.0) const {
    const std::size_t nexprs = o->exprs.size();

    if (o->is_det()) {
      state const* const s = this->get_first();
      cstate const* const cs = da->hs_cs(s->hs);
      void const* const ls_d = da->get_lstate_data(s);
      double const* const y = da->get_lstate_da()->y(ls_d);
      double const* const x = o->v_px(y);

      for (std::size_t i = 0; i < o->nvars; i++) {
        X[i] = o->v_px(y)[i];

        if (o->do_sa() && !o->sa_no_der()) {
          for (std::size_t j = 0, jk = 0; j < o->nparams; j++) {
            X[(o->mlt_du_dc + j) * o->nvars + i] = o->v_dpx_dc(y, j)[i];

            if (!o->sa_no_2der()) {
              for (std::size_t k = 0; k <= j; k++, jk++) {
                X[(o->mlt_d2u_dc2 + jk) * o->nvars + i] =
                    o->v_d2px_dc2(y, jk)[i];
              }
            }
          }
        }

        Xstd[i] = std::sqrt(o->is_det_centered()
                                ? x[o->nvars + o->CI(i, i)]
                                : x[o->nvars + o->CI(i, i)] - x[i] * x[i]);
      }

      for (std::size_t i = 0; i < nexprs; i++) {
        F[i] = o->exprs[i](cs, x, o->param_vals.data());
        Fstd[i] = 0.0;
      }
    } else {
      std::fill(X, X + o->mlt_u_len * o->nvars, 0.0);
      std::fill(F, F + o->mlt_u_len * nexprs, 0.0);
      std::fill(Xstd, Xstd + o->mlt_u_len * o->nvars, 0.0);
      std::fill(Fstd, Fstd + o->mlt_u_len * nexprs, 0.0);

      if (o->is_hybrid()) {
        for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
          void const* const ls_d = da->get_lstate_data(s);
          double const* const y = da->get_lstate_da()->y(ls_d);
          const double p = o->v_p(y) * p_scale;

          cstate const* const cs = da->hs_cs(s->hs);
          double const* const x = da->s_x(s, p);
          int const* const v = da->cs_get(cs);

          for (std::size_t i = 0; i < o->nvars; i++) {
            const double x_ = (i < o->npvars) ? x[i] : v[i];
            X[i] += p * x_;

            if (o->do_sa() && !o->sa_no_der()) {
              for (std::size_t j = 0, jk = 0; j < o->nparams; j++) {
                X[(o->mlt_du_dc + j) * o->nvars + i] += o->v_dp_dc(y, j) * x_;
                if (i < o->npvars) {
                  X[(o->mlt_du_dc + j) * o->nvars + i] +=
                      p * o->v_dpx_dc(y, j)[i];
                }

                if (!o->sa_no_2der()) {
                  for (std::size_t k = 0; k <= j; k++, jk++) {
                    X[(o->mlt_d2u_dc2 + jk) * o->nvars + i] +=
                        o->v_d2p_dc2(y, jk) * x_;
                    if (i < o->npvars) {
                      X[(o->mlt_d2u_dc2 + jk) * o->nvars + i] +=
                          o->v_dp_dc(y, j) * o->v_dpx_dc(y, k)[i] +
                          o->v_dp_dc(y, k) * o->v_dpx_dc(y, j)[i] +
                          p * o->v_d2px_dc2(y, jk)[i];
                    }
                  }
                }
              }
            }

            Xstd[i] += p * x_ * x_;
          }
          for (std::size_t i = 0; i < nexprs; i++) {
            const double f_ = o->exprs[i](cs, x, o->param_vals.data());
            F[i] += p * f_;
            Fstd[i] += p * f_ * f_;
          }
        }
      } else {
        for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
          void const* const ls_d = da->get_lstate_data(s);
          double const* const y = da->get_lstate_da()->y(ls_d);
          const double p = o->v_p(y) * p_scale;

          cstate const* const cs = da->hs_cs(s->hs);
          double const* const x = nullptr;
          int const* const v = da->cs_get(cs);

          for (std::size_t i = 0; i < o->nvars; i++) {
            X[i] += p * v[i];

            if (o->do_sa() && !o->sa_no_der()) {
              for (std::size_t j = 0, jk = 0; j < o->nparams; j++) {
                X[(o->mlt_du_dc + j) * o->nvars + i] += o->v_dp_dc(y, j) * v[i];

                if (!o->sa_no_2der()) {
                  for (std::size_t k = 0; k <= j; k++, jk++) {
                    X[(o->mlt_d2u_dc2 + jk) * o->nvars + i] +=
                        o->v_d2p_dc2(y, jk) * v[i];
                  }
                }
              }
            }

            Xstd[i] += p * v[i] * v[i];
          }
          for (std::size_t i = 0; i < nexprs; i++) {
            const double f_ = o->exprs[i](cs, x, o->param_vals.data());
            F[i] += p * f_;
            Fstd[i] += p * f_ * f_;
          }
        }
      }

      for (std::size_t i = 0; i < o->nvars; i++) {
        Xstd[i] = std::sqrt(Xstd[i] - X[i] * X[i]);
      }
      for (std::size_t i = 0; i < nexprs; i++) {
        Fstd[i] = std::sqrt(Fstd[i] - F[i] * F[i]);
      }
    }
  }

  void comp_moments(double* const X,
                    solver_loader::task_info::dump_moments const* const ti,
                    const double p_scale = 1.0) const {
    const std::vector<solver_loader::model_info::var const*>& vars_ =
        ti->get_vars();
    const std::size_t nvars_ = vars_.size();
    const std::size_t nmoments_ = ti->get_nmoments();
    const std::size_t x_len_ = o->xlen(nvars_, nmoments_);

    const int cond_index = ti->get_cond_index();

    if (o->is_det()) {
      state const* const s = this->get_first();

      for (std::size_t i = 0; i < nvars_; i++) {
        const std::size_t vi = vars_[i]->get_index();

        X[i] = da->s_px(s)[vi];

        if (nmoments_ > 1) {
          sI[0] = i;
          sI_[i]++;

          assert(sI_[i] == 1);

          if (o->is_det_centered()) {
            if (ti->dump_raw()) {
              raw_moments_(da->s_px(s), X, nmoments_, vars_, sI.data(),
                           sI_.data(), 1, o->nvars, nvars_);
            } else {
              copy_moments_(da->s_px(s), X, nmoments_, vars_, sI.data(),
                            sI_.data(), 1, o->nvars, nvars_);
            }
          } else {
            if (ti->dump_raw()) {
              copy_moments_(da->s_px(s), X, nmoments_, vars_, sI.data(),
                            sI_.data(), 1, o->nvars, nvars_);
            } else {
              center_moments_(da->s_px(s), X, nmoments_, nvars_, sI.data(),
                              sI_.data(), 1, nvars_);
            }
          }

          sI_[i]--;

          assert(sI_[i] == 0);
        }
      }
    } else {
      std::fill(X, X + x_len_, 0.0);

      std::vector<double> COV(x_len_, 0.0);

      double p_sum = 0.0;

      if (o->is_hybrid()) {
        std::vector<solver_loader::model_info::var const*> dvars;

        for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
          void const* const ls_d = da->get_lstate_data(s);
          double const* const y = this->lstate_da->y(ls_d);
          const double p = o->v_p(y) * p_scale;
          cstate const* const cs = da->hs_cs(s->hs);
          double const* const x = da->s_x_cov(y);

          if (cond_index < 0 ||
              o->exprs[cond_index](cs, x, o->param_vals.data())) {
            int const* const v = da->cs_get(cs);

            for (std::size_t i = 0; i < nvars_; i++) {
              sI[0] = i;
              sI_[i]++;

              assert(sI_[i] == 1);

              double vv = 1.0;
              const std::size_t vi = vars_[i]->get_index();

              if (vi < o->npvars) {
                X[i] += p * x[vi];
                dvars.push_back(vars_[i]);
              } else {
                X[i] += p * v[vi];
                vv *= v[vi];
              }

              if (nmoments_ > 1 && std::fabs(vv) > 0.0) {
                if (o->is_det_centered()) {
                  s_raw_moments_hybrid_(v, x, cov_.data(), nmoments_, vars_,
                                        sI.data(), sI_.data(), 1, vv, dvars,
                                        nvars_);
                } else {
                  s_copy_moments_hybrid_(v, x, cov_.data(), nmoments_, vars_,
                                         sI.data(), sI_.data(), 1, vv, dvars,
                                         nvars_);
                }
              }

              if (vi < o->npvars) {
                dvars.pop_back();
              }
              sI_[i]--;

              assert(sI_[i] == 0);
              assert(dvars.empty());
            }

            for (std::size_t i = nvars_; i < x_len_; i++) {
              COV[i] += p * cov_[i];
            }

            p_sum += p;
          }
        }
      } else {
        for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
          cstate const* const cs = da->hs_cs(s->hs);
          double const* const x = nullptr;

          if (cond_index < 0 ||
              o->exprs[cond_index](cs, x, o->param_vals.data())) {
            const double p = da->s_p(s) * p_scale;
            int const* const v = da->cs_get(cs);

            for (std::size_t i = 0; i < nvars_; i++) {
              const std::size_t vi = vars_[i]->get_index();
              const double px_ = p * v[vi];

              X[i] += px_;

              if (nmoments_ > 1) {
                sI[0] = i;
                sI_[i]++;

                assert(sI_[i] == 1);

                s_raw_moments_stoch_(v, px_, COV.data(), nmoments_, vars_,
                                     sI.data(), sI_.data(), 1, nvars_);

                sI_[i]--;

                assert(sI_[i] == 0);
              }
            }

            p_sum += p;
          }
        }
      }

      if (cond_index >= 0 && p_sum > std::numeric_limits<double>::epsilon()) {
        const double _1_p_sum = 1.0 / p_sum;
        for (std::size_t i = 0; i < nvars_; i++) {
          X[i] *= _1_p_sum;
        }
        for (std::size_t i = nvars_; i < x_len_; i++) {
          COV[i] *= _1_p_sum;
        }
      }

      if (ti->dump_raw()) {
        for (std::size_t i = nvars_; i < x_len_; i++) {
          X[i] = COV[i];
        }
      } else {
        if (nmoments_ > 1) {
          for (std::size_t i = 0; i < nvars_; i++) {
            sI[0] = i;
            sI_[i]++;

            assert(sI_[i] == 1);

            center_moments_(COV.data(), X, nmoments_, nvars_, sI.data(),
                            sI_.data(), 1, nvars_);

            sI_[i]--;

            assert(sI_[i] == 0);
          }
        }
      }
    }
  }

  void dump_distr(std::ofstream& os,
                  solver_loader::task_info::dump_distr const* const ti,
                  const double p_scale = 1.0) const {
    double p_sum = 0.0;
    const int cond_index = ti->get_cond_index();
    if (cond_index >= 0) {
      for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
        const double p = da->s_p(s) * p_scale;
        cstate const* const cs = da->hs_cs(s->hs);
        double const* const x = da->s_x(s, p);

        if (o->exprs[cond_index](cs, x, o->param_vals.data())) {
          p_sum += p;
        }
      }
    }

    if (cond_index < 0 || std::numeric_limits<double>::epsilon()) {
      const double _1_p_sum = (cond_index < 0) ? 1.0 : 1.0 / p_sum;

      for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
        void* const ls_d = da->get_lstate_data(s);
        double const* const y = this->lstate_da->y(ls_d);

        const double p = o->v_p(y) * p_scale;
        cstate const* const cs = da->hs_cs(s->hs);
        double const* const x = da->s_x(s, p);

        if (cond_index < 0 ||
            o->exprs[cond_index](cs, x, o->param_vals.data())) {
          int const* v = da->cs_get(cs);

          for (std::size_t i = 0; i < o->npvars; i++) {
            os << x[i] << o->csvSep;
          }
          for (std::size_t i = o->npvars; i < o->nvars; i++) {
            os << v[i] << o->csvSep;
          }

          os << p* _1_p_sum;

          os << std::endl;
        }
      }
    }
  }

  bool dump_plot_distr_data(
      solver_loader::task_info::plot_distr const* const ti,
      dump const* const dmp, double p_scale = 1.0) const {
    const int cond_index = ti->get_cond_index();
    if (cond_index >= 0) {
      double p_sum = 0.0;

      for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
        const double p = da->s_p(s);
        cstate const* const cs = da->hs_cs(s->hs);
        double const* const x = da->s_x(s, p);

        if (o->exprs[cond_index](cs, x, o->param_vals.data())) {
          p_sum += p;
        }
      }

      p_scale /= p_sum;
    }

    std::ofstream os(dmp->get_local_path());
    os.precision(16);
    os << "count";

    std::vector<std::map<float, double>> V;

    const std::size_t nexprs = ti->get_exprs().size();

    int vmin = 0, vmax = 0;

    if (!nexprs) {
      V.resize(o->nvars);

      for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
        const double p = da->s_p(s);
        cstate const* const cs = da->hs_cs(s->hs);
        double const* const x = da->s_x(s, p);
        int const* const v = da->cs_get(cs);

        if (cond_index < 0 ||
            o->exprs[cond_index](cs, x, o->param_vals.data())) {
          for (std::size_t i = 0; i < o->nvars; i++) {
            const float v_ = (i < o->npvars) ? x[i] : v[i];

            V[i][v_] += p;

            if (s != this->get_first() || i) {
              if (v_ < vmin) {
                vmin = (int)(v_ - 0.5);
              }
              if (v_ > vmax) {
                vmax = (int)(v_ + 0.5);
              }
            } else {
              vmin = (int)(v_ - 0.5);
              vmax = (int)(v_ + 0.5);
            }
          }
        }
      }

      for (std::size_t i = 0; i < o->nvars; i++) {
        os << o->csvSep << o->sl->get_vars()[i]->get_name();
      }
    } else {
      V.resize(nexprs);

      for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
        void const* const ls_d = da->get_lstate_data(s);
        double const* const y = da->get_lstate_da()->y(ls_d);
        const double p = o->v_p(y);

        cstate const* const cs = da->hs_cs(s->hs);
        double const* const x = da->s_x(s, p);

        if (cond_index < 0 ||
            o->exprs[cond_index](cs, x, o->param_vals.data())) {
          for (std::size_t i = 0; i < nexprs; i++) {
            auto const& e = ti->get_exprs()[i];

            const float v_ =
                o->exprs[e.expr_index](cs, x, o->param_vals.data());

            const double val =
                (e.dparami == nullptr)
                    ? p
                    : ((e.dparamj == nullptr))
                          ? o->v_dp_dc(y, e.dparami->get_index_p())
                          : o->v_d2p_dc2(y, o->CI(e.dparami->get_index_p(),
                                                  e.dparamj->get_index_p()));

            V[i][v_] += val;

            if (s != this->get_first() || i) {
              if (v_ < vmin) {
                vmin = (int)(v_ - 0.5);
              }
              if (v_ > vmax) {
                vmax = (int)(v_ + 0.5);
              }
            } else {
              vmin = (int)(v_ - 0.5);
              vmax = (int)(v_ + 0.5);
            }
          }
        }
      }

      std::string tmp_;
      for (std::size_t i = 0; i < nexprs; i++) {
        auto const& e = ti->get_exprs()[i];
        if (e.line_props.get_title(tmp_)) {
          os << o->csvSep << tmp_;
        } else {
          os << o->csvSep << e.expr;
        }
      }
    }

    os << std::endl;
    for (int v = vmin; v <= vmax; v++) {
      os << v;
      for (std::size_t i = 0; i < V.size(); i++) {
        auto const& p = V[i].find(v);
        if (p != V[i].end()) {
          os << o->csvSep << p_scale * p->second;
        } else {
          os << o->csvSep << '0';
        }
      }
      os << std::endl;
    }

    os.close();

    if (ti->do_cmp()) {
      const std::size_t timepoint_index = dmp->get_timepoint_index();

      static char tmp[1024];
      sprintf(tmp, "%s/%s", ti->get_cmp_path().c_str(), dmp->get_local_fname());
      std::ifstream is_cmp(tmp);
      if (is_cmp.is_open()) {
        sprintf(tmp, "%s/nrm.csv", dmp->get_task()->get_local_dir());
        std::ofstream os_nrm_total(tmp, std::ios::out | std::ios::app);
        os_nrm_total.precision(16);
        if (os_nrm_total.is_open()) {
          sprintf(tmp, "%s/nrm_%lu.csv", dmp->get_task()->get_local_dir(),
                  timepoint_index);
          std::ofstream os_nrm(tmp);
          os_nrm.precision(16);
          if (os_nrm.is_open()) {
            double d = 0.0, d2 = 0.0, d_max = 0.0;
            double r = 0.0, r2 = 0.0, r_max = 0.0;
            double d_r_max = 0.0;
            std::string line;

            getline(is_cmp, line);
            while (getline(is_cmp, line)) {
              std::istringstream iss(line.c_str());
              std::string s;
              getline(iss, s, ',');
              const int v = atoi(s.c_str());
              os_nrm << v;
              std::size_t vi = 0;
              while (getline(iss, s, ',') && vi < V.size()) {
                const double p = atof(s.c_str());
                const double p2 = (V[vi].count(v)) ? p_scale * V[vi][v] : 0.0;
                os_nrm << o->csvSep << std::fabs(p - p2);
                o->sl->update_err(p, p2, d, d2, d_max, r, r2, r_max, d_r_max);
                vi++;
              }
              os_nrm << std::endl;
            }
            os_nrm.close();

            os_nrm_total << timepoint_index << o->csvSep;
            o->sl->write_nrm(os_nrm_total, d, d2, d_max, r, r2, r_max, d_r_max)
                << std::endl;
          } else {
            std::cerr << "Couldn't create file '" << tmp << "'" << std::endl;
          }
          os_nrm_total.close();
        } else {
          std::cerr << "Couldn't create file '" << tmp << "'" << std::endl;
        }
        is_cmp.close();
      } else {
        std::cerr << "Couldn't open file '" << tmp << "'" << std::endl;
      }
    }

    return true;
  }

  bool dump_plot_distr_2d_data(
      solver_loader::task_info::plot_distr_2d const* const ti,
      dump const* const dmp, double p_scale = 1.0) const {
    solver_loader::model_info::var const* const v1 =
        static_cast<solver_loader::parser::varAST const* const>(ti->get_et1())
            ->get_var();
    solver_loader::model_info::var const* const v2 =
        static_cast<solver_loader::parser::varAST const* const>(ti->get_et2())
            ->get_var();

    const std::size_t v1i = v1->get_index();
    const std::size_t v2i = v2->get_index();

    int v1min = 0, v1max = 0, v2min = 0, v2max = 0;

    const int cond_index = ti->get_cond_index();
    double p_sum = 0.0;
    for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
      const double p = da->s_p(s);
      cstate const* const cs = da->hs_cs(s->hs);
      double const* const x = da->s_x(s, p);

      if (cond_index < 0 || o->exprs[cond_index](cs, x, o->param_vals.data())) {
        p_sum += p;

        int const* const v = da->cs_get(cs);

        const int v1_ = (v1i < o->npvars) ? x[v1i] : v[v1i];
        const int v2_ = (v2i < o->npvars) ? x[v2i] : v[v2i];

        if (s != this->get_first()) {
          if (v1_ < v1min) {
            v1min = v1_;
          }
          if (v1_ > v1max) {
            v1max = v1_;
          }
          if (v2_ < v2min) {
            v2min = v2_;
          }
          if (v2_ > v2max) {
            v2max = v2_;
          }
        } else {
          v1min = v1_;
          v1max = v1_;
          v2min = v2_;
          v2max = v2_;
        }
      }
    }

    if (v1max == v1min) {
      v1min = std::max(0, v1min - 20);
      v1max += 20;
    }
    if (v2max == v2min) {
      v2min = std::max(0, v2min - 20);
      v2max += 20;
    }
    if (v1max > 1000) {
      v1max = 1000;
    }
    if (v2max > 1000) {
      v2max = 1000;
    }

    v1min = 0;
    v2min = 0;

    p_scale /= p_sum;

    std::vector<std::vector<double>> sp(
        v1max - v1min + 1, std::vector<double>(v2max - v2min + 1, 0.0));
    for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
      void const* const ls_d = da->get_lstate_data(s);
      double const* const y = da->get_lstate_da()->y(ls_d);
      const double p = o->v_p(y);

      cstate const* const cs = da->hs_cs(s->hs);
      double const* const x = da->s_x(s, p);

      if (cond_index < 0 || o->exprs[cond_index](cs, x, o->param_vals.data())) {
        int const* const v = da->cs_get(cs);

        const int v1_ = (v1i < o->npvars) ? x[v1i] : v[v1i];
        const int v2_ = (v2i < o->npvars) ? x[v2i] : v[v2i];

        if (v1_ >= v1min && v1_ <= v1max && v2_ >= v2min && v2_ <= v2max) {
          const double val =
              (ti->get_dparami() == nullptr)
                  ? p
                  : ((ti->get_dparamj() == nullptr))
                        ? o->v_dp_dc(y, ti->get_dparami()->get_index_p())
                        : o->v_d2p_dc2(y,
                                       o->CI(ti->get_dparami()->get_index_p(),
                                             ti->get_dparamj()->get_index_p()));

          sp[v1_ - v1min][v2_ - v2min] += val;
        }
      }
    }

    std::ofstream os(dmp->get_local_path());
    for (int j = v2min; j <= v2max; j++) {
      os << sp[0][j - v2min] * p_scale;
      for (int i = v1min + 1; i <= v1max; i++) {
        os << o->csvSep << sp[i - v1min][j - v2min] * p_scale;
      }
      os << std::endl;
    }
    os.close();

    return true;
  }

  void distr_err_precomp() {
    for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
      da->s_p(s) = -da->s_p(s);
    }
  }

  void distr_err_update(int const* const v, const double p, double& d,
                        double& d2, double& d_max, double& r, double& r2,
                        double& r_max, double& d_r_max) {
    hstate* const hs = this->get_hstate_allocator()->get_last();
    o->sl->get_model()->cs_set(da->hs_cs(hs), v);
    da->hs_update_hash(hs);
    state* const s = this->find(hs);
    if (s != nullptr) {
      const double sp = -da->s_p(s);
      if (sp > 0.0) {
        da->s_p(s) = sp;
        o->sl->update_err(p, sp, d, d2, d_max, r, r2, r_max, d_r_max);
      }
    } else {
      o->sl->update_err(p, 0.0, d, d2, d_max, r, r2, r_max, d_r_max);
    }
  }

  void distr_err_postcomp(double& d, double& d2, double& d_max, double& r,
                          double& r2, double& r_max, double& d_r_max) {
    for (state* s = this->get_first(); s != nullptr; s = s->lnext) {
      const double p = -da->s_p(s);
      if (p > 0.0) {
        da->s_p(s) = p;
        o->sl->update_err(p, 0.0, d, d2, d_max, r, r2, r_max, d_r_max);
      }
    }
  }

 private:
  options const* const o;

  lstate_data_accessor_t const* const lstate_da;
  hstate_data_accessor_t const* const hstate_da;
  hstate_succ_data_accessor_t const* const hstate_succ_da;

  container::allocator<lstate_t>* const lstate_allocator;
  container::allocator<hstate_t>* const hstate_allocator;
  container::allocator<hstate_succ_t>* const hstate_succ_allocator;

  dataAccessor const* const da;

  mutable std::vector<double> cov_;

  mutable std::vector<int> iv_;

  mutable std::vector<std::size_t> sI, sK, sIk;
  mutable std::vector<std::size_t> sI_, sJ_, sK_, sIk_;
};
}

#endif
