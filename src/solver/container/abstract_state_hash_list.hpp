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

#ifndef SOLVER_CONTAINER_ABSTRACT_STATE_HASH_LIST_HPP_
#define SOLVER_CONTAINER_ABSTRACT_STATE_HASH_LIST_HPP_

#include <cstdint>
#include "allocator.hpp"
#include "../../solver_loader/base.hpp"

namespace solver {
namespace container {

template <typename lstate_t, typename hstate_t, typename hstate_succ_t,
          std::size_t ht_len = 1 << 20>
class abstractStateHashListT {
 public:
  abstractStateHashListT(allocator<lstate_t>* const lstate_allocator,
                         allocator<hstate_t>* const hstate_allocator,
                         allocator<hstate_succ_t>* const hstate_succ_allocator,
                         hs_data_accessor const* const hstate_da)
      : lstate_allocator(lstate_allocator),
        hstate_allocator(hstate_allocator),
        hstate_succ_allocator(hstate_succ_allocator),
        hstate_da(hstate_da),
        first(nullptr),
        last(nullptr),
        first_erased(nullptr),
        first_temp(nullptr),
        nactive(0),
        ntotal(0),
        ht_len1(ht_len - 1) {
    memset(HT, 0, sizeof(HT));
  }

  ~abstractStateHashListT() {}

  hstate_t* get_hs_first() const { return hs_first; }

  lstate_t* get_first() const { return first; }

  lstate_t* get_last() const { return last; }

  bool empty() const { return first == nullptr; }

  lstate_t* get_first_erased() const { return first_erased; }

  lstate_t* get_first_temp() const { return first_temp; }

  std::size_t get_nactive() const { return nactive; }

  std::size_t get_ntotal() const { return ntotal; }

  std::size_t get_ht_len() const { return ht_len; }

  allocator<lstate_t>* get_lstate_allocator() const { return lstate_allocator; }

  allocator<hstate_t>* get_hstate_allocator() const { return hstate_allocator; }

  allocator<hstate_succ_t>* get_hstate_succ_allocator() const {
    return hstate_succ_allocator;
  }

  lstate_t* add(hstate_t* const hs) {
    lstate_t* s;

    if (first_erased != nullptr) {
      s = first_erased;
      first_erased = first_erased->lnext;
      if (first_erased != nullptr) {
        first_erased->lprev = nullptr;
      }
    } else {
      s = lstate_allocator->malloc();

      ntotal++;
    }

    hs->ls = s;
    s->hs = hs;

    add_last_(s);

    return s;
  }

  lstate_t* add(hstate_t* const hs, bool& fnew) {
    hstate_t*& ht = HT[hs->hash & ht_len1];
    hstate_t* hs_ = ht;
    while (hs_ != nullptr) {
      if (hs_equal(hs_, hs)) {
        fnew = false;
        if (hs_->ls != nullptr) {
          return hs_->ls;
        }
        return add(hs_);
      }
      hs_ = hs_->hnext;
    }

    fnew = true;
    hs->hnext = ht;
    ht = hs;
    hs->lnext = hs_first;
    hs_first = hs;
    return add(hs);
  }

  lstate_t* find(hstate_t* const hs) {
    hstate_t*& ht = HT[hs->hash & ht_len1];
    hstate_t* hs_ = ht;
    while (hs_ != nullptr) {
      if (hs_equal(hs_, hs)) {
        if (hs_->ls != nullptr) {
          return hs_->ls;
        }
        return nullptr;
      }
      hs_ = hs_->hnext;
    }
    return nullptr;
  }

  void erase(lstate_t* const s) {
    s->hs->ls = nullptr;
    s->hs = nullptr;

    erase_(s);
    add_(s, first_erased);
  }

  void erase_all() {
    if (first != nullptr) {
      if (first_erased != nullptr) {
        last->lnext = first_erased;
        first_erased->lprev = last;
        first_erased = first;
      } else {
        first_erased = first;
      }

      first = nullptr;
      last = nullptr;

      nactive = 0;
    }
  }

  void make_temp(lstate_t* const s) {
    erase_(s);
    add_(s, first_temp);
  }

  void restore_temp(lstate_t* const s) {
    erase_(s, first_temp);
    add_last_(s);
  }

  void erase_temp(lstate_t* const s) {
    s->hs->ls = nullptr;
    s->hs = nullptr;

    erase_(s, first_temp);
    add_(s, first_erased);
  }

  void sort(
      std::function<bool(lstate_t const* const, lstate_t const* const)> cmp) {
    sort_(first, last, cmp);
  }

 private:
  void const* get_hstate_data(hstate_t const* const hs) const {
    return hstate_allocator->get_data(hs);
  }

  bool hs_equal(hstate_t const* const hs1, hstate_t const* const hs2) const {
    return hs1->hash == hs2->hash &&
           hstate_da->cs_equal(get_hstate_data(hs1), get_hstate_data(hs2));
  }

  void add_(lstate_t* const s, lstate_t*& first_) {
    s->lprev = nullptr;
    s->lnext = first_;
    if (first_ != nullptr) {
      first_->lprev = s;
    }
    first_ = s;
  }

  void erase_(lstate_t* const s, lstate_t*& first_) {
    if (s->lprev != nullptr) {
      s->lprev->lnext = s->lnext;
    } else {
      first_ = s->lnext;
    }
    if (s->lnext != nullptr) {
      s->lnext->lprev = s->lprev;
    }
  }

  void erase_(lstate_t* const s) {
    if (s->lprev != nullptr) {
      s->lprev->lnext = s->lnext;
    } else {
      first = s->lnext;
    }
    if (s->lnext != nullptr) {
      s->lnext->lprev = s->lprev;
    } else {
      last = s->lprev;
    }
    nactive--;
  }

  void add_last_(lstate_t* const s) {
    s->lnext = nullptr;
    s->lprev = last;
    if (last != nullptr) {
      last->lnext = s;
    } else {
      first = s;
    }
    last = s;
    nactive++;
  }

  void sort_(
      lstate_t*& head, lstate_t*& tail,
      std::function<bool(lstate_t const* const, lstate_t const* const)> cmp) {
    assert(head != nullptr);
    assert(tail != nullptr);

    if (head != tail) {
      lstate_t* p = head;
      lstate_t* q = head;
      lstate_t* prev = nullptr;

      while (q != nullptr) {
        prev = p;
        p = p->lnext;
        q = q->lnext;
        if (q != nullptr) {
          q = q->lnext;
        }
      }

      assert(prev != nullptr);

      prev->lnext = nullptr;
      p->lprev = nullptr;

      sort_(p, tail, cmp);
      sort_(head, prev, cmp);

      q = head;
      head = nullptr;
      prev = nullptr;

      while (p != nullptr && q != nullptr) {
        lstate_t* t;
        if (cmp(p, q)) {
          t = p;
          p = p->lnext;
        } else {
          t = q;
          q = q->lnext;
        }

        t->lprev = prev;
        if (prev != nullptr) {
          prev->lnext = t;
        } else {
          head = t;
        }
        prev = t;
      }

      while (p != nullptr) {
        p->lprev = prev;
        if (prev != nullptr) {
          prev->lnext = p;
        } else {
          head = p;
        }
        prev = p;
        p = p->lnext;
      }

      while (q != nullptr) {
        q->lprev = prev;
        if (prev != nullptr) {
          prev->lnext = q;
        } else {
          head = q;
        }
        prev = q;
        q = q->lnext;
      }

      tail = prev;
    }

    head->lprev = nullptr;
    tail->lnext = nullptr;
  }

  hstate_t* HT[ht_len];

  allocator<lstate_t>* const lstate_allocator;
  allocator<hstate_t>* const hstate_allocator;
  allocator<hstate_succ_t>* const hstate_succ_allocator;

  hs_data_accessor const* const hstate_da;

  hstate_t* hs_first;

  lstate_t* first;
  lstate_t* last;

  lstate_t* first_erased;
  lstate_t* first_temp;

  std::size_t nactive, ntotal;

  const std::size_t ht_len1;
};
}
}

#endif
