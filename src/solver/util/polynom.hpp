/*
 *  polynom.hpp
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

#ifndef SOLVER_UTIL_POLYNOM_HPP_
#define SOLVER_UTIL_POLYNOM_HPP_

#include <iostream>
#include <ostream>
#include <list>
#include "../model/var.hpp"

namespace solver {
namespace util {

struct monomel {
  model::var const* v;
  unsigned int d;
};
typedef std::list<monomel> monomel_list;
typedef monomel_list::iterator monomel_list_iter;
typedef monomel_list::const_iterator monomel_list_citer;

struct monom {
  double c;
  monomel_list e;
};
typedef std::list<monom> monom_list;
typedef monom_list::iterator monom_list_iter;
typedef monom_list::const_iterator monom_list_citer;

class polynom {
 private:
  monom_list M;

  monom* getmonom(const monom& m) const {
    for (monom_list_citer m_ = M.begin(); m_ != M.end(); m_++) {
      bool f_ = true;
      if (!m.e.empty()) {
        for (monomel_list_citer e_ = m.e.begin(); e_ != m.e.end(); e_++) {
          bool f__ = true;
          for (monomel_list_citer e2 = m_->e.begin(); e2 != m_->e.end(); e2++) {
            if (e_->v == e2->v && e_->d == e2->d) {
              f__ = false;
              break;
            }
          }
          if (f__) {
            f_ = false;
            break;
          }
        }
      } else {
        if (!m_->e.empty()) {
          f_ = false;
        }
      }

      if (f_) {
        return const_cast<monom*>(&(*m_));
      }
    }

    return nullptr;
  }

 public:
  polynom() {}

  polynom(const double c) {
    monom m_ = {c};
    M.push_back(m_);
  }

  ~polynom() { M.clear(); }

  const monom_list& get_monoms() const { return M; }

  void addmonom(const monom& m2) {
    if (m2.c != 0.0) {
      monom* m = getmonom(m2);
      if (m != nullptr) {
        m->c += m2.c;
      } else {
        M.push_back(m2);
      }
    }
  }

  void add(const polynom& p2) {
    for (monom_list_citer m2 = p2.get_monoms().begin();
         m2 != p2.get_monoms().end(); m2++) {
      addmonom(*m2);
    }
  }

  void subtract(const polynom& p2) {
    for (monom_list_citer m2 = p2.get_monoms().begin();
         m2 != p2.get_monoms().end(); m2++) {
      monom* m = getmonom(*m2);
      if (m != nullptr) {
        m->c -= m2->c;
      } else {
        monom m(*m2);
        m.c = -m.c;
        M.push_back(m);
      }
    }
  }

  void multiply(const double c) {
    if (c != 0.0) {
      for (monom_list_iter m_ = M.begin(); m_ != M.end(); m_++) {
        m_->c *= c;
      }
    } else {
      M.clear();
      monom m_ = {c};
      M.push_back(m_);
    }
  }

  void multiply(const polynom& p2) {
    monom_list M_(M);
    M.clear();
    for (monom_list_citer m1 = M_.begin(); m1 != M_.end(); m1++) {
      for (monom_list_citer m2 = p2.get_monoms().begin();
           m2 != p2.get_monoms().end(); m2++) {
        monom m_(*m1);
        m_.c *= m2->c;
        for (monomel_list_citer e2 = m2->e.begin(); e2 != m2->e.end(); e2++) {
          bool f_ = true;
          for (monomel_list_iter e_ = m_.e.begin(); e_ != m_.e.end(); e_++) {
            if (e_->v == e2->v) {
              e_->d += e2->d;
              f_ = false;
              break;
            }
          }
          if (f_) {
            m_.e.push_back(*e2);
          }
        }

        addmonom(m_);
      }
    }
  }

  std::ostream& write(std::ostream& os) const {
    if (!M.empty()) {
      bool f_ = false;
      for (monom_list_citer m_ = M.begin(); m_ != M.end(); m_++) {
        if (m_->c != 0.0) {
          if (f_ && m_->c > 0.0) {
            os << '+';
          }
          os << m_->c;
          for (monomel_list_citer e_ = m_->e.begin(); e_ != m_->e.end(); e_++) {
            if (e_->d > 0) {
              os << "*x" << e_->v->get_index() + 1;
              if (e_->d > 1) {
                os << "^" << e_->d;
              }
            }
          }
          f_ = true;
        }
      }
    } else {
      os << "0";
    }

    return os;
  }

  std::ostream& writeD(std::ostream& os, model::var const* const v) const {
    bool f_ = true;
    if (!M.empty()) {
      bool f___ = false;
      for (monom_list_citer m_ = M.begin(); m_ != M.end(); m_++) {
        if (m_->c != 0.0) {
          bool f__ = false;
          for (monomel_list_citer e_ = m_->e.begin(); e_ != m_->e.end(); e_++) {
            if (e_->v == v && e_->d > 0) {
              f__ = true;
              break;
            }
          }

          if (f__) {
            f_ = false;

            if (f___ && m_->c > 0.0) {
              os << '+';
            }
            os << m_->c;
            for (monomel_list_citer e_ = m_->e.begin(); e_ != m_->e.end();
                 e_++) {
              if (e_->v != v) {
                if (e_->d > 0) {
                  os << "*x" << e_->v->get_index() + 1;
                  if (e_->d > 1) {
                    os << "^" << e_->d;
                  }
                }
              } else {
                if (e_->d > 1) {
                  os << "*x" << e_->v->get_index() + 1;
                  if (e_->d > 2) {
                    os << "^" << e_->d - 1;
                  }
                }
              }
            }
            f___ = true;
          }
        }
      }
    }

    if (f_) {
      os << "0";
    }

    return os;
  }
};
}
}

#endif
