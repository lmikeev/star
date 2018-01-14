/*
 *  allocator.hpp
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

#ifndef SOLVER_CONTAINER_ALLOCATOR_HPP_
#define SOLVER_CONTAINER_ALLOCATOR_HPP_

#include <boost/pool/pool.hpp>
#include "../options.hpp"

namespace solver {
namespace container {

class data_accessor {
 protected:
  options const* const o;

  const std::size_t size;

 public:
  data_accessor(options const* const o, const std::size_t size)
      : o(o), size(size) {}
  virtual ~data_accessor() {}

  std::size_t get_size() const { return size; }
};

class hs_data_accessor : public data_accessor {
 public:
  hs_data_accessor(options const* const o, const std::size_t size)
      : data_accessor(o, size) {}
  virtual ~hs_data_accessor() {}

  virtual std::size_t get_cs_size() const = 0;
  virtual bool cs_equal(void const* const d1, void const* const d2) const = 0;
};

template <typename T>
class allocator {
 private:
  mutable boost::pool<>* p;

  mutable T* last;

  void malloc_() const {
    last = static_cast<T*>(p->malloc());
    memset(last, 0, p->get_requested_size());
  }

 public:
  allocator(data_accessor const* const da)
      : p(new boost::pool<>(sizeof(T) + da->get_size())) {
    malloc_();
  }

  allocator() : p(new boost::pool<>(sizeof(T))) { malloc_(); }

  ~allocator() {}

  std::size_t get_requested_size() const { return p->get_requested_size(); }

  T* get_last() const { return last; }

  T* malloc() const {
    T* const ret = last;
    malloc_();
    return ret;
  }

  static void* get_data(T* const h) { return static_cast<void*>(h + 1); }

  static void const* get_data(T const* const h) {
    return static_cast<void const*>(h + 1);
  }
};
}
}

#endif
