/*
 *  type.hpp
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

#ifndef SOLVER_LOADER_MODEL_INFO_TYPE_HPP_
#define SOLVER_LOADER_MODEL_INFO_TYPE_HPP_

#include "obj.hpp"
#include "val.hpp"

namespace solver_loader {

namespace model_info {

enum type_t { T_AUTO, T_BOOLEAN, T_INTEGER, T_REAL };

class type : public obj {
 private:
  const int length;

 public:
  type(const std::string& name, trsys const* const ts, const int length = -1)
      : obj(name, ts), length(length) {}
  virtual ~type() {}

  int get_length() const { return length; }

  bool is_array() const { return length > 0; }

  type_t get_type() const {
    if (is_boolean()) {
      return T_BOOLEAN;
    }
    if (is_integer()) {
      return T_INTEGER;
    }
    if (is_real()) {
      return T_REAL;
    }
    return T_AUTO;
  }

  virtual bool is_alias() const { return false; }

  virtual bool is_boolean() const { return false; }

  virtual bool is_integer() const { return false; }

  virtual bool is_species() const { return false; }

  virtual bool is_subrange() const { return false; }

  virtual bool is_real() const { return false; }

  virtual bool is_ordinal() const { return is_boolean() || is_integer(); }
};

class type_alias : public type {
 private:
  type const* const typ;

 public:
  type_alias(const std::string& name, trsys const* const ts,
             type const* const typ)
      : type(name, ts), typ(typ) {}

  bool is_alias() const { return true; }

  type const* get_type() const { return typ; }

  bool is_ordinal() const { return typ->is_ordinal(); }
};

class type_boolean : public type {
 public:
  type_boolean(const std::string& name, trsys const* const ts)
      : type(name, ts) {}

  bool is_boolean() const { return true; }
};

class type_integer : public type {
 public:
  type_integer(const std::string& name, trsys const* const ts)
      : type(name, ts) {}

  bool is_integer() const { return true; }
};

class type_species : public type_integer {
 public:
  type_species(const std::string& name, trsys const* const ts)
      : type_integer(name, ts) {}

  bool is_species() const { return true; }
};

class type_subrange : public type_integer {
 private:
  val* const minvalue;
  val* const maxvalue;

 public:
  type_subrange(const std::string& name, trsys const* const ts,
                val* const minvalue = nullptr, val* const maxvalue = nullptr)
      : type_integer(name, ts), minvalue(minvalue), maxvalue(maxvalue) {}

  val* get_minvalue() const { return minvalue; }

  val* get_maxvalue() const { return maxvalue; }

  bool is_subrange() const { return true; }
};

class type_real : public type {
 public:
  type_real(const std::string& name, trsys const* const ts) : type(name, ts) {}

  bool is_real() const { return true; }
};
}
}

#endif
