/*
 *  ast.hpp
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

#ifndef SOLVER_LOADER_PARSER_AST_HPP_
#define SOLVER_LOADER_PARSER_AST_HPP_

#include <vector>
#include <iostream>
#include <map>
#include <cassert>
#include "../model_info/fnc.hpp"
#include "../model_info/type.hpp"

#if HAVE_LIBSBML
#include <sbml/SBMLTypes.h>
#endif

#ifdef LLVM_JIT
#define _GNU_SOURCE
#define __STDC_CONSTANT_MACROS
#define __STDC_FORMAT_MACROS
#define __STDC_LIMIT_MACROS
#include <llvm/Analysis/Passes.h>
#include <llvm/ExecutionEngine/ExecutionEngine.h>
#include <llvm/ExecutionEngine/JIT.h>
#include <llvm/IR/DataLayout.h>
#include <llvm/IR/DerivedTypes.h>
#include <llvm/IR/IRBuilder.h>
#include <llvm/IR/LLVMContext.h>
#include <llvm/IR/Module.h>
#include <llvm/IR/Verifier.h>
#include <llvm/PassManager.h>
#include <llvm/Support/TargetSelect.h>
#include <llvm/Transforms/Scalar.h>
#endif

namespace solver_loader {

class base;

namespace model_info {
class cnst;
class var;
class clock;
class transition;
struct fnc_arg;
}

namespace writer {
class base;
}

namespace parser {

enum token_t {
  TOK_UNKNOWN,

  TOK_ID,
  TOK_ID_PRIME,

  TOK_TYPE,
  TOK_CONST,
  TOK_VAR,
  TOK_MOLECULE,
  TOK_CLOCK,
  TOK_FUNCTION,

  TOK_STDTYPE,
  TOK_STDFNC,

  TOK_GUARDEDCMDS,
  TOK_CHEMREACTIONS,
  TOK_PTA,

  TOK_INIT,

  TOK_END,

  TOK_BOOLNUM,
  TOK_INTNUM,
  TOK_REALNUM,

  TOK_NOT,

  TOK_MULT,
  TOK_DIV,

  TOK_PLUS,
  TOK_MINUS,

  TOK_LESS,
  TOK_LEQ,
  TOK_GR,
  TOK_GREQ,

  TOK_EQ,
  TOK_INEQ,

  TOK_AND,

  TOK_OR,

  TOK_R_OP,
  TOK_R_CL,
  TOK_S_OP,
  TOK_S_CL,
  TOK_C_OP,
  TOK_C_CL,

  TOK_QUESTION,
  TOK_EXCL,
  TOK_COMMA,
  TOK_DOT,
  TOK_DDOT,
  TOK_COLON,
  TOK_SEMICOLON,
  TOK_RARROW,
  TOK_DARROW,
  TOK_AT,
  TOK_TILDE,

  TOK_EOF
};

enum stdfnc_id {
  STDF_TIME = 0,
  STDF_MASS_ACTION,
  STDF_ABS,
  STDF_MIN,
  STDF_MAX,
  STDF_FLOOR,
  STDF_CEIL,
  STDF_MOD,
  STDF_EXP,
  STDF_SQRT,
  STDF_POW,
  STDF_LN,
  STDF_LOG10,
  STDF_LOG,
  STDF_SIN,
  STDF_COS,
  STDF_TAN,
  STDF_ASIN,
  STDF_ACOS,
  STDF_ATAN,
  STDF_ATAN2,
  STDF_SINH,
  STDF_COSH,
  STDF_TANH,
  STDF_ASINH,
  STDF_ACOSH,
  STDF_ATANH
};

struct stdfnc {
  stdfnc_id const id;
  char const* name;
  int const narg;
  model_info::type_t const rtyp;

  stdfnc_id get_id() const { return id; }

  char const* get_name() const { return name; }

  int get_narg() const { return narg; }

  model_info::type_t get_rtype() const { return rtyp; }
};

extern writer::base const* const std_writer;
extern writer::base const* const cpp_writer;
extern writer::base const* const matlab_writer;

class AST {
 public:
  virtual ~AST() {}

#ifdef LLVM_JIT
  virtual llvm::Value* codegen(llvm::IRBuilder<>& builder,
                               std::stringstream& last_err) const = 0;
#endif

  bool is_leaf() const { return is_number() || is_cnst() || is_var(); }

  bool is_integral() const { return is_boolnum() || is_intnum(); }

  bool is_number() const { return is_integral() || is_realnum(); }

  bool is_fcall() const { return is_stdfcall() || is_userfcall(); }

  virtual bool is_boolnum() const { return false; }

  virtual bool is_intnum() const { return false; }

  virtual bool is_realnum() const { return false; }

  virtual bool is_cnst() const { return false; }

  virtual bool is_var() const { return false; }

  virtual bool is_clock() const { return false; }

  virtual bool is_unary() const { return false; }

  virtual bool is_binary() const { return false; }

  virtual bool is_polynom() const { return false; }

  virtual bool is_stdfcall() const { return false; }

  virtual bool is_userfcall() const { return false; }

  virtual bool is_userfarg() const { return false; }

  virtual bool is_ifelse() const { return false; }

  virtual std::ostream& write(
      std::ostream& os, writer::base const* const w = std_writer,
      solver_loader::base const* const sl = nullptr) const = 0;

  static AST const* const et_0;
  static AST const* const et_1;

  static stdfnc const* get_stdf(const stdfnc_id id);
};

class unaryAST : public AST {
 private:
  const token_t op_tok;
  AST const* const term;

 public:
  unaryAST(const token_t op_tok, AST const* const term)
      : op_tok(op_tok), term(term) {}

  unaryAST(unaryAST const* const et)
      : unaryAST(et->get_op_tok(), et->get_term()) {}

  int get_op() const {
    switch (op_tok) {
      case TOK_PLUS:
        return 0;

      case TOK_MINUS:
        return 1;

      case TOK_NOT:
        return 2;

      default:
        break;
    }
    assert(false);
    return -1;
  }

  token_t get_op_tok() const { return op_tok; }

  AST const* get_term() const { return term; }

  bool is_unary() const { return true; }

  std::ostream& write(std::ostream& os, writer::base const* const w,
                      solver_loader::base const* const sl = nullptr) const;

#ifdef LLVM_JIT
  llvm::Value* codegen(llvm::IRBuilder<>& builder,
                       std::stringstream& last_err) const {
    llvm::Value* const v = term->codegen(builder, last_err);
    if (v == nullptr) {
      return nullptr;
    }

    switch (get_op_tok()) {
      case TOK_PLUS:
        return v;

      case TOK_MINUS:
        return builder.CreateFNeg(v);

      case TOK_NOT:
        return builder.CreateNot(v);

      default:
        last_err << "not supported unary operation (" << get_op_tok() << ")";
    }

    return nullptr;
  }
#endif
};

class binaryAST : public AST {
 private:
  const token_t op_tok;
  AST const* const lterm;
  AST const* const rterm;

 public:
  binaryAST(const token_t op_tok, AST const* const lterm,
            AST const* const rterm)
      : op_tok(op_tok), lterm(lterm), rterm(rterm) {}

  binaryAST(binaryAST const* const et)
      : binaryAST(et->get_op_tok(), et->get_lterm(), et->get_rterm()) {}

  int get_op() const { return op_tok - TOK_MULT; }

  token_t get_op_tok() const { return op_tok; }

  AST const* get_lterm() const { return lterm; }

  AST const* get_rterm() const { return rterm; }

  bool is_binary() const { return true; }

  virtual std::ostream& write(
      std::ostream& os, writer::base const* const w,
      solver_loader::base const* const sl = nullptr) const;

#ifdef LLVM_JIT
  llvm::Value* codegen(llvm::IRBuilder<>& builder,
                       std::stringstream& last_err) const {
    llvm::Value* const lv = lterm->codegen(builder, last_err);
    if (lv == nullptr) {
      return nullptr;
    }

    llvm::Value* const rv = rterm->codegen(builder, last_err);
    if (rv == nullptr) {
      return nullptr;
    }

    switch (get_op_tok()) {
      case TOK_MULT:
        return builder.CreateFMul(lv, rv);

      case TOK_DIV:
        return builder.CreateFDiv(lv, rv);

      case TOK_PLUS:
        return builder.CreateFAdd(lv, rv);

      case TOK_MINUS:
        return builder.CreateFSub(lv, rv);

      case TOK_AND:
        return builder.CreateAnd(lv, rv);

      case TOK_OR:
        return builder.CreateOr(lv, rv);

      default:
        last_err << "not supported binary operation (" << get_op_tok() << ")";
    }

    return nullptr;
  }
#endif
};

class boolnumAST : public AST {
 private:
  const bool value;

 public:
  boolnumAST(const bool value) : value(value) {}

  boolnumAST(boolnumAST const* const et) : boolnumAST(et->get_value()) {}

  bool get_value() const { return value; }

  bool is_boolnum() const { return true; }

  std::ostream& write(std::ostream& os, writer::base const* const w,
                      solver_loader::base const* const sl = nullptr) const;

#ifdef LLVM_JIT
  llvm::Value* codegen(llvm::IRBuilder<>& builder,
                       std::stringstream& last_err) const {
    return llvm::ConstantInt::get(llvm::getGlobalContext(),
                                  llvm::APInt(32, value, true));
  }
#endif
};

class intnumAST : public AST {
 private:
  const int value;

 public:
  intnumAST(const int value) : value(value) {}

  intnumAST(intnumAST const* const et) : intnumAST(et->get_value()) {}

  int get_value() const { return value; }

  bool is_intnum() const { return true; }

  std::ostream& write(std::ostream& os, writer::base const* const w,
                      solver_loader::base const* const sl = nullptr) const;

#ifdef LLVM_JIT
  llvm::Value* codegen(llvm::IRBuilder<>& builder,
                       std::stringstream& last_err) const {
    return llvm::ConstantInt::get(llvm::getGlobalContext(),
                                  llvm::APInt(32, value, true));
  }
#endif
};

class realnumAST : public AST {
 private:
  const double value;

 public:
  realnumAST(const double value) : value(value) {}

  realnumAST(realnumAST const* const et) : realnumAST(et->get_value()) {}

  double get_value() const { return value; }

  bool is_realnum() const { return true; }

  std::ostream& write(std::ostream& os, writer::base const* const w,
                      solver_loader::base const* const sl = nullptr) const;

#ifdef LLVM_JIT
  llvm::Value* codegen(llvm::IRBuilder<>& builder,
                       std::stringstream& last_err) const {
    return llvm::ConstantFP::get(llvm::getGlobalContext(),
                                 llvm::APFloat(value));
  }
#endif
};

class cnstAST : public AST {
 private:
  model_info::cnst const* const c;

 public:
  cnstAST(model_info::cnst const* const c) : c(c) {}

  cnstAST(cnstAST const* const et) : cnstAST(et->get_cnst()) {}

  bool is_cnst() const { return true; }

  model_info::cnst const* get_cnst() const { return c; }

  std::ostream& write(std::ostream& os, writer::base const* const w,
                      solver_loader::base const* const sl = nullptr) const;

#ifdef LLVM_JIT
  llvm::Value* codegen(llvm::IRBuilder<>& builder,
                       std::stringstream& last_err) const {}
#endif
};

class varAST : public AST {
 private:
  model_info::var const* const v;
  const bool isprime;

 public:
  varAST(model_info::var const* const v, const bool isprime = false)
      : v(v), isprime(isprime) {}

  varAST(varAST const* const et) : varAST(et->get_var(), et->is_prime()) {}

  bool is_var() const { return true; }

  model_info::var const* get_var() const { return v; }

  bool is_prime() const { return isprime; }

  std::ostream& write(std::ostream& os, writer::base const* const w,
                      solver_loader::base const* const sl = nullptr) const;

#ifdef LLVM_JIT
  llvm::Value* codegen(llvm::IRBuilder<>& builder,
                       std::stringstream& last_err) const {}
#endif
};

class clockAST : public AST {
 private:
  model_info::clock const* const cl;

 public:
  clockAST(model_info::clock const* const cl) : cl(cl) {}

  clockAST(clockAST const* const et) : clockAST(et->get_clock()) {}

  bool is_clock() const { return true; }

  model_info::clock const* get_clock() const { return cl; }

  std::ostream& write(std::ostream& os, writer::base const* const w,
                      solver_loader::base const* const sl = nullptr) const;

#ifdef LLVM_JIT
  llvm::Value* codegen(llvm::IRBuilder<>& builder,
                       std::stringstream& last_err) const {}
#endif
};

struct monom {
  double k;
  std::map<model_info::cnst const*, int> c;
  std::map<model_info::var const*, int> v;

  monom(const double k = 1.0) : k(k) {}

  monom(const monom& a) : k(a.k), c(a.c), v(a.v) {}

  monom(const monom& a, const monom& b) : k(a.k * b.k) {
    for (auto const& ci : a.c) {
      c[ci.first] += ci.second;
    }
    for (auto const& vi : a.v) {
      v[vi.first] += vi.second;
    }
    for (auto const& ci : b.c) {
      c[ci.first] += ci.second;
    }
    for (auto const& vi : b.v) {
      v[vi.first] += vi.second;
    }
  }

  bool cmp(const monom& b) const {
    if (c.size() != b.c.size() || v.size() != b.v.size()) {
      return false;
    }

    auto aci = c.begin();
    auto bci = b.c.begin();
    while (aci != c.end()) {
      if (aci->first != bci->first || aci->second != bci->second) {
        return false;
      }
      aci++;
      bci++;
    }

    assert(bci == b.c.end());

    auto avi = v.begin();
    auto bvi = b.v.begin();
    while (avi != v.end()) {
      if (avi->first != bvi->first || avi->second != bvi->second) {
        return false;
      }
      avi++;
      bvi++;
    }

    assert(bvi == b.v.end());

    return true;
  }
};

class polynomAST : public AST {
 protected:
  std::vector<monom> terms;

 public:
  polynomAST() {}

  polynomAST(const double c) {
    monom a(c);
    terms.push_back(a);
  }

  polynomAST(const std::vector<monom>& terms) : terms(terms) {}

  polynomAST(polynomAST const* const et) : polynomAST(et->get_terms()) {}

  bool is_polynom() const { return true; }

  std::vector<monom>& get_terms() { return terms; }

  const std::vector<monom>& get_terms() const { return terms; }

  std::ostream& write(std::ostream& os, writer::base const* const w,
                      solver_loader::base const* const sl = nullptr) const;

  bool is_const() const {
    return terms.size() == 1 && terms.front().c.empty() &&
           terms.front().v.empty();
  }

#ifdef LLVM_JIT
  llvm::Value* codegen(llvm::IRBuilder<>& builder,
                       std::stringstream& last_err) {}
#endif
};

class callfAST : public AST {
 protected:
  const std::vector<AST const*> args;

 public:
  callfAST(const std::vector<AST const*>& args) : args(args) {}

  callfAST(callfAST const* const et) : callfAST(et->get_args()) {}

  const std::vector<AST const*>& get_args() const { return args; }

#ifdef LLVM_JIT
  llvm::Value* codegen(llvm::IRBuilder<>& builder,
                       std::stringstream& last_err) {}
#endif
};

class userfcallAST : public callfAST {
 private:
  model_info::fnc const* const f;
  const std::vector<AST const*> args;

 public:
  userfcallAST(model_info::fnc const* const f,
               const std::vector<AST const*>& args)
      : callfAST(args), f(f), args(args) {}

  userfcallAST(userfcallAST const* const et)
      : userfcallAST(et->get_f(), et->get_args()) {}

  bool is_userfcall() const { return true; }

  model_info::fnc const* get_f() const { return f; }

  const std::vector<AST const*>& get_args() const { return args; }

  std::ostream& write(std::ostream& os, writer::base const* const w,
                      solver_loader::base const* const sl = nullptr) const;

#ifdef LLVM_JIT
  llvm::Value* codegen(llvm::IRBuilder<>& builder,
                       std::stringstream& last_err) {}
#endif
};

class userfargAST : public AST {
 private:
  std::vector<model_info::fnc_arg> const* const fargs;
  const std::size_t index;

 public:
  userfargAST(std::vector<model_info::fnc_arg> const* const fargs,
              const std::size_t index)
      : fargs(fargs), index(index) {}

  userfargAST(userfargAST const* const et)
      : userfargAST(et->get_fargs(), et->get_index()) {}

  bool is_userfarg() const { return true; }

  std::vector<model_info::fnc_arg> const* get_fargs() const { return fargs; }

  std::size_t get_index() const { return index; }

  std::ostream& write(std::ostream& os, writer::base const* const w,
                      solver_loader::base const* const sl = nullptr) const;
};

class stdfcallAST : public callfAST {
 private:
  stdfnc const* const stdf;
  const std::vector<AST const*> args;

 public:
  stdfcallAST(stdfnc const* const stdf, const std::vector<AST const*>& args)
      : callfAST(args), stdf(stdf), args(args) {}

  stdfcallAST(stdfcallAST const* const et)
      : stdfcallAST(et->get_stdf(), et->get_args()) {}

  bool is_stdfcall() const { return true; }

  stdfnc const* get_stdf() const { return stdf; }

  const std::vector<AST const*>& get_args() const { return args; }

  std::ostream& write(std::ostream& os, writer::base const* const w,
                      solver_loader::base const* const sl = nullptr) const;

#ifdef LLVM_JIT
  llvm::Value* codegen(llvm::IRBuilder<>& builder,
                       std::stringstream& last_err) {}
#endif
};

class ifelseAST : public AST {
 private:
  AST const* const cond;
  AST const* const then_term;
  AST const* const else_term;

 public:
  ifelseAST(AST const* const cond, AST const* const then_term,
            AST const* const else_term)
      : cond(cond), then_term(then_term), else_term(else_term) {}

  ifelseAST(ifelseAST const* const et)
      : ifelseAST(et->get_cond(), et->get_then_term(), et->get_else_term()) {}

  bool is_ifelse() const { return true; }

  AST const* get_cond() const { return cond; }

  AST const* get_then_term() const { return then_term; }

  AST const* get_else_term() const { return else_term; }

  std::ostream& write(std::ostream& os, writer::base const* const w,
                      solver_loader::base const* const sl = nullptr) const;

#ifdef LLVM_JIT
  llvm::Value* codegen(llvm::IRBuilder<>& builder,
                       std::stringstream& last_err) {}
#endif
};

AST const* build_expr(AST const* const et,
                      model_info::transition const* const tr = nullptr,
                      userfcallAST const* const etf = nullptr);
AST const* build_product(AST const* const et, model_info::var const* const v,
                         model_info::transition const* const tr = nullptr,
                         userfcallAST const* const etf = nullptr);
AST const* build_product(AST const* const etl, AST const* const etr,
                         model_info::transition const* const tr = nullptr,
                         userfcallAST const* const etf = nullptr);
AST const* build_fraction(AST const* const etn, AST const* const etd,
                          model_info::transition const* const tr = nullptr,
                          userfcallAST const* const etf = nullptr);
AST const* build_sum(AST const* const etl, AST const* const etr,
                     model_info::transition const* const tr = nullptr,
                     userfcallAST const* const etf = nullptr);
AST const* build_difference(AST const* const etl, parser::AST const* const etr,
                            model_info::transition const* const tr = nullptr,
                            userfcallAST const* const etf = nullptr);
AST const* build_negative(AST const* const et,
                          model_info::transition const* const tr = nullptr,
                          userfcallAST const* const etf = nullptr);
AST const* build_power(AST const* const et, AST const* const etp,
                       model_info::transition const* const tr = nullptr,
                       userfcallAST const* const etf = nullptr);

#if HAVE_LIBSBML

struct lparam_s {
  model_info::cnst const* c;
  std::string name;
};

bool sbml_read_expr(
    ASTNode const* const anode, model_info::trsys const* const context,
    const int lr, AST const*& et, std::stringstream& last_err,
    std::vector<model_info::fnc_arg> const* const fargs = nullptr,
    std::vector<lparam_s> const* const lparams = nullptr);

#endif
}

void compute_value(model_info::val* const value, const bool comp_init = false);
}

#endif
