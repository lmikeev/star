/*
 *  sbml.cpp
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

#include <algorithm>
#include <cmath>
#include "sbml.hpp"
#include "../model_info/trsys.hpp"
#include "../parser/parser.hpp"

#undef isnan

namespace solver_loader {
namespace parser {

#if HAVE_LIBSBML

bool sbml_read_compartments(Model const* const model,
                            model_info::trsys* const ts, std::stringstream&) {
  ListOfCompartments const* const compartments = model->getListOfCompartments();
  for (std::size_t i = 0; i < compartments->size(); ++i) {
    Compartment const* const c = compartments->get(i);

    ts->add_cnst(new model_info::cnst(
        c->getId(), ts, model_info::trsys::typeReal,
        new model_info::val(
            new realnumAST(!std::isnan(c->getSize()) ? c->getSize() : 1.0))));
  }
  return true;
}

bool sbml_read_parameters(Model const* const model, model_info::trsys* const ts,
                          std::stringstream&) {
  ListOfParameters const* const parameters = model->getListOfParameters();
  for (std::size_t i = 0; i < parameters->size(); ++i) {
    Parameter const* const p = parameters->get(i);

    ts->add_cnst(new model_info::cnst(
        p->getId(), ts, model_info::trsys::typeReal,
        new model_info::val(new realnumAST(p->getValue()))));
  }
  return true;
}

bool sbml_read_species(Model const* const model, model_info::trsys* const ts,
                       model_info::ic_s& is, std::stringstream& last_err) {
  ListOfSpecies const* const species = model->getListOfSpecies();
  for (std::size_t i = 0; i < species->size(); ++i) {
    Species const* const s = species->get(i);
    const std::string& sid = s->getId();

    {
      model_info::cnst const* compartment = nullptr;

      if (!s->getHasOnlySubstanceUnits()) {
        compartment = ts->findcnst(s->getCompartment());
        if (compartment == nullptr) {
          last_err << "can't find compartment '" << s->getCompartment() << "'";
          return false;
        }
      }

      ts->add_var(new model_info::var(sid, ts, model_info::trsys::typeSpecies,
                                      compartment));

      parser::AST const* const value_et = new realnumAST(s->getInitialAmount());

      model_info::ic_s_i isi;
      isi.v = ts->get_vars().back();
      isi.value = new model_info::val(value_et);
      is.li.push_back(isi);
      isi.v->set_init(value_et);
    }
  }
  return true;
}

bool sbml_read_function_definitions(Model const* const model,
                                    model_info::trsys* const ts,
                                    std::stringstream& last_err) {
  ListOfFunctionDefinitions const* const functions =
      model->getListOfFunctionDefinitions();
  for (std::size_t i = 0; i < functions->size(); ++i) {
    FunctionDefinition const* const f = functions->get(i);

    const std::size_t narg = f->getNumArguments();
    std::vector<model_info::fnc_arg>* const args =
        new std::vector<model_info::fnc_arg>;

    for (std::size_t j = 0; j < narg; ++j) {
      model_info::fnc_arg a;
      a.name = f->getArgument(j)->getName();
      a.typ = model_info::trsys::typeReal;
      args->push_back(a);
    }

    parser::AST const* body;
    if (!sbml_read_expr(f->getBody(), ts, 0, body, last_err, args)) {
      return false;
    }

    ts->add_fnc(new model_info::fnc(f->getId(), ts, model_info::trsys::typeReal,
                                    args, body));
  }
  return true;
}

bool sbml_read_reactions(Model const* const model, model_info::trsys* const ts,
                         std::stringstream& last_err,
                         base const* const sl = nullptr) {
  ListOfSpecies const* const species = model->getListOfSpecies();
  ListOfReactions const* const reactions = model->getListOfReactions();
  for (std::size_t i = 0; i < reactions->size(); ++i) {
    Reaction const* const r = reactions->get(i);
    ListOfSpeciesReferences const* const s_in = r->getListOfReactants();
    ListOfSpeciesReferences const* const s_out = r->getListOfProducts();

    std::vector<model_info::chemreaction_item> reactants;

    for (std::size_t j = 0; j < s_in->size(); ++j) {
      SpeciesReference const* s_ =
          static_cast<SpeciesReference const*>(s_in->get(j));
      const std::string& sid = s_->getSpecies();
      Species const* const s = species->get(sid);

      if (!s->getBoundaryCondition()) {
        model_info::var const* const v = ts->findvar(sid);
        if (v == nullptr) {
          last_err << "unknown identifier '" << sid << "'";
          return false;
        }

        const double c = s_->getStoichiometry();

        auto const& crit = std::find_if(
            reactants.begin(), reactants.end(),
            [&](const model_info::chemreaction_item& ri) { return ri.v == v; });

        if (crit == reactants.end()) {
          const model_info::chemreaction_item cri = {v, c};
          reactants.push_back(cri);
        } else {
          crit->c += c;
        }
      }
    }

    std::vector<model_info::chemreaction_item> stoichiometry;
    stoichiometry.assign(reactants.begin(), reactants.end());
    for (auto& ri : stoichiometry) {
      ri.c = -ri.c;
    }

    std::vector<model_info::chemreaction_item> products;

    for (std::size_t j = 0; j < s_out->size(); ++j) {
      SpeciesReference const* const s_ =
          static_cast<SpeciesReference const*>(s_out->get(j));
      const std::string& sid = s_->getSpecies();

      {
        Species const* const s = species->get(sid);
        if (!s->getBoundaryCondition()) {
          model_info::var const* const v = ts->findvar(sid);
          if (v == nullptr) {
            last_err << "unknown identifier '" << sid << "'";
            return false;
          }

          const double c = s_->getStoichiometry();

          auto const& crit =
              std::find_if(products.begin(), products.end(),
                           [&](const model_info::chemreaction_item& ri) {
                return ri.v == v;
              });

          if (crit == products.end()) {
            const model_info::chemreaction_item cri = {v, c};
            products.push_back(cri);
          } else {
            crit->c += c;
          }

          auto const& crit2 =
              std::find_if(stoichiometry.begin(), stoichiometry.end(),
                           [&](const model_info::chemreaction_item& ri) {
                return ri.v == v;
              });

          if (crit2 == stoichiometry.end()) {
            const model_info::chemreaction_item cri = {v, c};
            stoichiometry.push_back(cri);
          } else {
            crit2->c += c;
          }
        }
      }
    }

    stoichiometry.erase(
        std::remove_if(stoichiometry.begin(), stoichiometry.end(),
                       [](const model_info::chemreaction_item& ri) {
          return std::fabs(ri.c) < std::numeric_limits<double>::epsilon();
        }),
        stoichiometry.end());

    std::vector<lparam_s> lparams;

    ListOfParameters const* const lparameters =
        r->getKineticLaw()->getListOfParameters();
    for (std::size_t j = 0; j < lparameters->size(); ++j) {
      Parameter const* const lp = lparameters->get(j);

      static char tmp[1024];
      sprintf(tmp, "r%lu_%s", i + 1, lp->getId().c_str());

      model_info::cnst const* const c = new model_info::cnst(
          tmp, ts, model_info::trsys::typeReal,
          new model_info::val(new realnumAST(lp->getValue())));

      ts->add_cnst(c);

      lparam_s p;
      p.c = c;
      p.name = lp->getId();
      lparams.push_back(p);
    }

    AST const* opt_guard = nullptr;

    AST const* guard;
    std::vector<model_info::transition_update_item> updates;

    if (!set_chemreaction_guard_updates(opt_guard, reactants, stoichiometry,
                                        guard, updates, last_err, sl)) {
      return false;
    }

    if (r->getReversible()) {
      parser::AST const* rate;
      if (!sbml_read_expr(r->getKineticLaw()->getMath(), ts, -1, rate, last_err,
                          nullptr, &lparams)) {
        return false;
      }

      ts->add_transition(new model_info::chemreaction(
          "", ts, guard, rate, updates, reactants, products, stoichiometry));

      parser::AST const* reverse_rate;
      if (!sbml_read_expr(r->getKineticLaw()->getMath(), ts, 1, reverse_rate,
                          last_err, nullptr, &lparams)) {
        return false;
      }

      for (auto& ri : stoichiometry) {
        ri.c = -ri.c;
      }

      AST const* reverse_guard;
      std::vector<model_info::transition_update_item> reverse_updates;

      if (!set_chemreaction_guard_updates(opt_guard, products, stoichiometry,
                                          reverse_guard, reverse_updates,
                                          last_err, sl)) {
        return false;
      }

      ts->add_transition(new model_info::chemreaction(
          "", ts, reverse_guard, reverse_rate, reverse_updates, products,
          reactants, stoichiometry, opt_guard));
    } else {
      parser::AST const* rate;
      if (!sbml_read_expr(r->getKineticLaw()->getMath(), ts, 0, rate, last_err,
                          nullptr, &lparams)) {
        return false;
      }

      ts->add_transition(new model_info::chemreaction(
          "", ts, guard, rate, updates, reactants, products, stoichiometry));
    }
  }
  return true;
}

bool sbml_read_rules(Model const* const model, model_info::trsys* const ts,
                     model_info::ic_s& is, std::stringstream& last_err) {
  ListOfInitialAssignments const* const initassignments =
      model->getListOfInitialAssignments();
  for (std::size_t i = 0; i < initassignments->size(); ++i) {
    InitialAssignment const* const ia = initassignments->get(i);

    model_info::var const* const v = ts->findvar(ia->getSymbol());
    if (v == nullptr) {
      last_err << "can't find variable '" << ia->getSymbol() << "'";
      return false;
    }

    parser::AST const* value_et;
    if (!sbml_read_expr(ia->getMath(), ts, 0, value_et, last_err)) {
      return false;
    }

    auto ili =
        std::find_if(is.li.begin(), is.li.end(),
                     [=](const model_info::ic_s_i& is) { return is.v == v; });
    if (ili != is.li.end()) {
      ili->value = new model_info::val(value_et);
    } else {
      model_info::ic_s_i isi;
      isi.v = v;
      isi.value = new model_info::val(value_et);
      is.li.push_back(isi);
    }

    v->set_init(value_et);
  }

  return true;
}

bool build_ASTNode_from_expr(
    ASTNode*& anode, AST const* const et,
    model_info::chemreaction const* const cr = nullptr) {
  if (et->is_number()) {
  } else if (et->is_cnst()) {
    anode = new ASTNode(AST_NAME);
    anode->setName(
        static_cast<cnstAST const*>(et)->get_cnst()->get_name().c_str());
  } else if (et->is_var()) {
    anode = new ASTNode(AST_NAME);
    anode->setName(
        static_cast<varAST const*>(et)->get_var()->get_name().c_str());
  } else if (et->is_stdfcall()) {
    stdfcallAST const* const etu = static_cast<stdfcallAST const*>(et);
    if (etu->get_stdf()->get_id() == STDF_MASS_ACTION) {
      assert(cr != nullptr);
    } else {
    }
  } else {
  }

  return true;
}

bool sbml_write_reactions(Model* const model,
                          model_info::trsys const* const ts) {
  for (auto t : ts->get_transitions()) {
    model_info::chemreaction const* const cr =
        static_cast<model_info::chemreaction const*>(t);
    if (!cr->is_reversible() || cr->get_reverse()) {
      Reaction* const reaction = model->createReaction();

      static unsigned int r_id;
      r_id++;

      if (cr->get_name() != "") {
        reaction->setId(cr->get_name());
      } else {
        reaction->setId(std::to_string(r_id));
      }

      reaction->setFast(false);
      reaction->setReversible(cr->is_reversible());

      const std::vector<model_info::chemreaction_item>& reactants =
          cr->get_reactants();
      const std::vector<model_info::chemreaction_item>& products =
          cr->get_products();

      for (auto const& re : reactants) {
        SpeciesReference* spr = reaction->createReactant();
        spr->setSpecies(re.v->get_name());
        spr->setConstant(false);
      }

      for (auto const& re : products) {
        SpeciesReference* spr = reaction->createProduct();
        spr->setSpecies(re.v->get_name());
        spr->setConstant(false);
      }

      KineticLaw* const kl = reaction->createKineticLaw();

      ASTNode* astRate = nullptr;
      if (!build_ASTNode_from_expr(astRate, cr->get_rate(), cr)) {
        return false;
      }

      if (cr->get_reverse()) {
        ASTNode* astReverseRate = nullptr;
        if (!build_ASTNode_from_expr(astReverseRate,
                                     cr->get_reverse()->get_rate(), cr)) {
          return false;
        }

        ASTNode* const astMinus = new ASTNode(AST_MINUS);
        astMinus->addChild(astRate);
        astMinus->addChild(astReverseRate);
        astRate = astMinus;
      }

      ASTNode* const astComp = new ASTNode(AST_NAME);
      char const* const compname = "comp";
      astComp->setName(compname);
      ASTNode* const astMath = new ASTNode(AST_TIMES);
      astMath->addChild(astComp);
      astMath->addChild(astRate);

      kl->setMath(astMath);
    }
  }

  for (auto& tss : ts->get_trsyss()) {
    if (!sbml_write_reactions(model, tss)) {
      return false;
    }
  }

  return true;
}

#endif

bool sbml_read(char const* const file_name, model_info::trsys*& ts,
               std::stringstream& last_err, base const* const sl) {
#if HAVE_LIBSBML
  SBMLDocument const* const document = readSBML(file_name);

  const unsigned int errors = document->getNumErrors();
  if (errors > 0) {
    return false;
  }

  ts = new model_info::trsys("", nullptr);

  Model const* const model = document->getModel();

  model_info::ic_s is;

  if (!sbml_read_compartments(model, ts, last_err) ||
      !sbml_read_parameters(model, ts, last_err) ||
      !sbml_read_species(model, ts, is, last_err) ||
      !sbml_read_function_definitions(model, ts, last_err) ||
      !sbml_read_reactions(model, ts, last_err, sl) ||
      !sbml_read_rules(model, ts, is, last_err)) {
    ts = nullptr;
    return false;
  }

  std::vector<model_info::ic_s> ls;
  is.p = 1.0;
  ls.push_back(is);

  ts->add_ic(new model_info::ic("", ts, ls));

  return true;
#else
  last_err << "please recompile using -Dsbml-support=ON";
  return false;
#endif
}

bool sbml_write(char const* const file_name, model_info::trsys const* const ts,
                std::stringstream& last_err) {
#if HAVE_LIBSBML
  if (ts->get_type() != model_info::TS_CHEMREACTIONS) {
    last_err << "model should contain chemical reactions only";
    return false;
  }

  SBMLDocument* const document = new SBMLDocument();

  Model* const model = document->createModel();
  model->setId("exported_model");

  Compartment* const comp = model->createCompartment();
  char const* const compname = "comp";
  comp->setId(compname);
  comp->setSize(1.0);
  comp->setConstant(true);

  for (auto const c : ts->get_cnsts()) {
    compute_value(c->get_value());

    Parameter* const p = model->createParameter();
    p->setId(c->get_name().c_str());
    p->setConstant(true);
    if (c->get_type()->is_boolean()) {
      p->setValue(static_cast<bool>(c->get_value()->get()));
    } else if (c->get_type()->is_integer()) {
      p->setValue(static_cast<int>(c->get_value()->get()));
    } else {
      p->setValue(static_cast<double>(c->get_value()->get()));
    }
  }

  for (auto const v : ts->get_vars()) {
    Species* const sp = model->createSpecies();
    sp->setCompartment(compname);
    sp->setId(v->get_name());

    double v0 = 0.0;

    if (ts->get_ics().size() == 1 &&
        ts->get_ics().front()->get_states().size() == 1) {
      auto const& vi =
          std::find_if(ts->get_ics().front()->get_states().front().li.begin(),
                       ts->get_ics().front()->get_states().front().li.end(),
                       [&](model_info::ic_s_i const& i) { return v == i.v; });
      if (vi != ts->get_ics().front()->get_states().front().li.end()) {
        compute_value(vi->value, true);
        v0 = static_cast<double>(vi->value->get());
      }
    }

    sp->setInitialAmount(v0);
    sp->setHasOnlySubstanceUnits(false);
    sp->setBoundaryCondition(false);
    sp->setConstant(false);
  }

  if (!sbml_write_reactions(model, ts)) {
    return false;
  }

  unsigned int numConsistencyErrors = 0;
  unsigned int numConsistencyWarnings = 0;
  unsigned int numValidationErrors = 0;
  unsigned int numValidationWarnings = 0;
  unsigned int numCheckFailures = document->checkInternalConsistency();
  std::string consistencyMessages;
  std::string validationMessages;

  if (numCheckFailures > 0) {
    for (unsigned int i = 0; i < numCheckFailures; i++) {
      const SBMLError* err = document->getError(i);
      if (err->isFatal() || err->isError()) {
        ++numConsistencyErrors;
      } else {
        ++numConsistencyWarnings;
      }
    }
    std::ostringstream oss;
    document->printErrors(oss);
    consistencyMessages = oss.str();
  }

  if (numConsistencyErrors) {
    consistencyMessages += "Further validation aborted.";
  } else {
    numCheckFailures = document->checkConsistency();
    if (numCheckFailures > 0) {
      for (unsigned int i = 0; i < numCheckFailures; i++) {
        const SBMLError* err = document->getError(i);
        if (err->isFatal() || err->isError()) {
          ++numValidationErrors;
        } else {
          ++numValidationWarnings;
        }
      }
      std::ostringstream oss;
      document->printErrors(oss);
      validationMessages = oss.str();
    }
  }

  if (numConsistencyErrors || numValidationErrors) {
    if (numConsistencyErrors) {
      last_err << "ERROR: encountered " << numConsistencyErrors
               << " consistency error" << (numConsistencyErrors == 1 ? "" : "s")
               << " in model '" << document->getModel()->getId() << "'."
               << std::endl;
    }
    if (numConsistencyWarnings) {
      last_err << "Notice: encountered " << numConsistencyWarnings
               << " consistency warning"
               << (numConsistencyWarnings == 1 ? "" : "s") << " in model '"
               << document->getModel()->getId() << "'." << std::endl;
    }
    last_err << std::endl << consistencyMessages;

    if (numValidationErrors) {
      last_err << "ERROR: encountered " << numValidationErrors
               << " validation error" << (numValidationErrors == 1 ? "" : "s")
               << " in model '" << document->getModel()->getId() << "'."
               << std::endl;
    }
    if (numValidationWarnings) {
      last_err << "Notice: encountered " << numValidationWarnings
               << " validation warning"
               << (numValidationWarnings == 1 ? "" : "s") << " in model '"
               << document->getModel()->getId() << "'." << std::endl;
    }
    last_err << std::endl << validationMessages;

    return false;
  }

  if (!writeSBMLToFile(document, file_name)) {
    last_err << "failed to write '" << file_name << "'";
    return false;
  }

  return true;
#else
  last_err << "please recompile using -Dno-sbml=OFF";
  return false;
#endif
}
}
}
