// Copyright (c) 2021 Johannes Delport
// This code is licensed under MIT license (see LICENSE for details)
//
// aether_sims D11: cuprate-memristor primitive. See Memristor.hpp for the
// model equations and provenance (Guenkel et al., Adv. Electron. Mater.
// 2026 e00861, doi 10.1002/aelm.202500861).

#include "JoSIM/Memristor.hpp"

#include <cmath>
#include <cstdio>
#include <utility>

#include "JoSIM/Constants.hpp"
#include "JoSIM/Errors.hpp"
#include "JoSIM/Misc.hpp"

using namespace JoSIM;

Memristor::Memristor(const std::pair<tokens_t, string_o>& s,
                     const NodeConfig& ncon, const nodemap& nm,
                     std::unordered_set<std::string>& lm, nodeconnections& nc,
                     Input& iObj, Spread& spread, int64_t& bi) {
  at_ = iObj.argAnal;
  // Check if the label has already been defined
  if (lm.count(s.first.at(0)) != 0) {
    Errors::invalid_component_errors(ComponentErrors::DUPLICATE_LABEL,
                                     s.first.at(0));
  }
  // Set the label
  netlistInfo.label_ = s.first.at(0);
  // Add the label to the known labels list
  lm.emplace(s.first.at(0));
  // Set the model for this memristor instance (4th token)
  set_model(s.first, iObj.netlist.models_new, s.second);
  // Companion resistance default: the calibration LRS resistance (max
  // device conductance) -- the one-step-lag fixed point then contracts
  // for any resistive surrounding network. Initial state/temperature.
  r0_ = model_.mrRmin();
  tLoc_ = model_.mrTcal();
  double spr = 1.0;
  for (int64_t i = 4; i < s.first.size(); ++i) {
    const auto& t = s.first.at(i);
    if (t.rfind("R0=", 0) == 0) {
      r0_ = parse_param(t.substr(3), iObj.parameters, s.second);
    } else if (t.rfind("G0=", 0) == 0) {
      g_ = parse_param(t.substr(3), iObj.parameters, s.second);
    } else if (t.rfind("TLOC=", 0) == 0) {
      tLoc_ = parse_param(t.substr(5), iObj.parameters, s.second);
    } else if (t.rfind("SPREAD=", 0) == 0) {
      spr = parse_param(t.substr(7), iObj.parameters, s.second);
    }
  }
  // SPREAD= jitters the calibration resistances (device-to-device spread)
  if (spr != 1.0) {
    model_.mrRmax(spread.spread_value(model_.mrRmax(), Spread::RES, spr));
    model_.mrRmin(spread.spread_value(model_.mrRmin(), Spread::RES, spr));
  }
  if (g_ < 0.0) g_ = 0.0;
  if (g_ > 1.0) g_ = 1.0;
  gPrevAccept_ = gPrev2_ = g_;
  netlistInfo.value_ = r0_;
  // Set the node configuration type
  indexInfo.nodeConfig_ = ncon;
  // Set current index and increment it
  indexInfo.currentIndex_ = bi++;
  // Set the node indices, using token 2 and 3
  set_node_indices(tokens_t(s.first.begin() + 1, s.first.begin() + 3), nm, nc);
  // Set the non zero, column index and row pointer vectors
  set_matrix_info();
  // Append the branch-row value (identical shape to a Resistor of R0)
  if (at_ == AnalysisType::Voltage) {
    matrixInfo.nonZeros_.emplace_back(-r0_);
  } else if (at_ == AnalysisType::Phase) {
    matrixInfo.nonZeros_.emplace_back(-((2.0 * iObj.transSim.tstep()) / 3.0) *
                                      (r0_ / Constants::SIGMA));
  }
  // Check the initial temperature against the calibration range
  set_temperature(tLoc_);
}

void Memristor::set_model(const tokens_t& t,
                          const vector_pair_t<Model, string_o>& models,
                          const string_o& subc) {
  bool found = false;
  for (auto& m : models) {
    if (m.first.modelName() == t.at(3) && m.first.mtype() == 1) {
      if ((m.second && subc && m.second.value() == subc.value()) ||
          (!m.second && !subc)) {
        model_ = m.first;
        found = true;
        break;
      }
    }
  }
  if (!found) {
    // Fall back to a global (non-subcircuit) memristor model
    for (auto& m : models) {
      if (m.first.modelName() == t.at(3) && m.first.mtype() == 1 &&
          !m.second) {
        model_ = m.first;
        found = true;
        break;
      }
    }
  }
  if (!found) {
    Errors::invalid_component_errors(ComponentErrors::MODEL_NOT_DEFINED,
                                     Misc::vector_to_string(t));
  }
}

double Memristor::rmax_t() const {
  return model_.mrRmax() * std::exp(-model_.mrAh() * (tLoc_ - model_.mrTcal()));
}

double Memristor::rmin_t() const {
  return model_.mrRmin() * std::exp(model_.mrAl() * (tLoc_ - model_.mrTcal()));
}

double Memristor::mmax_t() const {
  // SCLC shallow-trap law m = 1 + Tt/T, pinned to the calibration value:
  // m(T) = 1 + (mmax_cal - 1) * Tcal / T.
  return 1.0 + (model_.mrMmax() - 1.0) * model_.mrTcal() / tLoc_;
}

double Memristor::r_of_g(const double& g) const {
  const double rmax = rmax_t(), rmin = rmin_t();
  return (rmax * rmin) / (rmax * g + rmin * (1.0 - g));
}

double Memristor::m_of_g(const double& g) const {
  return model_.mrMmin() * g + mmax_t() * (1.0 - g);
}

double Memristor::i_model(const double& v, const double& g) const {
  const double av = std::abs(v);
  if (av == 0.0) return 0.0;
  const double i = std::pow(av, m_of_g(g)) / r_of_g(g);
  return v < 0.0 ? -i : i;
}

void Memristor::update_state(const double& v, const double& dt) {
  // Exact exponential update of the memory equation (paper Eqs (9)/(10));
  // unconditionally stable for any dt, exact for piecewise-constant v.
  const double s = model_.mrTscale();
  if (v < 0.0) {
    // SET / potentiation, driven by the SCLC current magnitude
    const double iscl = std::abs(i_model(v, g_));
    const double kp = s * model_.mrKp0() * std::pow(iscl, model_.mrAlpha());
    if (kp > 0.0) {
      const double e = std::exp(-kp * dt);
      g_ = 1.0 + (g_ - 1.0) * e;
    }
  } else if (v > 0.0) {
    // RESET / depression toward the voltage-dependent floor gmin(V)
    const double gmin =
        1.0 / (1.0 + std::exp((v - model_.mrVlog()) / model_.mrXi()));
    if (g_ > gmin) {
      const double kd = s * model_.mrKd0() * std::exp(model_.mrEtad() * v);
      const double e = std::exp(-kd * dt);
      g_ = gmin + (g_ - gmin) * e;
    }
  }
  if (g_ < 0.0) g_ = 0.0;
  if (g_ > 1.0) g_ = 1.0;
  // Refresh the Thevenin RHS correction for the next solve
  eN_ = v - r0_ * i_model(v, g_);
}

void Memristor::set_temperature(const double& t) {
  tLoc_ = t;
  if (tLoc_ < model_.mrTrmin() || tLoc_ > model_.mrTrmax()) {
    extrapolated_ = true;
    if (!extrapWarned_) {
      std::fprintf(stderr,
                   "WARNING: memristor %s at T_loc = %.6g K is outside the "
                   "model calibration range [%.6g, %.6g] K -- EXTRAPOLATED "
                   "(D11 flag)\n",
                   netlistInfo.label_.c_str(), tLoc_, model_.mrTrmin(),
                   model_.mrTrmax());
      extrapWarned_ = true;
    }
  }
}

// Update timestep based on a scalar factor (mirrors Resistor)
void Memristor::update_timestep(const double& factor) {
  if (at_ == AnalysisType::Phase) {
    matrixInfo.nonZeros_.back() = factor * factor * matrixInfo.nonZeros_.back();
  }
}
