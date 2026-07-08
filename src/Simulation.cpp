// Copyright (c) 2021 Johannes Delport
// This code is licensed under MIT license (see LICENSE for details)

#include "JoSIM/Simulation.hpp"

#include <cmath>
#include <iostream>

#include "JoSIM/AnalysisType.hpp"
#include "JoSIM/Components.hpp"
#include "JoSIM/Constants.hpp"
#include "JoSIM/Matrix.hpp"
#include "JoSIM/Model.hpp"
#include "JoSIM/ProgressBar.hpp"

using namespace JoSIM;

Simulation::Simulation(Input &iObj, Matrix &mObj) {
  // Batch: setup + factor + run the whole transient + free, retrying with a
  // smaller step if requested. Behaviour-identical to the original loop; the
  // body is now factored into prepare()/trans_sim()/finish() so the stepped
  // path can reuse it (aether_sims D4).
  while (needsTR_) {
    prepare(iObj, mObj);
    trans_sim(mObj);
    if (needsTR_) {
      reduce_step(iObj, mObj);
    }
    finish();
  }
}

Simulation::Simulation(Input &iObj, Matrix &mObj, bool deferRun) {
  // Stepped / co-simulation entry (aether_sims D4): do setup + factor and the
  // pre-t=0 stabilization, then leave the solver factored so the caller can
  // drive step(i, mObj) in a loop and finish() (or let the destructor free).
  // No reduce_step retry here -- keep dt small enough (the global restart is a
  // P3 concern).
  (void)deferRun;
  prepare(iObj, mObj);
  run_startup(mObj);
}

Simulation::~Simulation() { finish(); }

void Simulation::prepare(Input &iObj, Matrix &mObj) {
  // Generic setup for the given step size
  setup(iObj, mObj);
  // Solver setup
#ifdef SLU
  lu.create_matrix(mObj.rp.size() - 1, mObj.nz, mObj.ci, mObj.rp);
  lu.factorize();
#else
  simOK_ = klu_l_defaults(&Common_);
  assert(simOK_);
  Symbolic_ = klu_l_analyze(mObj.rp.size() - 1, &mObj.rp.front(),
                            &mObj.ci.front(), &Common_);
  Numeric_ = klu_l_factor(&mObj.rp.front(), &mObj.ci.front(), &mObj.nz.front(),
                          Symbolic_, &Common_);
#endif
  kluReady_ = true;
  // b_ and the time axis are initialized here (moved out of trans_sim) so
  // step() works whether called from trans_sim or from a stepped driver.
  results.timeAxis.clear();
  b_.resize(mObj.rp.size(), 0.0);
}

void Simulation::finish() {
  if (!kluReady_) return;
#ifdef SLU
  lu.free();
#else
  klu_l_free_symbolic(&Symbolic_, &Common_);
  klu_l_free_numeric(&Numeric_, &Common_);
#endif
  kluReady_ = false;
}

namespace {
// Find a junction device by its JoSIM label, or nullptr.
JJ *find_jj(Matrix &mObj, const std::string &label) {
  for (const auto &j : mObj.components.junctionIndices) {
    JJ &jj = std::get<JJ>(mObj.components.devices.at(j));
    if (jj.netlistInfo.label_ == label) return &jj;
  }
  return nullptr;
}
}  // namespace

void Simulation::set_jj_temperature(Matrix &mObj, const std::string &label,
                                    double T) {
  if (JJ *jj = find_jj(mObj, label)) jj->update_temperature(T);
}

void Simulation::set_all_temperatures(Matrix &mObj, double T) {
  for (const auto &j : mObj.components.junctionIndices) {
    std::get<JJ>(mObj.components.devices.at(j)).update_temperature(T);
  }
}

void Simulation::set_jj_control_current(Matrix &mObj, const std::string &label,
                                        double ictrl) {
  // aether_sims D8: drive a junction's Ic(Ictrl) law directly (the static
  // per-junction trim of OQ-8 / the Python-hook co-sim path). For junctions
  // with a CTRL= binding the engine overwrites this after every solve.
  if (JJ *jj = find_jj(mObj, label)) jj->update_control_current(ictrl);
}

double Simulation::jj_ic(Matrix &mObj, const std::string &label) {
  JJ *jj = find_jj(mObj, label);
  return jj ? jj->model_.ic() : 0.0;
}

double Simulation::jj_phase(Matrix &mObj, const std::string &label) {
  JJ *jj = find_jj(mObj, label);
  if (!jj) return 0.0;
  if (atyp_ == AnalysisType::Voltage) return x_.at(jj->variableIndex_);
  // Phase mode: the phase is the node difference.
  double pp = jj->indexInfo.posIndex_ ? x_.at(jj->indexInfo.posIndex_.value()) : 0.0;
  double pn = jj->indexInfo.negIndex_ ? x_.at(jj->indexInfo.negIndex_.value()) : 0.0;
  return pp - pn;
}

double Simulation::jj_voltage(Matrix &mObj, const std::string &label) {
  JJ *jj = find_jj(mObj, label);
  if (!jj) return 0.0;
  if (atyp_ == AnalysisType::Phase) return x_.at(jj->variableIndex_);
  // Voltage mode: the voltage is the node difference.
  double vp = jj->indexInfo.posIndex_ ? x_.at(jj->indexInfo.posIndex_.value()) : 0.0;
  double vn = jj->indexInfo.negIndex_ ? x_.at(jj->indexInfo.negIndex_.value()) : 0.0;
  return vp - vn;
}

double Simulation::jj_current(Matrix &mObj, const std::string &label) {
  JJ *jj = find_jj(mObj, label);
  if (!jj || !jj->indexInfo.currentIndex_) return 0.0;
  return x_.at(jj->indexInfo.currentIndex_.value());
}

double Simulation::jj_power(Matrix &mObj, const std::string &label) {
  JJ *jj = find_jj(mObj, label);
  return jj ? jj->power() : 0.0;
}

void Simulation::setup(Input &iObj, Matrix &mObj) {
  // Simulation setup
  simSize_ = iObj.transSim.simsize();
  atyp_ = iObj.argAnal;
  minOut_ = iObj.argMin;
  needsLU_ = false;
  needsTR_ = false;
  stepSize_ = iObj.transSim.tstep();
  prstep_ = iObj.transSim.prstep();
  prstart_ = iObj.transSim.prstart();
  startup_ = iObj.transSim.startup();
  x_.clear();
  x_.resize(mObj.branchIndex, 0.0);
  if (!mObj.relevantTraces.empty()) {
    results.xVector.resize(mObj.branchIndex);
    for (const auto &i : mObj.relevantIndices) {
      results.xVector.at(i).emplace();
    }
  } else {
    results.xVector.resize(mObj.branchIndex, std::vector<double>(0));
  }
}

bool Simulation::solve_only(int64_t i, Matrix &mObj) {
  // Setup the b matrix, then solve Ax=b into x_. No result storage.
  setup_b(mObj, i, i * stepSize_);
  if (needsTR_) return true;
  x_ = b_;
#ifdef SLU
  lu.solve(x_);
#else
  simOK_ =
      klu_l_tsolve(Symbolic_, Numeric_, mObj.rp.size() - 1, 1, &x_.front(),
                   &Common_);
  if (!simOK_) Errors::simulation_errors(SimulationErrors::MATRIX_SINGULAR);
#endif
  // aether_sims D8: feed each CTRL=-coupled junction the control branch
  // current from this solve; its Ic(Ictrl) applies from the next step.
  // Memoryless one-step-lag coupling -- the engine-native equivalent of the
  // spec's post_step (read Ictrl) -> pre_step (push Ic) hook pair, and it runs
  // on every path (batch, run_main, stepped and startup).
  for (const auto& cp : mObj.components.ctrlCouplings) {
    std::get<JJ>(mObj.components.devices.at(cp.first))
        .update_control_current(x_.at(cp.second));
  }
  return false;
}

bool Simulation::step(int64_t i, Matrix &mObj) {
  // One main-loop transient step: solve, then store the requested traces.
  if (solve_only(i, mObj)) return true;
  // Store results (only requested, to prevent massive memory usage)
  for (int64_t j = 0; j < results.xVector.size(); ++j) {
    if (results.xVector.at(j)) {
      results.xVector.at(j).value().emplace_back(x_.at(j));
    }
  }
  // Store the time step
  results.timeAxis.emplace_back(i * stepSize_);
  return false;
}

void Simulation::run_startup(Matrix &mObj) {
  if (!startup_) return;
  // Stabilize the simulation before starting at t=0
  int64_t startup = 2 * pow(10, (abs(log10(stepSize_)) - 12) * 2 + 1);
  if (startup > 1000) startup = 1000;
  for (int64_t i = -startup; i < 0; ++i) {
    if (solve_only(i, mObj)) return;
  }
}

void Simulation::trans_sim(Matrix &mObj) {
  ProgressBar bar;
  if (!minOut_) {
    bar.create_thread();
    bar.set_bar_width(30);
    bar.fill_bar_progress_with("O");
    bar.fill_bar_remainder_with(" ");
    bar.set_status_text("Simulating");
    bar.set_total((float)simSize_);
  }
  run_startup(mObj);
  if (needsTR_) return;
  // Start the simulation loop
  for (int64_t i = 0; i < simSize_; ++i) {
    if (!minOut_) {
      bar.update(static_cast<float>(i));
    }
    double t = i * stepSize_;
    if (pre_step_hook_) pre_step_hook_(i, t);
    if (step(i, mObj)) return;
    if (post_step_hook_) post_step_hook_(i, t);
  }
  if (!minOut_) {
    bar.complete();
    std::cout << "\n";
  }
}

void Simulation::reduce_and_restart(Input &iObj, Matrix &mObj) {
  // Free the old factorization, halve the timestep + rebuild (reduce_step),
  // then re-prepare and re-run the startup. The caller restarts step() at i=0.
  finish();
  reduce_step(iObj, mObj);
  prepare(iObj, mObj);
  run_startup(mObj);
}

void Simulation::run_main(Matrix &mObj) {
  // Main loop with hooks, no progress bar -- for the deferred (stepped)
  // C++-callback co-sim path. Startup is assumed already run by the deferred
  // constructor.
  for (int64_t i = 0; i < simSize_; ++i) {
    double t = i * stepSize_;
    if (pre_step_hook_) pre_step_hook_(i, t);
    if (step(i, mObj)) return;
    if (post_step_hook_) post_step_hook_(i, t);
  }
}

void Simulation::reduce_step(Input &iObj, Matrix &mObj) {
  iObj.transSim.tstep(iObj.transSim.tstep() / 2);
  bool tempMinOut = iObj.argMin;
  if (!iObj.argMin) iObj.argMin = true;
  Matrix newmObj;
  mObj = newmObj;
  // Create the matrix in csr format
  for (auto &i : mObj.sourcegen) {
    i.clearMisc();
  }
  mObj.create_matrix(iObj);
  find_relevant_traces(iObj, mObj);
  //// Dump expanded Netlist since it is no longer needed
  // iObj.netlist.expNetlist.clear();
  // iObj.netlist.expNetlist.shrink_to_fit();
  if (!tempMinOut) iObj.argMin = tempMinOut;
  results.xVector.clear();
  results.timeAxis.clear();
}

void Simulation::setup_b(Matrix &mObj, int64_t i, double step, double factor) {
  // Clear b matrix and reset
  b_.clear();
  b_.resize(mObj.rp.size(), 0.0);
  // Handle jj
  handle_jj(mObj, i, step, factor);
  if (needsTR_) return;
  // Re-factorize the LU if any jj transitions
  if (needsLU_) {
    mObj.create_nz();
#ifdef SLU
      lu.factorize(true);
#else
      klu_l_free_numeric(&Numeric_, &Common_);
      Numeric_ = klu_l_factor(&mObj.rp.front(), &mObj.ci.front(),
                              &mObj.nz.front(), Symbolic_, &Common_);
#endif
    needsLU_ = false;
  }
  // Handle current sources
  handle_cs(mObj, step, i);
  // Handle resistors
  handle_resistors(mObj, step);
  // Handle inductors
  handle_inductors(mObj, factor);
  // Handle capacitors
  handle_capacitors(mObj);
  // Handle voltage sources
  handle_vs(mObj, i, step, factor);
  // Handle phase sources
  handle_ps(mObj, i, step, factor);
  // Handle ccvs
  handle_ccvs(mObj);
  // Handle vccs
  handle_vccs(mObj);
  // Handle transmission lines
  handle_tx(mObj, i, step);
}

void Simulation::handle_cs(Matrix &mObj, double &step, const int64_t &i) {
  for (const auto &j : mObj.components.currentsources) {
    if (j.indexInfo.nodeConfig_ == NodeConfig::POSGND) {
      b_.at(j.indexInfo.posIndex_.value()) -=
          (mObj.sourcegen.at(j.sourceIndex_).value(step));
    } else if (j.indexInfo.nodeConfig_ == NodeConfig::GNDNEG) {
      b_.at(j.indexInfo.negIndex_.value()) +=
          (mObj.sourcegen.at(j.sourceIndex_).value(step));
    } else {
      b_.at(j.indexInfo.posIndex_.value()) -=
          (mObj.sourcegen.at(j.sourceIndex_).value(step));
      b_.at(j.indexInfo.negIndex_.value()) +=
          (mObj.sourcegen.at(j.sourceIndex_).value(step));
    }
  }
}

void Simulation::handle_resistors(Matrix &mObj, double &step) {
  for (const auto &j : mObj.components.resistorIndices) {
    auto &temp = std::get<Resistor>(mObj.components.devices.at(j));
    NodeConfig &nc = temp.indexInfo.nodeConfig_;
    if (nc == NodeConfig::POSGND) {
      if (temp.thermalNoise) {
        b_.at(temp.indexInfo.posIndex_.value()) -=
            temp.thermalNoise.value().value(step);
      }
      temp.pn1_ = (x_.at(temp.indexInfo.posIndex_.value()));
    } else if (nc == NodeConfig::GNDNEG) {
      if (temp.thermalNoise) {
        b_.at(temp.indexInfo.negIndex_.value()) +=
            temp.thermalNoise.value().value(step);
      }
      temp.pn1_ = (-x_.at(temp.indexInfo.negIndex_.value()));
    } else {
      if (temp.thermalNoise) {
        b_.at(temp.indexInfo.posIndex_.value()) -=
            temp.thermalNoise.value().value(step);
        b_.at(temp.indexInfo.negIndex_.value()) +=
            temp.thermalNoise.value().value(step);
      }
      temp.pn1_ = (x_.at(temp.indexInfo.posIndex_.value()) -
                   x_.at(temp.indexInfo.negIndex_.value()));
    }
    if (atyp_ == AnalysisType::Phase) {
      // 4/3 φp1 - 1/3 φp2
      b_.at(temp.indexInfo.currentIndex_.value()) =
          (4.0 / 3.0) * temp.pn1_ - (1.0 / 3.0) * temp.pn2_;
      temp.pn4_ = temp.pn3_;
      temp.pn3_ = temp.pn2_;
      temp.pn2_ = temp.pn1_;
    }
  }
}

void Simulation::handle_inductors(Matrix &mObj, double factor) {
  for (const auto &j : mObj.components.inductorIndices) {
    auto &temp = std::get<Inductor>(mObj.components.devices.at(j));
    if (atyp_ == AnalysisType::Voltage) {
      // -2L/h Ip + L/2h Ip2
      b_.at(temp.indexInfo.currentIndex_.value()) =
          -(2.0 * temp.netlistInfo.value_ / (stepSize_ * factor)) *
              x_.at(temp.indexInfo.currentIndex_.value()) +
          ((temp.netlistInfo.value_ / (2.0 * (stepSize_ * factor))) *
           temp.In2_);
      // -2M/h Im + M/2h Im2
      for (const auto &m : temp.get_mutualInductance()) {
        Inductor &mi = std::get<Inductor>(mObj.components.devices.at(m.first));
        b_.at(temp.indexInfo.currentIndex_.value()) +=
            (-((2 * m.second) / (stepSize_ * factor)) *
                 x_.at(mi.indexInfo.currentIndex_.value()) +
             (m.second / (2.0 * (stepSize_ * factor))) * mi.In2_);
      }
      temp.In4_ = temp.In3_;
      temp.In3_ = temp.In2_;
      temp.In2_ = x_.at(temp.indexInfo.currentIndex_.value());
    }
  }
}

void Simulation::handle_capacitors(Matrix &mObj) {
  for (const auto &j : mObj.components.capacitorIndices) {
    auto &temp = std::get<Capacitor>(mObj.components.devices.at(j));
    if (temp.indexInfo.posIndex_ && !temp.indexInfo.negIndex_) {
      temp.pn1_ = (x_.at(temp.indexInfo.posIndex_.value()));
    } else if (!temp.indexInfo.posIndex_ && temp.indexInfo.negIndex_) {
      temp.pn1_ = (-x_.at(temp.indexInfo.negIndex_.value()));
    } else {
      temp.pn1_ = (x_.at(temp.indexInfo.posIndex_.value()) -
                   x_.at(temp.indexInfo.negIndex_.value()));
    }
    if (atyp_ == AnalysisType::Voltage) {
      // 4/3 Vp1 - 1/3 Vp2
      b_.at(temp.indexInfo.currentIndex_.value()) =
          (4.0 / 3.0) * temp.pn1_ - (1.0 / 3.0) * temp.pn2_;
    } else if (atyp_ == AnalysisType::Phase) {
      // (8/3)φn-1 - (22/9)φn-2 + (8/9)φn-3 - (1/9)φn-4
      b_.at(temp.indexInfo.currentIndex_.value()) =
          (8.0 / 3.0) * temp.pn1_ - (22.0 / 9.0) * temp.pn2_ +
          (8.0 / 9.0) * temp.pn3_ - (1.0 / 9.0) * temp.pn4_;
      temp.pn7_ = temp.pn6_;
      temp.pn6_ = temp.pn5_;
      temp.pn5_ = temp.pn4_;
      temp.pn4_ = temp.pn3_;
      temp.pn3_ = temp.pn2_;
    }
    temp.pn2_ = temp.pn1_;
  }
}

void Simulation::handle_jj(Matrix &mObj, int64_t &i, double &step,
                           double factor) {
  for (const auto &j : mObj.components.junctionIndices) {
    auto &temp = std::get<JJ>(mObj.components.devices.at(j));
    const auto &model = temp.model_;
    if (temp.indexInfo.posIndex_ && !temp.indexInfo.negIndex_) {
      if (temp.thermalNoise) {
        b_.at(temp.indexInfo.posIndex_.value()) -=
            temp.thermalNoise.value().value(step);
      }
      temp.pn1_ = (x_.at(temp.indexInfo.posIndex_.value()));
    } else if (!temp.indexInfo.posIndex_ && temp.indexInfo.negIndex_) {
      if (temp.thermalNoise) {
        b_.at(temp.indexInfo.negIndex_.value()) +=
            temp.thermalNoise.value().value(step);
      }
      temp.pn1_ = (-x_.at(temp.indexInfo.negIndex_.value()));
    } else {
      if (temp.thermalNoise) {
        b_.at(temp.indexInfo.posIndex_.value()) -=
            temp.thermalNoise.value().value(step);
        b_.at(temp.indexInfo.negIndex_.value()) +=
            temp.thermalNoise.value().value(step);
      }
      temp.pn1_ = (x_.at(temp.indexInfo.posIndex_.value()) -
                   x_.at(temp.indexInfo.negIndex_.value()));
    }
    if (i > 0) {
      if (atyp_ == AnalysisType::Voltage) {
        temp.vn1_ = temp.pn1_;
        temp.pn1_ = x_.at(temp.variableIndex_);
      } else if (atyp_ == AnalysisType::Phase) {
        temp.vn1_ = x_.at(temp.variableIndex_);
        temp.pn1_ = temp.pn1_;
      }
    }
    // Guess voltage (V0)
    double v0 =
        (5.0 / 2.0) * temp.vn1_ - 2.0 * temp.vn2_ + (1.0 / 2.0) * temp.vn3_;
    // Phase guess (P0)
    temp.phi0_ = (4.0 / 3.0) * temp.pn1_ - (1.0 / 3.0) * temp.pn2_ +
                 ((1.0 / Constants::SIGMA) * ((2.0 * (stepSize_)) / 3.0)) * v0;
    // Ensure timestep is not too large
    if ((double)i / (double)simSize_ > 0.01) {
      if (abs(temp.phi0_ - temp.pn1_) > (0.20 * 2 * Constants::PI)) {
        needsTR_ = true;
        return;
        /*Errors::simulation_errors(
          SimulationErrors::PHASEGUESS_TOO_LARGE, temp.netlistInfo.label_);*/
      }
    }
    // (hbar / 2 * e) ( -(2 / h) φp1 + (1 / 2h) φp2 )
    if (atyp_ == AnalysisType::Voltage) {
      b_.at(temp.variableIndex_) =
          (Constants::SIGMA) * (-(2.0 / (stepSize_)) * temp.pn1_ +
                                (1.0 / (2.0 * (stepSize_))) * temp.pn2_);
      // (4 / 3) φp1 - (1/3) φp2
    } else if (atyp_ == AnalysisType::Phase) {
      b_.at(temp.variableIndex_) =
          (4.0 / 3.0) * temp.pn1_ - (1.0 / 3.0) * temp.pn2_;
    }
    temp.pn4_ = temp.pn3_;
    temp.pn3_ = temp.pn2_;
    temp.pn2_ = temp.pn1_;
    temp.vn6_ = temp.vn5_;
    temp.vn5_ = temp.vn4_;
    temp.vn4_ = temp.vn3_;
    temp.vn3_ = temp.vn2_;
    // Update junction transition
    if (model.rtype() == 1) {
      auto testLU = temp.update_value(v0);
      if (testLU && !needsLU_) {
        needsLU_ = true;
      }
    }
    // Ic * sin (phi * (φ0 - φ))
    double ic_sin_phi = 0.0;
    for (int harm = 0; harm < temp.model_.cpr().size(); ++harm) {
      ic_sin_phi += temp.model_.ic() *
                    (temp.model_.cpr().at(harm) *
                     sin((harm + 1) * (temp.phi0_ - temp.model_.phiOff())));
    }
    // aether_sims OQ-2: the plain Ic*sin branch is used for non-temperature
    // models AND for the YBCO weak link (ictemp == 1), whose supercurrent is
    // Ic(T)*sin(phi) with Ic(T) from update_temperature(). Only the BCS tDep
    // model uses the Ambegaokar-Baratoff / Haberkorn branch below.
    if (!temp.model_.tDep() || temp.model_.ictemp() == 1) {
      // -(hR / h + 2RC) * (Ic sin (φ0) - 2C / h Vp1 + C/2h Vp2 + It)
      b_.at(temp.indexInfo.currentIndex_.value()) =
          (temp.matrixInfo.nonZeros_.back()) *
          (ic_sin_phi - (((2 * model.c()) / (stepSize_)) * temp.vn1_) +
           ((model.c() / (2.0 * (stepSize_))) * temp.vn2_) + temp.it_);
    } else {
      double sin2_half_phi = 0.0;
      for (int harm = 0; harm < temp.model_.cpr().size(); ++harm) {
        sin2_half_phi +=
            temp.model_.cpr().at(harm) *
            sin((harm + 1) * (temp.phi0_ - temp.model_.phiOff()) / 2);
      }
      sin2_half_phi = sin2_half_phi * sin2_half_phi;
      double sin_phi = 0.0;
      for (int harm = 0; harm < temp.model_.cpr().size(); ++harm) {
        sin_phi += temp.model_.cpr().at(harm) *
                   sin((harm + 1) * (temp.phi0_ - temp.model_.phiOff()));
      }
      double sqrt_part = sqrt(1 - model.d() * sin2_half_phi);
      b_.at(temp.indexInfo.currentIndex_.value()) =
          // -(hR / h + 2RC) *(
          (temp.matrixInfo.nonZeros_.back()) *
          ((
               // (π * Δ / 2 * e * Rn) * g(Ictrl)  [ctrlScale_: aether_sims D8
               // -- this branch's supercurrent amplitude comes from del_/Rn,
               // not model_.ic(), so the control trim applies here directly]
               temp.ctrlScale_ *
               ((Constants::PI * temp.del_) /
                (2 * Constants::EV * temp.model_.rn()))
               // * (sin(φ0 - φ) / √(1 - D * sin²((φ0 - φ) / 2))
               * (sin_phi / sqrt_part)
               // * tanh(Δ / (2 * kB * T) * √(1 - D * sin²((φ0 - φ) / 2)))
               * tanh(temp.del_ / (2 * Constants::BOLTZMANN * model.t()) *
                      sqrt_part))
           // - 2C / h Vp1
           - (((2 * model.c()) / stepSize_) * temp.vn1_)
           // + C/2h Vp2
           + ((model.c() / (2.0 * stepSize_)) * temp.vn2_)
           // + It)
           + temp.it_);
    }
    // aether_sims D3: resistive dissipation from the junction voltage.
    temp.update_power(temp.vn1_);
    temp.vn2_ = temp.vn1_;
  }
}

void Simulation::handle_vs(Matrix &mObj, const int64_t &i, double &step,
                           double factor) {
  for (const auto &j : mObj.components.vsIndices) {
    auto &temp = std::get<VoltageSource>(mObj.components.devices.at(j));
    JoSIM::NodeConfig &nc = temp.indexInfo.nodeConfig_;
    if (atyp_ == AnalysisType::Voltage) {
      // Vn
      b_.at(temp.indexInfo.currentIndex_.value()) =
          (mObj.sourcegen.at(temp.sourceIndex_).value(step));
    } else if (atyp_ == AnalysisType::Phase) {
      if (nc == NodeConfig::POSGND) {
        temp.pn1_ = (x_.at(temp.indexInfo.posIndex_.value()));
      } else if (nc == NodeConfig::GNDNEG) {
        temp.pn1_ = (-x_.at(temp.indexInfo.negIndex_.value()));
      } else {
        temp.pn1_ = (x_.at(temp.indexInfo.posIndex_.value()) -
                     x_.at(temp.indexInfo.negIndex_.value()));
      }
      // (2e/hbar)(2h/3)Vn + (4/3)φn-1 - (1/3)φn-2
      b_.at(temp.indexInfo.currentIndex_.value()) =
          ((2 * stepSize_ * factor) / (3 * Constants::SIGMA)) *
              (mObj.sourcegen.at(temp.sourceIndex_).value(step)) +
          (4.0 / 3.0) * temp.pn1_ - (1.0 / 3.0) * temp.pn2_;
      temp.pn4_ = temp.pn3_;
      temp.pn3_ = temp.pn2_;
      temp.pn2_ = temp.pn1_;
    }
  }
}

void Simulation::handle_ps(Matrix &mObj, const int64_t &i, double &step,
                           double factor) {
  for (const auto &j : mObj.components.psIndices) {
    auto &temp = std::get<PhaseSource>(mObj.components.devices.at(j));
    if (atyp_ == AnalysisType::Phase) {
      // φn
      b_.at(temp.indexInfo.currentIndex_.value()) =
          (mObj.sourcegen.at(temp.sourceIndex_).value(step));
    } else if (atyp_ == AnalysisType::Voltage) {
      if (i == 0) {
        b_.at(temp.indexInfo.currentIndex_.value()) =
            (Constants::SIGMA / (stepSize_ * factor)) *
            ((3.0 / 2.0) * (mObj.sourcegen.at(temp.sourceIndex_).value(step)));
      } else if (i == 1) {
        b_.at(temp.indexInfo.currentIndex_.value()) =
            (Constants::SIGMA / (stepSize_ * factor)) *
            ((3.0 / 2.0) * (mObj.sourcegen.at(temp.sourceIndex_).value(step)) -
             2.0 * (mObj.sourcegen.at(temp.sourceIndex_)
                        .value(step - (stepSize_ * factor))));
      } else {
        b_.at(temp.indexInfo.currentIndex_.value()) =
            (Constants::SIGMA / (stepSize_ * factor)) *
            ((3.0 / 2.0) * (mObj.sourcegen.at(temp.sourceIndex_).value(step)) -
             2.0 * (mObj.sourcegen.at(temp.sourceIndex_)
                        .value(step - (stepSize_ * factor))) +
             0.5 * (mObj.sourcegen.at(temp.sourceIndex_)
                        .value(step - (2 * stepSize_ * factor))));
      }
    }
  }
}

void Simulation::handle_ccvs(Matrix &mObj) {
  for (const auto &j : mObj.components.ccvsIndices) {
    auto &temp = std::get<CCVS>(mObj.components.devices.at(j));
    if (atyp_ == AnalysisType::Phase) {
      if (temp.indexInfo.posIndex_ && !temp.indexInfo.negIndex_) {
        temp.pn1_ = (x_.at(temp.indexInfo.posIndex_.value()));
      } else if (!temp.indexInfo.posIndex_ && temp.indexInfo.negIndex_) {
        temp.pn1_ = (-x_.at(temp.indexInfo.negIndex_.value()));
      } else {
        temp.pn1_ = (x_.at(temp.indexInfo.posIndex_.value()) -
                     x_.at(temp.indexInfo.negIndex_.value()));
      }
      b_.at(temp.indexInfo.currentIndex_.value()) =
          (4.0 / 3.0) * temp.pn1_ - (1.0 / 3.0) * temp.pn2_;
      temp.pn4_ = temp.pn3_;
      temp.pn3_ = temp.pn2_;
      temp.pn2_ = temp.pn1_;
    }
  }
}

void Simulation::handle_vccs(Matrix &mObj) {
  for (const auto &j : mObj.components.vccsIndices) {
    auto &temp = std::get<VCCS>(mObj.components.devices.at(j));
    if (atyp_ == AnalysisType::Phase) {
      if (temp.posIndex2_ && !temp.negIndex2_) {
        temp.pn1_ = (x_.at(temp.posIndex2_.value()));
      } else if (!temp.posIndex2_ && temp.negIndex2_) {
        temp.pn1_ = (-x_.at(temp.negIndex2_.value()));
      } else {
        temp.pn1_ =
            (x_.at(temp.posIndex2_.value()) - x_.at(temp.negIndex2_.value()));
      }
      b_.at(temp.indexInfo.currentIndex_.value()) =
          (4.0 / 3.0) * temp.pn1_ - (1.0 / 3.0) * temp.pn2_;
      temp.pn4_ = temp.pn3_;
      temp.pn3_ = temp.pn2_;
      temp.pn2_ = temp.pn1_;
    }
  }
}

void Simulation::handle_tx(Matrix &mObj, const int64_t &i, double &step,
                           double factor) {
  for (const auto &j : mObj.components.txIndices) {
    auto &temp = std::get<TransmissionLine>(mObj.components.devices.at(j));
    // Z0
    double &Z = temp.netlistInfo.value_;
    // Td == k
    int64_t &k = temp.timestepDelay_;
    // Shorthands
    JoSIM::NodeConfig &nc = temp.indexInfo.nodeConfig_,
                      &nc2 = temp.nodeConfig2_;
    int_o &posInd = temp.indexInfo.posIndex_,
          &negInd = temp.indexInfo.negIndex_, &posInd2 = temp.posIndex2_,
          &negInd2 = temp.negIndex2_;
    int64_t &curInd = temp.indexInfo.currentIndex_.value(),
            &curInd2 = temp.currentIndex2_;
    if (atyp_ == AnalysisType::Voltage) {
      if (i >= k) {
        // φ1n-k
        if (nc == NodeConfig::POSGND) {
          temp.nk_1_ = results.xVector.at(posInd.value()).value().at(i - k);
        } else if (nc == NodeConfig::GNDNEG) {
          temp.nk_1_ = -results.xVector.at(negInd.value()).value().at(i - k);
        } else {
          temp.nk_1_ = results.xVector.at(posInd.value()).value().at(i - k) -
                       results.xVector.at(negInd.value()).value().at(i - k);
        }
        // φ2n-k
        if (nc2 == NodeConfig::POSGND) {
          temp.nk_2_ = results.xVector.at(posInd2.value()).value().at(i - k);
        } else if (nc2 == NodeConfig::GNDNEG) {
          temp.nk_2_ = -results.xVector.at(negInd2.value()).value().at(i - k);
        } else {
          temp.nk_2_ = results.xVector.at(posInd2.value()).value().at(i - k) -
                       results.xVector.at(negInd2.value()).value().at(i - k);
        }
        // I1n-k
        double &I1nk = results.xVector.at(curInd).value().at(i - k);
        // I2n-k
        double &I2nk = results.xVector.at(curInd2).value().at(i - k);
        // I1 = ZI2n-k + V2n-k
        b_.at(curInd) = Z * I2nk + temp.nk_2_;
        // I2 = ZI1n-k + V1n-k
        b_.at(curInd2) = Z * I1nk + temp.nk_1_;
      }
    } else if (atyp_ == AnalysisType::Phase) {
      if (i > 0) {
        temp.n2_1_ = temp.n1_1_;
        temp.n2_2_ = temp.n1_2_;
        // φ1n-1
        if (nc == NodeConfig::POSGND) {
          temp.n1_1_ = (x_.at(posInd.value()));
        } else if (nc == NodeConfig::GNDNEG) {
          temp.n1_1_ = (-x_.at(negInd.value()));
        } else {
          temp.n1_1_ = (x_.at(posInd.value()) - x_.at(negInd.value()));
        }
        // φ2n-1
        if (nc2 == NodeConfig::POSGND) {
          temp.n1_2_ = (x_.at(posInd2.value()));
        } else if (nc2 == NodeConfig::GNDNEG) {
          temp.n1_2_ = (-x_.at(negInd2.value()));
        } else {
          temp.n1_2_ = (x_.at(posInd2.value()) - x_.at(negInd2.value()));
        }
      }
      if (i >= k) {
        // φ1n-k
        if (nc == NodeConfig::POSGND) {
          temp.nk_1_ = results.xVector.at(posInd.value()).value().at(i - k);
        } else if (nc == NodeConfig::GNDNEG) {
          temp.nk_1_ = -results.xVector.at(negInd.value()).value().at(i - k);
        } else {
          temp.nk_1_ = results.xVector.at(posInd.value()).value().at(i - k) -
                       results.xVector.at(negInd.value()).value().at(i - k);
        }
        // φ2n-k
        if (nc2 == NodeConfig::POSGND) {
          temp.nk_2_ = results.xVector.at(posInd2.value()).value().at(i - k);
        } else if (nc2 == NodeConfig::GNDNEG) {
          temp.nk_2_ = -results.xVector.at(negInd2.value()).value().at(i - k);
        } else {
          temp.nk_2_ = results.xVector.at(posInd2.value()).value().at(i - k) -
                       results.xVector.at(negInd2.value()).value().at(i - k);
        }
        // I1n-k
        double &I1nk = results.xVector.at(curInd).value().at(i - k);
        // I2n-k
        double &I2nk = results.xVector.at(curInd2).value().at(i - k);
        if (i == k) {
          // I1 = Z(2e/hbar)(2h/3)I2n-k + (4/3)φ1n-1 - (1/3)φ1n-2 + φ2n-k
          b_.at(curInd) = (Z / Constants::SIGMA) *
                              ((2.0 * stepSize_ * factor) / 3.0) * I2nk +
                          (4.0 / 3.0) * temp.n1_1_ - (1.0 / 3.0) * temp.n2_1_ +
                          temp.nk_2_;
          // I2 = Z(2e/hbar)(2h/3)I1n-k + (4/3)φ2n-1 - (1/3)φ2n-2 + φ1n-k
          b_.at(curInd2) = (Z / Constants::SIGMA) *
                               ((2.0 * stepSize_ * factor) / 3.0) * I1nk +
                           (4.0 / 3.0) * temp.n1_2_ - (1.0 / 3.0) * temp.n2_2_ +
                           temp.nk_1_;
        } else if (i == k + 1) {
          // φ1n-k-1
          if (nc == NodeConfig::POSGND) {
            temp.nk1_1_ =
                results.xVector.at(posInd.value()).value().at(i - k - 1);
          } else if (nc == NodeConfig::GNDNEG) {
            temp.nk1_1_ =
                -results.xVector.at(negInd.value()).value().at(i - k - 1);
          } else {
            temp.nk1_1_ =
                results.xVector.at(posInd.value()).value().at(i - k - 1) -
                results.xVector.at(negInd.value()).value().at(i - k - 1);
          }
          // φ2n-k-1
          if (nc2 == NodeConfig::POSGND) {
            temp.nk1_2_ =
                results.xVector.at(posInd2.value()).value().at(i - k - 1);
          } else if (nc2 == NodeConfig::GNDNEG) {
            temp.nk1_2_ =
                -results.xVector.at(negInd2.value()).value().at(i - k - 1);
          } else {
            temp.nk1_2_ =
                results.xVector.at(posInd2.value()).value().at(i - k - 1) -
                results.xVector.at(negInd2.value()).value().at(i - k - 1);
          }
          // I1 = Z(2e/hbar)(2h/3)I2n-k + (4/3)φ1n-1 - (1/3)φ1n-2 +
          //      φ2n-k - (4/3)φ2n-k-1
          b_.at(curInd) = (Z / Constants::SIGMA) *
                              ((2.0 * stepSize_ * factor) / 3.0) * I2nk +
                          (4.0 / 3.0) * temp.n1_1_ - (1.0 / 3.0) * temp.n2_1_ +
                          temp.nk_2_ - (4.0 / 3.0) * temp.nk1_2_;
          // I2 = Z(2e/hbar)(2h/3)I1n-k + (4/3)φ2n-1 - (1/3)φ2n-2 +
          //      φ1n-k - (4/3)φ1n-k-1
          b_.at(curInd2) = (Z / Constants::SIGMA) *
                               ((2.0 * stepSize_ * factor) / 3.0) * I1nk +
                           (4.0 / 3.0) * temp.n1_2_ - (1.0 / 3.0) * temp.n2_2_ +
                           temp.nk_1_ - (4.0 / 3.0) * temp.nk1_1_;
        } else if (i > k + 1) {
          // φ1n-k-1
          if (nc == NodeConfig::POSGND) {
            temp.nk1_1_ =
                results.xVector.at(posInd.value()).value().at(i - k - 1);
          } else if (nc == NodeConfig::GNDNEG) {
            temp.nk1_1_ =
                -results.xVector.at(negInd.value()).value().at(i - k - 1);
          } else {
            temp.nk1_1_ =
                results.xVector.at(posInd.value()).value().at(i - k - 1) -
                results.xVector.at(negInd.value()).value().at(i - k - 1);
          }
          // φ2n-k-1
          if (nc2 == NodeConfig::POSGND) {
            temp.nk1_2_ =
                results.xVector.at(posInd2.value()).value().at(i - k - 1);
          } else if (nc2 == NodeConfig::GNDNEG) {
            temp.nk1_2_ =
                -results.xVector.at(negInd2.value()).value().at(i - k - 1);
          } else {
            temp.nk1_2_ =
                results.xVector.at(posInd2.value()).value().at(i - k - 1) -
                results.xVector.at(negInd2.value()).value().at(i - k - 1);
          }
          // φ1n-k-2
          if (nc == NodeConfig::POSGND) {
            temp.nk2_1_ =
                results.xVector.at(posInd.value()).value().at(i - k - 2);
          } else if (nc == NodeConfig::GNDNEG) {
            temp.nk2_1_ =
                -results.xVector.at(negInd.value()).value().at(i - k - 2);
          } else {
            temp.nk2_1_ =
                results.xVector.at(posInd.value()).value().at(i - k - 2) -
                results.xVector.at(negInd.value()).value().at(i - k - 2);
          }
          // φ2n-k-2
          if (nc2 == NodeConfig::POSGND) {
            temp.nk2_2_ =
                results.xVector.at(posInd2.value()).value().at(i - k - 2);
          } else if (nc2 == NodeConfig::GNDNEG) {
            temp.nk2_2_ =
                -results.xVector.at(negInd2.value()).value().at(i - k - 2);
          } else {
            temp.nk2_2_ =
                results.xVector.at(posInd2.value()).value().at(i - k - 2) -
                results.xVector.at(negInd2.value()).value().at(i - k - 2);
          }
          // I1 = Z(2e/hbar)(2h/3)I2n-k + (4/3)φ1n-1 - (1/3)φ1n-2 +
          //      φ2n-k - (4/3)φ2n-k-1 + (1/3)φ2n-k-2
          b_.at(curInd) = (Z / Constants::SIGMA) *
                              ((2.0 * stepSize_ * factor) / 3.0) * I2nk +
                          (4.0 / 3.0) * temp.n1_1_ - (1.0 / 3.0) * temp.n2_1_ +
                          temp.nk_2_ - (4.0 / 3.0) * temp.nk1_2_ +
                          (1.0 / 3.0) * temp.nk2_2_;
          // I2 = Z(2e/hbar)(2h/3)I1n-k + (4/3)φ2n-1 - (1/3)φ2n-2 +
          //      φ1n-k - (4/3)φ1n-k-1 + (1/3)φ1n-k-2
          b_.at(curInd2) = (Z / Constants::SIGMA) *
                               ((2.0 * stepSize_ * factor) / 3.0) * I1nk +
                           (4.0 / 3.0) * temp.n1_2_ - (1.0 / 3.0) * temp.n2_2_ +
                           temp.nk_1_ - (4.0 / 3.0) * temp.nk1_1_ +
                           (1.0 / 3.0) * temp.nk2_1_;
        }
      } else {
        // I1 = (4/3)φ1n-1 - (1/3)φ1n-2
        b_.at(curInd) = (4.0 / 3.0) * temp.n1_1_ - (1.0 / 3.0) * temp.n2_1_;
        // I2 = (4/3)φ2n-1 - (1/3)φ2n-2
        b_.at(curInd2) = (4.0 / 3.0) * temp.n1_2_ - (1.0 / 3.0) * temp.n2_2_;
      }
    }
  }
}