// Copyright (c) 2021 Johannes Delport
// This code is licensed under MIT license (see LICENSE for details)
#ifndef JOSIM_SIMULATION_H
#define JOSIM_SIMULATION_H

#include <suitesparse/klu.h>

#include <cassert>
#include <functional>

#include "JoSIM/Errors.hpp"
#include "JoSIM/LUSolve.hpp"
#include "JoSIM/Matrix.hpp"
#include "JoSIM/Misc.hpp"

namespace JoSIM {

#define TRANSIENT 0
#define DC 1
#define AC 2
#define PHASE 3
#define NONE_SPECIFIED 4

class Results {
 public:
  std::vector<std::optional<std::vector<double>>> xVector;
  std::vector<double> timeAxis;
};

class Simulation {
 private:
  bool SLU = false;
  std::vector<double> x_, b_;
  int64_t simSize_;
  JoSIM::AnalysisType atyp_;
  bool minOut_;
  bool needsLU_;
  bool needsTR_ = true;
  bool startup_;
  bool kluReady_ = false;  // KLU factored and not yet freed (aether_sims D4)
  double stepSize_, prstep_, prstart_;
#ifdef SLU
  LUSolve lu;
#else
  int64_t simOK_;
  klu_l_symbolic* Symbolic_;
  klu_l_common Common_;
  klu_l_numeric* Numeric_;
#endif

  void setup(Input& iObj, Matrix& mObj);
  void trans_sim(Matrix& mObj);
  void setup_b(Matrix& mObj, int64_t i, double step, double factor = 1);
  void reduce_step(Input& iObj, Matrix& mObj);

  // aether_sims D4 helpers (shared by the batch and stepped paths).
  void prepare(Input& iObj, Matrix& mObj);  // setup + KLU analyze/factor
  void run_startup(Matrix& mObj);           // pre-t=0 stabilization loop
  bool solve_only(int64_t i, Matrix& mObj); // setup_b + solve, no store

  void handle_cs(Matrix& mObj, double& step, const int64_t& i);
  void handle_resistors(Matrix& mObj, double& step);
  void handle_inductors(Matrix& mObj, double factor = 1);
  void handle_capacitors(Matrix& mObj);
  void handle_jj(Matrix& mObj, int64_t& i, double& step, double factor = 1);
  void handle_vs(Matrix& mObj, const int64_t& i, double& step,
                 double factor = 1);
  void handle_ps(Matrix& mObj, const int64_t& i, double& step,
                 double factor = 1);
  void handle_ccvs(Matrix& mObj);
  void handle_vccs(Matrix& mObj);
  void handle_tx(Matrix& mObj, const int64_t& i, double& step,
                 double factor = 1);

 public:
  Results results;

  // Co-simulation hooks (aether_sims D4). If set, trans_sim() invokes them
  // around each main-loop step with (index, time); aether_sims uses pre_ to
  // push temperature down and post_ to read dissipated power up.
  std::function<void(int64_t, double)> pre_step_hook_;
  std::function<void(int64_t, double)> post_step_hook_;

  // Batch: runs the whole transient in the constructor (unchanged behaviour).
  Simulation(Input& iObj, Matrix& mObj);
  // Deferred: setup + factor (+ startup stabilization) but NOT the main loop.
  // Drive it with step() in a loop, then finish() (or let the destructor free).
  Simulation(Input& iObj, Matrix& mObj, bool deferRun);
  ~Simulation();

  // One main-loop transient step (setup_b + solve + store). Returns true if a
  // timestep reduction is requested (needsTR_); keep dt small enough when
  // driving stepped so this stays false (the global restart is a P3 concern).
  bool step(int64_t i, Matrix& mObj);

  // Run the whole main loop (post-startup) invoking the hooks each step. For
  // the C++-callback co-sim on a deferred Simulation: set the hooks, then call
  // run_main(). (The Python-driven alternative is a step() loop.)
  void run_main(Matrix& mObj);

  // Free the factored solver (idempotent). Called by the stepped driver when
  // done; also run by the destructor.
  void finish();

  // D7 / reduce_step support for the stepped path. needs_reduction() is true
  // when the last step() requested a timestep reduction (a JJ phase guess
  // jumped too far). reduce_and_restart() halves the timestep, rebuilds the
  // matrix and re-prepares; the caller must then restart its step() loop from
  // i=0 (sim_size() changes). This mirrors what the batch constructor does
  // internally, so a stepped run can match batch on stiff circuits.
  bool needs_reduction() const { return needsTR_; }
  void reduce_and_restart(Input& iObj, Matrix& mObj);

  // Per-junction electro-thermal access by JoSIM label (aether_sims D1/D2).
  // set_*_temperature drive JJ::update_temperature(); jj_ic reads the current
  // (temperature-dependent) critical current.
  void set_jj_temperature(Matrix& mObj, const std::string& label, double T);
  void set_all_temperatures(Matrix& mObj, double T);
  // aether_sims D8: drive a junction's Ic(Ictrl) law with an explicit control
  // current (static trim / Python-hook path). Junctions with a CTRL= netlist
  // binding are instead fed by the engine after every solve.
  void set_jj_control_current(Matrix& mObj, const std::string& label,
                              double ictrl);
  double jj_ic(Matrix& mObj, const std::string& label);
  // Live per-junction state read from the current solution (call after step()).
  double jj_phase(Matrix& mObj, const std::string& label);
  double jj_voltage(Matrix& mObj, const std::string& label);
  double jj_current(Matrix& mObj, const std::string& label);
  double jj_power(Matrix& mObj, const std::string& label);

  int64_t sim_size() const { return simSize_; }
  double step_size() const { return stepSize_; }
};
}  // namespace JoSIM
#endif