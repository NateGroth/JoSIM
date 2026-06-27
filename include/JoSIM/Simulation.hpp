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

  int64_t sim_size() const { return simSize_; }
  double step_size() const { return stepSize_; }
};
}  // namespace JoSIM
#endif