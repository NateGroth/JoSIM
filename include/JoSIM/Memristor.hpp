// Copyright (c) 2021 Johannes Delport
// This code is licensed under MIT license (see LICENSE for details)
//
// aether_sims D11: cuprate-memristor primitive (trap-based compact model).
//
// Model source: Guenkel, Barrera, Balcells, Mestres, Miranda, Palau, Sune,
// "SPICE-Compatible Compact Modeling of Cuprate-Based Memristors Across a
// Wide Temperature Range", Adv. Electron. Mater. 2026, e00861
// (doi 10.1002/aelm.202500861). State variable g in [0,1] is the fraction
// of the YBCO layer occupied by the low-resistance layer (LRL):
//
//   memory:   dg/dt = Kp(1-g) - H(g - gmin) Kd (g - gmin)
//   gmin(V) = 1/(1 + exp((V - Vlog)/xi))                      (logistic)
//   Kp(V)   = H(-V) Kp0 [ |V|^m(g) / R(g) ]^alpha             (SCLC-driven)
//   Kd(V)   = H(V)  Kd0 exp(etad V)                           (detrapping)
//   current:  I(g,V) = sgn(V) |V|^m(g) / R(g)                 (SCLC)
//   R(g)    = Rmax Rmin / (Rmax g + Rmin (1-g))
//   m(g)    = mmin g + mmax (1-g)
//
// Temperature (T_loc, a device state coupled to the Python lumped thermal
// network exactly like the JJ's D1 path):
//   mmax(T) = 1 + (mmax_cal - 1) Tcal / T     (SCLC trap law m = 1 + Tt/T)
//   Rmax(T) = Rmax_cal exp(-AH (T - Tcal))    (HRS shrinks exponentially)
//   Rmin(T) = Rmin_cal exp(+AL (T - Tcal))    (LRS grows exponentially)
//   mmin, Kp0, alpha, Kd0, etad, Vlog, xi     held T-independent (paper).
// Calibration range of record: Table 1 fits at 150-300 K (measured device
// range 80-300 K). T_loc outside [TRMIN,TRMAX] flags EXTRAPOLATED and
// warns once (D11 requirement -- the 65 K operating point IS extrapolated).
//
// MNA formulation: fixed companion resistance R0 stamped exactly like a
// Resistor (constant matrix, no per-step refactor); the nonlinearity and
// state ride the branch-row RHS as a per-step Thevenin correction
//   V - R0 I = e_n,   e_n = V_prev - R0 I_model(V_prev, g_n).
// With R0 = Rmin(Tcal) (max device conductance) the one-step-lag fixed
// point contracts for any resistive surrounding network (G_dev <= G0).
// The device law is explicit in V_prev -- same one-step-lag posture as the
// D8 CTRL coupling; the ms-scale device dynamics are >> any engine step.
// The g update uses the paper's exact exponential (recursive) solution,
// Eqs (9)/(10), unconditionally stable for any timestep.
#ifndef JOSIM_MEMRISTOR_HPP
#define JOSIM_MEMRISTOR_HPP

#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "JoSIM/AnalysisType.hpp"
#include "JoSIM/BasicComponent.hpp"
#include "JoSIM/Input.hpp"
#include "JoSIM/Model.hpp"
#include "JoSIM/ParameterName.hpp"
#include "JoSIM/Parameters.hpp"
#include "JoSIM/Spread.hpp"

namespace JoSIM {

/*
 YMlabel V+ V- model [R0=ohm] [G0=initial g] [TLOC=K]

 Branch relation (voltage mode), R0 fixed:
 ⎡ 0  0   1⎤ ⎡V+⎤   ⎡  0⎤
 ⎜ 0  0  -1⎟ ⎜V-⎟ = ⎜  0⎟
 ⎣ 1 -1 -R0⎦ ⎣Io⎦   ⎣ e_n⎦
 (PHASE mode: same nonzero pattern as a Resistor of value R0, and the
 branch RHS additionally carries (2h/3)(1/sigma) e_n.)
*/
class Memristor : public BasicComponent {
 private:
  AnalysisType at_;

 public:
  Model model_;
  // BDF2 phase history (phase mode), mirroring Resistor
  double pn1_ = 0.0, pn2_ = 0.0, pn3_ = 0.0, pn4_ = 0.0;
  double r0_ = 0.0;          //: fixed companion resistance (ohm)
  double g_ = 0.0;           //: state variable (LRL fraction), [0,1]
  double gPrevAccept_ = 0.0; //: g history for step_back
  double gPrev2_ = 0.0;
  double tLoc_;              //: local element temperature (K)
  double eN_ = 0.0;          //: current RHS correction (V)
  double vPrev_ = 0.0;       //: device voltage at previous step (V)
  double iPrev_ = 0.0;       //: branch current at previous step (A)
  double power_ = 0.0;       //: instantaneous V*I (W)
  bool extrapolated_ = false;    //: T_loc left the calibration range
  bool extrapWarned_ = false;    //: warn-once latch

  Memristor(const std::pair<tokens_t, string_o>& s, const NodeConfig& ncon,
            const nodemap& nm, std::unordered_set<std::string>& lm,
            nodeconnections& nc, Input& iObj, Spread& spread, int64_t& bi);

  void set_model(const tokens_t& t,
                 const vector_pair_t<Model, string_o>& models,
                 const string_o& subc);

  // -- temperature-dependent parameter laws (evaluated at tLoc_) --------
  double rmax_t() const;
  double rmin_t() const;
  double mmax_t() const;
  double r_of_g(const double& g) const;
  double m_of_g(const double& g) const;
  double i_model(const double& v, const double& g) const;

  // Advance g by dt with the device voltage held at v (exact exponential
  // update, paper Eqs (9)/(10)), then refresh e_n for the next solve.
  void update_state(const double& v, const double& dt);

  void set_temperature(const double& t);

  void update_timestep(const double& factor) override;

  void step_back() override {
    pn2_ = pn4_;
    g_ = gPrev2_;
  }
};  // class Memristor

}  // namespace JoSIM
#endif
