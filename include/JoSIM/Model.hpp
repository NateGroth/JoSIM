// Copyright (c) 2021 Johannes Delport
// This code is licensed under MIT license (see LICENSE for details)
#ifndef JOSIM_MODEL_HPP
#define JOSIM_MODEL_HPP

#include "JoSIM/Constants.hpp"
#include "JoSIM/Parameters.hpp"
#include "JoSIM/TypeDefines.hpp"

namespace JoSIM {
class Model {
 private:
  std::string modelName_;
  double vg_;
  double ic_;
  std::vector<double> cpr_;
  int64_t rtype_;
  double rn_;
  double r0_;
  double c_;
  double t_;
  double tc_;
  double deltaV_;
  double d_;
  double icFct_;
  double phiOff_;
  bool tDep_;
  // aether_sims D2: selectable Ic(T) law. ictemp_ = 0 -> BCS / Ambegaokar-
  // Baratoff (default, preserves upstream behaviour); 1 -> YBCO weak link
  // Ic(T) = Ic0 (1 - T/Tc)^wlpow_.
  int64_t ictemp_;
  double wlpow_;
  // aether_sims D8: phenomenological multiterminal coupling law
  // Ic(Ictrl) = Ic0 * g(x), x = Ictrl / ctrlNorm_. The transfer function g is
  // a fitted *input* (digitized from measurement, e.g. Wisne & Chandrasekhar
  // arXiv:2507.14357), never an engine prediction. icctrl_ = 0 -> off
  // (upstream parity); 1 -> POLY, g = sum_k ctrlCoef_[k] x^k; 2 -> TABLE,
  // piecewise-linear through ctrlTab_ = {x1,g1, x2,g2, ...}. x is clamped to
  // the fitted domain and g floored at 0. A table given only on x >= 0 is
  // treated as even, g(|x|).
  int64_t icctrl_;
  double ctrlNorm_;
  std::vector<double> ctrlCoef_;
  std::vector<double> ctrlTab_;
  // aether_sims T3: first-order lag (seconds) on the Ic(Ictrl) coupling --
  // the control current seen by the law is a low-pass of the branch current.
  // 0 (default) = memoryless D8 semantics, bit-exact legacy behavior.
  double ctrlLag_;
  // aether_sims T4-D9: odd fraction s of the signed coupling law
  // Ic_eff+- = Ic0*[(1-s)*f_even(|x|) + s*(1 -+ x)]; 0 = legacy even law.
  double ctrlOdd_;
  // aether_sims D11: model discriminator. 0 = JJ (upstream semantics,
  // untouched); 1 = MEMRISTOR (cuprate trap-based compact model,
  // Guenkel et al. Adv. Electron. Mater. 2026 e00861). Memristor models
  // reuse this class so the single .MODEL parse/lookup path serves both.
  int64_t mtype_;
  // -- D11 memristor parameters (valid when mtype_ == 1) ----------------
  double mrRmax_;    //: HRS resistance at Tcal (ohm)
  double mrRmin_;    //: LRS resistance at Tcal (ohm)
  double mrMmax_;    //: HRS SCLC exponent at Tcal
  double mrMmin_;    //: LRS SCLC exponent (T-independent, paper Table 1)
  double mrKp0_;     //: potentiation rate prefactor ((A)^-alpha / s)
  double mrAlpha_;   //: potentiation current exponent
  double mrKd0_;     //: depression rate prefactor (1/s)
  double mrEtad_;    //: depression field-acceleration (1/V)
  double mrVlog_;    //: gmin logistic center (V)
  double mrXi_;      //: gmin logistic width (V)
  double mrTcal_;    //: calibration temperature of the R/m values (K)
  double mrAh_;      //: HRS resistance T-slope, Rmax ~ exp(-AH (T-Tcal))
  double mrAl_;      //: LRS resistance T-slope, Rmin ~ exp(+AL (T-Tcal))
  double mrTrmin_;   //: calibration-range lower bound (K)
  double mrTrmax_;   //: calibration-range upper bound (K)
  double mrTscale_;  //: rate time-scale factor (test economy; law-invariant)

 public:
  Model()
      : vg_(2.8E-3),
        ic_(1E-3),
        cpr_({1.0}),
        rtype_(1),
        rn_(5),
        r0_(30),
        c_(2.5E-12),
        t_(4.2),
        tc_(9.1),
        deltaV_(0.1E-3),
        d_(0),
        icFct_(Constants::PI / 4),
        phiOff_(0),
        tDep_(false),
        ictemp_(0),
        wlpow_(1.0),
        icctrl_(0),
        ctrlNorm_(1.0),
        ctrlCoef_({1.0}),
        ctrlTab_({}),
        ctrlLag_(0.0),
        ctrlOdd_(0.0),
        mtype_(0),
        // D11 defaults: paper Table 1 at 150 K with Kd0/etad/Vlog/xi held
        // T-independent; AH/AL from the Table-1 150->300 K endpoints;
        // calibration range = the fitted 150-300 K.
        mrRmax_(1.5E6),
        mrRmin_(1.25E4),
        mrMmax_(4.2),
        mrMmin_(3.0),
        mrKp0_(1E23),
        mrAlpha_(5.0),
        mrKd0_(1.5),
        mrEtad_(0.75),
        mrVlog_(1.6),
        mrXi_(0.25),
        mrTcal_(150.0),
        mrAh_(0.012077),
        mrAl_(0.0068635),
        mrTrmin_(150.0),
        mrTrmax_(300.0),
        mrTscale_(1.0){};

  std::string modelName() const { return modelName_; }
  void modelName(const std::string& n) { modelName_ = n; }
  double vg() const { return vg_; }
  void vg(const double& v) { vg_ = v; }
  double ic() const { return ic_; }
  void ic(const double& i) { ic_ = i; }
  std::vector<double> cpr() const { return cpr_; }
  void cpr(const std::vector<double>& i) { cpr_ = i; }
  int64_t rtype() const { return rtype_; }
  void rtype(const int64_t& r) { rtype_ = r; }
  double rn() const { return rn_; }
  void rn(const double& r) { rn_ = r; }
  double r0() const { return r0_; }
  void r0(const double& r) { r0_ = r; }
  double c() const { return c_; }
  void c(const double& c) { c_ = c; }
  double t() const { return t_; }
  void t(const double& t) { t_ = t; }
  double tc() const { return tc_; }
  void tc(const double& t) { tc_ = t; }
  double deltaV() const { return deltaV_; }
  void deltaV(const double& d) { deltaV_ = d; }
  double d() const { return d_; }
  void d(const double& t) { d_ = t; }
  double icFct() const { return icFct_; }
  void icFct(const double& r) { icFct_ = r; }
  double phiOff() const { return phiOff_; }
  void phiOff(const double& o) { phiOff_ = o; }
  bool tDep() { return tDep_; }
  void tDep(bool b) { tDep_ = b; }
  int64_t ictemp() const { return ictemp_; }
  void ictemp(const int64_t& i) { ictemp_ = i; }
  double wlpow() const { return wlpow_; }
  void wlpow(const double& n) { wlpow_ = n; }
  int64_t icctrl() const { return icctrl_; }
  void icctrl(const int64_t& i) { icctrl_ = i; }
  double ctrlNorm() const { return ctrlNorm_; }
  void ctrlNorm(const double& n) { ctrlNorm_ = n; }
  double ctrlLag() const { return ctrlLag_; }
  void ctrlLag(const double& l) { ctrlLag_ = l; }
  double ctrlOdd() const { return ctrlOdd_; }
  void ctrlOdd(const double& s) { ctrlOdd_ = s; }
  const std::vector<double>& ctrlCoef() const { return ctrlCoef_; }
  void ctrlCoef(const std::vector<double>& c) { ctrlCoef_ = c; }
  const std::vector<double>& ctrlTab() const { return ctrlTab_; }
  void ctrlTab(const std::vector<double>& t) { ctrlTab_ = t; }
  // -- D11 memristor accessors ------------------------------------------
  int64_t mtype() const { return mtype_; }
  void mtype(const int64_t& m) { mtype_ = m; }
  double mrRmax() const { return mrRmax_; }
  void mrRmax(const double& v) { mrRmax_ = v; }
  double mrRmin() const { return mrRmin_; }
  void mrRmin(const double& v) { mrRmin_ = v; }
  double mrMmax() const { return mrMmax_; }
  void mrMmax(const double& v) { mrMmax_ = v; }
  double mrMmin() const { return mrMmin_; }
  void mrMmin(const double& v) { mrMmin_ = v; }
  double mrKp0() const { return mrKp0_; }
  void mrKp0(const double& v) { mrKp0_ = v; }
  double mrAlpha() const { return mrAlpha_; }
  void mrAlpha(const double& v) { mrAlpha_ = v; }
  double mrKd0() const { return mrKd0_; }
  void mrKd0(const double& v) { mrKd0_ = v; }
  double mrEtad() const { return mrEtad_; }
  void mrEtad(const double& v) { mrEtad_ = v; }
  double mrVlog() const { return mrVlog_; }
  void mrVlog(const double& v) { mrVlog_ = v; }
  double mrXi() const { return mrXi_; }
  void mrXi(const double& v) { mrXi_ = v; }
  double mrTcal() const { return mrTcal_; }
  void mrTcal(const double& v) { mrTcal_ = v; }
  double mrAh() const { return mrAh_; }
  void mrAh(const double& v) { mrAh_ = v; }
  double mrAl() const { return mrAl_; }
  void mrAl(const double& v) { mrAl_ = v; }
  double mrTrmin() const { return mrTrmin_; }
  void mrTrmin(const double& v) { mrTrmin_ = v; }
  double mrTrmax() const { return mrTrmax_; }
  void mrTrmax(const double& v) { mrTrmax_ = v; }
  double mrTscale() const { return mrTscale_; }
  void mrTscale(const double& v) { mrTscale_ = v; }
  // aether_sims D8: evaluate the dimensionless suppression factor g(Ictrl).
  // Memoryless by construction (pure function of the instantaneous control
  // current -- required until D7 per-step LTE control lands, since a
  // reduce_step restart discards accumulated history).
  double ctrl_scale(double ictrl) const {
    if (icctrl_ == 0) return 1.0;
    double x = (ctrlNorm_ != 0.0) ? ictrl / ctrlNorm_ : 0.0;
    double g = 1.0;
    if (icctrl_ == 1) {  // POLY over signed x, clamped to the unit domain
      if (x > 1.0) x = 1.0;
      if (x < -1.0) x = -1.0;
      g = 0.0;
      double xk = 1.0;
      for (const auto& c : ctrlCoef_) {
        g += c * xk;
        xk *= x;
      }
    } else if (icctrl_ == 2 && ctrlTab_.size() >= 4) {  // TABLE (pw-linear)
      // A table defined only on x >= 0 is even: evaluate at |x|.
      if (ctrlTab_.front() >= 0.0 && x < 0.0) x = -x;
      const auto& t = ctrlTab_;
      size_t n = t.size() / 2;
      if (x <= t[0]) {
        g = t[1];
      } else if (x >= t[2 * (n - 1)]) {
        g = t[2 * (n - 1) + 1];
      } else {
        for (size_t k = 0; k + 1 < n; ++k) {
          double x0 = t[2 * k], x1 = t[2 * (k + 1)];
          if (x >= x0 && x <= x1) {
            double f = (x1 > x0) ? (x - x0) / (x1 - x0) : 0.0;
            g = t[2 * k + 1] + f * (t[2 * (k + 1) + 1] - t[2 * k + 1]);
            break;
          }
        }
      }
    }
    return g < 0.0 ? 0.0 : g;
  }
  static void parse_model(const std::pair<tokens_t, string_o>& s,
                          vector_pair_t<Model, string_o>& models,
                          const param_map& p);
};
}  // namespace JoSIM

#endif