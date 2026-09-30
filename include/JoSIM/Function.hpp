// Copyright (c) 2021 Johannes Delport
// This code is licensed under MIT license (see LICENSE for details)
#ifndef JOSIM_FUNCTION_HPP
#define JOSIM_FUNCTION_HPP

#include <cstdint>
#include <limits>
#include <vector>

#include "Input.hpp"

namespace JoSIM {

enum class FunctionType {
  PWL = 0,
  PULSE = 1,
  SINUSOID = 2,
  CUS = 3,
  NOISE = 4,
  PWS = 5,
  DC = 6,
  EXP = 7
};

class Function {
 private:
  FunctionType fType_ = FunctionType::PWL;
  std::vector<double> timeValues_;
  std::vector<double> ampValues_;
  std::vector<double> miscValues_;
  // aether_sims T33-R1 (Johnson fix): NOISE(VA TD TSTEP) is a HELD process on
  // the absolute time grid TD + k*TSTEP -- sample k is N(0, VA^2 / (2 TSTEP)),
  // a pure function of (noiseKey_, k) -- and the engine sees its exact charge
  // in every step: value(t) = the average of the held process over
  // (t - h, t], h the engine step of the pass. Delivered power is then
  // independent of the step the engine takes (a halved-step restart included),
  // of TSTEP >= or < h, and of how many times a step evaluates the source (a
  // two-node element stamps it twice).
  uint64_t noiseKey_ = 0;
  double noiseSigma_ = 0.0;  // per-sample standard deviation VA / sqrt(2 TSTEP)
  double noiseH_ = 0.0;      // the engine step of the pass (parse time)
  int64_t noiseCk_[2] = {-1, -1};
  double noiseCs_[2] = {0.0, 0.0};
  int noiseCnext_ = 0;
  double noiseLastX_ = std::numeric_limits<double>::quiet_NaN();
  double noiseLastV_ = 0.0;
  double noise_sample(int64_t k);
  void parse_pwl(const tokens_t& t, const Input& iObj, const string_o& s);
  void parse_pulse(const tokens_t& t, const Input& iObj, const string_o& s);
  void parse_sin(const tokens_t& t, const Input& iObj, const string_o& s);
  void parse_cus(const tokens_t& t, const Input& iObj, const string_o& s);
  void parse_noise(const tokens_t& t, const Input& iObj, const string_o& s);
  void parse_dc(const tokens_t& t, const Input& iObj, const string_o& s);
  void parse_exp(const tokens_t& t, const Input& iObj, const string_o& s);
  double return_pwl(double& x);
  double return_pulse(double& x);
  double return_sin(double& x);
  double return_cus(double& x);
  double return_noise(double& x);
  double return_pws(double& x);
  double return_dc();
  double return_exp(double& x);

 public:
  Function(){};
  void parse_function(const std::string& str, const Input& iObj,
                      const string_o& subckt);
  double value(double x);
  void ampValues(std::vector<double> values);
  std::vector<double> ampValues() { return ampValues_; }
  void clearMisc() { miscValues_.clear(); }

};  // class Function

}  // namespace JoSIM

#endif  // JOSIM_FUNCTION_HPP
