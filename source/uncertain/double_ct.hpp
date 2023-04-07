// License Terms
//
// Copyright (c) 2019, California Institute of Technology ("Caltech").
// U.S. Government sponsorship acknowledged.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright notice,
//       this list of conditions and the following disclaimer.
//     * Redistributions must reproduce the above copyright notice, this list
//       of conditions and the following disclaimer in the documentation and/or
//       other materials provided with the distribution.
//     * Neither the name of Caltech nor its operating division, the Jet
//       Propulsion Laboratory, nor the names of its contributors may be used
//       to endorse or promote products derived from this software without
//       specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// udouble.h: This file includes classes for propagation of uncertainties
// and supporting classes and functions.
// By Evan Manning (manning@alumni.caltech.edu).

#pragma once

#include <functional>
#include <uncertain/source_set.hpp>

namespace uncertain {

// Correlation tracking class keeps an array of uncertainty
// components from various sources.  (Array implementation
// is specified by template parameter.)
template <class T>
class UDoubleCT {
 private:
  double value;
  T unc_components;
  size_t epoch;
  static SourceSet sources;

 public:
  // default constructor creates a new independent uncertainty element
  UDoubleCT(double val = 0.0, double unc = 0.0, const std::string &name = {}) : value(val) {
    epoch = sources.get_epoch();
    if (unc < 0.0) {
      throw std::runtime_error("Error: negative uncertainty: " + std::to_string(unc));
    }
    if (unc != 0.0) {
      std::string source_name;
      if (!name.empty()) {
        source_name = name;
      } else {
        std::stringstream os;
        os << "anon: ";
        uncertain_print(val, unc, os);
        source_name = os.str();
      }
      size_t new_source_num = sources.get_new_source(source_name);
      unc_components.set_element(new_source_num, unc);
    }
  }

  // copy constructor does not create a new independent uncertainty element
  UDoubleCT(const UDoubleCT &ud) = default;

  // \todo operator= ?

  ~UDoubleCT() = default;

  double mean() const { return value; }

  double deviation() const { return unc_components.norm(); }

  static void new_epoch() { sources.new_epoch(); }

  void print_uncertain_sources(std::ostream &os = std::cout) const {
    double total_uncertainty = this->deviation();
    if (total_uncertainty == 0.0)
      os << "No uncertainty";
    else
      for (unsigned int i = 0; i < sources.get_num_sources(); i++) {
        double unc_portion = this->unc_components[i] / total_uncertainty;
        unc_portion *= unc_portion;
        os << "[" << i << "] " << sources.get_source_name(i) << ": " << int_percent(unc_portion)
           << "% (" << this->unc_components[i] << ")\n";
      }
    os << std::endl;
  }

  UDoubleCT operator+() const { return *this; }

  UDoubleCT operator-() const {
    UDoubleCT retval;

    retval.value = -value;
    retval.unc_components = -unc_components;
    return retval;
  }

  UDoubleCT &operator+=(const UDoubleCT &b) {
    sources.check_epoch(epoch);
    sources.check_epoch(b.epoch);
    unc_components += b.unc_components;
    value += b.value;
    return *this;
  }

  UDoubleCT &operator+=(double b) {
    value += b;
    return *this;
  }

  friend UDoubleCT operator+(UDoubleCT a, const UDoubleCT &b) { return a += b; }

  friend UDoubleCT operator+(UDoubleCT a, double b) { return a += b; }

  friend UDoubleCT operator+(double b, UDoubleCT a) { return a += b; }

  UDoubleCT &operator-=(const UDoubleCT &b) {
    sources.check_epoch(epoch);
    sources.check_epoch(b.epoch);
    unc_components -= b.unc_components;
    value -= b.value;
    return *this;
  }

  UDoubleCT &operator-=(double b) {
    value -= b;
    return *this;
  }

  friend UDoubleCT operator-(UDoubleCT a, const UDoubleCT &b) { return a -= b; }

  friend UDoubleCT operator-(UDoubleCT a, double b) { return a -= b; }

  friend UDoubleCT operator-(double b, UDoubleCT a) {
    a -= b;
    return -a;
  }

  UDoubleCT &operator*=(const UDoubleCT &b) {
    sources.check_epoch(epoch);
    sources.check_epoch(b.epoch);
    unc_components *= b.value;
    unc_components += b.unc_components * value;
    value *= b.value;
    return *this;
  }

  UDoubleCT &operator*=(double b) {
    unc_components *= b;
    value *= b;
    return *this;
  }

  UDoubleCT operator++() { return (*this += 1.0); }

  UDoubleCT operator--() { return (*this -= 1.0); }

  UDoubleCT operator++(int) {
    UDoubleCT retval(*this);
    *this += 1.0;
    return retval;
  }

  UDoubleCT operator--(int) {
    UDoubleCT retval(*this);
    *this -= 1.0;
    return retval;
  }

  friend UDoubleCT operator*(UDoubleCT a, const UDoubleCT &b) { return a *= b; }

  friend UDoubleCT operator*(UDoubleCT a, double b) { return a *= b; }

  friend UDoubleCT operator*(double b, UDoubleCT a) { return a *= b; }

  UDoubleCT &operator/=(const UDoubleCT &b) {
    sources.check_epoch(epoch);
    sources.check_epoch(b.epoch);
    unc_components /= b.value;
    unc_components -= b.unc_components * (value / (b.value * b.value));
    value /= b.value;
    return *this;
  }

  UDoubleCT &operator/=(double b) {
    unc_components /= b;
    value /= b;
    return *this;
  }

  friend UDoubleCT operator/(UDoubleCT a, const UDoubleCT &b) { return a /= b; }

  friend UDoubleCT operator/(UDoubleCT a, double b) { return a /= b; }

  friend UDoubleCT operator/(const double a, const UDoubleCT &b) {
    UDoubleCT retval(0.0);

    retval.unc_components = b.unc_components * (-a / (b.value * b.value));
    retval.value = a / b.value;
    return retval;
  }

  friend std::ostream &operator<<(std::ostream &os, const UDoubleCT &ud) {
    uncertain_print(ud.mean(), ud.deviation(), os);
    return os;
  }

  friend std::istream &operator>>(std::istream &is, UDoubleCT &ud) {
    double mean, sigma;
    std::string source_name;

    uncertain_read(mean, sigma, is);
    std::stringstream os;
    os << "input: ";
    uncertain_print(mean, sigma, os);
    source_name = os.str();
    ud = UDoubleCT<T>(mean, sigma, source_name);
    return is;
  }

  static UDoubleCT func1(std::function<one_arg_ret(double)> func_w_moments, UDoubleCT arg) {
    one_arg_ret funcret = func_w_moments(arg.value);
    arg.value = funcret.value;
    arg.unc_components *= funcret.arg.slope;
    return arg;
  }

  static UDoubleCT func2(std::function<two_arg_ret(double, double)> func_w_moments,
                         const UDoubleCT &arg1, const UDoubleCT &arg2) {
    UDoubleCT<T>::sources.check_epoch(arg1.epoch);
    UDoubleCT<T>::sources.check_epoch(arg2.epoch);
    UDoubleCT retval(arg1);
    two_arg_ret funcret = func_w_moments(arg1.value, arg2.value);
    retval.value = funcret.value;
    retval.unc_components =
        arg1.unc_components * funcret.arg1.slope + arg2.unc_components * funcret.arg2.slope;
    return retval;
  }

  friend UDoubleCT sqrt(UDoubleCT arg) { return func1(&sqrt_w_moments, arg); }

  friend UDoubleCT sin(UDoubleCT arg) { return func1(&sin_w_moments, arg); }

  friend UDoubleCT cos(UDoubleCT arg) { return func1(&cos_w_moments, arg); }

  friend UDoubleCT tan(UDoubleCT arg) { return func1(&tan_w_moments, arg); }

  friend UDoubleCT asin(UDoubleCT arg) { return func1(&asin_w_moments, arg); }

  friend UDoubleCT acos(UDoubleCT arg) { return func1(&acos_w_moments, arg); }

  friend UDoubleCT atan(UDoubleCT arg) { return func1(&atan_w_moments, arg); }

  friend UDoubleCT ceil(UDoubleCT arg) { return func1(&ceil_w_moments, arg); }

  friend UDoubleCT floor(UDoubleCT arg) { return func1(&floor_w_moments, arg); }

  friend UDoubleCT fabs(UDoubleCT arg) { return func1(&fabs_w_moments, arg); }

  friend UDoubleCT exp(UDoubleCT arg) { return func1(&exp_w_moments, arg); }

  friend UDoubleCT log(UDoubleCT arg) { return func1(&log_w_moments, arg); }

  friend UDoubleCT log10(UDoubleCT arg) { return func1(&log10_w_moments, arg); }

  friend UDoubleCT sinh(UDoubleCT arg) { return func1(&sinh_w_moments, arg); }

  friend UDoubleCT cosh(UDoubleCT arg) { return func1(&cosh_w_moments, arg); }

  friend UDoubleCT tanh(UDoubleCT arg) { return func1(&tanh_w_moments, arg); }

  friend UDoubleCT fmod(const UDoubleCT &arg1, const UDoubleCT &arg2) {
    return func2(&fmod_w_moments, arg1, arg2);
  }

  friend UDoubleCT atan2(const UDoubleCT &arg1, const UDoubleCT &arg2) {
    return func2(&atan2_w_moments, arg1, arg2);
  }

  friend UDoubleCT pow(const UDoubleCT &arg1, const UDoubleCT &arg2) {
    return func2(&pow_w_moments, arg1, arg2);
  }

  friend UDoubleCT ldexp(UDoubleCT arg, const int intarg) {
    one_arg_ret funcret = ldexp_w_moments(arg.value, intarg);
    arg.value = funcret.value;
    arg.unc_components *= funcret.arg.slope;
    return arg;
  }

  friend UDoubleCT frexp(UDoubleCT arg, int *intarg) {
    one_arg_ret funcret = frexp_w_moments(arg.value, *intarg);
    arg.value = funcret.value;
    arg.unc_components *= funcret.arg.slope;
    return arg;
  }

  friend UDoubleCT modf(UDoubleCT arg, double *dblarg) {
    one_arg_ret funcret = modf_w_moments(arg.value, *dblarg);
    arg.value = funcret.value;
    arg.unc_components *= funcret.arg.slope;
    return arg;
  }
};

}  // namespace uncertain
