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

// udms.h: This file includes a class for simple propagation of Uncertainties
// according to a pure Gaussian model.

// By Evan Manning (manning@alumni.caltech.edu)
// Be warned: this is a particularly simple  model of uncertainties.
// It is designed to accompany an article to be published in the
// C/C++ Users Journal in March 1996.  A fuller collection of classes
// is available at the same ftp site (airs1.jpl.nasa.gov/pub/evan/c++)
// in udouble.h.

#pragma once

#include <uncertain/functions.hpp>

namespace uncertain {

// model uncertain number using only mean and sigma (pure Gaussian)
// This is the simplest possible model of uncertainty.  It ignores
// all second-order and higher-order effects, and when two
// uncertain numbers interact this model assumes either that they
// are 100% correlated or that they are 100% uncorrelated, depending
// on the template parameter.
template <bool is_correlated>
class UDoubleMS {
 private:
  double value;        // the central (expected) value
  double uncertainty;  // the uncertainty (standard deviation)

 public:
  // This is the default conversion from type double
  UDoubleMS(double val = 0.0, double unc = 0.0) : value(val), uncertainty(unc) {
    if ((unc < 0.0) && !is_correlated) {
      throw std::runtime_error("Error: negative uncertainty: " + std::to_string(unc));
    }
  }

  UDoubleMS(const UDoubleMS &ud) : value(ud.value), uncertainty(ud.uncertainty) {}

  ~UDoubleMS() = default;

  // read-only access to data members
  double mean() const { return value; }

  // \todo should this be fabs even when correlated?
  double deviation() const { return std::fabs(uncertainty); }

  UDoubleMS<is_correlated> operator+() const { return *this; }

  UDoubleMS<is_correlated> operator-() const {
    if (is_correlated)
      return UDoubleMS<is_correlated>(-value, -uncertainty);
    else
      return UDoubleMS<is_correlated>(-value, uncertainty);
  }

  friend UDoubleMS<is_correlated> operator+(UDoubleMS<is_correlated> a,
                                            const UDoubleMS<is_correlated> &b) {
    return a += b;
  }

  friend UDoubleMS<is_correlated> operator+(UDoubleMS<is_correlated> a, double b) { return a += b; }

  friend UDoubleMS<is_correlated> operator+(double b, UDoubleMS<is_correlated> a) { return a += b; }

  friend UDoubleMS<is_correlated> operator-(UDoubleMS<is_correlated> a,
                                            const UDoubleMS<is_correlated> &b) {
    return a -= b;
  }

  friend UDoubleMS<is_correlated> operator-(UDoubleMS<is_correlated> a, double b) { return a -= b; }

  friend UDoubleMS<is_correlated> operator-(double b, UDoubleMS<is_correlated> a) {
    a -= b;
    return -a;
  }

  UDoubleMS<is_correlated> operator++() { return (*this += 1.0); }

  UDoubleMS<is_correlated> operator--() { return (*this -= 1.0); }

  UDoubleMS<is_correlated> operator++(int) {
    UDoubleMS<is_correlated> retval(*this);
    *this += 1.0;
    return retval;
  }

  UDoubleMS<is_correlated> operator--(int) {
    UDoubleMS<is_correlated> retval(*this);
    *this -= 1.0;
    return retval;
  }

  friend UDoubleMS<is_correlated> operator*(UDoubleMS<is_correlated> a,
                                            const UDoubleMS<is_correlated> &b) {
    return a *= b;
  }

  friend UDoubleMS<is_correlated> operator*(UDoubleMS<is_correlated> a, double b) { return a *= b; }

  friend UDoubleMS<is_correlated> operator*(double b, UDoubleMS<is_correlated> a) { return a *= b; }

  friend UDoubleMS<is_correlated> operator/(UDoubleMS<is_correlated> a,
                                            const UDoubleMS<is_correlated> &b) {
    return a /= b;
  }

  friend UDoubleMS<is_correlated> operator/(UDoubleMS<is_correlated> a, double b) { return a /= b; }

  friend UDoubleMS<is_correlated> operator/(double b, UDoubleMS<is_correlated> a) {
    UDoubleMS<is_correlated> retval;
    retval.uncertainty = -b * a.uncertainty / (a.value * a.value);
    retval.value = b / a.value;
    return retval;
  }

  UDoubleMS<is_correlated> &operator+=(const UDoubleMS<is_correlated> &ud) {
    if (is_correlated)
      uncertainty += ud.uncertainty;
    else
      uncertainty = std::hypot(uncertainty, ud.uncertainty);
    value += ud.value;
    return *this;
  }

  UDoubleMS<is_correlated> &operator+=(double a) {
    value += a;
    return *this;
  }

  UDoubleMS<is_correlated> &operator-=(const UDoubleMS<is_correlated> &ud) {
    if (is_correlated)
      uncertainty -= ud.uncertainty;
    else
      uncertainty = std::hypot(uncertainty, ud.uncertainty);
    value -= ud.value;
    return *this;
  }

  UDoubleMS<is_correlated> &operator-=(double a) {
    value -= a;
    return *this;
  }

  UDoubleMS<is_correlated> &operator*=(const UDoubleMS<is_correlated> &ud) {
    if (is_correlated)
      uncertainty = uncertainty * ud.value + ud.uncertainty * value;
    else
      uncertainty = std::hypot(uncertainty * ud.value, ud.uncertainty * value);
    value *= ud.value;
    return *this;
  }

  UDoubleMS<is_correlated> &operator*=(double a) {
    value *= a;
    uncertainty *= a;
    return *this;
  }

  UDoubleMS<is_correlated> &operator/=(const UDoubleMS<is_correlated> &ud) {
    if (is_correlated)
      uncertainty = uncertainty / ud.value - (ud.uncertainty * value) / (ud.value * ud.value);
    else
      uncertainty =
          std::hypot(uncertainty / ud.value, (ud.uncertainty * value) / (ud.value * ud.value));
    value /= ud.value;
    return *this;
  }

  UDoubleMS<is_correlated> &operator/=(double a) {
    value /= a;
    uncertainty /= a;
    return *this;
  }

  friend std::ostream &operator<<(std::ostream &os, const UDoubleMS<is_correlated> &ud) {
    uncertain_print(ud.mean(), ud.deviation(), os);
    return os;
  }

  friend std::istream &operator>>(std::istream &is, UDoubleMS<is_correlated> &ud) {
    double mean, sigma;
    uncertain_read(mean, sigma, is);
    ud = UDoubleMS<is_correlated>(mean, sigma);
    return is;
  }

  // math library functions.  These functions multiply the uncertainty
  // by the derivative of the function at this value.
  friend UDoubleMS<is_correlated> ceil(UDoubleMS<is_correlated> arg) {
    arg.value = std::ceil(arg.value);
    arg.uncertainty = 0.0;
    return arg;
  }

  friend UDoubleMS<is_correlated> floor(UDoubleMS<is_correlated> arg) {
    arg.value = std::floor(arg.value);
    arg.uncertainty = 0.0;
    return arg;
  }

  friend UDoubleMS<is_correlated> fabs(UDoubleMS<is_correlated> arg) {
    if (is_correlated && (arg.value < 0.0)) arg.uncertainty *= -1.0;
    arg.value = std::fabs(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> ldexp(UDoubleMS<is_correlated> arg, int intarg) {
    if (is_correlated)
      arg.uncertainty *= std::ldexp(1.0, intarg);
    else
      arg.uncertainty *= std::fabs(std::ldexp(1.0, intarg));
    arg.value = std::ldexp(arg.value, intarg);
    return arg;
  }

  friend UDoubleMS<is_correlated> modf(UDoubleMS<is_correlated> arg, double *intpart) {
    arg.value = std::modf(arg.value, intpart);
    return arg;
  }

  friend UDoubleMS<is_correlated> frexp(UDoubleMS<is_correlated> arg, int *intarg) {
    arg.uncertainty *= std::pow(2.0, double(-*intarg));
    arg.value = std::frexp(arg.value, intarg);
    return arg;
  }

  friend UDoubleMS<is_correlated> fmod(const UDoubleMS<is_correlated> &arg1,
                                       const UDoubleMS<is_correlated> &arg2) {
    UDoubleMS<is_correlated> retval;
    double slope1, slope2;

    slope1 = 1.0 / arg2.value;
    if ((arg1.value / arg2.value) > 0.0)
      slope2 = -std::floor(arg1.value / arg2.value);
    else
      slope2 = std::floor(-arg1.value / arg2.value);
    if (is_correlated)
      retval.uncertainty = slope1 * arg1.uncertainty + slope2 * arg2.uncertainty;
    else
      retval.uncertainty = std::hypot(slope1 * arg1.uncertainty, slope2 * arg2.uncertainty);
    retval.value = std::fmod(arg1.value, arg2.value);
    return retval;
  }

  friend UDoubleMS<is_correlated> sqrt(UDoubleMS<is_correlated> arg) {
    arg.value = std::sqrt(arg.value);
    if (is_correlated)
      arg.uncertainty /= 2.0 * arg.value;
    else
      arg.uncertainty /= std::fabs(2.0 * arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> sin(UDoubleMS<is_correlated> arg) {
    if (is_correlated)
      arg.uncertainty *= std::cos(arg.value);
    else
      arg.uncertainty *= std::fabs(std::cos(arg.value));
    arg.value = std::sin(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> cos(UDoubleMS<is_correlated> arg) {
    if (is_correlated)
      arg.uncertainty *= -std::sin(arg.value);
    else
      arg.uncertainty *= std::fabs(std::sin(arg.value));
    arg.value = std::cos(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> tan(UDoubleMS<is_correlated> arg) {
    double costemp = std::cos(arg.value);
    arg.uncertainty /= costemp * costemp;
    arg.value = std::tan(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> asin(UDoubleMS<is_correlated> arg) {
    arg.uncertainty /= std::sqrt(1.0 - arg.value * arg.value);
    arg.value = std::asin(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> acos(UDoubleMS<is_correlated> arg) {
    if (is_correlated)
      arg.uncertainty /= -std::sqrt(1.0 - arg.value * arg.value);
    else
      arg.uncertainty /= std::sqrt(1.0 - arg.value * arg.value);
    arg.value = std::acos(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> atan(UDoubleMS<is_correlated> arg) {
    arg.uncertainty /= 1.0 + arg.value * arg.value;
    arg.value = std::atan(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> atan2(const UDoubleMS<is_correlated> &arg1,
                                        const UDoubleMS<is_correlated> &arg2) {
    UDoubleMS<is_correlated> retval;
    double slope1 = 1.0, slope2 = 1.0;
    double sum2 = arg2.value * arg2.value + arg1.value * arg1.value;

    if (sum2 != 0.0) {
      slope1 = arg2.value / sum2;
      slope2 = -arg1.value / sum2;
    }
    if (is_correlated)
      retval.uncertainty = slope1 * arg1.uncertainty + slope2 * arg2.uncertainty;
    else
      retval.uncertainty = std::hypot(slope1 * arg1.uncertainty, slope2 * arg2.uncertainty);
    retval.value = std::atan2(arg1.value, arg2.value);
    return retval;
  }

  friend UDoubleMS<is_correlated> exp(UDoubleMS<is_correlated> arg) {
    arg.value = std::exp(arg.value);
    if (is_correlated)
      arg.uncertainty *= arg.value;
    else
      arg.uncertainty *= std::fabs(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> log(UDoubleMS<is_correlated> arg) {
    if (is_correlated)
      arg.uncertainty /= arg.value;
    else
      arg.uncertainty /= std::fabs(arg.value);
    arg.value = std::log(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> log10(UDoubleMS<is_correlated> arg) {
    if (is_correlated)
      arg.uncertainty *= kLog10e / arg.value;
    else
      arg.uncertainty *= kLog10e / std::fabs(arg.value);
    arg.value = std::log10(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> sinh(UDoubleMS<is_correlated> arg) {
    arg.uncertainty *= std::cosh(arg.value);
    arg.value = std::sinh(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> cosh(UDoubleMS<is_correlated> arg) {
    if (is_correlated)
      arg.uncertainty *= std::sinh(arg.value);
    else
      arg.uncertainty *= std::fabs(std::sinh(arg.value));
    arg.value = std::cosh(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> tanh(UDoubleMS<is_correlated> arg) {
    double coshtemp = std::cosh(arg.value);
    arg.uncertainty /= coshtemp * coshtemp;
    arg.value = std::tanh(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> pow(const UDoubleMS<is_correlated> &arg1,
                                      const UDoubleMS<is_correlated> &arg2) {
    UDoubleMS<is_correlated> retval;
    double slope1, slope2;

    retval.value = std::pow(arg1.value, arg2.value);
    if (arg1.value == 0.0) {
      slope2 = 0.0;
      slope1 = 0.0;
      if (arg2.value == 1.0) slope1 = 1.0;
    } else if (arg1.value < 0.0) {
      // pow(arg1, arg2) for arg1 < 0.0 is only defined for integer arg2
      slope1 = arg2.value * retval.value / arg1.value;
      slope2 = 0.0;
    } else {
      slope1 = arg2.value * retval.value / arg1.value;
      slope2 = std::log(arg1.value) * retval.value;
    }
    if (is_correlated)
      retval.uncertainty = slope1 * arg1.uncertainty + slope2 * arg2.uncertainty;
    else
      retval.uncertainty = std::hypot(slope1 * arg1.uncertainty, slope2 * arg2.uncertainty);
    return retval;
  }

  // To propogate an uncertainty through a function for which the slope
  // is not known, we estimate the slope by comparing values for
  // f(mean + sigma) and f(mean - sigma).
  friend UDoubleMS<is_correlated> PropagateUncertaintiesBySlope(
      double (*certain_func)(double), const UDoubleMS<is_correlated> &arg) {
    UDoubleMS<is_correlated> retval;
    double sigma_up_value, sigma_down_value;

    retval.value = certain_func(arg.value);
    sigma_up_value = certain_func(arg.value + arg.uncertainty);
    sigma_down_value = certain_func(arg.value - arg.uncertainty);
    retval.uncertainty = (sigma_up_value - sigma_down_value) * 0.5;
    if (!is_correlated) retval.uncertainty = std::fabs(retval.uncertainty);

    return retval;
  }

  friend UDoubleMS<is_correlated> PropagateUncertaintiesBySlope(
      double (*certain_func)(double, double), const UDoubleMS<is_correlated> &arg1,
      const UDoubleMS<is_correlated> &arg2) {
    UDoubleMS<is_correlated> retval;

    retval.value = certain_func(arg1.value, arg2.value);
    if (is_correlated) {
      double up_val = certain_func(arg1.value + arg1.uncertainty, arg2.value + arg2.uncertainty);
      double down_val = certain_func(arg1.value - arg1.uncertainty, arg2.value - arg2.uncertainty);
      retval.uncertainty = 0.5 * (up_val - down_val);
    } else {
      double up_val1 = certain_func(arg1.value + arg1.uncertainty, arg2.value);
      double down_val1 = certain_func(arg1.value - arg1.uncertainty, arg2.value);
      double up_val2 = certain_func(arg1.value, arg2.value + arg2.uncertainty);
      double down_val2 = certain_func(arg1.value, arg2.value - arg2.uncertainty);
      retval.uncertainty = 0.5 * std::hypot(up_val1 - down_val1, up_val2 - down_val2);
    }
    return retval;
  }
};

// typedefs to hide the use of templates in the implementation
using UDoubleMSUncorr = UDoubleMS<false>;
using UDoubleMSCorr = UDoubleMS<true>;

}  // namespace uncertain
