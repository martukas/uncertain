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

#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>

// future direction: restrict namespace
// prints uncertainty to 2 digits and value to same precision
inline void uncertain_print(double mean, double sigma, std::ostream& os = std::cout)
{
  auto original_precision = os.precision();
  auto original_format = os.flags(std::ios::showpoint);

  // std::cerr << "<" << mean << " +/- " << sigma << "> " << std::endl;
  int precision;
  // special cases for zero, NaN, and Infinities (positive & negative)
  if ((sigma == 0.0) || (sigma != sigma) || (1.0 / sigma == 0.0))
  {
    precision = 0;
  }
  else
  {
    // round sigma to 2 digits
    int sigma_digits = 1 - int(floor(log10(fabs(sigma))));
    double round_10_pow = pow(10.0, sigma_digits);
    sigma = floor(sigma * round_10_pow + 0.5) / round_10_pow;

    // round mean to same # of digits
    mean = floor(mean * round_10_pow + 0.5) / round_10_pow;
    if (mean == 0.0)
    {
      if (sigma_digits > 0)
        precision = sigma_digits + 1;
      else
        precision = 1;
    }
    else
    {
      precision = int(floor(log10(fabs(mean)))) + sigma_digits + 1;
      if (precision < 1)
      {
        mean = 0.0;
        if (sigma_digits > 0)
          precision = sigma_digits + 1;
        else
          precision = 1;
      }
    }
  }
  os << std::setprecision(precision)
     << mean << " +/- "
     << std::setprecision(2)
     << sigma
     << std::setprecision(original_precision);
  os.flags(original_format);
}

// future direction: restrict namespace
// reads uncertainty as mean +/- sigma
inline void uncertain_read(double& mean, double& sigma, std::istream& is = std::cin)
{
  char plus, slash, minus;
  is >> mean >> plus >> slash >> minus >> sigma;
  if ((plus != '+') || (slash != '/') || (minus != '-'))
  {
    std::cerr << "Error: illegal characters encountered in reading "
                 "mean +/- sigma" << std::endl;
    exit(EXIT_FAILURE);
  }
}

// model uncertain number using only mean and sigma (pure Gaussian)
// This is the simplest possible model of uncertainty.  It ignores
// all second-order and higher-order effects, and when two
// uncertain numbers interact this model assumes either that they
// are 100% correlated or that they are 100% uncorrelated, depending
// on the template parameter. (The template parameter is conceptual
// a bool, but g++ 2.6.3 chokes on that.)
template<int is_correlated>
class UDoubleMS
{
 private:
  double value;          // the central (expected) value
  double uncertainty;    // the uncertainty (standard deviation)

 public:
  // This is the default conversion from type double
  UDoubleMS(const double val = 0.0, const double unc = 0.0)
      : value(val), uncertainty(unc)
  {
    if ((unc < 0.0) && !is_correlated)
    {
      std::cerr << "Error: negative uncertainty: " << unc << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  UDoubleMS(const UDoubleMS& ud)
      : value(ud.value), uncertainty(ud.uncertainty) {}

  ~UDoubleMS() {}

  UDoubleMS<is_correlated> operator+() const { return *this; }

  UDoubleMS<is_correlated> operator-() const
  {
    if (is_correlated)
      return UDoubleMS<is_correlated>(-value, -uncertainty);
    else
      return UDoubleMS<is_correlated>(-value, uncertainty);
  }

  friend UDoubleMS<is_correlated> operator+(UDoubleMS<is_correlated> a,
                                            const UDoubleMS<is_correlated>& b) { return a += b; }

  friend UDoubleMS<is_correlated> operator+(UDoubleMS<is_correlated> a,
                                            const double& b) { return a += b; }

  friend UDoubleMS<is_correlated> operator+(const double& b,
                                            UDoubleMS<is_correlated> a) { return a += b; }

  friend UDoubleMS<is_correlated> operator-(UDoubleMS<is_correlated> a,
                                            const UDoubleMS<is_correlated>& b) { return a -= b; }

  friend UDoubleMS<is_correlated> operator-(UDoubleMS<is_correlated> a,
                                            const double& b) { return a -= b; }

  friend UDoubleMS<is_correlated> operator-(const double& b,
                                            UDoubleMS<is_correlated> a)
  {
    a -= b;
    return -a;
  }

  UDoubleMS<is_correlated> operator++() { return (*this += 1.0); }

  UDoubleMS<is_correlated> operator--() { return (*this -= 1.0); }

  UDoubleMS<is_correlated> operator++(int)
  {
    UDoubleMS<is_correlated> retval(*this);
    *this += 1.0;
    return retval;
  }

  UDoubleMS<is_correlated> operator--(int)
  {
    UDoubleMS<is_correlated> retval(*this);
    *this -= 1.0;
    return retval;
  }

  friend UDoubleMS<is_correlated> operator*(UDoubleMS<is_correlated> a,
                                            const UDoubleMS<is_correlated>& b) { return a *= b; }

  friend UDoubleMS<is_correlated> operator*(UDoubleMS<is_correlated> a,
                                            const double& b) { return a *= b; }

  friend UDoubleMS<is_correlated> operator*(const double& b,
                                            UDoubleMS<is_correlated> a) { return a *= b; }

  friend UDoubleMS<is_correlated> operator/(UDoubleMS<is_correlated> a,
                                            const UDoubleMS<is_correlated>& b) { return a /= b; }

  friend UDoubleMS<is_correlated> operator/(UDoubleMS<is_correlated> a,
                                            const double& b) { return a /= b; }

  friend UDoubleMS<is_correlated> operator/(const double& b,
                                            UDoubleMS<is_correlated> a)
  {
    UDoubleMS<is_correlated> retval;
    retval.uncertainty = -b * a.uncertainty / (a.value * a.value);
    retval.value = b / a.value;
    return retval;
  }

  UDoubleMS<is_correlated>& operator+=(const UDoubleMS<is_correlated>& ud)
  {
    if (is_correlated)
      uncertainty += ud.uncertainty;
    else
      uncertainty = hypot(uncertainty, ud.uncertainty);
    value += ud.value;
    return *this;
  }

  UDoubleMS<is_correlated>& operator+=(const double& a)
  {
    value += a;
    return *this;
  }

  UDoubleMS<is_correlated>& operator-=(const UDoubleMS<is_correlated>& ud)
  {
    if (is_correlated)
      uncertainty -= ud.uncertainty;
    else
      uncertainty = hypot(uncertainty, ud.uncertainty);
    value -= ud.value;
    return *this;
  }

  UDoubleMS<is_correlated>& operator-=(const double& a)
  {
    value -= a;
    return *this;
  }

  UDoubleMS<is_correlated>& operator*=(const UDoubleMS<is_correlated>& ud)
  {
    if (is_correlated)
      uncertainty = uncertainty * ud.value + ud.uncertainty * value;
    else
      uncertainty = hypot(uncertainty * ud.value, ud.uncertainty * value);
    value *= ud.value;
    return *this;
  }

  UDoubleMS<is_correlated>& operator*=(const double& a)
  {
    value *= a;
    uncertainty *= a;
    return *this;
  }

  UDoubleMS<is_correlated>& operator/=(const UDoubleMS<is_correlated>& ud)
  {
    if (is_correlated)
      uncertainty = uncertainty / ud.value
          - (ud.uncertainty * value) / (ud.value * ud.value);
    else
      uncertainty = hypot(uncertainty / ud.value,
                          (ud.uncertainty * value) / (ud.value * ud.value));
    value /= ud.value;
    return *this;
  }

  UDoubleMS<is_correlated>& operator/=(const double& a)
  {
    value /= a;
    uncertainty /= a;
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& os, const UDoubleMS<is_correlated>& ud)
  {
    uncertain_print(ud.mean(), ud.deviation(), os);
    return os;
  }

  friend std::istream& operator>>(std::istream& is, UDoubleMS<is_correlated>& ud)
  {
    double mean, sigma;
    uncertain_read(mean, sigma, is);
    ud = UDoubleMS<is_correlated>(mean, sigma);
    return is;
  }

  // math library functions.  These functions multiply the uncertainty
  // by the derivative of the function at this value.
  friend UDoubleMS<is_correlated> ceil(UDoubleMS<is_correlated> arg)
  {
    arg.value = ceil(arg.value);
    arg.uncertainty = 0.0;
    return arg;
  }

  friend UDoubleMS<is_correlated> floor(UDoubleMS<is_correlated> arg)
  {
    arg.value = floor(arg.value);
    arg.uncertainty = 0.0;
    return arg;
  }

  friend UDoubleMS<is_correlated> fabs(UDoubleMS<is_correlated> arg)
  {
    if (is_correlated && (arg.value < 0.0))
      arg.uncertainty *= -1.0;
    arg.value = fabs(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> ldexp(UDoubleMS<is_correlated> arg,
                                        int intarg)
  {
    if (is_correlated)
      arg.uncertainty *= ldexp(1.0, intarg);
    else
      arg.uncertainty *= fabs(ldexp(1.0, intarg));
    arg.value = ldexp(arg.value, intarg);
    return arg;
  }

  friend UDoubleMS<is_correlated> modf(UDoubleMS<is_correlated> arg,
                                       double* intpart)
  {
    arg.value = modf(arg.value, intpart);
    return arg;
  }

  friend UDoubleMS<is_correlated> frexp(UDoubleMS<is_correlated> arg,
                                        int* intarg)
  {
    arg.uncertainty *= pow(2.0, double(-*intarg));
    arg.value = frexp(arg.value, intarg);
    return arg;
  }

  friend UDoubleMS<is_correlated> fmod(const UDoubleMS<is_correlated>& arg1,
                                       const UDoubleMS<is_correlated>& arg2)
  {
    UDoubleMS<is_correlated> retval;
    double slope1, slope2;

    slope1 = 1.0 / arg2.value;
    if ((arg1.value / arg2.value) > 0.0)
      slope2 = -floor(arg1.value / arg2.value);
    else
      slope2 = floor(-arg1.value / arg2.value);
    if (is_correlated)
      retval.uncertainty = slope1 * arg1.uncertainty
          + slope2 * arg2.uncertainty;
    else
      retval.uncertainty = hypot(slope1 * arg1.uncertainty,
                                 slope2 * arg2.uncertainty);
    retval.value = fmod(arg1.value, arg2.value);
    return retval;
  }

  friend UDoubleMS<is_correlated> sqrt(UDoubleMS<is_correlated> arg)
  {
    arg.value = sqrt(arg.value);
    if (is_correlated)
      arg.uncertainty /= 2.0 * arg.value;
    else
      arg.uncertainty /= fabs(2.0 * arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> sin(UDoubleMS<is_correlated> arg)
  {
    if (is_correlated)
      arg.uncertainty *= cos(arg.value);
    else
      arg.uncertainty *= fabs(cos(arg.value));
    arg.value = sin(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> cos(UDoubleMS<is_correlated> arg)
  {
    if (is_correlated)
      arg.uncertainty *= -sin(arg.value);
    else
      arg.uncertainty *= fabs(sin(arg.value));
    arg.value = cos(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> tan(UDoubleMS<is_correlated> arg)
  {
    double costemp = cos(arg.value);
    arg.uncertainty /= costemp * costemp;
    arg.value = tan(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> asin(UDoubleMS<is_correlated> arg)
  {
    arg.uncertainty /= sqrt(1.0 - arg.value * arg.value);
    arg.value = asin(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> acos(UDoubleMS<is_correlated> arg)
  {
    if (is_correlated)
      arg.uncertainty /= -sqrt(1.0 - arg.value * arg.value);
    else
      arg.uncertainty /= sqrt(1.0 - arg.value * arg.value);
    arg.value = acos(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> atan(UDoubleMS<is_correlated> arg)
  {
    arg.uncertainty /= 1.0 + arg.value * arg.value;
    arg.value = atan(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> atan2(const UDoubleMS<is_correlated>& arg1,
                                        const UDoubleMS<is_correlated>& arg2)
  {
    UDoubleMS<is_correlated> retval;
    double slope1 = 1.0, slope2 = 1.0;
    double sum2 = arg2.value * arg2.value + arg1.value * arg1.value;

    if (sum2 != 0.0)
    {
      slope1 = arg2.value / sum2;
      slope2 = -arg1.value / sum2;
    }
    if (is_correlated)
      retval.uncertainty = slope1 * arg1.uncertainty
          + slope2 * arg2.uncertainty;
    else
      retval.uncertainty = hypot(slope1 * arg1.uncertainty,
                                 slope2 * arg2.uncertainty);
    retval.value = atan2(arg1.value, arg2.value);
    return retval;
  }

  friend UDoubleMS<is_correlated> exp(UDoubleMS<is_correlated> arg)
  {
    arg.value = exp(arg.value);
    if (is_correlated)
      arg.uncertainty *= arg.value;
    else
      arg.uncertainty *= fabs(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> log(UDoubleMS<is_correlated> arg)
  {
    if (is_correlated)
      arg.uncertainty /= arg.value;
    else
      arg.uncertainty /= fabs(arg.value);
    arg.value = log(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> log10(UDoubleMS<is_correlated> arg)
  {
    if (is_correlated)
      arg.uncertainty *= 0.43429448189 / arg.value;
    else
      arg.uncertainty *= 0.43429448189 / fabs(arg.value);
    arg.value = log10(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> sinh(UDoubleMS<is_correlated> arg)
  {
    arg.uncertainty *= cosh(arg.value);
    arg.value = sinh(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> cosh(UDoubleMS<is_correlated> arg)
  {
    if (is_correlated)
      arg.uncertainty *= sinh(arg.value);
    else
      arg.uncertainty *= fabs(sinh(arg.value));
    arg.value = cosh(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> tanh(UDoubleMS<is_correlated> arg)
  {
    double coshtemp = cosh(arg.value);
    arg.uncertainty /= coshtemp * coshtemp;
    arg.value = tanh(arg.value);
    return arg;
  }

  friend UDoubleMS<is_correlated> pow(const UDoubleMS<is_correlated>& arg1,
                                      const UDoubleMS<is_correlated>& arg2)
  {
    UDoubleMS<is_correlated> retval;
    double slope1, slope2;

    retval.value = pow(arg1.value, arg2.value);
    if (arg1.value == 0.0)
    {
      slope2 = 0.0;
      slope1 = 0.0;
      if (arg2.value == 1.0)
        slope1 = 1.0;
    }
    else if (arg1.value < 0.0)
    {
      // pow(arg1, arg2) for arg1 < 0.0 is only defined for integer arg2
      slope1 = arg2.value * retval.value / arg1.value;
      slope2 = 0.0;
    }
    else
    {
      slope1 = arg2.value * retval.value / arg1.value;
      slope2 = log(arg1.value) * retval.value;
    }
    if (is_correlated)
      retval.uncertainty = slope1 * arg1.uncertainty
          + slope2 * arg2.uncertainty;
    else
      retval.uncertainty = hypot(slope1 * arg1.uncertainty,
                                 slope2 * arg2.uncertainty);
    return retval;
  }

  // read-only access to data members
  double mean() const { return value; }

  double deviation() const { return fabs(uncertainty); }

  // To propogate an uncertainty through a function for which the slope
  // is not known, we estimate the slope by comparing values for
  // f(mean + sigma) and f(mean - sigma).
  friend UDoubleMS<is_correlated> PropagateUncertaintiesBySlope(
      double (* certain_func)(double),
      const UDoubleMS<is_correlated>& arg)
  {
    UDoubleMS<is_correlated> retval;
    double sigma_up_value, sigma_down_value;

    retval.value = certain_func(arg.value);
    sigma_up_value = certain_func(arg.value + arg.uncertainty);
    sigma_down_value = certain_func(arg.value - arg.uncertainty);
    retval.uncertainty = (sigma_up_value - sigma_down_value) * 0.5;
    if (!is_correlated)
      retval.uncertainty = fabs(retval.uncertainty);

    return retval;
  }

  friend UDoubleMS<is_correlated> PropagateUncertaintiesBySlope(
      double (* certain_func)(double, double),
      const UDoubleMS<is_correlated>& arg1,
      const UDoubleMS<is_correlated>& arg2)
  {
    UDoubleMS<is_correlated> retval;

    retval.value = certain_func(arg1.value, arg2.value);
    if (is_correlated)
    {
      double up_val = certain_func(arg1.value + arg1.uncertainty,
                                   arg2.value + arg2.uncertainty);
      double down_val = certain_func(arg1.value - arg1.uncertainty,
                                     arg2.value - arg2.uncertainty);
      retval.uncertainty = 0.5 * (up_val - down_val);
    }
    else
    {
      double up_val1 = certain_func(arg1.value + arg1.uncertainty,
                                    arg2.value);
      double down_val1 = certain_func(arg1.value - arg1.uncertainty,
                                      arg2.value);
      double up_val2 = certain_func(arg1.value,
                                    arg2.value + arg2.uncertainty);
      double down_val2 = certain_func(arg1.value,
                                      arg2.value - arg2.uncertainty);
      retval.uncertainty = 0.5 * hypot(up_val1 - down_val1,
                                       up_val2 - down_val2);
    }
    return retval;
  }
};

// typedefs to hide the use of templates in the implementation
typedef UDoubleMS<0> UDoubleMSUncorr;
typedef UDoubleMS<1> UDoubleMSCorr;

