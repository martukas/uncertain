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

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace uncertain {

#ifndef M_PI
constexpr double kPi = 3.14159265358979323846;
constexpr double kHalfPi = 1.57079632679489661923;
constexpr double kLog10e = 0.434294481903251827651;
constexpr double kSqrt2 = 1.41421356237309504880;
constexpr double k1Sqrt2 = 0.707106781186547524401;
#else
// \todo could be different if you use (IEEE 128-bit) long double
constexpr double kPi = M_PI;
constexpr double kHalfPi = M_PI_2;
constexpr double kLog10e = M_LOG10E;
constexpr double kSqrt2 = M_SQRT2;
constexpr double k1Sqrt2 = M_SQRT1_2;
#endif

// square: just a notational convenience
inline double sqr(double a) { return a * a; }

// \todo may be unneeded since C++17
// This function takes the square root of the sum of the squares of
// numbers, which is equal to the length of the hypotenuse of a right
// triangle if the two arguments are the lengths of the legs.
inline double hypot(double a, double b, double c) { return std::sqrt(sqr(a) + sqr(b) + sqr(c)); }

// \todo this may be no longer true, since C++11
// std::hypot() could have problems with overflow.  A possible
// reformulation of the 2-argument version is:
// inline double std::hypot(const double &a, const double &b)
// {
//   double fa = fabs(a);
//   double fb = fabs(b);
//   if (fb > fa)
//   {
//     fa /= fb;
//     return fb * sqrt(1.0 + fa * fa);
//   }
//   else if (fb == fa) // making this a special case makes std::hypot(Inf,Inf)
//   work
//   {
//     return fb;
//   }
//   else
//   {
//     fb /= fa;
//     return fa * sqrt(1.0 + fb * fb);
//   }
// }
// but this would be less likely to be inlined and so might slow the
// UDoubleMS class down unacceptably.

// This function translates floating point ratios to percentages
inline int int_percent(double in) { return int(std::floor(in * 100.0 + 0.5)); }

enum class discontinuity_type {
  none,                 // eg sin(x) is continuous everywhere
  step,                 // eg floor(x) x-> 1.0
  infinite_wrap,        // eg 1/x x-> 0.0
  infinite_then_undef,  // eg log(x) x-> 0.0
  slope_only,           // eg fabs(x) x-> 0.0
  undefined_beyond      // eg asin(x) x -> 1.0
};

// prints uncertainty to 2 digits and value to same precision
void uncertain_print(double mean, double sigma, std::ostream &os = std::cout);

// reads uncertainty as mean +/- sigma
void uncertain_read(double &mean, double &sigma, std::istream &is = std::cin);

// \todo include skewing of distribution
void gauss_loss(double uncertainty, double disc_dist, const discontinuity_type &disc_type,
                std::string id_string, std::string func_str, double disc_thresh);

// This object tells all about the effects of an argument on a function
// return value.
struct arg_effect {
  double slope;
  double curve;
  double disc_dist;
  discontinuity_type disc_type;
};

// These structs are enhanced library returns for functions of one
// and two arguments:
struct one_arg_ret {
  double value;
  arg_effect arg;
};

struct two_arg_ret {
  double value;
  arg_effect arg1;
  arg_effect arg2;
};

// The *_w_moments() functions correspond to math function in
// the standard C library.  But these versions return the slope, curve,
// and distance to the nearest discontinuity as well as the value.
inline one_arg_ret ceil_w_moments(double arg) {
  one_arg_ret retval;
  retval.arg.slope = 0.0;
  retval.arg.curve = 0.0;
  retval.arg.disc_type = discontinuity_type::step;  // discontinuity at each integer
  double temp;
  retval.arg.disc_dist = std::modf(arg, &temp);
  if (retval.arg.disc_dist > 0.5) retval.arg.disc_dist = 1.0 - retval.arg.disc_dist;
  retval.value = std::ceil(arg);
  return retval;
}

inline one_arg_ret floor_w_moments(double arg) {
  one_arg_ret retval;
  retval.arg.slope = 0.0;
  retval.arg.curve = 0.0;
  retval.arg.disc_type = discontinuity_type::step;  // discontinuity at each integer
  double temp;
  retval.arg.disc_dist = std::modf(arg, &temp);
  if (retval.arg.disc_dist > 0.5) retval.arg.disc_dist = 1.0 - retval.arg.disc_dist;
  retval.value = std::floor(arg);
  return retval;
}

inline one_arg_ret fabs_w_moments(double arg) {
  one_arg_ret retval;
  if (arg > 0.0)
    retval.arg.slope = 1.0;
  else
    retval.arg.slope = -1.0;
  retval.value = std::fabs(arg);
  retval.arg.curve = 0.0;
  retval.arg.disc_type = discontinuity_type::slope_only;
  retval.arg.disc_dist = retval.value;  // discontinuity at 0.0
  return retval;
}

inline one_arg_ret ldexp_w_moments(double arg, const int intarg) {
  one_arg_ret retval;
  retval.value = std::ldexp(arg, intarg);
  retval.arg.slope = std::ldexp(1.0, intarg);
  retval.arg.curve = 0.0;
  retval.arg.disc_type = discontinuity_type::none;
  return retval;
}

inline one_arg_ret modf_w_moments(double arg, double &intpart) {
  one_arg_ret retval;
  retval.value = std::modf(arg, &intpart);
  retval.arg.slope = 1.0;
  retval.arg.curve = 0.0;
  retval.arg.disc_type = discontinuity_type::step;
  retval.arg.disc_dist = std::fabs(retval.value);
  if (retval.arg.disc_dist > 0.5) retval.arg.disc_dist = 1.0 - retval.arg.disc_dist;
  return retval;
}

inline one_arg_ret frexp_w_moments(double arg, int &intexp) {
  one_arg_ret retval;
  retval.value = std::frexp(arg, &intexp);
  retval.arg.slope = std::pow(2.0, double(-intexp));
  retval.arg.curve = 0.0;
  retval.arg.disc_type = discontinuity_type::step;
  double disc_loc = std::pow(2.0, double(intexp));
  retval.arg.disc_dist = std::fabs(disc_loc - std::fabs(arg));
  double alt_dist = std::fabs(0.5 * disc_loc - std::fabs(arg));
  if (retval.arg.disc_dist > alt_dist) retval.arg.disc_dist = alt_dist;
  return retval;
}

inline two_arg_ret fmod_w_moments(double arg1, double arg2) {
  two_arg_ret retval;
  retval.value = std::fmod(arg1, arg2);
  retval.arg1.slope = 1.0 / arg2;
  if ((arg1 / arg2) > 0.0)
    retval.arg2.slope = -std::floor(arg1 / arg2);
  else
    retval.arg2.slope = std::floor(-arg1 / arg2);
  retval.arg1.curve = 0.0;
  retval.arg2.curve = 0.0;
  retval.arg1.disc_type = discontinuity_type::step;
  retval.arg1.disc_dist = std::fabs(retval.value);
  if (retval.arg1.disc_dist > std::fabs(arg2) * 0.5)
    retval.arg1.disc_dist = std::fabs(arg2) - retval.arg1.disc_dist;
  retval.arg2.disc_type = discontinuity_type::step;
  double rat = std::fabs(arg1 / arg2);
  double a2floortarg = std::fabs(arg1) / std::floor(rat);
  double a2ceiltarg = std::fabs(arg1) / std::ceil(rat);
  if (std::fabs(a2floortarg - std::fabs(arg2)) < std::fabs(a2ceiltarg - std::fabs(arg2)))
    retval.arg2.disc_dist = std::fabs(a2floortarg - std::fabs(arg2));
  else
    retval.arg2.disc_dist = std::fabs(a2ceiltarg - std::fabs(arg2));
  return retval;
}

inline one_arg_ret sqrt_w_moments(double arg) {
  one_arg_ret retval;
  retval.value = std::sqrt(arg);
  retval.arg.slope = 1.0 / (2.0 * retval.value);
  retval.arg.curve = -0.25 / (retval.value * retval.value * retval.value);
  retval.arg.disc_dist = arg;  // discontinuity at 0.0
  retval.arg.disc_type = discontinuity_type::undefined_beyond;
  return retval;
}

inline one_arg_ret sin_w_moments(double arg) {
  one_arg_ret retval;
  retval.arg.disc_type = discontinuity_type::none;
  retval.arg.slope = std::cos(arg);
  retval.arg.curve = -std::sin(arg);
  retval.value = -retval.arg.curve;
  return retval;
}

inline one_arg_ret cos_w_moments(double arg) {
  one_arg_ret retval;
  retval.arg.disc_type = discontinuity_type::none;
  retval.arg.slope = -std::sin(arg);
  retval.arg.curve = -std::cos(arg);
  retval.value = -retval.arg.curve;
  return retval;
}

inline one_arg_ret tan_w_moments(double arg) {
  one_arg_ret retval;
  double costemp = std::cos(arg);
  retval.arg.slope = 1.0 / (costemp * costemp);
  retval.arg.curve = -2.0 * retval.arg.slope / costemp;
  retval.arg.disc_type = discontinuity_type::infinite_wrap;
  retval.arg.disc_dist = std::fmod(arg - kHalfPi, kPi);
  if (retval.arg.disc_dist > kHalfPi) retval.arg.disc_dist = kPi - retval.arg.disc_dist;
  retval.value = std::tan(arg);
  return retval;
}

inline one_arg_ret asin_w_moments(double arg) {
  one_arg_ret retval;
  retval.arg.slope = 1.0 / std::sqrt(1.0 - arg * arg);
  retval.arg.curve = arg * retval.arg.slope * retval.arg.slope * retval.arg.slope;
  retval.arg.disc_type = discontinuity_type::undefined_beyond;
  if (arg > 0.0)
    retval.arg.disc_dist = 1.0 - arg;
  else
    retval.arg.disc_dist = arg + 1.0;
  retval.value = std::asin(arg);
  return retval;
}

inline one_arg_ret acos_w_moments(double arg) {
  one_arg_ret retval;
  retval.arg.slope = -1.0 / std::sqrt(1.0 - arg * arg);
  retval.arg.curve = arg * retval.arg.slope * retval.arg.slope * retval.arg.slope;
  retval.arg.disc_type = discontinuity_type::undefined_beyond;
  if (arg > 0.0)
    retval.arg.disc_dist = 1.0 - arg;
  else
    retval.arg.disc_dist = arg + 1.0;
  retval.value = std::acos(arg);
  return retval;
}

inline one_arg_ret atan_w_moments(double arg) {
  one_arg_ret retval;
  retval.arg.slope = 1.0 / (1.0 + arg * arg);
  retval.arg.curve = 2.0 * arg * retval.arg.slope * retval.arg.slope;
  retval.arg.disc_type = discontinuity_type::none;
  retval.value = std::atan(arg);
  return retval;
}

inline two_arg_ret atan2_w_moments(double arg1, double arg2) {
  two_arg_ret retval;
  double sum2 = arg2 * arg2 + arg1 * arg1;

  if (sum2 == 0.0) {
    retval.arg1.slope = 1.0;
    retval.arg2.slope = 1.0;
    retval.arg1.curve = 0.0;
    retval.arg2.curve = 0.0;
    retval.arg1.disc_type = discontinuity_type::none;
    retval.arg2.disc_type = discontinuity_type::step;
    retval.arg2.disc_dist = 0.0;
  } else {
    retval.arg1.slope = arg2 / sum2;
    retval.arg2.slope = -arg1 / sum2;
    retval.arg1.curve = -2.0 * arg1 * arg2 / (sum2 * sum2);
    retval.arg2.curve = -retval.arg1.curve;
    if (arg1 == 0.0) {
      retval.arg2.disc_type = discontinuity_type::step;
      retval.arg2.disc_dist = std::fabs(arg2);
    } else
      retval.arg2.disc_type = discontinuity_type::none;

    if (arg2 >= 0.0)
      retval.arg1.disc_type = discontinuity_type::none;
    else {
      retval.arg1.disc_type = discontinuity_type::step;
      retval.arg1.disc_dist = std::fabs(arg1);
    }
  }
  retval.value = std::atan2(arg1, arg2);
  return retval;
}

inline one_arg_ret exp_w_moments(double arg) {
  one_arg_ret retval;
  retval.value = std::exp(arg);
  retval.arg.slope = retval.value;
  retval.arg.curve = retval.arg.slope;
  retval.arg.disc_type = discontinuity_type::none;
  return retval;
}

inline one_arg_ret log_w_moments(double arg) {
  one_arg_ret retval;
  retval.arg.slope = 1.0 / arg;
  retval.arg.curve = -retval.arg.slope * retval.arg.slope;
  retval.arg.disc_dist = arg;
  retval.arg.disc_type = discontinuity_type::undefined_beyond;
  retval.value = std::log(arg);
  return retval;
}

inline one_arg_ret log10_w_moments(double arg) {
  one_arg_ret retval;
  retval.arg.slope = kLog10e / arg;
  retval.arg.curve = -retval.arg.slope / arg;
  retval.arg.disc_dist = arg;
  retval.arg.disc_type = discontinuity_type::undefined_beyond;
  retval.value = std::log10(arg);
  return retval;
}

inline one_arg_ret sinh_w_moments(double arg) {
  one_arg_ret retval;
  retval.arg.slope = std::cosh(arg);
  retval.value = std::sinh(arg);
  retval.arg.curve = retval.value;
  retval.arg.disc_type = discontinuity_type::none;
  return retval;
}

inline one_arg_ret cosh_w_moments(double arg) {
  one_arg_ret retval;
  retval.arg.slope = std::sinh(arg);
  retval.value = std::cosh(arg);
  retval.arg.curve = retval.value;
  retval.arg.disc_type = discontinuity_type::none;
  return retval;
}

inline one_arg_ret tanh_w_moments(double arg) {
  one_arg_ret retval;
  double coshtemp = std::cosh(arg);
  retval.arg.slope = 1.0 / (coshtemp * coshtemp);
  retval.arg.curve = -2.0 * std::sinh(arg) / (coshtemp * coshtemp * coshtemp);
  retval.arg.disc_type = discontinuity_type::none;
  retval.value = std::tanh(arg);
  return retval;
}

inline two_arg_ret pow_w_moments(double arg1, double arg2) {
  two_arg_ret retval;
  retval.value = std::pow(arg1, arg2);
  if (arg1 == 0.0) {
    retval.arg1.slope = 0.0;
    retval.arg1.curve = 0.0;
    if (arg2 == 1.0)
      retval.arg1.slope = 1.0;
    else if (arg2 == 2.0)
      retval.arg1.curve = 2.0;
    retval.arg1.disc_type = discontinuity_type::none;  // pow(arg1, 0.0) is 1.0 for any arg1
    // pow(0.0, arg2) is undefined for every negative arg2
    retval.arg2.disc_type = discontinuity_type::undefined_beyond;
    // pow(0.0, arg2) is 0.0 for every positive arg2, undef for neg arg2
    retval.arg2.slope = 0.0;
    retval.arg2.curve = 0.0;
    if (arg2 > 0.0)
      retval.arg2.disc_dist = arg2;
    else
      retval.arg2.disc_dist = 0.0;
  } else if (arg1 < 0.0) {
    // pow(arg1, arg2) for arg1 < 0.0 is only defined for integer arg2
    retval.arg1.slope = arg2 * retval.value / arg1;
    retval.arg1.curve = arg2 * (arg2 - 1.0) * retval.value / (arg1 * arg1);
    retval.arg1.disc_type = discontinuity_type::none;
    retval.arg2.slope = 0.0;
    retval.arg2.curve = 0.0;
    retval.arg2.disc_type = discontinuity_type::undefined_beyond;
    retval.arg2.disc_dist = 0.0;
  } else {
    retval.arg1.slope = arg2 * retval.value / arg1;
    retval.arg1.curve = arg2 * (arg2 - 1.0) * retval.value / (arg1 * arg1);
    if (arg2 == std::floor(arg2))  // if arg2 is an integer
    {
      retval.arg1.disc_type = discontinuity_type::none;
    } else {
      retval.arg1.disc_type = discontinuity_type::undefined_beyond;
      retval.arg1.disc_dist = arg1;
    }
    retval.arg2.slope = std::log(arg1) * retval.value;
    retval.arg2.curve = std::log(arg1) * retval.arg2.slope;
    retval.arg2.disc_type = discontinuity_type::none;
  }
  return retval;
}

// This function returns an approximation of the inverse Gaussian denstity
// function to within 4.5e-4.  From Abromowitz & Stegun's _Handbook_of_
// _Mathematical_Functions_ formula 26.2.23
inline double inverse_gaussian_density(double p) {
  if (p <= 0.0) {
    throw std::runtime_error("inverse_gaussian_density() called for negative value: " +
                             std::to_string(p));
  } else if (p > 0.5) {
    throw std::runtime_error("inverse_gaussian_density() called for too large value: " +
                             std::to_string(p));
  }

  const double c0 = 2.515517, c1 = 0.802853, c2 = 0.010328;
  const double d1 = 1.432788, d2 = 0.189269, d3 = 0.001308;
  double t = std::sqrt(std::log(1.0 / (p * p)));

  return t - (c0 + t * (c1 + c2 * t)) / (1.0 + t * (d1 + t * (d2 + t * d3)));
}

// Compare function to sort by absolute value
inline int abs_double_compare(const void *a, const void *b) {
  double fa = std::fabs(*(double *)a);
  double fb = std::fabs(*(double *)b);
  if (fa == fb)
    return 0;
  else if (fa < fb)
    return -1;
  else
    return 1;
}

}  // namespace uncertain
