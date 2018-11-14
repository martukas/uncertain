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

// Warning: this file contains an object (UDoubleInit) to insure that
// srand() gets called exactly once to seed the random number
// generator.  But this may cause problems if used with source
// code that already calls srand() and/or uses that family of
// random number generators in any files that don't #include this file.

// Warning: This header has not yet been corrected to work with more
// than one source file.


#pragma once

#include <uncertain/functions.hpp>
#include <sstream>
#include <functional>

namespace uncertain
{

// model uncertain number using only mean and sigma, like UDoubleMS,
// but also use some knowledge of curve & skew & discontinuities.
// (The "C" is for "Curve").  For simple cases this gives better
// answers than the simple UDoubleMS classes, but when multiple
// functions with significant curve are applied after each other
// it starts giving much worse answers.
//
// Ultimately this work will be more useful when it is incorporated
// into a correlation tracking class.
template<bool is_correlated>
class UDoubleMSC
{
 private:
  double value;
  double uncertainty;

  // Warn whenever discontinuity is
  // closer than discontinuity_thresh
  // sigmas from value
  static double discontinuity_thresh;

 public:
  static void set_disc_thresh(double new_thresh)
  {
    discontinuity_thresh = new_thresh;
  }

  // This is the default conversion from type double
  UDoubleMSC(const double val = 0.0, const double unc = 0.0)
      : value(val), uncertainty(unc)
  {
    if ((unc < 0.0) && !is_correlated)
    {
      throw std::runtime_error("Error: negative uncertainty: " + std::to_string(unc));
    }
  }

  UDoubleMSC(const UDoubleMSC& ud)
      : value(ud.value), uncertainty(ud.uncertainty)
  {}

  ~UDoubleMSC() = default;

  // read-only access to data members
  double mean() const
  { return value; }

  double deviation() const
  { return fabs(uncertainty); }

  UDoubleMSC<is_correlated> operator+() const
  { return *this; }

  UDoubleMSC<is_correlated> operator-() const
  {
    if (is_correlated)
      return UDoubleMSC<is_correlated>(-value, -uncertainty);
    else
      return UDoubleMSC<is_correlated>(-value, uncertainty);
  }

  friend UDoubleMSC<is_correlated> operator+(UDoubleMSC<is_correlated> a,
                                             const UDoubleMSC<is_correlated>& b)
  { return a += b; }

  friend UDoubleMSC<is_correlated> operator-(UDoubleMSC<is_correlated> a,
                                             const UDoubleMSC<is_correlated>& b)
  { return a -= b; }

  UDoubleMSC<is_correlated> operator++()
  { return (*this += 1.0); }

  UDoubleMSC<is_correlated> operator--()
  { return (*this -= 1.0); }

  UDoubleMSC<is_correlated> operator++(int)
  {
    UDoubleMSC<is_correlated> retval(*this);
    *this += 1.0;
    return retval;
  }

  UDoubleMSC<is_correlated> operator--(int)
  {
    UDoubleMSC<is_correlated> retval(*this);
    *this -= 1.0;
    return retval;
  }

  friend UDoubleMSC<is_correlated> operator*(UDoubleMSC<is_correlated> a,
                                             const UDoubleMSC<is_correlated>& b)
  { return a *= b; }

  friend UDoubleMSC<is_correlated> operator/(UDoubleMSC<is_correlated> a,
                                             const UDoubleMSC<is_correlated>& b)
  { return a /= b; }

  UDoubleMSC<is_correlated>& operator+=(const UDoubleMSC<is_correlated>& ud)
  {
    if (is_correlated)
      uncertainty += ud.uncertainty;
    else
      uncertainty = std::hypot(uncertainty, ud.uncertainty);
    value += ud.value;
    return *this;
  }

  UDoubleMSC<is_correlated>& operator-=(const UDoubleMSC<is_correlated>& ud)
  {
    if (is_correlated)
      uncertainty -= ud.uncertainty;
    else
      uncertainty = std::hypot(uncertainty, ud.uncertainty);
    value -= ud.value;
    return *this;
  }

  UDoubleMSC<is_correlated>& operator*=(const UDoubleMSC<is_correlated>& ud)
  {
    if (is_correlated)
    {
      double second_order_correlated_adjust = uncertainty * ud.uncertainty;
      double unctemp = uncertainty * ud.value + ud.uncertainty * value;
      if (fabs(unctemp) <= fabs(uncertainty * ud.uncertainty) * 1e-10)
        uncertainty *= ud.uncertainty * kSqrt2;
      else
        uncertainty = unctemp * sqrt(1.0 + 2.0 * sqr(uncertainty
                                                         * ud.uncertainty
                                                         / unctemp));

      value *= ud.value;
      value += second_order_correlated_adjust;
    }
    else
    {
      uncertainty = hypot(uncertainty * ud.value, ud.uncertainty * value,
                          uncertainty * ud.uncertainty);
      value *= ud.value;
    }
    return *this;
  }

  UDoubleMSC<is_correlated>& operator/=(const UDoubleMSC<is_correlated>& ud)
  {
    double second_order_correlated_adjust;
    if (is_correlated)
    {
      second_order_correlated_adjust = (ud.uncertainty / (ud.value * ud.value))
          * (value * ud.uncertainty / ud.value - uncertainty);
      uncertainty = (uncertainty / ud.value
          - (ud.uncertainty * value) / (ud.value * ud.value))
          * sqrt(1.0 + 2.0 * sqr(ud.uncertainty / ud.value));

      if (ud.uncertainty != 0.0)
      {
        double disc_dist = fabs(ud.value / ud.uncertainty);
        if (disc_dist < discontinuity_thresh)
        {
          std::cerr << "correlated division by " << ud << " is ";
          std::cerr << disc_dist
                    << " sigmas from an infinite wrap discontinuity" << std::endl;
        }
      }
    }
    else
    {
      second_order_correlated_adjust = ud.uncertainty * ud.uncertainty
          * value
          / (ud.value * ud.value * ud.value);
      double inverted_sigma = -(ud.uncertainty / sqr(ud.value))
          * sqrt(1.0 + 2.0 * sqr(ud.uncertainty
                                     / ud.value));
      uncertainty = hypot(uncertainty / ud.value, inverted_sigma * value,
                          uncertainty * inverted_sigma);
    }
    value /= ud.value;
    value += second_order_correlated_adjust;
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& os, const UDoubleMSC<is_correlated>& ud)
  {
    uncertain_print(ud.mean(), ud.deviation(), os);
    return os;
  }

  friend std::istream& operator>>(std::istream& is, UDoubleMSC<is_correlated>& ud)
  {
    double mean, sigma;
    uncertain_read(mean, sigma, is);
    ud = UDoubleMSC<is_correlated>(mean, sigma);
    return is;
  }

  static UDoubleMSC<is_correlated> func1(std::function<one_arg_ret(double)> func_w_moments,
                                         UDoubleMSC<is_correlated> arg,
                                         std::string funcname)
  {
    std::stringstream os;
    os << funcname << "(" << arg << ") ";
    one_arg_ret funcret = func_w_moments(arg.value);
    arg.value = funcret.value + sqr(arg.uncertainty)
        * funcret.arg.curve / 2.0;
    gauss_loss(arg.uncertainty, funcret.arg.disc_dist,
               funcret.arg.disc_type, "", os.str(),
               UDoubleMSC<is_correlated>::discontinuity_thresh);
    if (funcret.arg.slope == 0.0)
      arg.uncertainty *= arg.uncertainty * funcret.arg.curve
          * k1Sqrt2;
    else
    {
      arg.uncertainty *= funcret.arg.slope
          * sqrt(1.0 + 0.5 * sqr(funcret.arg.curve
                                     * arg.uncertainty
                                     / funcret.arg.slope));
      if (!is_correlated)
        arg.uncertainty = fabs(arg.uncertainty);
    }
    return arg;
  }

  // \todo enhance the below to account for moments with terms of each
  //  argument (e.g. for atan2(x,y), we now ignore d/dx(d/dy(atan2(x,y)))
  static UDoubleMSC<is_correlated> func2(std::function<two_arg_ret(double, double)> func_w_moments,
                                         UDoubleMSC<is_correlated> arg1,
                                         UDoubleMSC<is_correlated> arg2,
                                         std::string funcname)
  {
    UDoubleMSC<is_correlated> retval;
    double unc1, unc2;
    std::stringstream os;
    os << funcname << "(" << arg1 << ", " << arg2 << ") ";
    std::string str = os.str();
    two_arg_ret funcret = func_w_moments(arg1.value, arg2.value);
    retval.value = funcret.value
        + 0.5 * (funcret.arg1.curve * sqr(arg1.uncertainty)
            + funcret.arg2.curve * sqr(arg2.uncertainty));\
         gauss_loss(arg1.uncertainty, funcret.arg1.disc_dist,
                    funcret.arg1.disc_type, " on 1st argument",
                    str, UDoubleMSC<is_correlated>::discontinuity_thresh);
    gauss_loss(arg2.uncertainty, funcret.arg2.disc_dist,
               funcret.arg2.disc_type, " on 2nd argument",
               str, UDoubleMSC<is_correlated>::discontinuity_thresh);
    if (funcret.arg1.slope == 0.0)
      unc1 = sqr(arg1.uncertainty) * funcret.arg1.curve * k1Sqrt2;
    else
      unc1 = arg1.uncertainty * funcret.arg1.slope
          * sqrt(1.0 + 0.5 * sqr(funcret.arg1.curve * arg1.uncertainty
                                     / funcret.arg1.slope));
    if (funcret.arg2.slope == 0.0)
      unc2 = sqr(arg2.uncertainty) * funcret.arg2.curve * k1Sqrt2;
    else
      unc2 = arg2.uncertainty * funcret.arg2.slope
          * sqrt(1.0 + 0.5 * sqr(funcret.arg2.curve * arg2.uncertainty
                                     / funcret.arg2.slope));
    if (is_correlated)
      retval.uncertainty = unc1 + unc2;
    else
      retval.uncertainty = std::hypot(unc1, unc2);
    return retval;
  }

  friend UDoubleMSC<is_correlated> sqrt(UDoubleMSC<is_correlated> arg)
  {
    return func1(&sqrt_w_moments, arg, "sqrt");
  }

  friend UDoubleMSC<is_correlated> sin(UDoubleMSC<is_correlated> arg)
  {
    return func1(&sin_w_moments, arg, "sin");
  }

  friend UDoubleMSC<is_correlated> cos(UDoubleMSC<is_correlated> arg)
  {
    return func1(&cos_w_moments, arg, "cos");
  }

  friend UDoubleMSC<is_correlated> tan(UDoubleMSC<is_correlated> arg)
  {
    return func1(&tan_w_moments, arg, "tan");
  }

  friend UDoubleMSC<is_correlated> asin(UDoubleMSC<is_correlated> arg)
  {
    return func1(&asin_w_moments, arg, "asin");
  }

  friend UDoubleMSC<is_correlated> acos(UDoubleMSC<is_correlated> arg)
  {
    return func1(&acos_w_moments, arg, "acos");
  }

  friend UDoubleMSC<is_correlated> atan(UDoubleMSC<is_correlated> arg)
  {
    return func1(&atan_w_moments, arg, "atan");
  }

  friend UDoubleMSC<is_correlated> ceil(UDoubleMSC<is_correlated> arg)
  {
    return func1(&ceil_w_moments, arg, "ceil");
  }

  friend UDoubleMSC<is_correlated> floor(UDoubleMSC<is_correlated> arg)
  {
    return func1(&floor_w_moments, arg, "floor");
  }

  friend UDoubleMSC<is_correlated> fabs(UDoubleMSC<is_correlated> arg)
  {
    return func1(&fabs_w_moments, arg, "fabs");
  }

  friend UDoubleMSC<is_correlated> exp(UDoubleMSC<is_correlated> arg)
  {
    return func1(&exp_w_moments, arg, "exp");
  }

  friend UDoubleMSC<is_correlated> log(UDoubleMSC<is_correlated> arg)
  {
    return func1(&log_w_moments, arg, "log");
  }

  friend UDoubleMSC<is_correlated> log10(UDoubleMSC<is_correlated> arg)
  {
    return func1(&log10_w_moments, arg, "sqrtlog10");
  }

  friend UDoubleMSC<is_correlated> sinh(UDoubleMSC<is_correlated> arg)
  {
    return func1(&sinh_w_moments, arg, "sinh");
  }

  friend UDoubleMSC<is_correlated> cosh(UDoubleMSC<is_correlated> arg)
  {
    return func1(&cosh_w_moments, arg, "cosh");
  }

  friend UDoubleMSC<is_correlated> tanh(UDoubleMSC<is_correlated> arg)
  {
    return func1(&tanh_w_moments, arg, "tanh");
  }

  friend UDoubleMSC<is_correlated> fmod(UDoubleMSC<is_correlated> arg1, UDoubleMSC<is_correlated> arg2)
  {
    return func2(&fmod_w_moments, arg1, arg2, "fmod");
  }

  friend UDoubleMSC<is_correlated> atan2(UDoubleMSC<is_correlated> arg1, UDoubleMSC<is_correlated> arg2)
  {
    return func2(&atan2_w_moments, arg1, arg2, "atan2");
  }

  friend UDoubleMSC<is_correlated> pow(UDoubleMSC<is_correlated> arg1, UDoubleMSC<is_correlated> arg2)
  {
    return func2(&pow_w_moments, arg1, arg2, "pow");
  }

  friend UDoubleMSC<is_correlated> ldexp(UDoubleMSC<is_correlated> arg,
                                         const int intarg)
  {
    std::stringstream os;
    os << "ldexp(" << arg << ", " << intarg << ") ";
    one_arg_ret funcret = ldexp_w_moments(arg.value, intarg);
    arg.value = funcret.value + sqr(arg.uncertainty) * funcret.arg.curve / 2.0;
    gauss_loss(arg.uncertainty, funcret.arg.disc_dist, funcret.arg.disc_type,
               "", os.str(), UDoubleMSC<is_correlated>::discontinuity_thresh);
    if (funcret.arg.slope == 0.0)
      arg.uncertainty *= arg.uncertainty * funcret.arg.curve * k1Sqrt2;
    else
    {
      arg.uncertainty *= funcret.arg.slope
          * sqrt(1.0 + 0.5 * sqr(funcret.arg.curve
                                     * arg.uncertainty
                                     / funcret.arg.slope));
      if (!is_correlated)
        arg.uncertainty = fabs(arg.uncertainty);
    }
    return arg;
  }

  friend UDoubleMSC<is_correlated> frexp(UDoubleMSC<is_correlated> arg,
                                         int* intarg)
  {
    std::stringstream os;
    os << "frexp(" << arg << ", " << *intarg << ") ";
    one_arg_ret funcret = frexp_w_moments(arg.value, *intarg);
    arg.value = funcret.value + sqr(arg.uncertainty) * funcret.arg.curve / 2.0;
    gauss_loss(arg.uncertainty, funcret.arg.disc_dist, funcret.arg.disc_type,
               "", os.str(), UDoubleMSC<is_correlated>::discontinuity_thresh);
    if (funcret.arg.slope == 0.0)
      arg.uncertainty *= arg.uncertainty * funcret.arg.curve * k1Sqrt2;
    else
    {
      arg.uncertainty *= funcret.arg.slope
          * sqrt(1.0 + 0.5 * sqr(funcret.arg.curve
                                     * arg.uncertainty
                                     / funcret.arg.slope));
      if (!is_correlated)
        arg.uncertainty = fabs(arg.uncertainty);
    }
    return arg;
  }

  friend UDoubleMSC<is_correlated> modf(UDoubleMSC<is_correlated> arg,
                                        double* dblarg)
  {
    std::stringstream os;
    os << "modf(" << arg << ", " << *dblarg << ") ";
    one_arg_ret funcret = modf_w_moments(arg.value, *dblarg);
    arg.value = funcret.value + sqr(arg.uncertainty) * funcret.arg.curve / 2.0;
    gauss_loss(arg.uncertainty, funcret.arg.disc_dist, funcret.arg.disc_type,
               "", os.str(), UDoubleMSC<is_correlated>::discontinuity_thresh);
    if (funcret.arg.slope == 0.0)
      arg.uncertainty *= arg.uncertainty * funcret.arg.curve * k1Sqrt2;
    else
    {
      arg.uncertainty *= funcret.arg.slope
          * sqrt(1.0 + 0.5 * sqr(funcret.arg.curve
                                     * arg.uncertainty
                                     / funcret.arg.slope));
      if (!is_correlated)
        arg.uncertainty = fabs(arg.uncertainty);
    }
    return arg;
  }

  friend UDoubleMSC<is_correlated> PropagateUncertaintiesBySlope(
      double (* certain_func)(double),
      const UDoubleMSC<is_correlated>& arg)
  {
    UDoubleMSC<is_correlated> retval;
    double core_value, slope, curve;
    double sigma_up_value, sigma_down_value;

    core_value = certain_func(arg.value);
    sigma_up_value = certain_func(arg.value + arg.uncertainty);
    sigma_down_value = certain_func(arg.value - arg.uncertainty);
    if (arg.uncertainty)
    {
      slope = (sigma_up_value - sigma_down_value) * 0.5 / arg.uncertainty;
      curve = (sigma_up_value + sigma_down_value - 2.0 * core_value)
          / sqr(arg.uncertainty);
    }
    else
      slope = curve = 0.0;

    retval.value = core_value + sqr(arg.uncertainty) * curve / 2.0;
    if (slope == 0.0)
      retval.uncertainty = sqr(arg.uncertainty) * curve * k1Sqrt2;
    else
    {
      retval.uncertainty = arg.uncertainty * slope
          * sqrt(1.0 + 0.5 * sqr(curve * arg.uncertainty
                                     / slope));
      if (!is_correlated)
        retval.uncertainty = fabs(retval.uncertainty);
    }
    return retval;
  }

  friend UDoubleMSC<is_correlated> PropagateUncertaintiesBySlope(
      double (* certain_func)(double, double),
      const UDoubleMSC<is_correlated>& arg1,
      const UDoubleMSC<is_correlated>& arg2)
  {
    UDoubleMSC<is_correlated> retval;
    double core_value, slope1, curve1, slope2, curve2;
    double up1, down1, up2, down2, unc1, unc2;

    core_value = certain_func(arg1.value, arg2.value);
    up1 = certain_func(arg1.value + arg1.uncertainty, arg2.value);
    down1 = certain_func(arg1.value - arg1.uncertainty, arg2.value);
    up2 = certain_func(arg1.value, arg2.value + arg2.uncertainty);
    down2 = certain_func(arg1.value, arg2.value - arg2.uncertainty);

    if (arg1.uncertainty)
    {
      slope1 = (up1 - down1) * 0.5 / arg1.uncertainty;
      curve1 = (up1 + down1 - 2.0 * core_value) / sqr(arg1.uncertainty);
    }
    else
      slope1 = curve1 = 0.0;

    if (arg2.uncertainty)
    {
      slope2 = (up2 - down2) * 0.5 / arg2.uncertainty;
      curve2 = (up2 + down2 - 2.0 * core_value) / sqr(arg2.uncertainty);
    }
    else
      slope2 = curve2 = 0.0;

    retval.value = core_value + 0.5 * (sqr(arg1.uncertainty) * curve1
        + sqr(arg2.uncertainty) * curve2);
    if (slope1 == 0.0)
      unc1 = sqr(arg1.uncertainty) * curve1 * k1Sqrt2;
    else
      unc1 = arg1.uncertainty * slope1 * sqrt(1.0 + 0.5 * sqr(curve1
                                                                  * arg1.uncertainty / slope1));
    if (slope2 == 0.0)
      unc2 = sqr(arg2.uncertainty) * curve2 * k1Sqrt2;
    else
      unc2 = arg2.uncertainty * slope2 * sqrt(1.0 + 0.5 * sqr(curve2
                                                                  * arg2.uncertainty / slope2));
    if (is_correlated)
      retval.uncertainty = unc1 + unc2;
    else
      retval.uncertainty = std::hypot(unc1, unc2);
    return retval;
  }
// \todo add Monte-Carlo propogation
};

// \todo do not hardcode this
template<>
double UDoubleMSC<false>::discontinuity_thresh = 3.0;

template<>
double UDoubleMSC<true>::discontinuity_thresh = 0.0;


// typedefs to hide the use of templates in the implementation
using UDoubleMSCUncorr = UDoubleMSC<false>;
using UDoubleMSCCorr = UDoubleMSC<true>;

}
