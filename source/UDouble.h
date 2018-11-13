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

#include "functions.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <time.h>

// \todo add exceptions

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

  static double discontinuity_thresh;  // Warn whenever discontinuity is
  // closer than discontinuity_thresh
  // sigmas from value

 public:
  static void set_disc_thresh(const double& new_thresh)
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
      : value(ud.value), uncertainty(ud.uncertainty) {}

  ~UDoubleMSC() {}

  UDoubleMSC<is_correlated> operator+() const { return *this; }

  UDoubleMSC<is_correlated> operator-() const
  {
    if (is_correlated)
      return UDoubleMSC<is_correlated>(-value, -uncertainty);
    else
      return UDoubleMSC<is_correlated>(-value, uncertainty);
  }

  friend UDoubleMSC<is_correlated> operator+(UDoubleMSC<is_correlated> a,
                                             const UDoubleMSC<is_correlated>& b) { return a += b; }

  friend UDoubleMSC<is_correlated> operator-(UDoubleMSC<is_correlated> a,
                                             const UDoubleMSC<is_correlated>& b) { return a -= b; }

  UDoubleMSC<is_correlated> operator++() { return (*this += 1.0); }

  UDoubleMSC<is_correlated> operator--() { return (*this -= 1.0); }

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
                                             const UDoubleMSC<is_correlated>& b) { return a *= b; }

  friend UDoubleMSC<is_correlated> operator/(UDoubleMSC<is_correlated> a,
                                             const UDoubleMSC<is_correlated>& b) { return a /= b; }

  UDoubleMSC<is_correlated>& operator+=(const UDoubleMSC<is_correlated>& ud)
  {
    if (is_correlated)
      uncertainty += ud.uncertainty;
    else
      uncertainty = hypot(uncertainty, ud.uncertainty);
    value += ud.value;
    return *this;
  }

  UDoubleMSC<is_correlated>& operator-=(const UDoubleMSC<is_correlated>& ud)
  {
    if (is_correlated)
      uncertainty -= ud.uncertainty;
    else
      uncertainty = hypot(uncertainty, ud.uncertainty);
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
        uncertainty *= ud.uncertainty * 1.4142135623;
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

#define UDoubleMSCfunc1(func) \
      UDoubleMSC<is_correlated> func(UDoubleMSC<is_correlated> arg) \
      { \
         std::stringstream os; \
         os << #func "(" << arg << ") "; \
         one_arg_ret funcret = func ## _w_moments(arg.value); \
         arg.value = funcret.value + sqr(arg.uncertainty) \
                      * funcret.arg.curve / 2.0; \
         gauss_loss(arg.uncertainty, funcret.arg.disc_dist, \
                    funcret.arg.disc_type, "", os.str(), \
                    UDoubleMSC<is_correlated>::discontinuity_thresh); \
         if (funcret.arg.slope == 0.0) \
            arg.uncertainty *= arg.uncertainty * funcret.arg.curve \
                               * 0.7071067781; \
         else \
         { \
            arg.uncertainty *= funcret.arg.slope \
                            * sqrt(1.0 + 0.5 * sqr(funcret.arg.curve \
                                                   * arg.uncertainty \
                                                   / funcret.arg.slope)); \
            if (!is_correlated) \
               arg.uncertainty = fabs(arg.uncertainty); \
         } \
         return arg; \
      }
// \todo enhance the below to account for moments with terms of
// each argument (e.g. for atan2(x,y), we now ignore d/dx(d/dy(atan2(x,y)))
#define UDoubleMSCfunc2(func) \
      UDoubleMSC<is_correlated> func(const UDoubleMSC<is_correlated>& arg1, \
                                    const UDoubleMSC<is_correlated>& arg2) \
      { \
         UDoubleMSC<is_correlated> retval; \
         double unc1, unc2; \
         std::stringstream os; \
         os << #func "(" << arg1 << ", " << arg2 << ") "; \
         std::string str = os.str(); \
         two_arg_ret funcret = func ## _w_moments(arg1.value, arg2.value); \
         retval.value = funcret.value \
                        + 0.5 * (funcret.arg1.curve * sqr(arg1.uncertainty) \
                                 + funcret.arg2.curve * sqr(arg2.uncertainty));\
         gauss_loss(arg1.uncertainty, funcret.arg1.disc_dist, \
                    funcret.arg1.disc_type, " on 1st argument", \
                    str, UDoubleMSC<is_correlated>::discontinuity_thresh); \
         gauss_loss(arg2.uncertainty, funcret.arg2.disc_dist, \
                    funcret.arg2.disc_type, " on 2nd argument", \
                    str, UDoubleMSC<is_correlated>::discontinuity_thresh); \
         if (funcret.arg1.slope == 0.0) \
            unc1 = sqr(arg1.uncertainty) * funcret.arg1.curve * 0.7071067781; \
         else \
            unc1 = arg1.uncertainty * funcret.arg1.slope \
                   * sqrt(1.0 + 0.5 * sqr(funcret.arg1.curve * arg1.uncertainty \
                                          / funcret.arg1.slope)); \
         if (funcret.arg2.slope == 0.0) \
            unc2 = sqr(arg2.uncertainty) * funcret.arg2.curve * 0.7071067781; \
         else \
            unc2 = arg2.uncertainty * funcret.arg2.slope \
                   * sqrt(1.0 + 0.5 * sqr(funcret.arg2.curve * arg2.uncertainty \
                                          / funcret.arg2.slope)); \
         if (is_correlated) \
            retval.uncertainty = unc1 + unc2; \
         else \
            retval.uncertainty = hypot(unc1, unc2); \
         return retval; \
      }

  friend UDoubleMSCfunc1(sqrt)

  friend UDoubleMSCfunc1(sin)

  friend UDoubleMSCfunc1(cos)

  friend UDoubleMSCfunc1(tan)

  friend UDoubleMSCfunc1(asin)

  friend UDoubleMSCfunc1(acos)

  friend UDoubleMSCfunc1(atan)

  friend UDoubleMSCfunc1(ceil)

  friend UDoubleMSCfunc1(floor)

  friend UDoubleMSCfunc1(fabs)

  friend UDoubleMSCfunc2(fmod)

  friend UDoubleMSCfunc2(atan2)

  friend UDoubleMSCfunc1(exp)

  friend UDoubleMSCfunc1(log)

  friend UDoubleMSCfunc1(log10)

  friend UDoubleMSCfunc1(sinh)

  friend UDoubleMSCfunc1(cosh)

  friend UDoubleMSCfunc1(tanh)

  friend UDoubleMSCfunc2(pow)

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
      arg.uncertainty *= arg.uncertainty * funcret.arg.curve * 0.7071067781;
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
      arg.uncertainty *= arg.uncertainty * funcret.arg.curve * 0.7071067781;
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
      arg.uncertainty *= arg.uncertainty * funcret.arg.curve * 0.7071067781;
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

  // read-only access to data members
  double mean() const { return value; }

  double deviation() const { return fabs(uncertainty); }

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
      retval.uncertainty = sqr(arg.uncertainty) * curve * 0.7071067781;
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
      unc1 = sqr(arg1.uncertainty) * curve1 * 0.7071067781;
    else
      unc1 = arg1.uncertainty * slope1 * sqrt(1.0 + 0.5 * sqr(curve1
                                                                  * arg1.uncertainty / slope1));
    if (slope2 == 0.0)
      unc2 = sqr(arg2.uncertainty) * curve2 * 0.7071067781;
    else
      unc2 = arg2.uncertainty * slope2 * sqrt(1.0 + 0.5 * sqr(curve2
                                                                  * arg2.uncertainty / slope2));
    if (is_correlated)
      retval.uncertainty = unc1 + unc2;
    else
      retval.uncertainty = hypot(unc1, unc2);
    return retval;
  }
// \todo add Monte-Carlo propogation
};

template<>
double UDoubleMSC<false>::discontinuity_thresh = 3.0;

template<>
double UDoubleMSC<true>::discontinuity_thresh = 0.0;


// \todo make template variables or class consts
#define MAX_UNC_ELEMENTS 5

// A class for sources of uncertainties.
class UncertainSourceSet
{
 private:
  unsigned long num_sources;
  unsigned long source_epoch;
  std::string source_name[MAX_UNC_ELEMENTS];
  std::string class_name;
 public:
  UncertainSourceSet(std::string cname = "") : num_sources(0), source_epoch(0)
  {
    for (unsigned int i = 0; i < MAX_UNC_ELEMENTS; i++)
      source_name[i][0] = 0;
    class_name = cname;
  }

  unsigned long get_epoch() const { return source_epoch; }

  void check_epoch(const unsigned long& epoch) const
  {
    if (epoch != source_epoch)
    {
      throw std::runtime_error("Bad epoch: "
                                   + std::to_string(epoch) + " expected: " + std::to_string(source_epoch)
                                   + " in class " + class_name);
    }
  }

  void new_epoch()
  {
    for (unsigned int i = 0; i < num_sources; i++)
      source_name[i][0] = 0;
    source_epoch++;
    num_sources = 0;
  }

  int can_get_new_source() const
  {
    if (num_sources >= MAX_UNC_ELEMENTS)
    {
      std::stringstream ss;
      ss << "Already at maximum number of permissible uncertainty elements: "
         << MAX_UNC_ELEMENTS << "(" << num_sources << ")" << "in class "
         << class_name << ".Change the value of MAX_UNC_ELEMENTS and recompile"
         << " or use new_epoch().";
      throw std::runtime_error(ss.str());
    }
    return 1;
  }

  unsigned long get_new_source(std::string name)
  {
    // std::cerr << "trying to get new source (" << num_sources << ") for "
    //      << name << std::endl;
    (void) can_get_new_source();
    source_name[num_sources] = name;
    return num_sources++;
  }

  unsigned long get_num_sources() const
  {
    return num_sources;
  }

  std::string get_source_name(const unsigned long i) const
  {
    if (i >= num_sources)
    {
      throw std::runtime_error("get_source_name called with illegal source number: "
                                   + std::to_string(i));
    }
    return source_name[i];
  }
};

// Specialized array class has only those members needed to be an array
// of uncertainty elements used as the template parameter in UDoubleCT<>.
// This is the simplest possible implementation.
template<int size>
class SimpleArray
{
 private:
  double element[size];

 public:
  SimpleArray(double initval = 0.0)
  {
    for (int i = 0; i < size; i++)
      element[i] = initval;
  }

  SimpleArray(const SimpleArray& a)
  {
    for (int i = 0; i < size; i++)
      element[i] = a.element[i];
  }

  ~SimpleArray() {}

  SimpleArray operator-() const
  {
    SimpleArray retval;

    for (int i = 0; i < size; i++)
      retval.element[i] = -element[i];
    return retval;
  }

  SimpleArray& operator+=(const SimpleArray& b)
  {
    for (int i = 0; i < size; i++)
      element[i] += b.element[i];
    return *this;
  }

  friend SimpleArray operator+(SimpleArray a, const SimpleArray& b) { return a += b; }

  SimpleArray& operator-=(const SimpleArray& b)
  {
    for (int i = 0; i < size; i++)
      element[i] -= b.element[i];
    return *this;
  }

  SimpleArray& operator*=(const double& b)
  {
    for (int i = 0; i < size; i++)
      element[i] *= b;
    return *this;
  }

  friend SimpleArray operator*(SimpleArray a, const double& b) { return a *= b; }

  SimpleArray& operator/=(const double& b)
  {
    for (int i = 0; i < size; i++)
      element[i] /= b;
    return *this;
  }

  double operator[](int subscript)
  {
    if (subscript < 0)
    {
      throw std::runtime_error("Error: negative subscript: "
                                   + std::to_string(subscript));
    }
    if (subscript >= size)
    {
      throw std::runtime_error("Error: oversize subscript: "
                                   + std::to_string(subscript)
                                   + " Greater than " + std::to_string(size - 1));
    }
    return element[subscript];
  }

  void setelement(int subscript, double value)
  {
    if (subscript < 0)
    {
      throw std::runtime_error("Error: negative subscript: "
                                   + std::to_string(subscript));
    }
    if (subscript >= size)
    {
      throw std::runtime_error("Error: oversize subscript: "
                                   + std::to_string(subscript) +
          " Greater than " + std::to_string(size - 1));
    }
    element[subscript] = value;
  }

  double norm() const
  {
    double tot = 0.0;
    for (int i = 0; i < size; i++)
      tot += element[i] * element[i];
    return sqrt(tot);
  }
};

// Specialized array class has only those members needed to be an array
// of uncertainty elements used as the template parameter in UDoubleCT<>.
// This class is like SimpleArray but adds an undistributed factor for
// greater efficiency in the common case when all members of an array
// are to be multiplied by some factor.
template<int size>
class ArrayWithScale
{
 private:
  double element[size];
  double scale;

 public:
  ArrayWithScale(double initval = 0.0)
  {
    for (int i = 0; i < size; i++)
      element[i] = initval;
    scale = 1.0;
  }

  ArrayWithScale(const ArrayWithScale& a)
  {
    for (int i = 0; i < size; i++)
      element[i] = a.element[i];
    scale = a.scale;
  }

  ~ArrayWithScale() {}

  ArrayWithScale operator-() const
  {
    ArrayWithScale retval;

    for (int i = 0; i < size; i++)
      retval.element[i] = element[i];
    retval.scale = -scale;
    return retval;
  }

  ArrayWithScale& operator+=(const ArrayWithScale& b)
  {
    if (scale)
    {
      double scale_factor = b.scale / scale;
      for (int i = 0; i < size; i++)
        element[i] += b.element[i] * scale_factor;
    }
    else
    {
      scale = b.scale;
      for (int i = 0; i < size; i++)
        element[i] = b.element[i];
    }
    return *this;
  }

  friend ArrayWithScale operator+(ArrayWithScale a, const ArrayWithScale& b) { return a += b; }

  ArrayWithScale& operator-=(const ArrayWithScale& b)
  {
    if (scale)
    {
      double scale_factor = b.scale / scale;
      for (int i = 0; i < size; i++)
        element[i] -= b.element[i] * scale_factor;
    }
    else
    {
      scale = -b.scale;
      for (int i = 0; i < size; i++)
        element[i] = b.element[i];
    }
    return *this;
  }

  ArrayWithScale& operator*=(const double& b)
  {
    scale *= b;
    return *this;
  }

  friend ArrayWithScale operator*(ArrayWithScale a, const double& b) { return a *= b; }

  ArrayWithScale& operator/=(const double& b)
  {
    scale /= b;
    return *this;
  }

  double operator[](int subscript)
  {
    if (subscript < 0)
    {
      throw std::runtime_error("Error: negative subscript: " + std::to_string(subscript));
    }
    if (subscript >= size)
    {
      throw std::runtime_error("Error: oversize subscript: " + std::to_string(subscript) +
                " Greater than " + std::to_string(size - 1));
    }
    return element[subscript] * scale;
  }

  void setelement(int subscript, double value)
  {
    if (subscript < 0)
    {
      throw std::runtime_error("Error: negative subscript: " + std::to_string(subscript));
    }
    if (subscript >= size)
    {
      throw std::runtime_error("Error: oversize subscript: " + std::to_string(subscript) +
                " Greater than " + std::to_string(size - 1));
    }
    if (scale != 0.0)
      element[subscript] = value / scale;
  }

  double norm() const
  {
    if (scale == 0.0)
      return 0.0;
    double tot = 0.0;
    for (int i = 0; i < size; i++)
      tot += element[i] * element[i];
    return sqrt(tot * scale * scale);
  }
};

// Correlation tracking class keeps an array of uncertainty
// components from various sources.  (Array implementation
// is specified by template parameter.)
template<class T>
class UDoubleCT
{
 private:
  double value;
  T unc_components;
  unsigned long epoch;
  static UncertainSourceSet sources;

 public:
  // default constructor creates a new independent uncertainty element
  UDoubleCT(const double val = 0.0, const double unc = 0.0,
            std::string name = "")
      : value(val), unc_components(0.0)
  {
    epoch = sources.get_epoch();
    if (unc < 0.0)
    {
      throw std::runtime_error("Error: negative uncertainty: " + std::to_string(unc));
    }
    if (unc != 0.0)
    {
      std::string source_name;
      if (!name.empty())
      {
        source_name = name;
      }
      else
      {
        std::stringstream os;
        os << "anon: ";
        uncertain_print(val, unc, os);
        source_name = os.str();
      }
      unsigned long new_source_num = sources.get_new_source(source_name);
      unc_components.setelement(new_source_num, unc);
    }
  }

  // copy constructor does not create a new independent uncertainty element
  UDoubleCT(const UDoubleCT& ud) : value(ud.value), epoch(ud.epoch)
  {
    unc_components = ud.unc_components;
  }
  // operator= ?

  ~UDoubleCT() {}

  static void new_epoch() { sources.new_epoch(); }

  void print_uncertain_sources(std::ostream& os = std::cout)
  {
    double total_uncertainty = this->deviation();
    if (total_uncertainty == 0.0)
      os << "No uncertainty";
    else
      for (unsigned int i = 0; i < sources.get_num_sources(); i++)
      {
        double unc_portion = this->unc_components[i] / total_uncertainty;
        unc_portion *= unc_portion;
        os << sources.get_source_name(i) << ": "
           << int_percent(unc_portion) << "%" << std::endl;
      }
    os << std::endl;
  }

  UDoubleCT operator+() const { return *this; }

  UDoubleCT operator-() const
  {
    UDoubleCT retval;

    retval.value = -value;
    retval.unc_components = -unc_components;
    return retval;
  }

  UDoubleCT& operator+=(const UDoubleCT& b)
  {
    sources.check_epoch(epoch);
    sources.check_epoch(b.epoch);
    unc_components += b.unc_components;
    value += b.value;
    return *this;
  }

  UDoubleCT& operator+=(const double& b)
  {
    value += b;
    return *this;
  }

  friend UDoubleCT operator+(UDoubleCT a, const UDoubleCT& b) { return a += b; }

  friend UDoubleCT operator+(UDoubleCT a, const double& b) { return a += b; }

  friend UDoubleCT operator+(const double& b, UDoubleCT a) { return a += b; }

  UDoubleCT& operator-=(const UDoubleCT& b)
  {
    sources.check_epoch(epoch);
    sources.check_epoch(b.epoch);
    unc_components -= b.unc_components;
    value -= b.value;
    return *this;
  }

  UDoubleCT& operator-=(const double& b)
  {
    value -= b;
    return *this;
  }

  friend UDoubleCT operator-(UDoubleCT a, const UDoubleCT& b) { return a -= b; }

  friend UDoubleCT operator-(UDoubleCT a, const double& b) { return a -= b; }

  friend UDoubleCT operator-(const double& b, UDoubleCT a)
  {
    a -= b;
    return -a;
  }

  UDoubleCT& operator*=(const UDoubleCT& b)
  {
    sources.check_epoch(epoch);
    sources.check_epoch(b.epoch);
    unc_components *= b.value;
    unc_components += b.unc_components * value;
    value *= b.value;
    return *this;
  }

  UDoubleCT& operator*=(const double& b)
  {
    unc_components *= b;
    value *= b;
    return *this;
  }

  UDoubleCT operator++() { return (*this += 1.0); }

  UDoubleCT operator--() { return (*this -= 1.0); }

  UDoubleCT operator++(int)
  {
    UDoubleCT retval(*this);
    *this += 1.0;
    return retval;
  }

  UDoubleCT operator--(int)
  {
    UDoubleCT retval(*this);
    *this -= 1.0;
    return retval;
  }

  friend UDoubleCT operator*(UDoubleCT a, const UDoubleCT& b) { return a *= b; }

  friend UDoubleCT operator*(UDoubleCT a, const double& b) { return a *= b; }

  friend UDoubleCT operator*(const double& b, UDoubleCT a) { return a *= b; }

  UDoubleCT& operator/=(const UDoubleCT& b)
  {
    sources.check_epoch(epoch);
    sources.check_epoch(b.epoch);
    unc_components /= b.value;
    unc_components -= b.unc_components * (value / (b.value * b.value));
    value /= b.value;
    return *this;
  }

  UDoubleCT& operator/=(const double& b)
  {
    unc_components /= b;
    value /= b;
    return *this;
  }

  friend UDoubleCT operator/(UDoubleCT a, const UDoubleCT& b) { return a /= b; }

  friend UDoubleCT operator/(UDoubleCT a, const double& b) { return a /= b; }

  friend UDoubleCT operator/(const double a, const UDoubleCT& b)
  {
    UDoubleCT retval(0.0);

    retval.unc_components = b.unc_components * (-a / (b.value * b.value));
    retval.value = a / b.value;
    return retval;
  }

  friend std::ostream& operator<<(std::ostream& os, const UDoubleCT& ud)
  {
    uncertain_print(ud.mean(), ud.deviation(), os);
    return os;
  }

  friend std::istream& operator>>(std::istream& is, UDoubleCT& ud)
  {
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

#define UDoubleCTfunc1(func) \
   UDoubleCT func(UDoubleCT arg) \
   { \
      one_arg_ret funcret = func ## _w_moments(arg.value); \
      arg.value = funcret.value; \
      arg.unc_components *= funcret.arg.slope; \
      return arg; \
   }
#define UDoubleCTfunc2(func) \
   UDoubleCT func(const UDoubleCT& arg1, const UDoubleCT& arg2) \
   { \
      UDoubleCT<T>::sources.check_epoch(arg1.epoch); \
      UDoubleCT<T>::sources.check_epoch(arg2.epoch); \
      UDoubleCT retval(arg1); \
      two_arg_ret funcret = func ## _w_moments(arg1.value, arg2.value); \
      retval.value = funcret.value; \
      retval.unc_components = arg1.unc_components * funcret.arg1.slope \
                              + arg2.unc_components * funcret.arg2.slope; \
      return retval; \
   }

  friend UDoubleCTfunc1(sqrt)

  friend UDoubleCTfunc1(sin)

  friend UDoubleCTfunc1(cos)

  friend UDoubleCTfunc1(tan)

  friend UDoubleCTfunc1(asin)

  friend UDoubleCTfunc1(acos)

  friend UDoubleCTfunc1(atan)

  friend UDoubleCTfunc2(atan2)

  friend UDoubleCTfunc1(ceil)

  friend UDoubleCTfunc1(floor)

  friend UDoubleCTfunc1(fabs)

  friend UDoubleCTfunc2(fmod)

  friend UDoubleCTfunc1(exp)

  friend UDoubleCTfunc1(log)

  friend UDoubleCTfunc1(log10)

  friend UDoubleCTfunc1(sinh)

  friend UDoubleCTfunc1(cosh)

  friend UDoubleCTfunc1(tanh)

  friend UDoubleCTfunc2(pow)

  friend UDoubleCT ldexp(UDoubleCT arg, const int intarg)
  {
    one_arg_ret funcret = ldexp_w_moments(arg.value, intarg);
    arg.value = funcret.value;
    arg.unc_components *= funcret.arg.slope;
    return arg;
  }

  friend UDoubleCT frexp(UDoubleCT arg, int* intarg)
  {
    one_arg_ret funcret = frexp_w_moments(arg.value, *intarg);
    arg.value = funcret.value;
    arg.unc_components *= funcret.arg.slope;
    return arg;
  }

  friend UDoubleCT modf(UDoubleCT arg, double* dblarg)
  {
    one_arg_ret funcret = modf_w_moments(arg.value, *dblarg);
    arg.value = funcret.value;
    arg.unc_components *= funcret.arg.slope;
    return arg;
  }

  double mean() const { return value; }

  double deviation() const { return unc_components.norm(); }

};

typedef SimpleArray<MAX_UNC_ELEMENTS> SizedSimpleArray;
typedef ArrayWithScale<MAX_UNC_ELEMENTS> SizedArrayWithScale;

template<>
UncertainSourceSet UDoubleCT<SizedSimpleArray>::sources("Simple Array");

template<>
UncertainSourceSet UDoubleCT<SizedArrayWithScale>::sources("Array with Scale");

typedef UDoubleCT<SizedSimpleArray> UDoubleCTSA;
typedef UDoubleCT<SizedArrayWithScale> UDoubleCTAA;
