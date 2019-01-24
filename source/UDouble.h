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

#include "UDoubleMS.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <strstream>
#include <iomanip>
#include <time.h>

#define PI 3.14159265358979323846
#define HALF_PI (PI / 2.0)

// future direction: add exceptions


// future direction: these functions should be static in implementation
//                   file or hidden in a namespace

// This function takes the square root of the sum of the squares of
// numbers, which is equal to the length of the hypotenuse of a right
// triangle if the two arguments are the lengths of the legs.
inline double hypot(const double &a, const double &b, const double &c)
   { return sqrt(a * a + b * b + c * c); }
// hypot() could have problems with overflow.  A possible
// reformulation of the 2-argument version is:
// inline double hypot(const double &a, const double &b)
// {
//   double fa = fabs(a);
//   double fb = fabs(b);
//   if (fb > fa)
//   {
//     fa /= fb;
//     return fb * sqrt(1.0 + fa * fa);
//   }
//   else if (fb == fa) // making this a special case makes hypot(Inf,Inf) work
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

// square: just a notational convenience
inline double sqr(const double &a)
   { return a * a; }


typedef enum
{
    none,                // eg sin(x) is continuous everywhere
    step,                // eg floor(x) x-> 1.0
    infinite_wrap,       // eg 1/x x-> 0.0
    infinite_then_undef, // eg log(x) x-> 0.0
    slope_only,          // eg fabs(x) x-> 0.0
    undefined_beyond     // eg asin(x) x -> 1.0
} discontinuity_type;


// This function translates floating point ratios to percentages
inline int int_percent(const double &in)
{ return int(floor(in * 100.0 + 0.5)); }

// This object tells all about the effects of an argument on a function
// return value.
typedef struct {
   double slope;
   double curve;
   double disc_dist;
   discontinuity_type disc_type;
} arg_effect;

// These structs are enhanced library returns for functions of one
// and two arguments:
typedef struct {
   double value;
   arg_effect arg;
} one_arg_ret;

typedef struct {
   double value;
   arg_effect arg1;
   arg_effect arg2;
} two_arg_ret;

// The *_w_moments() functions correspond to math function in
// the standard C library.  But these versions return the slope, curve,
// and distance to the nearest discontinuity as well as the value.
inline one_arg_ret ceil_w_moments(const double& arg)
{
   one_arg_ret retval;
   retval.arg.slope = 0.0;
   retval.arg.curve = 0.0;
   retval.arg.disc_type = step; // discontinuity at each integer
   double temp;
   retval.arg.disc_dist = modf(arg, &temp);
   if (retval.arg.disc_dist > 0.5)
      retval.arg.disc_dist = 1.0 - retval.arg.disc_dist;
   retval.value = ceil(arg);
   return retval;
}
inline one_arg_ret floor_w_moments(const double& arg)
{
   one_arg_ret retval;
   retval.arg.slope = 0.0;
   retval.arg.curve = 0.0;
   retval.arg.disc_type = step; // discontinuity at each integer
   double temp;
   retval.arg.disc_dist = modf(arg, &temp);
   if (retval.arg.disc_dist > 0.5)
      retval.arg.disc_dist = 1.0 - retval.arg.disc_dist;
   retval.value = floor(arg);
   return retval;
}
inline one_arg_ret fabs_w_moments(const double& arg)
{
   one_arg_ret retval;
   if (arg > 0.0)
      retval.arg.slope = 1.0;
   else
      retval.arg.slope = -1.0;
   retval.value = fabs(arg);
   retval.arg.curve = 0.0;
   retval.arg.disc_type = slope_only;
   retval.arg.disc_dist = retval.value; // discontinuity at 0.0
   return retval;
}
inline one_arg_ret ldexp_w_moments(const double& arg, const int intarg)
{
   one_arg_ret retval;
   retval.value = ldexp(arg, intarg);
   retval.arg.slope = ldexp(1.0, intarg);
   retval.arg.curve = 0.0;
   retval.arg.disc_type = none;
   return retval;
}
inline one_arg_ret modf_w_moments(const double& arg, double& intpart)
{
   one_arg_ret retval;
   retval.value = modf(arg, &intpart);
   retval.arg.slope = 1.0;
   retval.arg.curve = 0.0;
   retval.arg.disc_type = step;
   retval.arg.disc_dist = fabs(retval.value);
   if (retval.arg.disc_dist > 0.5)
      retval.arg.disc_dist = 1.0 - retval.arg.disc_dist;
   return retval;
}
inline one_arg_ret frexp_w_moments(const double& arg, int& intexp)
{
   one_arg_ret retval;
   retval.value = frexp(arg, &intexp);
   retval.arg.slope = pow(2.0, double(-intexp));
   retval.arg.curve = 0.0;
   retval.arg.disc_type = step;
   double disc_loc = pow(2.0, double(intexp));
   retval.arg.disc_dist = fabs(disc_loc - fabs(arg));
   double alt_dist = fabs(0.5 * disc_loc - fabs(arg));
   if (retval.arg.disc_dist > alt_dist)
      retval.arg.disc_dist = alt_dist;
   return retval;
}
inline two_arg_ret fmod_w_moments(const double& arg1, const double& arg2)
{
   two_arg_ret retval;
   retval.value = fmod(arg1, arg2);
   retval.arg1.slope = 1.0 / arg2;
   if ((arg1/arg2) > 0.0)
      retval.arg2.slope = -floor(arg1/arg2);
   else
      retval.arg2.slope = floor(-arg1/arg2);
   retval.arg1.curve = 0.0;
   retval.arg2.curve = 0.0;
   retval.arg1.disc_type = step;
   retval.arg1.disc_dist = fabs(retval.value);
   if (retval.arg1.disc_dist > fabs(arg2) * 0.5)
      retval.arg1.disc_dist = fabs(arg2) - retval.arg1.disc_dist;
   retval.arg2.disc_type = step;
   double rat = fabs(arg1 / arg2);
   double a2floortarg = fabs(arg1) / floor(rat);
   double a2ceiltarg = fabs(arg1) / ceil(rat);
   if (fabs(a2floortarg - fabs(arg2)) < fabs(a2ceiltarg - fabs(arg2)))
      retval.arg2.disc_dist = fabs(a2floortarg - fabs(arg2));
   else
      retval.arg2.disc_dist = fabs(a2ceiltarg - fabs(arg2));
   return retval;
}
inline one_arg_ret sqrt_w_moments(const double& arg)
{
   one_arg_ret retval;
   retval.value = sqrt(arg);
   retval.arg.slope = 1.0 / (2.0 * retval.value);
   retval.arg.curve = -0.25 / (retval.value * retval.value * retval.value);
   retval.arg.disc_dist = arg; // discontinuity at 0.0
   retval.arg.disc_type = undefined_beyond;
   return retval;
}
inline one_arg_ret sin_w_moments(const double& arg)
{
   one_arg_ret retval;
   retval.arg.disc_type = none;
   retval.arg.slope = cos(arg);
   retval.arg.curve = -sin(arg);
   retval.value = -retval.arg.curve;
   return retval;
}
inline one_arg_ret cos_w_moments(const double& arg)
{
   one_arg_ret retval;
   retval.arg.disc_type = none;
   retval.arg.slope = -sin(arg);
   retval.arg.curve = -cos(arg);
   retval.value = -retval.arg.curve;
   return retval;
}
inline one_arg_ret tan_w_moments(const double& arg)
{
   one_arg_ret retval;
   double costemp = cos(arg);
   retval.arg.slope = 1.0 / (costemp * costemp);
   retval.arg.curve = -2.0 * retval.arg.slope / costemp;
   retval.arg.disc_type = infinite_wrap;
   retval.arg.disc_dist = fmod(arg - HALF_PI, PI);
   if (retval.arg.disc_dist > HALF_PI)
      retval.arg.disc_dist = PI - retval.arg.disc_dist;
   retval.value = tan(arg);
   return retval;
}
inline one_arg_ret asin_w_moments(const double& arg)
{
   one_arg_ret retval;
   retval.arg.slope = 1.0 / sqrt(1.0 - arg * arg);
   retval.arg.curve = arg * retval.arg.slope * retval.arg.slope
                      * retval.arg.slope;
   retval.arg.disc_type = undefined_beyond;
   if (arg > 0.0)
      retval.arg.disc_dist = 1.0 - arg;
   else
      retval.arg.disc_dist = arg + 1.0;
   retval.value = asin(arg);
   return retval;
}
inline one_arg_ret acos_w_moments(const double& arg)
{
   one_arg_ret retval;
   retval.arg.slope = -1.0 / sqrt(1.0 - arg * arg);
   retval.arg.curve = arg * retval.arg.slope * retval.arg.slope
                      * retval.arg.slope;
   retval.arg.disc_type = undefined_beyond;
   if (arg > 0.0)
      retval.arg.disc_dist = 1.0 - arg;
   else
      retval.arg.disc_dist = arg + 1.0;
   retval.value = acos(arg);
   return retval;
}
inline one_arg_ret atan_w_moments(const double& arg)
{
   one_arg_ret retval;
   retval.arg.slope = 1.0 / (1.0 + arg * arg);
   retval.arg.curve = 2.0 * arg * retval.arg.slope * retval.arg.slope;
   retval.arg.disc_type = none;
   retval.value = atan(arg);
   return retval;
}
inline two_arg_ret atan2_w_moments(const double& arg1, const double& arg2)
{
   two_arg_ret retval;
   double sum2 = arg2 * arg2 + arg1 * arg1;

   if (sum2 == 0.0)
   {
      retval.arg1.slope = 1.0;
      retval.arg2.slope = 1.0;
      retval.arg1.curve = 0.0;
      retval.arg2.curve = 0.0;
      retval.arg1.disc_type = none;
      retval.arg2.disc_type = step;
      retval.arg2.disc_dist = 0.0;
   }
   else
   {
      retval.arg1.slope = arg2 / sum2;
      retval.arg2.slope = -arg1 / sum2;
      retval.arg1.curve = -2.0 * arg1 * arg2 / (sum2 * sum2);
      retval.arg2.curve = -retval.arg1.curve;
      if (arg1 == 0.0)
      {
         retval.arg2.disc_type = step;
         retval.arg2.disc_dist = fabs(arg2);
      }
      else
         retval.arg2.disc_type = none;

      if (arg2 >= 0.0)
         retval.arg1.disc_type = none;
      else
      {
         retval.arg1.disc_type = step;
         retval.arg1.disc_dist = fabs(arg1);
      }
   }
   retval.value = atan2(arg1, arg2);
   return retval;
}
inline one_arg_ret exp_w_moments(const double& arg)
{
   one_arg_ret retval;
   retval.value = exp(arg);
   retval.arg.slope = retval.value;
   retval.arg.curve = retval.arg.slope;
   retval.arg.disc_type = none;
   return retval;
}
inline one_arg_ret log_w_moments(const double& arg)
{
   one_arg_ret retval;
   retval.arg.slope = 1.0 / arg;
   retval.arg.curve = -retval.arg.slope * retval.arg.slope;
   retval.arg.disc_dist = arg;
   retval.arg.disc_type = undefined_beyond;
   retval.value = log(arg);
   return retval;
}
inline one_arg_ret log10_w_moments(const double& arg)
{
   one_arg_ret retval;
   retval.arg.slope = 0.43429448189 / arg;
   retval.arg.curve = -retval.arg.slope / arg;
   retval.arg.disc_dist = arg;
   retval.arg.disc_type = undefined_beyond;
   retval.value = log10(arg);
   return retval;
}
inline one_arg_ret sinh_w_moments(const double& arg)
{
   one_arg_ret retval;
   retval.arg.slope = cosh(arg);
   retval.value = sinh(arg);
   retval.arg.curve = retval.value;
   retval.arg.disc_type = none;
   return retval;
}
inline one_arg_ret cosh_w_moments(const double& arg)
{
   one_arg_ret retval;
   retval.arg.slope = sinh(arg);
   retval.value = cosh(arg);
   retval.arg.curve = retval.value;
   retval.arg.disc_type = none;
   return retval;
}
inline one_arg_ret tanh_w_moments(const double& arg)
{
   one_arg_ret retval;
   double coshtemp = cosh(arg);
   retval.arg.slope = 1.0 / (coshtemp * coshtemp);
   retval.arg.curve = -2.0 * sinh(arg) / (coshtemp * coshtemp *coshtemp);
   retval.arg.disc_type = none;
   retval.value = tanh(arg);
   return retval;
}
inline two_arg_ret pow_w_moments(const double& arg1, const double& arg2)
{
   two_arg_ret retval;
   retval.value = pow(arg1, arg2);
   if (arg1 == 0.0)
   {
      retval.arg1.slope = 0.0;
      retval.arg1.curve = 0.0;
      if (arg2 == 1.0)
         retval.arg1.slope = 1.0;
      else if (arg2 == 2.0)
         retval.arg1.curve = 2.0;
      retval.arg1.disc_type = none;  // pow(arg1, 0.0) is 1.0 for any arg1
      // pow(0.0, arg2) is undefined for every negative arg2
      retval.arg2.disc_type = undefined_beyond;
      // pow(0.0, arg2) is 0.0 for every positive arg2, undef for neg arg2
      retval.arg2.slope = 0.0;
      retval.arg2.curve = 0.0;
      if (arg2 > 0.0)
         retval.arg2.disc_dist = arg2;
      else
         retval.arg2.disc_dist = 0.0;
   }
   else if (arg1 < 0.0)
   {
      // pow(arg1, arg2) for arg1 < 0.0 is only defined for integer arg2
      retval.arg1.slope = arg2 * retval.value / arg1;
      retval.arg1.curve = arg2 * (arg2 - 1.0) * retval.value / (arg1 * arg1);
      retval.arg1.disc_type = none;
      retval.arg2.slope = 0.0;
      retval.arg2.curve = 0.0;
      retval.arg2.disc_type = undefined_beyond;
      retval.arg2.disc_dist = 0.0;
   }
   else
   {
      retval.arg1.slope = arg2 * retval.value / arg1;
      retval.arg1.curve = arg2 * (arg2 - 1.0) * retval.value / (arg1 * arg1);
      if (arg2 == floor(arg2)) // if arg2 is an integer
      {
         retval.arg1.disc_type = none;
      }
      else
      {
         retval.arg1.disc_type = undefined_beyond;
         retval.arg1.disc_dist = arg1;
      }
      retval.arg2.slope = log(arg1) * retval.value;
      retval.arg2.curve = log(arg1) * retval.arg2.slope;
      retval.arg2.disc_type = none;
   }
   return retval;
}

// future direction: restrict namespace
// This function returns an approximation of the inverse Gaussian denstity
// function to within 4.5e-4.  From Abromowitz & Stegun's _Handbook_of_
// _Mathematical_Functions_ formula 26.2.23
inline double inverse_gaussian_density(const double& p)
{
   if (p <= 0.0) {
      std::cerr << "inverse_gaussian_density() called for negative value: "
           << p << std::endl;
      exit(EXIT_FAILURE);
   } else if (p > 0.5) {
      std::cerr << "inverse_gaussian_density() called for too large value: "
           << p << std::endl;
      exit(EXIT_FAILURE);
   }

   const double c0 = 2.515517, c1 = 0.802853, c2 = 0.010328;
   const double d1 = 1.432788, d2 = 0.189269, d3 = 0.001308;
   double t = sqrt(log(1.0 / (p * p)));

   return t - (c0 + t * (c1 + c2 * t))
              / (1.0 + t * (d1 + t * (d2 + t * d3)));
}


// future direction: include skewing of distribution
inline void gauss_loss(const double &uncertainty, const double &disc_dist,
                       const discontinuity_type &disc_type,
                       const char * const id_string,
                       const char * const func_str,
                       const double &disc_thresh)
   {
      double scaled_disc_dist = fabs(disc_dist / uncertainty);
      if ((scaled_disc_dist < disc_thresh) && (disc_type != none))
      {
         int original_precision = std::cerr.precision();
         std::cerr << std::setprecision(2);
         std::cerr << func_str << "is " << scaled_disc_dist << " sigmas"
              << id_string;
         if (disc_type == step)
            std::cerr << " from a step discontinuity" << std::endl;
         else if (disc_type == infinite_wrap)
            std::cerr << " from an infinite wrap discontinuity" << std::endl;
         else if (disc_type == infinite_then_undef)
            std::cerr << " from an infinite "
                 << "discontinuity beyond which it is undefined" << std::endl;
         else if (disc_type == slope_only)
            std::cerr << " from a discontinuity in slope" << std::endl;
         else if (disc_type == undefined_beyond)
            std::cerr << " from a point beyond which it is undefined" << std::endl;
         else
            std::cerr << " from unknown discontinuity " << disc_type << std::endl;
         std::cerr << std::setprecision(original_precision);
      }
   }



// model uncertain number using only mean and sigma, like UDoubleMS,
// but also use some knowledge of curve & skew & discontinuities.
// (The "C" is for "Curve").  For simple cases this gives better
// answers than the simple UDoubleMS classes, but when multiple
// functions with significant curve are applied after each other
// it starts giving much worse answers.
//
// Ultimately this work will be more useful when it is incorporated
// into a correlation tracking class.
template <int is_correlated>
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
         std::cerr << "Error: negative uncertainty: " << unc << std::endl;
         exit(EXIT_FAILURE);
      }
   }
   UDoubleMSC(const UDoubleMSC& ud)
      : value(ud.value), uncertainty(ud.uncertainty) {}

   ~UDoubleMSC(void) {}

   UDoubleMSC<is_correlated> operator +(void) const
      {return *this;}
   UDoubleMSC<is_correlated> operator -(void) const
   {
      if (is_correlated)
         return UDoubleMSC<is_correlated>(-value, -uncertainty);
      else
         return UDoubleMSC<is_correlated>(-value, uncertainty);
   }
   friend UDoubleMSC<is_correlated> operator +(UDoubleMSC<is_correlated> a,
                                      const UDoubleMSC<is_correlated>& b)
      { return a += b; }
   friend UDoubleMSC<is_correlated> operator -(UDoubleMSC<is_correlated> a,
                                      const UDoubleMSC<is_correlated>& b)
      { return a -= b; }
   UDoubleMSC<is_correlated> operator ++(void)
      { return (*this += 1.0); }
   UDoubleMSC<is_correlated> operator --(void)
      { return (*this -= 1.0); }
   UDoubleMSC<is_correlated> operator ++(int)
   {
      UDoubleMSC<is_correlated> retval(*this);
      *this += 1.0;
      return retval;
   }
   UDoubleMSC<is_correlated> operator --(int)
   {
      UDoubleMSC<is_correlated> retval(*this);
      *this -= 1.0;
      return retval;
   }
   friend UDoubleMSC<is_correlated> operator *(UDoubleMSC<is_correlated> a,
                                      const UDoubleMSC<is_correlated>& b)
      { return a *= b; }
   friend UDoubleMSC<is_correlated> operator /(UDoubleMSC<is_correlated> a,
                                      const UDoubleMSC<is_correlated>& b)
      { return a /= b; }
   UDoubleMSC<is_correlated> &operator +=(const UDoubleMSC<is_correlated>& ud)
   {
      if (is_correlated)
         uncertainty += ud.uncertainty;
      else
         uncertainty = hypot(uncertainty, ud.uncertainty);
      value += ud.value;
      return *this;
   }
   UDoubleMSC<is_correlated> &operator -=(const UDoubleMSC<is_correlated>& ud)
   {
      if (is_correlated)
         uncertainty -= ud.uncertainty;
      else
         uncertainty = hypot(uncertainty, ud.uncertainty);
      value -= ud.value;
      return *this;
   }
   UDoubleMSC<is_correlated> &operator *=(const UDoubleMSC<is_correlated>& ud)
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
   UDoubleMSC<is_correlated> &operator /=(const UDoubleMSC<is_correlated>& ud)
   {
      double second_order_correlated_adjust;
      if (is_correlated)
      {
         second_order_correlated_adjust = (ud.uncertainty / (ud.value*ud.value))
                       * (value * ud.uncertainty / ud.value - uncertainty);
         uncertainty = (uncertainty / ud.value
                        - (ud.uncertainty * value) / (ud.value * ud.value))
                         * sqrt(1.0 + 2.0 * sqr(ud.uncertainty / ud.value));
 
         if (ud.uncertainty != 0.0)
         {
            double disc_dist = fabs(ud.value/ud.uncertainty);
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
   friend std::ostream& operator <<(std::ostream &os, const UDoubleMSC<is_correlated> &ud)
   {
      uncertain_print(ud.mean(), ud.deviation(), os);
      return os;
   }
   friend std::istream& operator >>(std::istream &is, UDoubleMSC<is_correlated> &ud)
   {
      double mean, sigma;
      uncertain_read(mean, sigma, is);
      ud = UDoubleMSC<is_correlated>(mean, sigma);
      return is;
   }
   #define UDoubleMSCfunc1(func) \
      UDoubleMSC<is_correlated> func(UDoubleMSC<is_correlated> arg) \
      { \
         std::ostrstream os; \
         os << #func "(" << arg << ") " << std::ends; \
         one_arg_ret funcret = func ## _w_moments(arg.value); \
         arg.value = funcret.value + sqr(arg.uncertainty) \
                      * funcret.arg.curve / 2.0; \
         gauss_loss(arg.uncertainty, funcret.arg.disc_dist, \
                    funcret.arg.disc_type, "", os.str(), \
                    UDoubleMSC<is_correlated>::discontinuity_thresh); \
         os.rdbuf()->freeze(0); \
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
// future direction: enhance the below to account for moments with terms of
// each argument (e.g. for atan2(x,y), we now ignore d/dx(d/dy(atan2(x,y)))
   #define UDoubleMSCfunc2(func) \
      UDoubleMSC<is_correlated> func(const UDoubleMSC<is_correlated>& arg1, \
                                    const UDoubleMSC<is_correlated>& arg2) \
      { \
         UDoubleMSC<is_correlated> retval; \
         double unc1, unc2; \
         std::ostrstream os; \
         os << #func "(" << arg1 << ", " << arg2 << ") " << std::ends; \
         char *str = os.str(); \
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
         os.rdbuf()->freeze(0); \
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
      std::ostrstream os;
      os << "ldexp(" << arg << ", " << intarg << ") " << std::ends;
      one_arg_ret funcret = ldexp_w_moments(arg.value, intarg);
      arg.value = funcret.value + sqr(arg.uncertainty) * funcret.arg.curve / 2.0;
      gauss_loss(arg.uncertainty, funcret.arg.disc_dist, funcret.arg.disc_type,
                 "", os.str(), UDoubleMSC<is_correlated>::discontinuity_thresh);
      os.rdbuf()->freeze(0);
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
                                          int *intarg)
   {
      std::ostrstream os;
      os << "frexp(" << arg << ", " << *intarg << ") " << std::ends;
      one_arg_ret funcret = frexp_w_moments(arg.value, *intarg);
      arg.value = funcret.value + sqr(arg.uncertainty) * funcret.arg.curve / 2.0;
      gauss_loss(arg.uncertainty, funcret.arg.disc_dist, funcret.arg.disc_type,
                 "", os.str(), UDoubleMSC<is_correlated>::discontinuity_thresh);
      os.rdbuf()->freeze(0);
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
                                          double *dblarg)
   {
      std::ostrstream os;
      os << "modf(" << arg << ", " << *dblarg << ") " << std::ends;
      one_arg_ret funcret = modf_w_moments(arg.value, *dblarg);
      arg.value = funcret.value + sqr(arg.uncertainty) * funcret.arg.curve / 2.0;
      gauss_loss(arg.uncertainty, funcret.arg.disc_dist, funcret.arg.disc_type,
                 "", os.str(), UDoubleMSC<is_correlated>::discontinuity_thresh);
      os.rdbuf()->freeze(0);
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
   double mean(void) const {return value;}
   double deviation(void) const {return fabs(uncertainty);}

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
                                   double (*certain_func)(double, double),
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
// future direction: add Monte-Carlo propogation
};

template<>
double UDoubleMSC<0>::discontinuity_thresh = 3.0;

template<>
double UDoubleMSC<1>::discontinuity_thresh = 0.0;


// future direction: make template variables or class consts
#define MAX_UNC_ELEMENTS 5
#define MAX_SRC_NAME 64

// A class for sources of uncertainties.
class UncertainSourceSet {
private:
   unsigned long num_sources;
   unsigned long source_epoch;
   char source_name[MAX_UNC_ELEMENTS][MAX_SRC_NAME + 1];
   char class_name[100];
public:
   UncertainSourceSet(char *cname = "") : num_sources(0), source_epoch(0)
   {
      for (unsigned int i = 0; i < MAX_UNC_ELEMENTS; i++)
         source_name[i][0] = 0;
      strcpy(class_name, cname);
   }
   unsigned long get_epoch(void) const { return source_epoch; } const
   void check_epoch(const unsigned long & epoch) const
   {
      if (epoch != source_epoch)
      {
         std::cerr << "Bad epoch: " << epoch << " expected: " << source_epoch
              << " in class " << class_name << std::endl;
         exit(EXIT_FAILURE);
      }
   }
   void new_epoch(void)
   {
      for (unsigned int i = 0; i < num_sources; i++)
         source_name[i][0] = 0;
      source_epoch++;
      num_sources = 0;
   }
   int can_get_new_source(void) const
   {
      if (num_sources >= MAX_UNC_ELEMENTS)
      {
         std::cerr << "Cannot make new UncertainSource" << std::endl;
         std::cerr << "Already have maximum number of permissible uncertainty "
              << "elements: " << MAX_UNC_ELEMENTS << "(" << num_sources
              << ")" << "in class " << class_name << std::endl;
         std::cerr << "Change the value of MAX_UNC_ELEMENTS and recompile"
              << " or use new_epoch()" << std::endl;
         exit(EXIT_FAILURE);
      }
      return 1;
   }
   unsigned long get_new_source(const char * const name)
   {
      // std::cerr << "trying to get new source (" << num_sources << ") for "
      //      << name << std::endl;
      (void)can_get_new_source();
      strncpy(source_name[num_sources], name, MAX_SRC_NAME);
      source_name[num_sources][MAX_SRC_NAME] = 0;
      return num_sources++;
   }
   unsigned long get_num_sources(void) const
   {
      return num_sources;
   }
   const char * const get_source_name(const unsigned long i) const
   {
      if (i >= num_sources)
      {
         std::cerr << "get_source_name called with illegal source number: "
              << i << std::endl;
         exit(EXIT_FAILURE);
      }
      return source_name[i];
   }
};


// Specialized array class has only those members needed to be an array
// of uncertainty elements used as the template parameter in UDoubleCT<>.
// This is the simplest possible implementation.
template <int size>
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
   ~SimpleArray(void) {}
   SimpleArray operator -(void) const
   {
      SimpleArray retval;

      for (int i = 0; i < size; i++)
         retval.element[i] = -element[i];
      return retval;
   }
   SimpleArray &operator +=(const SimpleArray& b)
   {
      for (int i = 0; i < size; i++)
         element[i] += b.element[i];
      return *this;
   }
   friend SimpleArray operator +(SimpleArray a, const SimpleArray& b)
      { return a += b; }
   SimpleArray &operator -=(const SimpleArray& b)
   {
      for (int i = 0; i < size; i++)
         element[i] -= b.element[i];
      return *this;
   }
   SimpleArray &operator *=(const double& b)
   {
      for (int i = 0; i < size; i++)
         element[i] *= b;
      return *this;
   }
   friend SimpleArray operator *(SimpleArray a, const double& b)
      { return a *= b; }
   SimpleArray &operator /=(const double& b)
   {
      for (int i = 0; i < size; i++)
         element[i] /= b;
      return *this;
   }
   double operator[](int subscript)
   {
      if (subscript < 0)
      {
         std::cerr << "Error: negative subscript: " << subscript << std::endl;
         exit(EXIT_FAILURE);
      }
      if (subscript >= size)
      {
         std::cerr << "Error: oversize subscript: " << subscript
              << " Greater than " << (size - 1) << std::endl;
         exit(EXIT_FAILURE);
      }
      return element[subscript];
   }
   void setelement(int subscript, double value)
   {
      if (subscript < 0)
      {
         std::cerr << "Error: negative subscript: " << subscript << std::endl;
         exit(EXIT_FAILURE);
      }
      if (subscript >= size)
      {
         std::cerr << "Error: oversize subscript: " << subscript
              << " Greater than " << (size - 1) << std::endl;
         exit(EXIT_FAILURE);
      }
      element[subscript] = value;
   }
   double norm(void) const
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
template <int size>
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
   ~ArrayWithScale(void) {}
   ArrayWithScale operator -(void) const
   {
      ArrayWithScale retval;

      for (int i = 0; i < size; i++)
         retval.element[i] = element[i];
      retval.scale = -scale;
      return retval;
   }
   ArrayWithScale &operator +=(const ArrayWithScale& b)
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
   friend ArrayWithScale operator +(ArrayWithScale a, const ArrayWithScale& b)
      { return a += b; }
   ArrayWithScale &operator -=(const ArrayWithScale& b)
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
   ArrayWithScale &operator *=(const double& b)
   {
      scale *= b;
      return *this;
   }
   friend ArrayWithScale operator *(ArrayWithScale a, const double& b)
      { return a *= b; }
   ArrayWithScale &operator /=(const double& b)
   {
      scale /= b;
      return *this;
   }
   double operator[](int subscript)
   {
      if (subscript < 0)
      {
         std::cerr << "Error: negative subscript: " << subscript << std::endl;
         exit(EXIT_FAILURE);
      }
      if (subscript >= size)
      {
         std::cerr << "Error: oversize subscript: " << subscript
              << " Greater than " << (size - 1) << std::endl;
         exit(EXIT_FAILURE);
      }
      return element[subscript] * scale;
   }
   void setelement(int subscript, double value)
   {
      if (subscript < 0)
      {
         std::cerr << "Error: negative subscript: " << subscript << std::endl;
         exit(EXIT_FAILURE);
      }
      if (subscript >= size)
      {
         std::cerr << "Error: oversize subscript: " << subscript
              << " Greater than " << (size - 1) << std::endl;
         exit(EXIT_FAILURE);
      }
      if (scale != 0.0)
         element[subscript] = value / scale;
   }
   double norm(void) const
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
template <class T>
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
              const char * const name = "")
   : value(val), unc_components(0.0)
   {
      epoch = sources.get_epoch();
      if (unc < 0.0)
      {
         std::cerr << "Error: negative uncertainty: " << unc << std::endl;
         exit(EXIT_FAILURE);
      }
      if (unc != 0.0)
      {
         char source_name[MAX_SRC_NAME + 1];
         if (name && name[0])
         {
            strncpy(source_name, name, MAX_SRC_NAME);
            source_name[MAX_SRC_NAME] = 0;
         }
         else
         {
            std::ostrstream os;
            os << "anon: ";
            uncertain_print(val, unc, os);
            os << std::ends;
            strncpy(source_name, os.str(), MAX_SRC_NAME);
            os.rdbuf()->freeze(0);
            source_name[MAX_SRC_NAME] = 0;
         }
         unsigned long new_source_num = sources.get_new_source(source_name);
         unc_components.setelement(new_source_num, unc);
      }
   }
   // copy constructor does not create a new independent uncertainty element
   UDoubleCT(const UDoubleCT& ud) : value(ud.value),
                                      epoch(ud.epoch)
   {
      unc_components = ud.unc_components;
   }
   // operator= ?

   ~UDoubleCT(void) {}

   static void new_epoch(void) { sources.new_epoch(); }
   void print_uncertain_sources(std::ostream &os = std::cout)
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

   UDoubleCT operator +(void) const
      {return *this;}
   UDoubleCT operator -(void) const
   {
      UDoubleCT retval;

      retval.value = -value;
      retval.unc_components = -unc_components;
      return retval;
   }
   UDoubleCT &operator +=(const UDoubleCT& b)
   {
      sources.check_epoch(epoch);
      sources.check_epoch(b.epoch);
      unc_components += b.unc_components;
      value += b.value;
      return *this;
   }
   UDoubleCT &operator +=(const double& b)
   {
      value += b;
      return *this;
   }
   friend UDoubleCT operator +(UDoubleCT a, const UDoubleCT& b)
      { return a += b; }
   friend UDoubleCT operator +(UDoubleCT a, const double& b)
      { return a += b; }
   friend UDoubleCT operator +(const double& b, UDoubleCT a)
      { return a += b; }
   UDoubleCT &operator -=(const UDoubleCT& b)
   {
      sources.check_epoch(epoch);
      sources.check_epoch(b.epoch);
      unc_components -= b.unc_components;
      value -= b.value;
      return *this;
   }
   UDoubleCT &operator -=(const double& b)
   {
      value -= b;
      return *this;
   }
   friend UDoubleCT operator -(UDoubleCT a, const UDoubleCT& b)
      { return a -= b; }
   friend UDoubleCT operator -(UDoubleCT a, const double& b)
      { return a -= b; }
   friend UDoubleCT operator -(const double& b, UDoubleCT a)
      { a -= b; return -a; }
   UDoubleCT &operator *=(const UDoubleCT& b)
   {
      sources.check_epoch(epoch);
      sources.check_epoch(b.epoch);
      unc_components *= b.value;
      unc_components += b.unc_components * value;
      value *= b.value;
      return *this;
   }
   UDoubleCT &operator *=(const double& b)
   {
      unc_components *= b;
      value *= b;
      return *this;
   }
   UDoubleCT operator ++(void)
      { return (*this += 1.0); }
   UDoubleCT operator --(void)
      { return (*this -= 1.0); }
   UDoubleCT operator ++(int)
   {
      UDoubleCT retval(*this);
      *this += 1.0;
      return retval;
   }
   UDoubleCT operator --(int)
   {
      UDoubleCT retval(*this);
      *this -= 1.0;
      return retval;
   }
   friend UDoubleCT operator *(UDoubleCT a, const UDoubleCT& b)
      { return a *= b; }
   friend UDoubleCT operator *(UDoubleCT a, const double& b)
      { return a *= b; }
   friend UDoubleCT operator *(const double& b, UDoubleCT a)
      { return a *= b; }
   UDoubleCT &operator /=(const UDoubleCT& b)
   {
      sources.check_epoch(epoch);
      sources.check_epoch(b.epoch);
      unc_components /= b.value;
      unc_components -= b.unc_components * (value / (b.value * b.value));
      value /= b.value;
      return *this;
   }
   UDoubleCT &operator /=(const double& b)
   {
      unc_components /= b;
      value /= b;
      return *this;
   }
   friend UDoubleCT operator /(UDoubleCT a, const UDoubleCT& b)
      { return a /= b; }
   friend UDoubleCT operator /(UDoubleCT a, const double& b)
      { return a /= b; }
   friend UDoubleCT operator /(const double a, const UDoubleCT& b)
   {
      UDoubleCT retval(0.0);

      retval.unc_components = b.unc_components * (-a / (b.value * b.value));
      retval.value = a / b.value;
      return retval;
   }
   friend std::ostream& operator <<(std::ostream &os, const UDoubleCT &ud)
   {
      uncertain_print(ud.mean(), ud.deviation(), os);
      return os;
   }
   friend std::istream& operator >>(std::istream &is, UDoubleCT &ud)
   {
      double mean, sigma;
      char source_name[MAX_SRC_NAME + 1];

      uncertain_read(mean, sigma, is);
      std::ostrstream os;
      os << "input: ";
      uncertain_print(mean, sigma, os);
      os << std::ends;
      strncpy(source_name, os.str(), MAX_SRC_NAME);
      os.rdbuf()->freeze(0);
      source_name[MAX_SRC_NAME] = 0;
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
   friend UDoubleCT frexp(UDoubleCT arg, int * intarg)
   {
      one_arg_ret funcret = frexp_w_moments(arg.value, *intarg);
      arg.value = funcret.value;
      arg.unc_components *= funcret.arg.slope;
      return arg;
   }
   friend UDoubleCT modf(UDoubleCT arg, double * dblarg)
   {
      one_arg_ret funcret = modf_w_moments(arg.value, *dblarg);
      arg.value = funcret.value;
      arg.unc_components *= funcret.arg.slope;
      return arg;
   }
   double mean(void) const {return value;}
   double deviation(void) const {return unc_components.norm();}

};

typedef SimpleArray<MAX_UNC_ELEMENTS> SizedSimpleArray;
typedef ArrayWithScale<MAX_UNC_ELEMENTS> SizedArrayWithScale;

template<>
UncertainSourceSet UDoubleCT<SizedSimpleArray> ::sources("Simple Array");

template<>
UncertainSourceSet UDoubleCT<SizedArrayWithScale> ::sources("Array with Scale");

typedef UDoubleCT<SizedSimpleArray> UDoubleCTSA;
typedef UDoubleCT<SizedArrayWithScale> UDoubleCTAA;


// Compare function to sort by absolute value
// future direction: restrict namespace
inline int abs_double_compare(const void * a, const void * b)
   {
      double fa = fabs(*(double *)a);
      double fb = fabs(*(double *)b);
      if (fa == fb)
         return 0;
      else if (fa < fb)
         return -1;
      else
         return 1;
   }
