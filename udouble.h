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

// These flags disable a variety of warnings on library header files
// for Gimpel Flexi-lint
//lint -save
//lint -e762 -e578 -e1917 -e1918 -e1510 -e1712 -e1730 -e655 -e1909 -e578
//lint -e1511 -e1907 -e36 -e1908 -e763 -e534
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream.h>
#include <strstream.h>
#include <iomanip.h>
#include <time.h>
//lint -restore

#define PI 3.14159265358979323846
#define HALF_PI (PI / 2.0)

// future direction: add exceptions


// future direction: these functions should be static in implementation
//                   file or hidden in a namespace (when namespaces come!)

// These functions take the square root of the sum of the squares of
// numbers, which is equal to the length of the hypotenuse of a right
// triangle if the two arguments are the lengths of the legs.
inline double hypot(const double &a, const double &b)
   { return sqrt(a * a + b * b); }
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


// I don't think sqr() is standard but I get a collision with a
// suitable sqr() in math.h on my hp using gcc
#if NEED_SQR
// square: just a notational convenience
inline double sqr(const double &a)
   { return a * a; }
#endif

// future direction: restrict namespace
// prints uncertainty to 2 digits and value to same precision
inline void uncertain_print(double mean, double sigma, ostream &os = cout)
   {
      int original_precision = os.precision();
      long original_format = os.flags(ios::showpoint);

      // cerr << "<" << mean << " +/- " << sigma << "> " << endl;
      int precision;
      // special cases for zero, NaN, and Infinities (positive & negative)
      if ((sigma == 0.0) || (sigma != sigma) || (1.0/sigma == 0.0))
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
      os << setprecision(precision)
         << mean << " +/- "
         << setprecision(2)
         << sigma
         << setprecision(original_precision);
      os.flags(original_format);
   }

// reads uncertainty as mean +/- sigma
inline void uncertain_read(double& mean, double& sigma, istream &is = cin)
   {
      char plus, slash, minus;
      is >> mean >> plus >> slash >> minus >> sigma;
      if ((plus != '+') || (slash != '/') || (minus != '-'))
      {
         cerr << "Error: illegal characters encountered in reading "
                 "mean +/- sigma" << endl;
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
template <int is_correlated>
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
         cerr << "Error: negative uncertainty: " << unc << endl;
         exit(EXIT_FAILURE);
      }
   }
   UDoubleMS(const UDoubleMS& ud)
      : value(ud.value), uncertainty(ud.uncertainty) {}

   ~UDoubleMS(void) {}

   UDoubleMS<is_correlated> operator +(void) const
      {return *this;}
   UDoubleMS<is_correlated> operator -(void) const
   {
      if (is_correlated)
         return UDoubleMS<is_correlated>(-value, -uncertainty);
      else
         return UDoubleMS<is_correlated>(-value, uncertainty);
   }
   friend UDoubleMS<is_correlated> operator +(UDoubleMS<is_correlated> a,
                                      const UDoubleMS<is_correlated>& b)
      { return a += b; }
   friend UDoubleMS<is_correlated> operator +(UDoubleMS<is_correlated> a,
                                      const double& b)
      { return a += b; }
   friend UDoubleMS<is_correlated> operator +(const double& b,
                                      UDoubleMS<is_correlated> a)
      { return a += b; }
   friend UDoubleMS<is_correlated> operator -(UDoubleMS<is_correlated> a,
                                      const UDoubleMS<is_correlated>& b)
      { return a -= b; }
   friend UDoubleMS<is_correlated> operator -(UDoubleMS<is_correlated> a,
                                      const double& b)
      { return a -= b; }
   friend UDoubleMS<is_correlated> operator -(const double& b,
                                      UDoubleMS<is_correlated> a)
      { a -= b; return -a; }
   UDoubleMS<is_correlated> operator ++(void)
      { return (*this += 1.0); }
   UDoubleMS<is_correlated> operator --(void)
      { return (*this -= 1.0); }
   UDoubleMS<is_correlated> operator ++(int)
   {
      UDoubleMS<is_correlated> retval(*this);
      *this += 1.0;
      return retval;
   }
   UDoubleMS<is_correlated> operator --(int)
   {
      UDoubleMS<is_correlated> retval(*this);
      *this -= 1.0;
      return retval;
   }
   friend UDoubleMS<is_correlated> operator *(UDoubleMS<is_correlated> a,
                                      const UDoubleMS<is_correlated>& b)
      { return a *= b; }
   friend UDoubleMS<is_correlated> operator *(UDoubleMS<is_correlated> a,
                                      const double& b)
      { return a *= b; }
   friend UDoubleMS<is_correlated> operator *(const double& b,
                                      UDoubleMS<is_correlated> a)
      { return a *= b; }
   friend UDoubleMS<is_correlated> operator /(UDoubleMS<is_correlated> a,
                                      const UDoubleMS<is_correlated>& b)
      { return a /= b; }
   friend UDoubleMS<is_correlated> operator /(UDoubleMS<is_correlated> a,
                                      const double& b)
      { return a /= b; }
   friend UDoubleMS<is_correlated> operator /(const double& b,
                                      UDoubleMS<is_correlated> a)
   {
      UDoubleMS<is_correlated> retval;
      retval.uncertainty = -b * a.uncertainty / (a.value * a.value);
      retval.value = b / a.value;
      return retval;
   }
   UDoubleMS<is_correlated> &operator +=(const UDoubleMS<is_correlated>& ud)
   {
      if (is_correlated)
         uncertainty += ud.uncertainty;
      else
         uncertainty = hypot(uncertainty, ud.uncertainty);
      value += ud.value;
      return *this;
   }
   UDoubleMS<is_correlated> &operator +=(const double& a)
   {
      value += a;
      return *this;
   }
   UDoubleMS<is_correlated> &operator -=(const UDoubleMS<is_correlated>& ud)
   {
      if (is_correlated)
         uncertainty -= ud.uncertainty;
      else
         uncertainty = hypot(uncertainty, ud.uncertainty);
      value -= ud.value;
      return *this;
   }
   UDoubleMS<is_correlated> &operator -=(const double& a)
   {
      value -= a;
      return *this;
   }
   UDoubleMS<is_correlated> &operator *=(const UDoubleMS<is_correlated>& ud)
   {
      if (is_correlated)
         uncertainty = uncertainty * ud.value + ud.uncertainty * value;
      else
         uncertainty = hypot(uncertainty * ud.value, ud.uncertainty * value);
      value *= ud.value;
      return *this;
   }
   UDoubleMS<is_correlated> &operator *=(const double& a)
   {
      value *= a;
      uncertainty *= a;
      return *this;
   }
   UDoubleMS<is_correlated> &operator /=(const UDoubleMS<is_correlated>& ud)
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
   UDoubleMS<is_correlated> &operator /=(const double& a)
   {
      value /= a;
      uncertainty /= a;
      return *this;
   }
   friend ostream& operator <<(ostream &os, const UDoubleMS<is_correlated> &ud)
   {
      uncertain_print(ud.mean(), ud.deviation(), os);
      return os;
   }
   friend istream& operator >>(istream &is, UDoubleMS<is_correlated> &ud)
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
                                         double *intpart)
   {
      arg.value = modf(arg.value, intpart);
      return arg;
   }
   friend UDoubleMS<is_correlated> frexp(UDoubleMS<is_correlated> arg,
                                         int *intarg)
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
      if ((arg1.value/arg2.value) > 0.0)
         slope2 = -floor(arg1.value/arg2.value);
      else
         slope2 = floor(-arg1.value/arg2.value);
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
   double mean(void) const {return value;}
   double deviation(void) const {return fabs(uncertainty);}

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
                                   double (*certain_func)(double, double),
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
      cerr << "inverse_gaussian_density() called for negative value: "
           << p << endl;
      exit(EXIT_FAILURE);
   } else if (p > 0.5) {
      cerr << "inverse_gaussian_density() called for too large value: "
           << p << endl;
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
         int original_precision = cerr.precision();
         cerr << setprecision(2);
         cerr << func_str << "is " << scaled_disc_dist << " sigmas"
              << id_string;
         if (disc_type == step)
            cerr << " from a step discontinuity" << endl;
         else if (disc_type == infinite_wrap)
            cerr << " from an infinite wrap discontinuity" << endl;
         else if (disc_type == infinite_then_undef)
            cerr << " from an infinite "
                 << "discontinuity beyond which it is undefined" << endl;
         else if (disc_type == slope_only)
            cerr << " from a discontinuity in slope" << endl;
         else if (disc_type == undefined_beyond)
            cerr << " from a point beyond which it is undefined" << endl;
         else
            cerr << " from unknown discontinuity " << disc_type << endl;
         cerr << setprecision(original_precision);
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
         cerr << "Error: negative uncertainty: " << unc << endl;
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
               cerr << "correlated division by " << ud << " is ";
               cerr << disc_dist
                    << " sigmas from an infinite wrap discontinuity" << endl;
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
   friend ostream& operator <<(ostream &os, const UDoubleMSC<is_correlated> &ud)
   {
      uncertain_print(ud.mean(), ud.deviation(), os);
      return os;
   }
   friend istream& operator >>(istream &is, UDoubleMSC<is_correlated> &ud)
   {
      double mean, sigma;
      uncertain_read(mean, sigma, is);
      ud = UDoubleMSC<is_correlated>(mean, sigma);
      return is;
   }
   #define UDoubleMSCfunc1(func) \
      UDoubleMSC<is_correlated> func(UDoubleMSC<is_correlated> arg) \
      { \
         ostrstream os; \
         os << #func "(" << arg << ") " << ends; \
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
         ostrstream os; \
         os << #func "(" << arg1 << ", " << arg2 << ") " << ends; \
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
      ostrstream os;
      os << "ldexp(" << arg << ", " << intarg << ") " << ends;
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
      ostrstream os;
      os << "frexp(" << arg << ", " << *intarg << ") " << ends;
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
      ostrstream os;
      os << "modf(" << arg << ", " << *dblarg << ") " << ends;
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

double UDoubleMSC<0>::discontinuity_thresh = 3.0;
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
         cerr << "Bad epoch: " << epoch << " expected: " << source_epoch
              << " in class " << class_name << endl;
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
         cerr << "Cannot make new UncertainSource" << endl;
         cerr << "Already have maximum number of permissible uncertainty "
              << "elements: " << MAX_UNC_ELEMENTS << "(" << num_sources
              << ")" << "in class " << class_name << endl;
         cerr << "Change the value of MAX_UNC_ELEMENTS and recompile"
              << " or use new_epoch()" << endl;
         exit(EXIT_FAILURE);
      }
      return 1;
   }
   unsigned long get_new_source(const char * const name)
   {
      // cerr << "trying to get new source (" << num_sources << ") for "
      //      << name << endl;
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
         cerr << "get_source_name called with illegal source number: "
              << i << endl;
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
         cerr << "Error: negative subscript: " << subscript << endl;
         exit(EXIT_FAILURE);
      }
      if (subscript >= size)
      {
         cerr << "Error: oversize subscript: " << subscript
              << " Greater than " << (size - 1) << endl;
         exit(EXIT_FAILURE);
      }
      return element[subscript];
   }
   void setelement(int subscript, double value)
   {
      if (subscript < 0)
      {
         cerr << "Error: negative subscript: " << subscript << endl;
         exit(EXIT_FAILURE);
      }
      if (subscript >= size)
      {
         cerr << "Error: oversize subscript: " << subscript
              << " Greater than " << (size - 1) << endl;
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
         cerr << "Error: negative subscript: " << subscript << endl;
         exit(EXIT_FAILURE);
      }
      if (subscript >= size)
      {
         cerr << "Error: oversize subscript: " << subscript
              << " Greater than " << (size - 1) << endl;
         exit(EXIT_FAILURE);
      }
      return element[subscript] * scale;
   }
   void setelement(int subscript, double value)
   {
      if (subscript < 0)
      {
         cerr << "Error: negative subscript: " << subscript << endl;
         exit(EXIT_FAILURE);
      }
      if (subscript >= size)
      {
         cerr << "Error: oversize subscript: " << subscript
              << " Greater than " << (size - 1) << endl;
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
         cerr << "Error: negative uncertainty: " << unc << endl;
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
            ostrstream os;
            os << "anon: ";
            uncertain_print(val, unc, os);
            os << ends;
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
   void print_uncertain_sources(ostream &os = cout)
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
               << int_percent(unc_portion) << "%" << endl;
         }
      os << endl;
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
   friend ostream& operator <<(ostream &os, const UDoubleCT &ud)
   {
      uncertain_print(ud.mean(), ud.deviation(), os);
      return os;
   }
   friend istream& operator >>(istream &is, UDoubleCT &ud)
   {
      double mean, sigma;
      char source_name[MAX_SRC_NAME + 1];

      uncertain_read(mean, sigma, is);
      ostrstream os;
      os << "input: ";
      uncertain_print(mean, sigma, os);
      os << ends;
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

UncertainSourceSet UDoubleCT<SizedSimpleArray> ::sources("Simple Array");
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

// Ensemble uncertainty class.  Represents a distribution by a
// set of n=esize possible values distributed at intervals of
// uniform probability throughout the
// distribution.  The order of the possible values is randomized
// so that all operations can be done element-by-element with
// a minimum of spurious correlations between independent sources
// of uncertainty.  Depending on esize this class can be anywhere
// from very expensive computationally to unusably expensive.
// But for small problems and big esizes it gives "perfect"
// answers.
template <unsigned int esize>
class UDoubleEnsemble
{
private:
   double ensemble[esize];
   static UncertainSourceSet sources;
   static double src_ensemble[MAX_UNC_ELEMENTS][esize];
   unsigned long epoch;

public:
   // The main constructor initializes a new source of uncertainty
   // (if there is uncertainty).
   UDoubleEnsemble(const double val = 0.0, const double unc = 0.0,
                   const char * const name = 0)
              : epoch(sources.get_epoch())
   {
      if (unc < 0.0)
      {
         cerr << "Error: negative uncertainty: " << unc << endl;
         exit(EXIT_FAILURE);
      }
      if (unc != 0.0)
      {
         static double gauss_ensemble[esize];
         static int gauss_ensemble_inited = 0;

         // The base ensemble of n=esize points needs be initialized only
         // once for each ensemble size.  Once it is initialized, each new
         // independent uncertainty element can be made by copying & shuffling
         // this array then scaling it to the appropriate uncertainty and
         // translating it to the appropriate mean.
         if (!gauss_ensemble_inited)
         {
            gauss_ensemble_inited = 1;
            if (esize & 1) // odd ensemble size
            {
               for (unsigned int i = 0; i < esize / 2; i++)
               {
                  double deviate = inverse_gaussian_density(
                                        (2.0 * (i + 1.0)) / (2.0 * esize));
                  gauss_ensemble[2 * i] = deviate;
                  gauss_ensemble[2 * i + 1] = -deviate;
               }
               gauss_ensemble[esize - 1] = 0.0;
            }
            else
            {
               for (unsigned int i = 0; i < esize / 2; i++)
               {
                  double deviate = 0.0;
                  double k = (2.0 * i + 1.0) / (2.0 * esize);
                  for (unsigned j = 0; j < 100; j++)
                     deviate += inverse_gaussian_density(k + (j - 49.5)
                                                          / (100.0 * esize));
                  deviate /= 100.0;
                  gauss_ensemble[2 * i] = deviate;
                  gauss_ensemble[2 * i + 1] = -deviate;
               }
            }
            // Move the points a little to make all the first 5 moments give
            // exact values.
            PerfectEnsemble(gauss_ensemble, esize);
         }
         for (unsigned int i = 0; i < esize; i++)
            ensemble[i] = val + gauss_ensemble[i] * unc;
         char source_name[MAX_SRC_NAME + 1];
         if (name && name[0])
         {
            strncpy(source_name, name, MAX_SRC_NAME);
            source_name[MAX_SRC_NAME] = 0;
         }
         else
         {
            ostrstream os;
            os << "anon: ";
            uncertain_print(val, unc, os);
            os << ends;
            strncpy(source_name, os.str(), MAX_SRC_NAME);
            os.rdbuf()->freeze(0);
            source_name[MAX_SRC_NAME] = 0;
         }
         unsigned long source_num = sources.get_num_sources();
         this->shuffle();
         if (sources.can_get_new_source())
         {
            memcpy(src_ensemble[source_num], ensemble, sizeof(ensemble));
            source_num = sources.get_new_source(source_name);
         }
      }
      else  // uncertainty is zero
         for (unsigned long i = 0; i < esize; i++)
            ensemble[i] = val;
   }
   // copy constructor does not introduce a new uncertainty element
   UDoubleEnsemble(const UDoubleEnsemble& ud)
              : epoch(ud.epoch)
   {
      for (unsigned int i = 0; i < esize; i++)
         ensemble[i] = ud.ensemble[i];
   }
   // constructor from an ensemble.
   // future direction: add similar function that shuffles its input
   UDoubleEnsemble(const double * const newensemble,
                   const char * const name = "")
              : epoch(sources.get_epoch())
   {
      for (unsigned int i = 0; i < esize; i++)
         ensemble[i] = newensemble[i];
      char source_name[MAX_SRC_NAME + 1];
      if (name && name[0])
      {
         strncpy(source_name, name, MAX_SRC_NAME);
         source_name[MAX_SRC_NAME] = 0;
      }
      else
      {
         ostrstream os;
         os << "anon from ensemble: " << ensemble[0] << ends;
         strncpy(source_name, os.str(), MAX_SRC_NAME);
         os.rdbuf()->freeze(0);
         source_name[MAX_SRC_NAME] = 0;
      }
      unsigned long source_num = sources.get_num_sources();
      memcpy(src_ensemble[source_num], ensemble, sizeof(ensemble));
      source_num = sources.get_new_source(source_name);
   }
// future direction: add constructors with other distributions.

   ~UDoubleEnsemble(void) {}


   UDoubleEnsemble<esize> operator +(void) const
   {
      return *this;
   }
   UDoubleEnsemble<esize> operator -(void) const
   {
      UDoubleEnsemble<esize> retval;
      unsigned int i;

      for (i = 0; i < esize; i++)
         retval.ensemble[i] = -ensemble[i];
      retval.epoch = epoch;
      return retval;
   }
   friend UDoubleEnsemble<esize> operator +(UDoubleEnsemble<esize> a,
                                     const UDoubleEnsemble<esize>& b)
      { return a += b; }
   friend UDoubleEnsemble<esize> operator +(UDoubleEnsemble<esize> a,
                                     const double& b)
      { return a += b; }
   friend UDoubleEnsemble<esize> operator +(const double& b,
                                     UDoubleEnsemble<esize> a)
      { return a += b; }
   friend UDoubleEnsemble<esize> operator -(UDoubleEnsemble<esize> a,
                                     const UDoubleEnsemble<esize>& b)
      { return a -= b; }
   friend UDoubleEnsemble<esize> operator -(UDoubleEnsemble<esize> a,
                                     const double& b)
      { return a -= b; }
   friend UDoubleEnsemble<esize> operator -(const double& b,
                                     UDoubleEnsemble<esize> a)
      { return -(a -= b); }
   UDoubleEnsemble<esize> operator ++(void)
      { return (*this += 1.0); }
   UDoubleEnsemble<esize> operator --(void)
      { return (*this -= 1.0); }
   UDoubleEnsemble<esize> operator ++(int)
   {
      UDoubleEnsemble<esize> retval(*this);
      *this += 1.0;
      return retval;
   }
   UDoubleEnsemble<esize> operator --(int)
   {
      UDoubleEnsemble<esize> retval(*this);
      *this -= 1.0;
      return retval;
   }
   friend UDoubleEnsemble<esize> operator *(UDoubleEnsemble<esize> a,
                                     const UDoubleEnsemble<esize>& b)
      { return a *= b; }
   friend UDoubleEnsemble<esize> operator *(UDoubleEnsemble<esize> a,
                                     const double& b)
      { return a *= b; }
   friend UDoubleEnsemble<esize> operator *(const double& b,
                                     UDoubleEnsemble<esize> a)
      { return a *= b; }
   friend UDoubleEnsemble<esize> operator /(UDoubleEnsemble<esize> a,
                                     const UDoubleEnsemble<esize>& b)
      { return a /= b; }
   friend UDoubleEnsemble<esize> operator /(UDoubleEnsemble<esize> a,
                                     const double& b)
      { return a /= b; }
   // this one promotes a to UDoubleEnsemble
   friend UDoubleEnsemble<esize> operator /(const double& a,
                                     const UDoubleEnsemble<esize>& b)
      { UDoubleEnsemble<esize> uda(a); return uda /= b; }
   UDoubleEnsemble<esize> &operator +=(const UDoubleEnsemble<esize>& ud)
   {
      sources.check_epoch(epoch);
      sources.check_epoch(ud.epoch);

      for (unsigned int i = 0; i < esize; i++)
         ensemble[i] += ud.ensemble[i];
      return *this;
   }
   UDoubleEnsemble<esize> &operator +=(const double& d)
   {
      for (unsigned int i = 0; i < esize; i++)
         ensemble[i] += d;
      return *this;
   }
   UDoubleEnsemble<esize> &operator -=(const UDoubleEnsemble<esize>& ud)
   {
      sources.check_epoch(epoch);
      sources.check_epoch(ud.epoch);

      for (unsigned int i = 0; i < esize; i++)
         ensemble[i] -= ud.ensemble[i];
      return *this;
   }
   UDoubleEnsemble<esize> &operator -=(const double& d)
   {
      for (unsigned int i = 0; i < esize; i++)
         ensemble[i] -= d;
      return *this;
   }
   UDoubleEnsemble<esize> &operator *=(const UDoubleEnsemble<esize>& ud)
   {
      sources.check_epoch(epoch);
      sources.check_epoch(ud.epoch);

      for (unsigned int i = 0; i < esize; i++)
         ensemble[i] *= ud.ensemble[i];
      return *this;
   }
   UDoubleEnsemble<esize> &operator *=(const double& d)
   {
      for (unsigned int i = 0; i < esize; i++)
         ensemble[i] *= d;
      return *this;
   }
   UDoubleEnsemble<esize> &operator /=(const UDoubleEnsemble<esize>& ud)
   {
      sources.check_epoch(epoch);
      sources.check_epoch(ud.epoch);

      for (unsigned int i = 0; i < esize; i++)
         ensemble[i] /= ud.ensemble[i];
      return *this;
   }
   UDoubleEnsemble<esize> &operator /=(const double& d)
   {
      for (unsigned int i = 0; i < esize; i++)
         ensemble[i] /= d;
      return *this;
   }
// future direction: add procedures to make persistant
   friend ostream& operator <<(ostream &os, const UDoubleEnsemble<esize> &ud)
   {
      double mean, sigma, skew, kurtosis, m5;
      moments(ud.ensemble, mean, sigma, skew, kurtosis, m5, esize);
      uncertain_print(mean, sigma, os);

      if (sigma != 0.0)
      {
         int original_precision = os.precision();
         long original_format = os.flags(ios::showpoint);
         os << setprecision(2)
            << " [" << skew << " : " << kurtosis << " : " << m5 << "]"
            << setprecision(original_precision);
         os.flags(original_format);
      }
      return os;
   }
   friend istream& operator >>(istream &is, UDoubleEnsemble<esize> &ud)
   {
      double mean, sigma;
      uncertain_read(mean, sigma, is);
      ud = UDoubleEnsemble<esize>(mean, sigma);
      return is;
   }
   #define UDoubleEnsemblefunc1(func) \
      UDoubleEnsemble<esize> func(UDoubleEnsemble<esize> arg) \
      { \
         for (unsigned int i = 0; i < esize; i++) \
            arg.ensemble[i] = func(arg.ensemble[i]); \
         return arg; \
      }
   #define UDoubleEnsemblefunc2(func) \
      UDoubleEnsemble<esize> func(const UDoubleEnsemble<esize>& arg1, \
                                  const UDoubleEnsemble<esize>& arg2) \
      { \
         UDoubleEnsemble<esize> retval(arg1); \
         for (unsigned int i = 0; i < esize; i++) \
            retval.ensemble[i] = func(arg1.ensemble[i], arg2.ensemble[i]); \
         return retval; \
      }
   friend UDoubleEnsemblefunc1(sqrt)
   friend UDoubleEnsemblefunc1(sin)
   friend UDoubleEnsemblefunc1(cos)
   friend UDoubleEnsemblefunc1(tan)
   friend UDoubleEnsemblefunc1(asin)
   friend UDoubleEnsemblefunc1(acos)
   friend UDoubleEnsemblefunc1(atan)
   friend UDoubleEnsemblefunc2(atan2)
   friend UDoubleEnsemblefunc1(ceil)
   friend UDoubleEnsemblefunc1(floor)
   friend UDoubleEnsemblefunc1(fabs)
   friend UDoubleEnsemblefunc2(fmod)
   friend UDoubleEnsemblefunc1(exp)
   friend UDoubleEnsemblefunc1(log)
   friend UDoubleEnsemblefunc1(log10)
   friend UDoubleEnsemblefunc1(sinh)
   friend UDoubleEnsemblefunc1(cosh)
   friend UDoubleEnsemblefunc1(tanh)
   friend UDoubleEnsemblefunc2(pow)
   friend UDoubleEnsemble<esize> ldexp(UDoubleEnsemble<esize> arg,
                                       const int intarg)
   {
      for (unsigned int i = 0; i < esize; i++)
         arg.ensemble[i] = ldexp(arg.ensemble[i], intarg);
      return arg;
   }
   friend UDoubleEnsemble<esize> frexp(UDoubleEnsemble<esize> arg, int * intarg)
   {
      // use library frexp on mean to get value of return in second arg
      frexp(arg.mean(), intarg);
      for (unsigned int i = 0; i < esize; i++)
      {
         int tempint;  // ignore return in second arg in loop
         arg.ensemble[i] = frexp(arg.ensemble[i], &tempint);
      }
      return arg;
   }
   friend UDoubleEnsemble<esize> modf(UDoubleEnsemble<esize> arg,
                                      double * dblarg)
   {
      // use library modf on mean to get value of return in second arg
      modf(arg.mean(), dblarg);
      for (unsigned int i = 0; i < esize; i++)
      {
         double tempdbl;  // ignore return in second arg in loop
         arg.ensemble[i] = modf(arg.ensemble[i], &tempdbl);
      }
      return arg;
   }


   double mean(void) const
   {
      double sum = 0.0;

      for (unsigned int i = 0; i < esize; i++)
         sum += ensemble[i];
      return sum / esize;
   }
   double deviation(void) const
   {
      unsigned int i;
      double diff, sum_2_diff = 0.0; // watch overflow!
      double value = this->mean();

      for (i = 0; i < esize; i++) {
         diff = ensemble[i] - value;
         sum_2_diff += diff * diff;
      }
      return sqrt(sum_2_diff / esize);
   }
   static void new_epoch(void) { sources.new_epoch(); }
   void print_uncertain_sources(ostream &os = cout)
   {
      if (deviation() == 0.0)
         os << "No uncertainty";
      else
      {
         double unaccounted_uncertainty = 1.0;
         for (unsigned int i = 0; i < sources.get_num_sources(); i++)
         {
            double unc_portion = this->correlation(src_ensemble[i]);
            unc_portion *= unc_portion;
            unaccounted_uncertainty -= unc_portion;
            os << sources.get_source_name(i) << ": "
               << int_percent(unc_portion) << "%" << endl;
         }
         os << "other: " << int_percent(unaccounted_uncertainty) << "%" << endl;
      }
      os << endl;
   }
   double correlation(const UDoubleEnsemble<esize>& ud,
                      const unsigned int offset = 0) const
   {
      unsigned int i;
      double diff, diff_ud;
      double value = this->mean();
      double ud_value = ud.mean();
      // watch overflow!
      double sum_2_diff = 0.0, sum_2_diff_ud = 0.0, sum_prod_diff = 0.0;

      for (i = 0; i < esize; i++) {
         diff = ensemble[i] - value;
         sum_2_diff += diff * diff;
         int j = (i + offset) % esize;
         diff_ud = ud.ensemble[j] - ud_value;
         sum_2_diff_ud += diff_ud * diff_ud;
         sum_prod_diff += diff * diff_ud;
      }
      if (!sum_2_diff)
         return 0.0;
      if (!sum_2_diff_ud)
         return 0.0;
      if (!sum_prod_diff)
         return 0.0;
      return sum_prod_diff / sqrt(sum_2_diff * sum_2_diff_ud);
   }
   double correlation(const double ens[esize],
                      const unsigned int offset = 0) const
   {
      unsigned int i;
      double diff, diff_ud;
      double value = this->mean();
      double ud_value = 0.0;
      // watch overflow!
      double sum_2_diff = 0, sum_2_diff_ud = 0, sum_prod_diff = 0;

      for (i = 0; i < esize; i++)
         ud_value += ens[i];
      ud_value /= esize;
      for (i = 0; i < esize; i++) {
         diff = ensemble[i] - value;
         sum_2_diff += diff * diff;
         int j = (i + offset) % esize;
         diff_ud = ens[j] - ud_value;
         sum_2_diff_ud += diff_ud * diff_ud;
         sum_prod_diff += diff * diff_ud;
      }
      if (!sum_2_diff)
         return 0.0;
      if (!sum_2_diff_ud)
         return 0.0;
      if (!sum_prod_diff)
         return 0.0;
      return sum_prod_diff / sqrt(sum_2_diff * sum_2_diff_ud);
   }
   // future direction: add function that gives a description
   // future direction: allow superimposing histograms
   void print_histogram(ostream& os = cout) const
   {
      // centered bins for each 0.5 sigmas from -4 sigmas to +4 sigmas
      // outliers go in the outer bins
      int bin[17] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      unsigned int i;
      double value = this->mean();
      double sigma = this->deviation();

      if (sigma == 0.0)
      {
         os << "No histogram when no uncertainty" << endl;
         return;
      }
      for (i = 0; i < esize; i++)
      {
         double normval = (ensemble[i] - value) / sigma;
         int intval = int(floor(2.0 * normval + 0.5) + 8.0);
         if (intval < 0)
            bin[0]++;
         else if (intval > 16)
            bin[16]++;
         else
            bin[intval]++;
      }
      int binmax = 1;
      for (i = 0; i < 17; i++)
        if (bin[i] >= binmax)
           binmax = bin[i];
      int scale_divisor = 1;
      while (binmax / scale_divisor > 74)
         scale_divisor++;
      os << "Histogram:  (each * represents ";
      if (scale_divisor == 1)
         os << "1 point)";
      else
         os << scale_divisor << " points)";
      os << setiosflags(ios::showpos) << endl;
      int bin_display[17]; // number of characters displayed in each bin
      for (i = 0; i < 17; i++)
      {
         bin_display[i] = int((bin[i] + 0.5) / scale_divisor);
      }
      int first_display_bin = 0, last_display_bin = 16;
      while (bin_display[first_display_bin] == 0)
         first_display_bin++;
      while (bin_display[last_display_bin] == 0)
         last_display_bin--;
      for (i = first_display_bin; i <= last_display_bin; i++)
      {
         if (i & 1)
            os << "    | ";
         else
            os << setw(3) << (i / 2 - 4) << " + ";
         for (int j = 0; j < bin_display[i]; j++)
            os << "*";
         os << endl;
      }
      os << resetiosflags(ios::showpos) << endl;
   }
   friend UDoubleEnsemble<esize> Invoke(double (*certainfunc)(double),
                          const UDoubleEnsemble<esize>& arg)
   {
      UDoubleEnsemble<esize> retval;

      for (unsigned i = 0; i < esize; i++)
         retval.ensemble[i] = certainfunc(arg.ensemble[i]);
      return retval;
   }
   friend UDoubleEnsemble<esize> Invoke(double (*certainfunc)(double, double),
                          const UDoubleEnsemble<esize>& arg1,
                          const UDoubleEnsemble<esize>& arg2)
   {
      UDoubleEnsemble<esize> retval;

      for (unsigned i = 0; i < esize; i++)
         retval.ensemble[i] = certainfunc(arg1.ensemble[i], arg2.ensemble[i]);
      return retval;
   }
   void shuffle(void)
   {
      for (unsigned i = 0; i < esize - 1; i++)
      {
         unsigned j = i + (unsigned long)rand() % (esize - i);
         if (j != i)
         {
            double temp = ensemble[i];
            ensemble[i] = ensemble[j];
            ensemble[j] = temp;
         }
      }
   }
   // This function moves points a little bit so the first 5 moments all
   // get measured at precisely the expected values.
   friend void PerfectEnsemble(double * ensemble, const unsigned int ens_size)
   {
         unsigned long i;
         double value, sigma, skew, kurtosis, m5;

         double * test_ensemble = new double[ens_size];
         for (int j = 0; j < 3; j++)
         {
            moments(ensemble, value, sigma, skew, kurtosis, m5, ens_size);

            for (i = 0; i < ens_size; i++)
               ensemble[i] -= value;
            for (i = 0; i < ens_size; i++)
               ensemble[i] /= sigma;
            // future work: improve kurtosis correction
            double kurtfact = 0.045;
            for (unsigned int k = 0; k < 5; k++)
            {
               for (i = 0; i < ens_size; i++)
                  test_ensemble[i] = ensemble[i]
                                     - kurtfact * kurtosis * ensemble[i]
                                       * ensemble[i] * ensemble[i];
               double test_value, test_sigma, test_skew, test_kurtosis;
               double test_m5;
               moments(test_ensemble, test_value, test_sigma, test_skew,
                       test_kurtosis, test_m5, ens_size);
               kurtfact /= 1 - test_kurtosis / kurtosis;
            }
            for (i = 0; i < ens_size; i++)
               ensemble[i] -= kurtosis * kurtfact * ensemble[i]
                                    * ensemble[i] * ensemble[i];

         }
         delete [] test_ensemble;
         moments(ensemble, value, sigma, skew, kurtosis, m5, ens_size);
         for (i = 0; i < ens_size; i++)
            ensemble[i] -= value;
         for (i = 0; i < ens_size; i++)
            ensemble[i] /= sigma;
   }
   // figure the moments (sigma, skew, kurtosis, & 5th moment) from an
   // ensemble given the mean.
   friend void moments_fixed_mean(double const * const ensemble,
                               const double& mean, double& sigma, double& skew,
                               double& kurtosis, double& m5,
                               const unsigned int ens_size)
   {
      unsigned int i;
      double gddiff2 = 0.0, gddiff3 = 0.0, gddiff4 = 0.0, gddiff5 = 0.0;
      double * gddiff = new double[ens_size];

      // Sorting ensemble by absolute value first increases accuracy
      // but at a cost to performance.  Sorting is O(nlog(n)), whereas
      // the rest of what's being done here is O(n) in ens_size.  Operationally
      // this routine should be called mostly only when there is output,
      // so the cost should be okay.
      for (i = 0; i < ens_size; i++)
         gddiff[i] = ensemble[i] - mean;
      qsort(gddiff, ens_size, sizeof(double), abs_double_compare);
      for (i = 0; i < ens_size; i++)
      {
         gddiff2 += gddiff[i] * gddiff[i];
         gddiff3 += gddiff[i] * gddiff[i] * gddiff[i];
         gddiff4 += gddiff[i] * gddiff[i] * gddiff[i] * gddiff[i];
         gddiff5 += gddiff[i] * gddiff[i] * gddiff[i] * gddiff[i] * gddiff[i];
      }
      double var = gddiff2 / ens_size;
      sigma = sqrt(var);
      skew = gddiff3 / (var * sigma * ens_size);
      kurtosis = gddiff4 / (var * var * ens_size) - 3;
      m5 = gddiff5 / (var * var * sigma * ens_size);
      // cerr << mean << " +/- " << sigma << " [s: " << skew << ", k: ";
      // cerr << kurtosis << " m5: " << m5 << "]" << endl;
      delete [] gddiff;
   }
   // figure the moments (mean, sigma, skew, kurtosis, & 5th moment) from an
   // ensemble
   friend void moments(double const * const ensemble,
                    double& mean, double& sigma, double& skew,
                    double& kurtosis, double& m5, const unsigned int ens_size)
   {
      unsigned int i;
      double gdsum = 0.0;
      double * gddiff = new double[ens_size];

      for (i = 0; i < ens_size; i++)
         gdsum += ensemble[i];
      mean = gdsum / ens_size;
      // Sorting ensemble by absolute value first increases accuracy.
      // See note in moments_fixed_mean().
      for (i = 0; i < ens_size; i++)
         gddiff[i] = ensemble[i] - mean;
      qsort(gddiff, ens_size, sizeof(double), abs_double_compare);
      gdsum = 0.0;
      for (i = 0; i < ens_size; i++)
         gdsum += gddiff[i];
      mean += gdsum / ens_size;

      moments_fixed_mean(ensemble, mean, sigma, skew, kurtosis, m5, ens_size);
      delete [] gddiff;
   }
};

// Gimpel lint doesn't like it with consts, so for lint use defines
#define ens_a_size 128u
#define ens_b_size 1024u
// static const unsigned int ens_a_size = 128u;
// static const unsigned int ens_b_size = 1024u;



// These global values are used to pass the extra argument to ldexp(),
// modf(), and frexp() around the interface for my_ldexp(), etc. so
// these functions can be used as arguments for
// PropagateUncertaintiesBySlope().
static int GlobalInt;
static double GlobalDouble;

// These function uses GlobalInt & GlobalDouble (above) to pass the
// extra arguments to ldexp(), etc. around the interface for my_ldexp(), etc.
double my_ldexp(double a) { return ldexp(a, GlobalInt); }
double my_frexp(double a) { return frexp(a, &GlobalInt); }
double my_modf(double a) { return modf(a, &GlobalDouble); }


// test class has members of all of above classes
// future direction: add timing of calls to subclasses
// future direction: add method to compare to values gotten with +=,
//                   propogateby...()
class UDoubleTest
{
private:
   UDoubleMSUncorr msu;
   UDoubleMSCorr msc;
   UDoubleMSC<0> mscu;
   UDoubleMSC<1> mscc;
   UDoubleCTSA ctsa;
   UDoubleCTAA ctaa;
   UDoubleEnsemble<ens_a_size> ens_a;
   UDoubleEnsemble<ens_b_size> ens_b;

public:
   UDoubleTest(const double val = 0.0, const double unc = 0.0,
              const char * const name = "")
   : msu(val, unc), msc(val, unc), mscu(val, unc), mscc(val, unc),
     ctsa(val, unc, name), ctaa(val, unc, name),
     ens_a(val, unc, name), ens_b(val, unc, name) {}
   UDoubleTest(const UDoubleTest& ud)
   : msu(ud.msu), msc(ud.msc), mscu(ud.mscu), mscc(ud.mscc),
     ctsa(ud.ctsa), ctaa(ud.ctaa),
     ens_a(ud.ens_a), ens_b(ud.ens_b) {}
   ~UDoubleTest(void) {}

   UDoubleTest operator +(void) const
   {
      UDoubleTest retval;

      retval.msu = +msu;
      retval.msc = +msc;
      retval.mscu = +mscu;
      retval.mscc = +mscc;
      retval.ctsa = +ctsa;
      retval.ctaa = +ctaa;
      retval.ens_a = +ens_a;
      retval.ens_b = +ens_b;
      return retval;
   }
   UDoubleTest operator -(void) const
   {
      UDoubleTest retval;

      retval.msu = -msu;
      retval.msc = -msc;
      retval.mscu = -mscu;
      retval.mscc = -mscc;
      retval.ctsa = -ctsa;
      retval.ctaa = -ctaa;
      retval.ens_a = -ens_a;
      retval.ens_b = -ens_b;
      return retval;
   }
   friend UDoubleTest operator +(const UDoubleTest& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu + b.msu;
      retval.msc = a.msc + b.msc;
      retval.mscu = a.mscu + b.mscu;
      retval.mscc = a.mscc + b.mscc;
      retval.ctsa = a.ctsa + b.ctsa;
      retval.ctaa = a.ctaa + b.ctaa;
      retval.ens_a = a.ens_a + b.ens_a;
      retval.ens_b = a.ens_b + b.ens_b;
      return retval;
   }
   friend UDoubleTest operator +(const double& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a + b.msu;
      retval.msc = a + b.msc;
      retval.mscu = a + b.mscu;
      retval.mscc = a + b.mscc;
      retval.ctsa = a + b.ctsa;
      retval.ctaa = a + b.ctaa;
      retval.ens_a = a + b.ens_a;
      retval.ens_b = a + b.ens_b;
      return retval;
   }
   friend UDoubleTest operator +(const UDoubleTest& a, const double& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu + b;
      retval.msc = a.msc + b;
      retval.mscu = a.mscu + b;
      retval.mscc = a.mscc + b;
      retval.ctsa = a.ctsa + b;
      retval.ctaa = a.ctaa + b;
      retval.ens_a = a.ens_a + b;
      retval.ens_b = a.ens_b + b;
      return retval;
   }
   friend UDoubleTest operator -(const UDoubleTest& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu - b.msu;
      retval.msc = a.msc - b.msc;
      retval.mscu = a.mscu - b.mscu;
      retval.mscc = a.mscc - b.mscc;
      retval.ctsa = a.ctsa - b.ctsa;
      retval.ctaa = a.ctaa - b.ctaa;
      retval.ens_a = a.ens_a - b.ens_a;
      retval.ens_b = a.ens_b - b.ens_b;
      return retval;
   }
   friend UDoubleTest operator -(const double& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a - b.msu;
      retval.msc = a - b.msc;
      retval.mscu = a - b.mscu;
      retval.mscc = a - b.mscc;
      retval.ctsa = a - b.ctsa;
      retval.ctaa = a - b.ctaa;
      retval.ens_a = a - b.ens_a;
      retval.ens_b = a - b.ens_b;
      return retval;
   }
   friend UDoubleTest operator -(const UDoubleTest& a, const double& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu - b;
      retval.msc = a.msc - b;
      retval.mscu = a.mscu - b;
      retval.mscc = a.mscc - b;
      retval.ctsa = a.ctsa - b;
      retval.ctaa = a.ctaa - b;
      retval.ens_a = a.ens_a - b;
      retval.ens_b = a.ens_b - b;
      return retval;
   }
   friend UDoubleTest operator *(const UDoubleTest& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu * b.msu;
      retval.msc = a.msc * b.msc;
      retval.mscu = a.mscu * b.mscu;
      retval.mscc = a.mscc * b.mscc;
      retval.ctsa = a.ctsa * b.ctsa;
      retval.ctaa = a.ctaa * b.ctaa;
      retval.ens_a = a.ens_a * b.ens_a;
      retval.ens_b = a.ens_b * b.ens_b;
      return retval;
   }
   friend UDoubleTest operator *(const double& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a * b.msu;
      retval.msc = a * b.msc;
      retval.mscu = a * b.mscu;
      retval.mscc = a * b.mscc;
      retval.ctsa = a * b.ctsa;
      retval.ctaa = a * b.ctaa;
      retval.ens_a = a * b.ens_a;
      retval.ens_b = a * b.ens_b;
      return retval;
   }
   friend UDoubleTest operator *(const UDoubleTest& a, const double& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu * b;
      retval.msc = a.msc * b;
      retval.mscu = a.mscu * b;
      retval.mscc = a.mscc * b;
      retval.ctsa = a.ctsa * b;
      retval.ctaa = a.ctaa * b;
      retval.ens_a = a.ens_a * b;
      retval.ens_b = a.ens_b * b;
      return retval;
   }
   friend UDoubleTest operator /(const UDoubleTest& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu / b.msu;
      retval.msc = a.msc / b.msc;
      retval.mscu = a.mscu / b.mscu;
      retval.mscc = a.mscc / b.mscc;
      retval.ctsa = a.ctsa / b.ctsa;
      retval.ctaa = a.ctaa / b.ctaa;
      retval.ens_a = a.ens_a / b.ens_a;
      retval.ens_b = a.ens_b / b.ens_b;
      return retval;
   }
   friend UDoubleTest operator /(const double& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a / b.msu;
      retval.msc = a / b.msc;
      retval.mscu = a / b.mscu;
      retval.mscc = a / b.mscc;
      retval.ctsa = a / b.ctsa;
      retval.ctaa = a / b.ctaa;
      retval.ens_a = a / b.ens_a;
      retval.ens_b = a / b.ens_b;
      return retval;
   }
   friend UDoubleTest operator /(const UDoubleTest& a, const double& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu / b;
      retval.msc = a.msc / b;
      retval.mscu = a.mscu / b;
      retval.mscc = a.mscc / b;
      retval.ctsa = a.ctsa / b;
      retval.ctaa = a.ctaa / b;
      retval.ens_a = a.ens_a / b;
      retval.ens_b = a.ens_b / b;
      return retval;
   }
   // preincrement and predecrement operators return value after changes
   UDoubleTest operator ++(void)
   {
      ++msu;
      ++msc;
      ++mscu;
      ++mscc;
      ++ctsa;
      ++ctaa;
      ++ens_a;
      ++ens_b;
      return *this;
   }
   UDoubleTest operator --(void)
   {
      --msu;
      --msc;
      --mscu;
      --mscc;
      --ctsa;
      --ctaa;
      --ens_a;
      --ens_b;
      return *this;
   }
   // postincrement and postdecrement operators return value before changes
   UDoubleTest operator ++(int)
   {
      UDoubleTest retval(*this);

      msu++;
      msc++;
      mscu++;
      mscc++;
      ctsa++;
      ctaa++;
      ens_a++;
      ens_b++;
      return retval;
   }
   UDoubleTest operator --(int)
   {
      UDoubleTest retval(*this);

      msu--;
      msc--;
      mscu--;
      mscc--;
      ctsa--;
      ctaa--;
      ens_a--;
      ens_b--;
      return retval;
   }
   UDoubleTest &operator +=(const UDoubleTest& ud)
   {
      msu += ud.msu;
      msc += ud.msc;
      mscu += ud.mscu;
      mscc += ud.mscc;
      ctsa += ud.ctsa;
      ctaa += ud.ctaa;
      ens_a += ud.ens_a;
      ens_b += ud.ens_b;
      return *this;
   }
   UDoubleTest &operator -=(const UDoubleTest& ud)
   {
      msu -= ud.msu;
      msc -= ud.msc;
      mscu -= ud.mscu;
      mscc -= ud.mscc;
      ctsa -= ud.ctsa;
      ctaa -= ud.ctaa;
      ens_a -= ud.ens_a;
      ens_b -= ud.ens_b;
      return *this;
   }
   UDoubleTest &operator *=(const UDoubleTest& ud)
   {
      msu *= ud.msu;
      msc *= ud.msc;
      mscu *= ud.mscu;
      mscc *= ud.mscc;
      ctsa *= ud.ctsa;
      ctaa *= ud.ctaa;
      ens_a *= ud.ens_a;
      ens_b *= ud.ens_b;
      return *this;
   }
   UDoubleTest &operator /=(const UDoubleTest& ud)
   {
      msu /= ud.msu;
      msc /= ud.msc;
      mscu /= ud.mscu;
      mscc /= ud.mscc;
      ctsa /= ud.ctsa;
      ctaa /= ud.ctaa;
      ens_a /= ud.ens_a;
      ens_b /= ud.ens_b;
      return *this;
   }
   UDoubleTest &operator +=(const double& d)
   {
      msu += d;
      msc += d;
      mscu += d;
      mscc += d;
      ctsa += d;
      ctaa += d;
      ens_a += d;
      ens_b += d;
      return *this;
   }
   UDoubleTest &operator -=(const double& d)
   {
      msu -= d;
      msc -= d;
      mscu -= d;
      mscc -= d;
      ctsa -= d;
      ctaa -= d;
      ens_a -= d;
      ens_b -= d;
      return *this;
   }
   UDoubleTest &operator *=(const double& d)
   {
      msu *= d;
      msc *= d;
      mscu *= d;
      mscc *= d;
      ctsa *= d;
      ctaa *= d;
      ens_a *= d;
      ens_b *= d;
      return *this;
   }
   UDoubleTest &operator /=(const double& d)
   {
      msu /= d;
      msc /= d;
      mscu /= d;
      mscc /= d;
      ctsa /= d;
      ctaa /= d;
      ens_a /= d;
      ens_b /= d;
      return *this;
   }
   friend istrstream& operator >>(istrstream &is, UDoubleTest &ud)
   {
      streampos init_pos = is.tellg();        // remember initial position
      ios::iostate init_state = is.rdstate(); // remember initial state
      is >> ud.msu;
      streampos final_pos = is.tellg();       // remember final position
      ios::iostate final_state = is.rdstate();// remember final state

      is.seekg(init_pos);                     // return to initial position
      is.clear(init_state);                   // return to initial state

      is >> ud.msc;
      if (final_pos != is.tellg())
         cerr << "msc does not leave the input stream in the same "
                 "position as msu" << endl;
      if (final_state != is.rdstate())
         cerr << "msc does not leave the input stream in the same "
                 "state as msu" << endl;

      is.seekg(init_pos);
      is.clear(init_state);

      is >> ud.mscu;
      if (final_pos != is.tellg())
         cerr << "mscu does not leave the input stream in the same "
                 "position as msu" << endl;
      if (final_state != is.rdstate())
         cerr << "mscu does not leave the input stream in the same "
                 "state as msu" << endl;

      is.seekg(init_pos);
      is.clear(init_state);

      is >> ud.mscc;
      if (final_pos != is.tellg())
         cerr << "mscc does not leave the input stream in the same "
                 "position as msu" << endl;
      if (final_state != is.rdstate())
         cerr << "mscc does not leave the input stream in the same "
                 "state as msu" << endl;

      is.seekg(init_pos);
      is.clear(init_state);

      is >> ud.ctsa;
      if (final_pos != is.tellg())
         cerr << "ctsa does not leave the input stream in the same "
                 "position as msu" << endl;
      if (final_state != is.rdstate())
         cerr << "ctsa does not leave the input stream in the same "
                 "state as msu" << endl;

      is.seekg(init_pos);
      is.clear(init_state);

      is >> ud.ctaa;
      if (final_pos != is.tellg())
         cerr << "ctaa does not leave the input stream in the same "
                 "position as msu" << endl;
      if (final_state != is.rdstate())
         cerr << "ctaa does not leave the input stream in the same "
                 "state as msu" << endl;

      is.seekg(init_pos);
      is.clear(init_state);

      is >> ud.ens_a;
      if (final_pos != is.tellg())
         cerr << "ens_a does not leave the input stream in the same "
                 "position as msu" << endl;
      if (final_state != is.rdstate())
         cerr << "ens_a does not leave the input stream in the same "
                 "state as msu" << endl;

      is.seekg(init_pos);
      is.clear(init_state);

      is >> ud.ens_b;
      if (final_pos != is.tellg())
         cerr << "ens_b does not leave the input stream in the same "
                 "position as msu" << endl;
      if (final_state != is.rdstate())
         cerr << "ens_b does not leave the input stream in the same "
                 "state as msu" << endl;

      return is;
   }
   friend istream& operator >>(istream &is, UDoubleTest &ud)
   {
      static bool warned = false;
      if (!warned)
      {
         warned = true;
         cerr << "warning: UDoubleTest istream extractor does not test"
                 " the extractors of the" << endl
              << "constituent classes.  Use the istrstream extractor for "
                 "those tests." << endl;
      }
      double mean, sigma;
      uncertain_read(mean, sigma, is);
      ud = UDoubleTest(mean, sigma);
      return is;
   }
   void print_nonrandom_part(ostream &os) const
   {
      os << "Uncorrelated:   " << msu << endl;
      os << "Correlated:     " << msc << endl;
      os << "Better Uncorr:  " << mscu << endl;
      os << "Better Corr:    " << mscc << endl;
      os << "Corr Track Simp:" << ctsa << endl;
      os << "Corr Track Adv: " << ctaa << endl;
   }
   // future direction: add shortprint() only prints what's different
   friend ostream& operator <<(ostream &os, const UDoubleTest &ud)
   {
      ostrstream os_msu, os_msc, os_mscu, os_mscc, os_ctsa, os_ctaa;
      ostrstream os_ens_a, os_ens_b;

      os_msu << ud.msu << ends; char *str_msu = os_msu.str();
      os_msc << ud.msc << ends; char *str_msc = os_msc.str();
      os_mscu << ud.mscu << ends; char *str_mscu = os_mscu.str();
      os_mscc << ud.mscc << ends; char *str_mscc = os_mscc.str();
      os_ctsa << ud.ctsa << ends; char *str_ctsa = os_ctsa.str();
      os_ctaa << ud.ctaa << ends; char *str_ctaa = os_ctaa.str();
      os_ens_a << ud.ens_a << ends; char *str_ens_a = os_ens_a.str();
      os_ens_b << ud.ens_b << ends; char *str_ens_b = os_ens_b.str();

      size_t len = strlen(str_msu);  // will ignore higher moments in ensembles
      if (strncmp(str_msu, str_msc, len)
         || strncmp(str_msu, str_mscu, len) || strncmp(str_msu, str_mscc, len)
         || strncmp(str_msu, str_ctsa, len) || strncmp(str_msu, str_ctaa, len)
         || strncmp(str_msu, str_ens_a, len) || strncmp(str_msu, str_ens_b, len))
      {
         ud.print_nonrandom_part(os);
         os << "Ensemble<" << ens_a_size << ">: " << ud.ens_a << endl;
         os << "Ensemble<" << ens_b_size << ">: " << ud.ens_b << endl;
      }
      else
      {  // if all the same print just one
         os << ud.msc;
      }
      os_msu.rdbuf()->freeze(0); os_msc.rdbuf()->freeze(0);
      os_mscu.rdbuf()->freeze(0); os_mscc.rdbuf()->freeze(0);
      os_ctsa.rdbuf()->freeze(0); os_ctaa.rdbuf()->freeze(0);
      os_ens_a.rdbuf()->freeze(0); os_ens_b.rdbuf()->freeze(0);

      return os;
   }
   // future direction: add check for equality of ud*.ms*, ud*.cta
   friend void print_multi(ostream &os,
                           const UDoubleTest &uda,
                           const UDoubleTest &udb,
                           const UDoubleTest &udc)
   {
      uda.print_nonrandom_part(os);
      os << "Ensemble<" << ens_a_size << ">: " << uda.ens_a << endl;
      os << "                " << udb.ens_a << endl;
      os << "                " << udc.ens_a << endl;
      os << "Ensemble<" << ens_b_size << ">: " << udb.ens_b << endl;
      os << "                " << udb.ens_b << endl;
      os << "                " << udc.ens_b << endl;
   }
   #define UDoubleTestfunc1(func) \
      UDoubleTest func(UDoubleTest arg) \
      { \
         ostrstream os, alt_os; \
         ostrstream osc, alt_osc; \
         ostrstream osa, alt_osa; \
         char *str, *slope_str, *a_str, *alt_a_str; \
         UDoubleMSUncorr alt_msu =PropagateUncertaintiesBySlope(func, arg.msu);\
         arg.msu = func(arg.msu); \
         os << arg.msu << ends; \
         alt_os << alt_msu << ends; \
         if (strcmp(slope_str = alt_os.str(), str = os.str())) \
            cerr << "Warning: different values for " << #func "(): " \
                 << str << " vs. " << slope_str << endl; \
         arg.msc = func(arg.msc); \
         os.rdbuf()->freeze(0); alt_os.rdbuf()->freeze(0); \
         \
         UDoubleMSC<0> alt_mscu = PropagateUncertaintiesBySlope(func, arg.mscu);\
         arg.mscu = func(arg.mscu); \
         osc << arg.mscu << ends; \
         alt_osc << alt_mscu << ends; \
         if (strcmp(slope_str = alt_osc.str(), str = osc.str())) \
            cerr << "Warning: different values for curved " << #func "(): " \
                 << str << " vs. " << slope_str << endl; \
         arg.mscc = func(arg.mscc); \
         osc.rdbuf()->freeze(0); alt_osc.rdbuf()->freeze(0); \
         \
         arg.ctsa = func(arg.ctsa); \
         arg.ctaa = func(arg.ctaa); \
         UDoubleEnsemble<ens_a_size> alt_ens_a = Invoke(func, arg.ens_a); \
         arg.ens_a = func(arg.ens_a); \
         osa << arg.ens_a << ends; \
         alt_osa << alt_ens_a << ends; \
         if (strcmp(a_str = alt_osa.str(), alt_a_str = osa.str())) \
            cerr << "Warning: different values for ensemble<" << ens_a_size \
                 << "> " << #func "(): " \
                 << a_str << " vs. " << alt_a_str << endl; \
         arg.ens_b = func(arg.ens_b); \
         return arg; \
      }
   #define UDoubleTestfunc2(func) \
      UDoubleTest func(const UDoubleTest& arg1, const UDoubleTest& arg2) \
      { \
         UDoubleTest retval; \
         ostrstream os, alt_os; \
         ostrstream osc, alt_osc; \
         char *str, *slope_str; \
         UDoubleMSUncorr alt_msu =PropagateUncertaintiesBySlope(func, arg1.msu, \
                                                                arg2.msu);\
         retval.msu = func(arg1.msu, arg2.msu); \
         os << retval.msu << ends; \
         alt_os << alt_msu << ends; \
         if (strcmp(slope_str = alt_os.str(), str = os.str())) \
            cerr << "Warning: different values for " << #func "(): " \
                 << str << " vs. " << slope_str << endl; \
         retval.msc = func(arg1.msc, arg2.msc); \
         os.rdbuf()->freeze(0); alt_os.rdbuf()->freeze(0); \
         \
         UDoubleMSC<0> alt_mscu =PropagateUncertaintiesBySlope(func, arg1.mscu, \
                                                               arg2.mscu);\
         retval.mscu = func(arg1.mscu, arg2.mscu); \
         osc << retval.mscu << ends; \
         alt_osc << alt_mscu << ends; \
         if (strcmp(slope_str = alt_osc.str(), str = osc.str())) \
            cerr << "Warning: different values for curved " << #func "(): " \
                 << str << " vs. " << slope_str << endl; \
         retval.mscc = func(arg1.mscc, arg2.mscc); \
         osc.rdbuf()->freeze(0); alt_osc.rdbuf()->freeze(0); \
         \
         retval.ctsa = func(arg1.ctsa, arg2.ctsa); \
         retval.ctaa = func(arg1.ctaa, arg2.ctaa); \
         retval.ens_a = func(arg1.ens_a, arg2.ens_a); \
         retval.ens_b = func(arg1.ens_b, arg2.ens_b); \
         return retval; \
      }
   friend UDoubleTestfunc1(sqrt)
   friend UDoubleTestfunc1(sin)
   friend UDoubleTestfunc1(cos)
   friend UDoubleTestfunc1(tan)
   friend UDoubleTestfunc1(asin)
   friend UDoubleTestfunc1(acos)
   friend UDoubleTestfunc1(atan)
   friend UDoubleTestfunc2(atan2)
   friend UDoubleTestfunc1(ceil)
   friend UDoubleTestfunc1(floor)
   friend UDoubleTestfunc1(fabs)
   friend UDoubleTestfunc2(fmod)
   friend UDoubleTestfunc1(exp)
   friend UDoubleTestfunc1(log)
   friend UDoubleTestfunc1(log10)
   friend UDoubleTestfunc1(sinh)
   friend UDoubleTestfunc1(cosh)
   friend UDoubleTestfunc1(tanh)
   friend UDoubleTestfunc2(pow)
   friend UDoubleTest ldexp(UDoubleTest arg, const int intarg)
   {
      ostrstream os, alt_os;
      ostrstream osc, alt_osc;
      char *str, *slope_str;

      GlobalInt = intarg;
      UDoubleMSUncorr alt_msu =PropagateUncertaintiesBySlope(my_ldexp, arg.msu);
      arg.msu = ldexp(arg.msu, intarg);
      os << arg.msu << ends;
      alt_os << alt_msu << ends;
      if (strcmp(slope_str = alt_os.str(), str = os.str()))
         cerr << "Warning: different values for ldexp(): "
              << str << " vs. " << slope_str << endl;
      arg.msc = ldexp(arg.msc, intarg);
      os.rdbuf()->freeze(0); alt_os.rdbuf()->freeze(0);

      UDoubleMSC<0> alt_mscu =PropagateUncertaintiesBySlope(my_ldexp, arg.mscu);
      arg.mscu = ldexp(arg.mscu, intarg);
      osc << arg.mscu << ends;
      alt_osc << alt_mscu << ends;
      if (strcmp(slope_str = alt_osc.str(), str = osc.str()))
         cerr << "Warning: different values for curved ldexp(): "
              << str << " vs. " << slope_str << endl;
      arg.mscc = ldexp(arg.mscc, intarg);
      osc.rdbuf()->freeze(0); alt_osc.rdbuf()->freeze(0);

      arg.ctsa = ldexp(arg.ctsa, intarg);
      arg.ctaa = ldexp(arg.ctaa, intarg);
      arg.ens_a = ldexp(arg.ens_a, intarg);
      arg.ens_b = ldexp(arg.ens_b, intarg);
      return arg;
   }
   friend UDoubleTest frexp(UDoubleTest arg, int * intarg)
   {
      ostrstream os, alt_os;
      ostrstream osc, alt_osc;
      char *str, *slope_str;

      UDoubleMSUncorr alt_msu =PropagateUncertaintiesBySlope(my_frexp, arg.msu);
      alt_os << alt_msu << " (" << GlobalInt << ")" << ends;
      arg.msu = frexp(arg.msu, intarg);
      os << arg.msu << " (" << *intarg << ")" << ends;
      if (strcmp(slope_str = alt_os.str(), str = os.str()))
         cerr << "Warning: different values for frexp(): "
              << str << " vs. " << slope_str << endl;
      arg.msc = frexp(arg.msc, intarg);
      os.rdbuf()->freeze(0); alt_os.rdbuf()->freeze(0);

      UDoubleMSC<0> alt_mscu =PropagateUncertaintiesBySlope(my_frexp, arg.mscu);
      alt_osc << alt_mscu << " (" << GlobalInt << ")" << ends;
      arg.mscu = frexp(arg.mscu, intarg);
      osc << arg.mscu << " (" << *intarg << ")" << ends;
      if (strcmp(slope_str = alt_osc.str(), str = osc.str()))
         cerr << "Warning: different values for curved frexp(): "
              << str << " vs. " << slope_str << endl;
      arg.mscc = frexp(arg.mscc, intarg);
      osc.rdbuf()->freeze(0); alt_osc.rdbuf()->freeze(0);

      arg.ctsa = frexp(arg.ctsa, intarg);
      arg.ctaa = frexp(arg.ctaa, intarg);
      arg.ens_a = frexp(arg.ens_a, intarg);
      arg.ens_b = frexp(arg.ens_b, intarg);
      return arg;
   }
   friend UDoubleTest modf(UDoubleTest arg, double * dblarg)
   {
      ostrstream os, alt_os;
      ostrstream osc, alt_osc;
      char *str, *slope_str;

      UDoubleMSUncorr alt_msu =PropagateUncertaintiesBySlope(my_modf, arg.msu);
      alt_os << alt_msu << " (" << GlobalDouble << ")" << ends;
      arg.msu = modf(arg.msu, dblarg);
      os << arg.msu << " (" << *dblarg << ")" << ends;
      if (strcmp(slope_str = alt_os.str(), str = os.str()))
         cerr << "Warning: different values for modf(): "
              << str << " vs. " << slope_str << endl;
      arg.msc = modf(arg.msc, dblarg);
      os.rdbuf()->freeze(0); alt_os.rdbuf()->freeze(0);

      UDoubleMSC<0> alt_mscu =PropagateUncertaintiesBySlope(my_modf, arg.mscu);
      alt_osc << alt_mscu << " (" << GlobalDouble << ")" << ends;
      arg.mscu = modf(arg.mscu, dblarg);
      osc << arg.mscu << " (" << *dblarg << ")" << ends;
      if (strcmp(slope_str = alt_osc.str(), str = osc.str()))
         cerr << "Warning: different values for curved modf(): "
              << str << " vs. " << slope_str << endl;
      arg.mscc = modf(arg.mscc, dblarg);
      osc.rdbuf()->freeze(0); alt_osc.rdbuf()->freeze(0);

      arg.ctsa = modf(arg.ctsa, dblarg);
      arg.ctaa = modf(arg.ctaa, dblarg);
      arg.ens_a = modf(arg.ens_a, dblarg);
      arg.ens_b = modf(arg.ens_b, dblarg);
      return arg;
   }

   static void new_epoch(void)
   {
      UDoubleCTSA::new_epoch();
      UDoubleCTAA::new_epoch();
      UDoubleEnsemble<ens_a_size>::new_epoch();
      UDoubleEnsemble<ens_b_size>::new_epoch();
   }
   void print_uncertain_sources(ostream &os = cout)
   {
      os << "Sources of uncertainty:" << endl;
      os << "Corr Tracking Simple:  " << endl;
      ctsa.print_uncertain_sources(os);
      os << "Corr Tracking Advanced:  " << endl;
      ctaa.print_uncertain_sources(os);
      os << "Ensemble<" << ens_a_size << ">: " << endl;
      ens_a.print_uncertain_sources(os);
      os << "Ensemble<" << ens_b_size << ">: " << endl;
      ens_b.print_uncertain_sources(os);
      os << endl;
   }
   void print_correlation(const UDoubleTest ud, ostream &os = cout) const
   {
      os << "Measurements of correlation:" << endl;
      os << "Ensemble<" << ens_a_size << ">: ";
      os << ens_a.correlation(ud.ens_a) << endl;
      os << "Ensemble<" << ens_b_size << ">: ";
      os << ens_b.correlation(ud.ens_b) << endl;
      os << endl;
   }
};
UncertainSourceSet UDoubleEnsemble<ens_a_size>::sources("Small Ensemble");
double UDoubleEnsemble<ens_a_size>::src_ensemble[MAX_UNC_ELEMENTS][ens_a_size];
UncertainSourceSet UDoubleEnsemble<ens_b_size>::sources("Large Ensemble");
double UDoubleEnsemble<ens_b_size>::src_ensemble[MAX_UNC_ELEMENTS][ens_b_size];

class UDoubleInit {
private:
   static unsigned short count; // # of UDoubleInit objects existing
public:
   UDoubleInit(void)
   {
      if (count++ == 0) {
         // if time() fails, it returns -1.  In that case we just go
         // ahead and seed the random # generator with (unsigned int)-1
         // for lack of anything better.
         srand((unsigned int)time(0));
      }
   }
   ~UDoubleInit(void)
   {
      if (--count == 0) {
         // nothing to be done here now
      }
   }
};

static UDoubleInit udi; // forces construction in every including file
unsigned short UDoubleInit::count;

