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

// These lint flags are for Gimpel Flexi-Lint
//lint -save
//lint -e762 -e578 -e1917 -e1918 -e1510 -e1712 -e1730 -e655 -e1909 -e578
//lint -e1511 -e1907 -e36 -e1908 -e763 -e534
#include <stdlib.h>
#include <math.h>
#include <iostream.h>
#include <strstream.h>
#include <iomanip.h>
//lint -restore

// This function takes the square root of the sum of the squares of
// numbers, which is equal to the length of the hypotenuse of a right
// triangle if the two arguments are the lengths of the legs.
inline double hypot(const double &a, const double &b)
   { return sqrt(a * a + b * b); }

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
template <int is_correlated>
class UDouble
{
private:
   double value;
   double uncertainty;

public:
   // This is the default conversion from type double
   UDouble(const double val = 0.0, const double unc = 0.0)
      : value(val), uncertainty(unc)
   {
      if ((unc < 0.0) && !is_correlated)
      {
         cerr << "Error: negative uncertainty: " << unc << endl;
         exit(EXIT_FAILURE);
      }
   }
   UDouble(const UDouble& ud)
      : value(ud.value), uncertainty(ud.uncertainty) {}

   ~UDouble(void) {}

   UDouble<is_correlated> operator +(void) const
      {return *this;}
   UDouble<is_correlated> operator -(void) const
   {
      if (is_correlated)
         return UDouble<is_correlated>(-value, -uncertainty);
      else
         return UDouble<is_correlated>(-value, uncertainty);
   }
   friend UDouble<is_correlated> operator +(UDouble<is_correlated> a,
                                      const UDouble<is_correlated>& b)
      { return a += b; }
   friend UDouble<is_correlated> operator +(UDouble<is_correlated> a,
                                      const double& b)
      { return a += b; }
   friend UDouble<is_correlated> operator +(const double& b,
                                      UDouble<is_correlated> a)
      { return a += b; }
   friend UDouble<is_correlated> operator -(UDouble<is_correlated> a,
                                      const UDouble<is_correlated>& b)
      { return a -= b; }
   friend UDouble<is_correlated> operator -(UDouble<is_correlated> a,
                                      const double& b)
      { return a -= b; }
   friend UDouble<is_correlated> operator -(const double& b,
                                      UDouble<is_correlated> a)
      { a -= b; return -a; }
   UDouble<is_correlated> operator ++(void)
      { return (*this += 1.0); }
   UDouble<is_correlated> operator --(void)
      { return (*this -= 1.0); }
   UDouble<is_correlated> operator ++(int)
   {
      UDouble<is_correlated> retval(*this);
      *this += 1.0;
      return retval;
   }
   UDouble<is_correlated> operator --(int)
   {
      UDouble<is_correlated> retval(*this);
      *this -= 1.0;
      return retval;
   }
   friend UDouble<is_correlated> operator *(UDouble<is_correlated> a,
                                      const UDouble<is_correlated>& b)
      { return a *= b; }
   friend UDouble<is_correlated> operator *(UDouble<is_correlated> a,
                                      const double& b)
      { return a *= b; }
   friend UDouble<is_correlated> operator *(const double& b,
                                      UDouble<is_correlated> a)
      { return a *= b; }
   friend UDouble<is_correlated> operator /(UDouble<is_correlated> a,
                                      const UDouble<is_correlated>& b)
      { return a /= b; }
   friend UDouble<is_correlated> operator /(UDouble<is_correlated> a,
                                      const double& b)
      { return a /= b; }
   friend UDouble<is_correlated> operator /(const double& b,
                                      UDouble<is_correlated> a)
   {
      UDouble<is_correlated> retval;
      retval.uncertainty = -b * a.uncertainty / (a.value * a.value);
      retval.value = b / a.value;
      return retval;
   }
   UDouble<is_correlated> &operator +=(const UDouble<is_correlated>& ud)
   {
      if (is_correlated)
         uncertainty += ud.uncertainty;
      else
         uncertainty = hypot(uncertainty, ud.uncertainty);
      value += ud.value;
      return *this;
   }
   UDouble<is_correlated> &operator +=(const double& a)
   {
      value += a;
      return *this;
   }
   UDouble<is_correlated> &operator -=(const UDouble<is_correlated>& ud)
   {
      if (is_correlated)
         uncertainty -= ud.uncertainty;
      else
         uncertainty = hypot(uncertainty, ud.uncertainty);
      value -= ud.value;
      return *this;
   }
   UDouble<is_correlated> &operator -=(const double& a)
   {
      value -= a;
      return *this;
   }
   UDouble<is_correlated> &operator *=(const UDouble<is_correlated>& ud)
   {
      if (is_correlated)
         uncertainty = uncertainty * ud.value + ud.uncertainty * value;
      else
         uncertainty = hypot(uncertainty * ud.value, ud.uncertainty * value);
      value *= ud.value;
      return *this;
   }
   UDouble<is_correlated> &operator *=(const double& a)
   {
      value *= a;
      uncertainty *= a;
      return *this;
   }
   UDouble<is_correlated> &operator /=(const UDouble<is_correlated>& ud)
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
   UDouble<is_correlated> &operator /=(const double& a)
   {
      value /= a;
      uncertainty /= a;
      return *this;
   }
   friend ostream& operator <<(ostream &os, const UDouble<is_correlated> &ud)
   {
      uncertain_print(ud.mean(), ud.deviation(), os);
      return os;
   }
   friend istream& operator >>(istream &is, UDouble<is_correlated> &ud)
   {
      double mean, sigma;
      uncertain_read(mean, sigma, is);
      ud = UDouble<is_correlated>(mean, sigma);
      return is;
   }

   // math library functions
   friend UDouble<is_correlated> ceil(UDouble<is_correlated> arg)
   {
      arg.value = ceil(arg.value);
      arg.uncertainty = 0.0;
      return arg;
   }
   friend UDouble<is_correlated> floor(UDouble<is_correlated> arg)
   {
      arg.value = floor(arg.value);
      arg.uncertainty = 0.0;
      return arg;
   }
   friend UDouble<is_correlated> fabs(UDouble<is_correlated> arg)
   {
      if (is_correlated && (arg.value < 0.0))
         arg.uncertainty *= -1.0;
      arg.value = fabs(arg.value);
      return arg;
   }
   friend UDouble<is_correlated> ldexp(UDouble<is_correlated> arg,
                                         int intarg)
   {
      if (is_correlated)
         arg.uncertainty *= ldexp(1.0, intarg);
      else
         arg.uncertainty *= fabs(ldexp(1.0, intarg));
      arg.value = ldexp(arg.value, intarg);
      return arg;
   }
   friend UDouble<is_correlated> modf(UDouble<is_correlated> arg,
                                         double *intpart)
   {
      arg.value = modf(arg.value, intpart);
      return arg;
   }
   friend UDouble<is_correlated> frexp(UDouble<is_correlated> arg,
                                         int *intarg)
   {
      arg.uncertainty *= pow(2.0, double(-*intarg));
      arg.value = frexp(arg.value, intarg);
      return arg;
   }
   friend UDouble<is_correlated> fmod(const UDouble<is_correlated>& arg1,
                                    const UDouble<is_correlated>& arg2)
   {
      UDouble<is_correlated> retval;
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
   friend UDouble<is_correlated> sqrt(UDouble<is_correlated> arg)
   {
      arg.value = sqrt(arg.value);
      if (is_correlated)
         arg.uncertainty /= 2.0 * arg.value;
      else
         arg.uncertainty /= fabs(2.0 * arg.value);
      return arg;
   }
   friend UDouble<is_correlated> sin(UDouble<is_correlated> arg)
   {
      if (is_correlated)
         arg.uncertainty *= cos(arg.value);
      else
         arg.uncertainty *= fabs(cos(arg.value));
      arg.value = sin(arg.value);
      return arg;
   }
   friend UDouble<is_correlated> cos(UDouble<is_correlated> arg)
   {
      if (is_correlated)
         arg.uncertainty *= -sin(arg.value);
      else
         arg.uncertainty *= fabs(sin(arg.value));
      arg.value = cos(arg.value);
      return arg;
   }
   friend UDouble<is_correlated> tan(UDouble<is_correlated> arg)
   {
      double costemp = cos(arg.value);
      arg.uncertainty /= costemp * costemp;
      arg.value = tan(arg.value);
      return arg;
   }
   friend UDouble<is_correlated> asin(UDouble<is_correlated> arg)
   {
      arg.uncertainty /= sqrt(1.0 - arg.value * arg.value);
      arg.value = asin(arg.value);
      return arg;
   }
   friend UDouble<is_correlated> acos(UDouble<is_correlated> arg)
   {
      if (is_correlated)
         arg.uncertainty /= -sqrt(1.0 - arg.value * arg.value);
      else
         arg.uncertainty /= sqrt(1.0 - arg.value * arg.value);
      arg.value = acos(arg.value);
      return arg;
   }
   friend UDouble<is_correlated> atan(UDouble<is_correlated> arg)
   {
      arg.uncertainty /= 1.0 + arg.value * arg.value;
      arg.value = atan(arg.value);
      return arg;
   }
   friend UDouble<is_correlated> atan2(const UDouble<is_correlated>& arg1,
                                    const UDouble<is_correlated>& arg2)
   {
      UDouble<is_correlated> retval;
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
   friend UDouble<is_correlated> exp(UDouble<is_correlated> arg)
   {
      arg.value = exp(arg.value);
      if (is_correlated)
         arg.uncertainty *= arg.value;
      else
         arg.uncertainty *= fabs(arg.value);
      return arg;
   }
   friend UDouble<is_correlated> log(UDouble<is_correlated> arg)
   {
      if (is_correlated)
        arg.uncertainty /= arg.value;
      else
        arg.uncertainty /= fabs(arg.value);
      arg.value = log(arg.value);
      return arg;
   }
   friend UDouble<is_correlated> log10(UDouble<is_correlated> arg)
   {
      if (is_correlated)
         arg.uncertainty *= 0.43429448189 / arg.value;
      else
         arg.uncertainty *= 0.43429448189 / fabs(arg.value);
      arg.value = log10(arg.value);
      return arg;
   }
   friend UDouble<is_correlated> sinh(UDouble<is_correlated> arg)
   {
      arg.uncertainty *= cosh(arg.value);
      arg.value = sinh(arg.value);
      return arg;
   }
   friend UDouble<is_correlated> cosh(UDouble<is_correlated> arg)
   {
      if (is_correlated)
         arg.uncertainty *= sinh(arg.value);
      else
         arg.uncertainty *= fabs(sinh(arg.value));
      arg.value = cosh(arg.value);
      return arg;
   }
   friend UDouble<is_correlated> tanh(UDouble<is_correlated> arg)
   {
      double coshtemp = cosh(arg.value);
      arg.uncertainty /= coshtemp * coshtemp;
      arg.value = tanh(arg.value);
      return arg;
   }
   friend UDouble<is_correlated> pow(const UDouble<is_correlated>& arg1,
                                    const UDouble<is_correlated>& arg2)
   {
      UDouble<is_correlated> retval;
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

   friend UDouble<is_correlated> PropagateUncertaintiesBySlope(
                                    double (* certain_func)(double),
                                    const UDouble<is_correlated>& arg)
   {
      UDouble<is_correlated> retval;
      double sigma_up_value, sigma_down_value;

      retval.value = certain_func(arg.value);
      sigma_up_value = certain_func(arg.value + arg.uncertainty);
      sigma_down_value = certain_func(arg.value - arg.uncertainty);
      retval.uncertainty = (sigma_up_value - sigma_down_value) * 0.5;
      if (!is_correlated)
         retval.uncertainty = fabs(retval.uncertainty);

      return retval;
   }
   friend UDouble<is_correlated> PropagateUncertaintiesBySlope(
                                   double (*certain_func)(double, double),
                                   const UDouble<is_correlated>& arg1,
                                   const UDouble<is_correlated>& arg2)
   {
      UDouble<is_correlated> retval;

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

typedef UDouble<0> UDoubleUncorr;
typedef UDouble<1> UDoubleCorr;
