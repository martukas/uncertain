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

// udmstest.cc: This file demonstrates the use of the classes in udms.h:
// a class for simple propagation of Uncertainties according to a pure
// Gaussian model.

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
#include <math.h>
#include <string.h>
#include <iostream.h>
//lint -restore
#include "udms.h"

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


// test class has UDoubleUncorr and UDoubleCorr members
class UDoubleTest
{
private:
   UDoubleUncorr msu;
   UDoubleCorr msc;

public:
   UDoubleTest(const double val = 0.0, const double unc = 0.0)
   : msu(val, unc), msc(val, unc) {}
   UDoubleTest(const UDoubleTest& ud)
   : msu(ud.msu), msc(ud.msc) {}

   UDoubleTest operator +(void) const
   {
      UDoubleTest retval;

      retval.msu = +msu;
      retval.msc = +msc;
      return retval;
   }
   UDoubleTest operator -(void) const
   {
      UDoubleTest retval;

      retval.msu = -msu;
      retval.msc = -msc;
      return retval;
   }
   friend UDoubleTest operator +(const UDoubleTest& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu + b.msu;
      retval.msc = a.msc + b.msc;
      return retval;
   }
   friend UDoubleTest operator +(const double& a, const UDoubleTest& b)
   {
      UDoubleTest retval;
      retval.msu = a + b.msu;
      retval.msc = a + b.msc;
      return retval;
   }
   friend UDoubleTest operator +(const UDoubleTest& a, const double& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu + b;
      retval.msc = a.msc + b;
      return retval;
   }
   friend UDoubleTest operator -(const UDoubleTest& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu - b.msu;
      retval.msc = a.msc - b.msc;
       return retval;
   }
   friend UDoubleTest operator -(const double& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a - b.msu;
      retval.msc = a - b.msc;
      return retval;
   }
   friend UDoubleTest operator -(const UDoubleTest& a, const double& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu - b;
      retval.msc = a.msc - b;
      return retval;
   }
   friend UDoubleTest operator *(const UDoubleTest& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu * b.msu;
      retval.msc = a.msc * b.msc;
       return retval;
   }
   friend UDoubleTest operator *(const double& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a * b.msu;
      retval.msc = a * b.msc;
      return retval;
   }
   friend UDoubleTest operator *(const UDoubleTest& a, const double& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu * b;
      retval.msc = a.msc * b;
      return retval;
   }
   friend UDoubleTest operator /(const UDoubleTest& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu / b.msu;
      retval.msc = a.msc / b.msc;
      return retval;
   }
   friend UDoubleTest operator /(const double& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a / b.msu;
      retval.msc = a / b.msc;

      return retval;
   }
   friend UDoubleTest operator /(const UDoubleTest& a, const double& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu / b;
      retval.msc = a.msc / b;
      return retval;
   }
   // preincrement and predecrement operators return value after changes
   UDoubleTest operator ++(void)
   {
      ++msu;
      ++msc;
      return *this;
   }
   UDoubleTest operator --(void)
   {
      --msu;
      --msc;
       return *this;
   }
   // postincrement and postdecrement operators return value before changes
   UDoubleTest operator ++(int)
   {
      UDoubleTest retval(*this);

      msu++;
      msc++;
      return retval;
   }
   UDoubleTest operator --(int)
   {
      UDoubleTest retval(*this);

      msu--;
      msc--;
      return retval;
   }
   UDoubleTest &operator +=(const UDoubleTest& ud)
   {
      msu += ud.msu;
      msc += ud.msc;
      return *this;
   }
   UDoubleTest &operator -=(const UDoubleTest& ud)
   {
      msu -= ud.msu;
      msc -= ud.msc;
       return *this;
   }
   UDoubleTest &operator *=(const UDoubleTest& ud)
   {
      msu *= ud.msu;
      msc *= ud.msc;
      return *this;
   }
   UDoubleTest &operator /=(const UDoubleTest& ud)
   {
      msu /= ud.msu;
      msc /= ud.msc;
      return *this;
   }
   UDoubleTest &operator +=(const double& d)
   {
      msu += d;
      msc += d;
      return *this;
   }
   UDoubleTest &operator -=(const double& d)
   {
      msu -= d;
      msc -= d;
      return *this;
   }
   UDoubleTest &operator *=(const double& d)
   {
      msu *= d;
      msc *= d;
      return *this;
   }
   UDoubleTest &operator /=(const double& d)
   {
      msu /= d;
      msc /= d;
      return *this;
   }
   friend ostream& operator <<(ostream &os, const UDoubleTest &ud)
   {
      ostrstream osu, osc;
      char *up, *cp;
      osu << ud.msu << ends;
      osc << ud.msc << ends;
      if (strcmp(up = osu.str(), cp = osc.str()))
         os << "Uncorrelated: " << up << "  Correlated: " << cp;
      else
         os << up;
      osu.freeze(0);
      osc.freeze(0);
      return os;
   }
   #define UDoubleTestfunc1(func) \
      UDoubleTest func(UDoubleTest arg) \
      { \
         ostrstream os, alt_os; \
         char *str, *slope_str; \
         UDoubleUncorr alt_msu =PropagateUncertaintiesBySlope(func, arg.msu);\
         arg.msu = func(arg.msu); \
         os << arg.msu << ends; \
         alt_os << alt_msu << ends; \
         if (strcmp(slope_str = alt_os.str(), str = os.str())) \
            cerr << "Warning: different values for " << #func "(): " \
                 << str << " vs. " << slope_str << endl; \
         arg.msc = func(arg.msc); \
         os.freeze(0); alt_os.freeze(0); \
         return arg; \
      }
   #define UDoubleTestfunc2(func) \
      UDoubleTest func(const UDoubleTest& arg1, const UDoubleTest& arg2) \
      { \
         UDoubleTest retval; \
         ostrstream os, alt_os; \
         char *str, *slope_str; \
         UDoubleUncorr alt_msu = PropagateUncertaintiesBySlope(func, \
                                                  arg1.msu, arg2.msu); \
         retval.msu = func(arg1.msu, arg2.msu); \
         os << retval.msu << ends; \
         alt_os << alt_msu << ends; \
         if (strcmp(slope_str = alt_os.str(), str = os.str())) \
            cerr << "Warning: different values for " << #func << "(): " \
                 << str << " vs. " << slope_str << endl; \
         retval.msc = func(arg1.msc, arg2.msc); \
         os.freeze(0); alt_os.freeze(0); \
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
      char *str, *slope_str;
      GlobalInt = intarg;
      UDoubleUncorr alt_msu =PropagateUncertaintiesBySlope(my_ldexp, arg.msu);
      arg.msu = ldexp(arg.msu, intarg);
      os << arg.msu << ends;
      alt_os << alt_msu << ends;
      if (strcmp(slope_str = alt_os.str(), str = os.str()))
         cerr << "Warning: different values for " << "ldexp(): "
              << str << " vs. " << slope_str << endl;
      arg.msc = ldexp(arg.msc, intarg);
      os.freeze(0); alt_os.freeze(0);
      return arg;
   }
   friend UDoubleTest frexp(UDoubleTest arg, int *intarg)
   {
      ostrstream os, alt_os;
      char *str, *slope_str;
      UDoubleUncorr alt_msu =PropagateUncertaintiesBySlope(my_frexp, arg.msu);
      arg.msu = frexp(arg.msu, intarg);
      os << arg.msu << ends;
      alt_os << alt_msu << ends;
      if (strcmp(slope_str = alt_os.str(), str = os.str()))
         cerr << "Warning: different values for " << "frexp(): "
              << str << " vs. " << slope_str << endl;
      arg.msc = frexp(arg.msc, intarg);
      os.freeze(0); alt_os.freeze(0);
      return arg;
   }
   friend UDoubleTest modf(UDoubleTest arg, double *dblarg)
   {
      ostrstream os, alt_os;
      char *str, *slope_str;
      UDoubleUncorr alt_msu = PropagateUncertaintiesBySlope(my_modf, arg.msu);
      arg.msu = modf(arg.msu, dblarg);
      os << arg.msu << ends;
      alt_os << alt_msu << ends;
      if (strcmp(slope_str = alt_os.str(), str = os.str()))
         cerr << "Warning: different values for " << "modf(): "
              << str << " vs. " << slope_str << endl;
      arg.msc = modf(arg.msc, dblarg);
      os.freeze(0); alt_os.freeze(0);
      return arg;
   }

};

#define LOGISTIC_ITERATIONS 100

template <class T>
T logistic(T r)
{
   unsigned int i;
   T x(0.6);

   for (i = 0; i < LOGISTIC_ITERATIONS; i++)
      x = r * x * (1.0 - x);
   return x;
}

// square: just a notational convenience
inline double sqr(const double &a)
   { return a * a; }
inline UDoubleTest sqr(const UDoubleTest &a)
   { return a * a; }

void infix_check (void) {
   UDoubleTest     a, b(0.0, 0.02);
   UDoubleTest     c(5.0, 0.03), t;

   a = UDoubleTest(1.0, 0.1);

   cout << endl << endl;
   cout << "a = " << a << endl;
   cout << "b = " << b << endl;
   cout << "c = " << c << endl;
   cout << endl;

   cout << "UNARY OPERATORS" << endl << endl;
   cout << "Unary operators always give the same result for correlated "
        << "& uncorrelated:" << endl;
   cout << "+a = " << +a << endl;
   cout << "-a = " << -a << endl;
   cout << "preincrement to " << ++a << endl;
   cout << "postincrement leaves at " << a++ << endl;
   cout << "a ends as " << a << endl << endl;
   cout << "-= 2.0 restores a to " << (a -= 2.0) << endl << endl;
   
   cout << "INFIX OPERATORS" << endl << endl;

   cout << "Operations where only one operand is uncertain give the same "
           "result for" << endl
        << "correlated & uncorrelated:" << endl;
   cout << "2 + a = " << (2.0 + a) << endl;
   cout << "a + 1 = " << (a + 1.0) << endl;
   cout << "a += 3 = " << (a += 3.0) << endl;

   cout << "2 - a = " << (2.0 - a) << endl;
   cout << "a - 1 = " << (a - 1.0) << endl;
   cout << "a -= 3 = " << (a -= 3.0) << endl;

   cout << "2 * a = " << (2.0 * a) << endl;
   cout << "a * 10 = " << (a * 10.0) << endl;
   cout << "a *= 3 = " << (a *= 3.0) << endl;

   cout << "2 / a = " << (2.0 / a) << endl;
   cout << "a / 10 = " << (a / 10.0) << endl;
   cout << "a /= 3 = " << (a /= 3.0) << endl;
   cout << endl;

   cout << "When both operands are the same, the correlated answer is right:"
        << endl;
   cout << "a + a = " << (t = a + a) << endl;
   cout << "a - a = " << (t = a - a) << endl;
   cout << "a * a = " << (t = a * a) << endl;
   cout << "a / a = " << (t = a / a) << endl << endl;

   cout << "When the operands are independent, the uncorrelated answer "
        << "is right:" << endl;
   cout << "a + b = " << (t = a + b) << endl;
   cout << "a - b = " << (t = a - b) << endl;
   cout << "a * b = " << (t = a * b) << endl;
   cout << "a * c = " << (t = a * c) << endl;
   cout << "c / a = " << (t = c / a) << endl;
   cout << endl;

   cout << "For these more complicated expressions, neither the uncorrelated"
        << endl
        << "nor the correlated uncertainty is right:" << endl;
   cout << "(c - a) / (c + a) = " << ((c - a) / (c + a)) << endl
        << " (should be 0.667 +/- 0.028)" << endl << endl;
   cout << "a + a + b = " << (a + a + b) << endl
        << " (should be 2.00 +/- 0.20)" << endl;
   UDoubleCorr tc = UDoubleCorr(1.0, 0.1) + UDoubleCorr(1.0, 0.1);
   cout << "But these classes can be coerced into giving the right uncertainty "
           "in this" << endl
        << "last case by separately performing the correlated and uncorrelated "
           "parts:" << endl
        << "first \"UDoubleCorr tc = UDoubleCorr(1.0, 0.1) + "
           "UDoubleCorr(1.0, 0.1);\"" << endl
        << "gives: " << tc << ".  Then the result of the correlated addition "
           "can" << endl
        << "be converted to the equivalent uncorrelated value and added with"
        << endl
        << "\"UDoubleUncorr(tc.mean(), tc.deviation()) + "
           "UDoubleUncorr(0.0, 0.02)\"," << endl
        << "giving the desired: "
        << (UDoubleUncorr(tc.mean(), tc.deviation())
            + UDoubleUncorr(0.0, 0.02))
        << endl;
}

void io_check(void) {
   const int size = 80;
   char input[size] = "3 +/- 4 1+/-2 2.+/-1 3+/-.3";
   istrstream ibuf(input, size);
   UDoubleUncorr a, b;
   UDoubleCorr x, y;

   cout << endl << endl << "EXTRACTOR" << endl;
   cout << "From \"" << input << "\" extracted: " << endl;
   ibuf >> a;
   cout << a << "," << endl;
   ibuf >> x;
   cout << x << "," << endl;
   ibuf >> b;
   cout << b << "," << endl;
   ibuf >> y;
   cout << y << "," << endl;
}

int main (void) {


   infix_check();
   io_check();

   UDoubleTest     a, b(0.0, 0.02);
   UDoubleTest     c(5.0, 0.03), d, t;

   a = UDoubleTest(1.0, 0.1);

   cout << endl << "MATH LIBRARY FUNCTIONS" << endl << endl;

   cout << "Math library functions with just one argument give the same "
           "result for" << endl
        << "correlated & uncorrelated models.  This answer is accurate "
           "to within" << endl
        << "the limits of the Gaussian model." << endl;
   cout << "sqrt(" << c << ") = " << (t = sqrt(c)) << endl;

   cout << "sin(" << c << ") = " << (t = sin(c)) << endl;

   cout << "sin(" << b << ") = " << (t = sin(b)) << endl;

   cout << "sin(" << a << ") = " << (t = sin(a)) << endl;

   cout << "cos(" << c << ") = " << (t = cos(c)) << endl;
   cout << "cos(" << b << ") = " << (t = cos(b)) << endl;
   cout << "cos(" << (0.02+b) << ") = " << (t = cos(0.02+b)) << endl;
   cout << "cos(" << (a) << ") = " << (t = cos(a)) << endl;
   t = ldexp(c+0.5, 5);
   cout << "ldexp(" << (c+0.5) << ", " << 5 << ") = " << t << endl;
   t = ldexp(a, 3);
   cout << "ldexp(" << a << ", " << 3 << ") = " << t << endl;
   cout << endl;

   cout << "sine squared + cosine squared equals one is a basic geometric "
           "truth.  Because" << endl
        << "there is only one source of uncertainty, the correlated model "
           "works better here." << endl;
   cout << "sin^2 + cos^2 (" << c << ") = " << endl
        << " " << (cos(c) * cos(c) + sin(c) * sin(c)) << endl;
   cout << "sin^2 + cos^2 (" << a << ") = " << endl
        << " " << (cos(a) * cos(a) + sin(a) * sin(a)) << endl << endl;

   cout << "Repeated application of trigonometric and exponential "
           "functions and" << endl
        << "their inverses should return the original value.  A final "
           "division" << endl
        << "by the original value should give 1.0 +/- 0.0 but the correlated "
           "division" << endl
        << "makes the uncorrelated class give the wrong uncertainty.  The "
           "second" << endl
        << "derivative of tan(1.0) is just big enough to cause a warning "
           "from the" << endl
        << "test class that propagating the uncertainty by slope gives a "
           "slightly" << endl
        << "different uncertainty." << endl;
   cout << "(asin(sin(atan(tan(acos(cos(log(exp(a)))))))) / a) = " << endl;
   cout << " " << (asin(sin(atan(tan(acos(cos(log(exp(a)))))))) / a) << endl
        << endl;

   cout << "This test uses the relationship of log() and log10() to check "
           "these functions." << endl;
   cout << "The result should be 1.0 +/- 0.0, but the correlated division "
        << "makes the" << endl
        << "uncorrelated class give the wrong uncertainty." << endl;
   cout << "(log10(a) * log(10.0) / log(a)) = " << endl
        << " " << (log10(a + 10.0) * log(10.0) / log(a + 10.0)) << endl
        << " " << (log10(b + 2.0) * log(10.0) / log(b + 2.0)) << endl
        << " " << (log10(c) * log(10.0) / log(c)) << endl << endl;

   cout << "checking hyperbolic trig functions: should be 1.0 +/- 0.0" << endl;
   cout << "(cosh(a)*cosh(a)/sqrt((1.0+sinh(a)*sinh(a)) / "
        << "(1.0 - tanh(a)*tanh(a)))) = " << endl;
   cout << (sqr(cosh(a)) / sqrt((1.0 + sqr(sinh(a))) / (1.0 - sqr(tanh(a)))))
        << "," << endl << " "
        << (sqr(cosh(b)) / sqrt((1.0 + sqr(sinh(b))) / (1.0 - sqr(tanh(b)))))
        << "," << endl << " "
        << (sqr(cosh(c)) / sqrt((1.0 + sqr(sinh(c))) / (1.0 - sqr(tanh(c)))))
        << endl << endl;

   cout << "Math library functions with two arguments give different "
           "results for" << endl
        << "correlated & uncorrelated models.  The correlated answer is "
           "better" << endl
        << "when the two arguments have the same source of uncertainty, "
           "otherwise" << endl
        << "the uncorrelated answer is better." << endl;
   cout << "pow(" << c << ", " << 2.0*a << ") = "
        << endl << " " << pow(c, 2.0*a) << endl;
   cout << "pow(" << c << ", " << c << ") = "
        << endl << " " << pow(c, c) << endl << endl;
   cout << "atan2(" << c << ", " << c << ") = "
        << endl << " " << atan2(c, c) << endl;
   cout << "atan2(" << a << ", " << c << ") = "
        << endl << " "  << atan2(a, c) << endl << endl;

   cout << "Checking pow() in terms of sqrt() and exp(): "
           "(should be 1.0 +/- 0.0)"
        << endl;
   cout << "(sqrt(c) / pow(c, 0.5)) =" << endl
        << (sqrt(c) / pow(c, 0.5)) << endl;
   cout << "(exp(c) / pow(exp(1.0), c)) =" << endl
        << (exp(c) / pow(exp(1.0), c)) << endl << endl;

   cout << "Checking atan2() in terms of atan(): (should be 1.0 +/- 0.0)"
        << endl;
   cout << "(atan(a) / atan2(a, 1.0)) = " << endl
        << " " << (atan(a) / atan2(a, 1.0)) << endl << endl;

   cout << "The library functions ceil(), floor(), fabs(), modf(), "
           "frexp(), and" << endl
        << "fmod() all exist primarily for their discontinuities.  All "
           "give bad" << endl
        << "results (and warnings from the test class) when used near a "
           "discontinuity." << endl
        << "Extreme care must be taken in using these functions with "
           "uncertain arguments." << endl;
   t = ceil(c);
   cout << "ceil(" << c << ") = " << endl << " " << t << endl;
   t = ceil(c + 0.5);
   cout << "ceil(" << (c + 0.5) << ") = " << t << endl;
   t = floor(c);
   cout << "floor(" << c << ") = " << t << endl;
   t = floor(c+0.5);
   cout << "floor(" << (c + 0.5) << ") = " << t << endl;
   t = fabs(c);
   cout << "fabs(" << c << ") = " << t << endl;
   t = fabs(b);
   cout << "fabs(" << b << ") = " << t << endl;
   t = fabs(b+0.02);
   cout << "fabs(" << (b + 0.02) << ") = " << t << endl;
   t = fmod(c+0.5, 1.0);
   cout << "fmod(" << (c+0.5) << ", " << 1.0 << ") = " << t << endl;
   t = fmod(5.5, 0.9 + 0.1*a);
   cout << "fmod(" << 5.5 << ", " << (0.9 + 0.1*a) << ") = " << t << endl;
   double dbl_temp;
   t = modf(c, &dbl_temp);
   cout << "modf(" << c << ", x) = " << t << endl;
   t = modf(c+0.5, &dbl_temp);
   cout << "modf(" << (c+0.5) << ", x) = " << t << endl;
   int int_temp;
   t = frexp(c, &int_temp);
   cout << "frexp(" << c << ", i) = " << t << endl;
   t = frexp(c+0.5, &int_temp);
   cout << "frexp(" << (c+0.5) << ", i) = " << t << endl;
   cout << endl;

   cout << "Repeated application of the logistic function "
           "(x = r * x * (1.0 - x))" << endl
        << "Can be chaotic or not depending on 'r'.  Because there is only "
           "one source" << endl
        << "of uncertainty, the uncorrelated model gives nonsense.  The "
           "correlated class" << endl
        << "does not accurately model the true distribution because the "
           "Gaussian model" << endl
        << "breaks down, but does show much greater uncertainty in the "
           "chaotic regime." << endl
        << "See _Chaos:_Making_a_New_Science_ by James Gleick pp. 70-78 "
           "for more on" << endl
        << "the logistic function." << endl;
   cout << LOGISTIC_ITERATIONS << " iterations of logistic function:" << endl;
   a = UDoubleTest(2.9, 0.000001);
   t = logistic(a);
   cout << "before 1st bifurcation:" << endl << t << endl << endl;

   b = UDoubleTest(3.1, 0.000001);
   t = logistic(b);
   cout << "after 1st bifurcation:" << endl << t << endl << endl;

   c = UDoubleTest(3.7, 0.000001);
   t = logistic(c);
   cout << "chaotic regime:" << endl << t << endl << endl;

   d = UDoubleTest(3.84, 0.000001);
   t = logistic(d);
   cout << "island of order in chaos:" << endl << t << endl << endl;

   return 0;
}
