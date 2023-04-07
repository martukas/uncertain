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

#include <cmath>
#include <cstring>
#include <iostream>
#include <string>

#include "UDoubleTest.hpp"

#define LOGISTIC_ITERATIONS 100

template <class T>
T logistic(T r) {
  unsigned int i;
  T x(0.6);

  for (i = 0; i < LOGISTIC_ITERATIONS; i++) x = r * x * (1.0 - x);
  return x;
}

// square: just a notational convenience
inline UDoubleTest sqr(const UDoubleTest& a) { return a * a; }

void infix_check() {
  UDoubleTest a, b(0.0, 0.02);
  UDoubleTest c(5.0, 0.03), t;

  a = UDoubleTest(1.0, 0.1);

  std::cout << std::endl << std::endl;
  std::cout << "a = " << a << std::endl;
  std::cout << "b = " << b << std::endl;
  std::cout << "c = " << c << std::endl;
  std::cout << std::endl;

  std::cout << "UNARY OPERATORS" << std::endl << std::endl;
  std::cout << "Unary operators always give the same result for correlated "
            << "& uncorrelated:" << std::endl;
  std::cout << "+a = " << +a << std::endl;
  std::cout << "-a = " << -a << std::endl;
  std::cout << "preincrement to " << ++a << std::endl;
  std::cout << "postincrement leaves at " << a++ << std::endl;
  std::cout << "a ends as " << a << std::endl << std::endl;
  std::cout << "-= 2.0 restores a to " << (a -= 2.0) << std::endl << std::endl;

  std::cout << "INFIX OPERATORS" << std::endl << std::endl;

  std::cout << "Operations where only one operand is uncertain give the same "
               "result for"
            << std::endl
            << "correlated & uncorrelated:" << std::endl;
  std::cout << "2 + a = " << (2.0 + a) << std::endl;
  std::cout << "a + 1 = " << (a + 1.0) << std::endl;
  std::cout << "a += 3 = " << (a += 3.0) << std::endl;

  std::cout << "2 - a = " << (2.0 - a) << std::endl;
  std::cout << "a - 1 = " << (a - 1.0) << std::endl;
  std::cout << "a -= 3 = " << (a -= 3.0) << std::endl;

  std::cout << "2 * a = " << (2.0 * a) << std::endl;
  std::cout << "a * 10 = " << (a * 10.0) << std::endl;
  std::cout << "a *= 3 = " << (a *= 3.0) << std::endl;

  std::cout << "2 / a = " << (2.0 / a) << std::endl;
  std::cout << "a / 10 = " << (a / 10.0) << std::endl;
  std::cout << "a /= 3 = " << (a /= 3.0) << std::endl;
  std::cout << std::endl;

  std::cout << "When both operands are the same, the correlated answer is right:" << std::endl;
  std::cout << "a + a = " << (t = a + a) << std::endl;
  std::cout << "a - a = " << (t = a - a) << std::endl;
  std::cout << "a * a = " << (t = a * a) << std::endl;
  std::cout << "a / a = " << (t = a / a) << std::endl << std::endl;

  std::cout << "When the operands are independent, the uncorrelated answer "
            << "is right:" << std::endl;
  std::cout << "a + b = " << (t = a + b) << std::endl;
  std::cout << "a - b = " << (t = a - b) << std::endl;
  std::cout << "a * b = " << (t = a * b) << std::endl;
  std::cout << "a * c = " << (t = a * c) << std::endl;
  std::cout << "c / a = " << (t = c / a) << std::endl;
  std::cout << std::endl;

  std::cout << "For these more complicated expressions, neither the uncorrelated" << std::endl
            << "nor the correlated uncertainty is right:" << std::endl;
  std::cout << "(c - a) / (c + a) = " << ((c - a) / (c + a)) << std::endl
            << " (should be 0.667 +/- 0.028)" << std::endl
            << std::endl;
  std::cout << "a + a + b = " << (a + a + b) << std::endl
            << " (should be 2.00 +/- 0.20)" << std::endl;
  uncertain::UDoubleMSCorr tc =
      uncertain::UDoubleMSCorr(1.0, 0.1) + uncertain::UDoubleMSCorr(1.0, 0.1);
  std::cout << "But these classes can be coerced into giving the right uncertainty "
               "in this"
            << std::endl
            << "last case by separately performing the correlated and uncorrelated "
               "parts:"
            << std::endl
            << "first \"UDoubleCorr tc = UDoubleCorr(1.0, 0.1) + "
               "UDoubleCorr(1.0, 0.1);\""
            << std::endl
            << "gives: " << tc
            << ".  Then the result of the correlated addition "
               "can"
            << std::endl
            << "be converted to the equivalent uncorrelated value and added with" << std::endl
            << "\"UDoubleUncorr(tc.mean(), tc.deviation()) + "
               "UDoubleUncorr(0.0, 0.02)\","
            << std::endl
            << "giving the desired: "
            << (uncertain::UDoubleMSUncorr(tc.mean(), tc.deviation()) +
                uncertain::UDoubleMSUncorr(0.0, 0.02))
            << std::endl;
}

void io_check() {
  std::string input = "3 +/- 4 1+/-2 2.+/-1 3+/-.3";
  std::stringstream ibuf;
  ibuf.str(input);
  uncertain::UDoubleMSUncorr a, b;
  uncertain::UDoubleMSCorr x, y;

  std::cout << std::endl << std::endl << "EXTRACTOR" << std::endl;
  std::cout << "From \"" << input << "\" extracted: " << std::endl;
  ibuf >> a;
  std::cout << a << "," << std::endl;
  ibuf >> x;
  std::cout << x << "," << std::endl;
  ibuf >> b;
  std::cout << b << "," << std::endl;
  ibuf >> y;
  std::cout << y << "," << std::endl;
}

int main() {
  infix_check();
  io_check();

  UDoubleTest a, b(0.0, 0.02);
  UDoubleTest c(5.0, 0.03), d, t;

  a = UDoubleTest(1.0, 0.1);

  std::cout << std::endl << "MATH LIBRARY FUNCTIONS" << std::endl << std::endl;

  std::cout << "Math library functions with just one argument give the same "
               "result for"
            << std::endl
            << "correlated & uncorrelated models.  This answer is accurate "
               "to within"
            << std::endl
            << "the limits of the Gaussian model." << std::endl;
  std::cout << "sqrt(" << c << ") = " << (t = sqrt(c)) << std::endl;

  std::cout << "sin(" << c << ") = " << (t = sin(c)) << std::endl;

  std::cout << "sin(" << b << ") = " << (t = sin(b)) << std::endl;

  std::cout << "sin(" << a << ") = " << (t = sin(a)) << std::endl;

  std::cout << "cos(" << c << ") = " << (t = cos(c)) << std::endl;
  std::cout << "cos(" << b << ") = " << (t = cos(b)) << std::endl;
  std::cout << "cos(" << (0.02 + b) << ") = " << (t = cos(0.02 + b)) << std::endl;
  std::cout << "cos(" << (a) << ") = " << (t = cos(a)) << std::endl;
  t = ldexp(c + 0.5, 5);
  std::cout << "ldexp(" << (c + 0.5) << ", " << 5 << ") = " << t << std::endl;
  t = ldexp(a, 3);
  std::cout << "ldexp(" << a << ", " << 3 << ") = " << t << std::endl;
  std::cout << std::endl;

  std::cout << "sine squared + cosine squared equals one is a basic geometric "
               "truth.  Because"
            << std::endl
            << "there is only one source of uncertainty, the correlated model "
               "works better here."
            << std::endl;
  std::cout << "sin^2 + cos^2 (" << c << ") = " << std::endl
            << " " << (cos(c) * cos(c) + sin(c) * sin(c)) << std::endl;
  std::cout << "sin^2 + cos^2 (" << a << ") = " << std::endl
            << " " << (cos(a) * cos(a) + sin(a) * sin(a)) << std::endl
            << std::endl;

  std::cout << "Repeated application of trigonometric and exponential "
               "functions and"
            << std::endl
            << "their inverses should return the original value.  A final "
               "division"
            << std::endl
            << "by the original value should give 1.0 +/- 0.0 but the correlated "
               "division"
            << std::endl
            << "makes the uncorrelated class give the wrong uncertainty.  The "
               "second"
            << std::endl
            << "derivative of tan(1.0) is just big enough to cause a warning "
               "from the"
            << std::endl
            << "test class that propagating the uncertainty by slope gives a "
               "slightly"
            << std::endl
            << "different uncertainty." << std::endl;
  std::cout << "(asin(sin(atan(tan(acos(cos(log(exp(a)))))))) / a) = " << std::endl;
  std::cout << " " << (asin(sin(atan(tan(acos(cos(log(exp(a)))))))) / a) << std::endl << std::endl;

  std::cout << "This test uses the relationship of log() and log10() to check "
               "these functions."
            << std::endl;
  std::cout << "The result should be 1.0 +/- 0.0, but the correlated division "
            << "makes the" << std::endl
            << "uncorrelated class give the wrong uncertainty." << std::endl;
  std::cout << "(log10(a) * log(10.0) / log(a)) = " << std::endl
            << " " << (log10(a + 10.0) * log(10.0) / log(a + 10.0)) << std::endl
            << " " << (log10(b + 2.0) * log(10.0) / log(b + 2.0)) << std::endl
            << " " << (log10(c) * log(10.0) / log(c)) << std::endl
            << std::endl;

  std::cout << "checking hyperbolic trig functions: should be 1.0 +/- 0.0" << std::endl;
  std::cout << "(cosh(a)*cosh(a)/sqrt((1.0+sinh(a)*sinh(a)) / "
            << "(1.0 - tanh(a)*tanh(a)))) = " << std::endl;
  std::cout << (sqr(cosh(a)) / sqrt((1.0 + sqr(sinh(a))) / (1.0 - sqr(tanh(a))))) << ","
            << std::endl
            << " " << (sqr(cosh(b)) / sqrt((1.0 + sqr(sinh(b))) / (1.0 - sqr(tanh(b))))) << ","
            << std::endl
            << " " << (sqr(cosh(c)) / sqrt((1.0 + sqr(sinh(c))) / (1.0 - sqr(tanh(c)))))
            << std::endl
            << std::endl;

  std::cout << "Math library functions with two arguments give different "
               "results for"
            << std::endl
            << "correlated & uncorrelated models.  The correlated answer is "
               "better"
            << std::endl
            << "when the two arguments have the same source of uncertainty, "
               "otherwise"
            << std::endl
            << "the uncorrelated answer is better." << std::endl;
  std::cout << "pow(" << c << ", " << 2.0 * a << ") = " << std::endl
            << " " << pow(c, 2.0 * a) << std::endl;
  std::cout << "pow(" << c << ", " << c << ") = " << std::endl
            << " " << pow(c, c) << std::endl
            << std::endl;
  std::cout << "atan2(" << c << ", " << c << ") = " << std::endl << " " << atan2(c, c) << std::endl;
  std::cout << "atan2(" << a << ", " << c << ") = " << std::endl
            << " " << atan2(a, c) << std::endl
            << std::endl;

  std::cout << "Checking pow() in terms of sqrt() and exp(): "
               "(should be 1.0 +/- 0.0)"
            << std::endl;
  std::cout << "(sqrt(c) / pow(c, 0.5)) =" << std::endl << (sqrt(c) / pow(c, 0.5)) << std::endl;
  std::cout << "(exp(c) / pow(exp(1.0), c)) =" << std::endl
            << (exp(c) / pow(exp(1.0), c)) << std::endl
            << std::endl;

  std::cout << "Checking atan2() in terms of atan(): (should be 1.0 +/- 0.0)" << std::endl;
  std::cout << "(atan(a) / atan2(a, 1.0)) = " << std::endl
            << " " << (atan(a) / atan2(a, 1.0)) << std::endl
            << std::endl;

  std::cout << "The library functions ceil(), floor(), fabs(), modf(), "
               "frexp(), and"
            << std::endl
            << "fmod() all exist primarily for their discontinuities.  All "
               "give bad"
            << std::endl
            << "results (and warnings from the test class) when used near a "
               "discontinuity."
            << std::endl
            << "Extreme care must be taken in using these functions with "
               "uncertain arguments."
            << std::endl;
  t = ceil(c);
  std::cout << "ceil(" << c << ") = " << std::endl << " " << t << std::endl;
  t = ceil(c + 0.5);
  std::cout << "ceil(" << (c + 0.5) << ") = " << t << std::endl;
  t = floor(c);
  std::cout << "floor(" << c << ") = " << t << std::endl;
  t = floor(c + 0.5);
  std::cout << "floor(" << (c + 0.5) << ") = " << t << std::endl;
  t = fabs(c);
  std::cout << "fabs(" << c << ") = " << t << std::endl;
  t = fabs(b);
  std::cout << "fabs(" << b << ") = " << t << std::endl;
  t = fabs(b + 0.02);
  std::cout << "fabs(" << (b + 0.02) << ") = " << t << std::endl;
  t = fmod(c + 0.5, 1.0);
  std::cout << "fmod(" << (c + 0.5) << ", " << 1.0 << ") = " << t << std::endl;
  t = fmod(5.5, 0.9 + 0.1 * a);
  std::cout << "fmod(" << 5.5 << ", " << (0.9 + 0.1 * a) << ") = " << t << std::endl;
  double dbl_temp;
  t = modf(c, &dbl_temp);
  std::cout << "modf(" << c << ", x) = " << t << std::endl;
  t = modf(c + 0.5, &dbl_temp);
  std::cout << "modf(" << (c + 0.5) << ", x) = " << t << std::endl;
  int int_temp;
  t = frexp(c, &int_temp);
  std::cout << "frexp(" << c << ", i) = " << t << std::endl;
  t = frexp(c + 0.5, &int_temp);
  std::cout << "frexp(" << (c + 0.5) << ", i) = " << t << std::endl;
  std::cout << std::endl;

  std::cout << "Repeated application of the logistic function "
               "(x = r * x * (1.0 - x))"
            << std::endl
            << "Can be chaotic or not depending on 'r'.  Because there is only "
               "one source"
            << std::endl
            << "of uncertainty, the uncorrelated model gives nonsense.  The "
               "correlated class"
            << std::endl
            << "does not accurately model the true distribution because the "
               "Gaussian model"
            << std::endl
            << "breaks down, but does show much greater uncertainty in the "
               "chaotic regime."
            << std::endl
            << "See _Chaos:_Making_a_New_Science_ by James Gleick pp. 70-78 "
               "for more on"
            << std::endl
            << "the logistic function." << std::endl;
  std::cout << LOGISTIC_ITERATIONS << " iterations of logistic function:" << std::endl;
  a = UDoubleTest(2.9, 0.000001);
  t = logistic(a);
  std::cout << "before 1st bifurcation:" << std::endl << t << std::endl << std::endl;

  b = UDoubleTest(3.1, 0.000001);
  t = logistic(b);
  std::cout << "after 1st bifurcation:" << std::endl << t << std::endl << std::endl;

  c = UDoubleTest(3.7, 0.000001);
  t = logistic(c);
  std::cout << "chaotic regime:" << std::endl << t << std::endl << std::endl;

  d = UDoubleTest(3.84, 0.000001);
  t = logistic(d);
  std::cout << "island of order in chaos:" << std::endl << t << std::endl << std::endl;

  return 0;
}
