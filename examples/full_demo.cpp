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


// utest.cc: This file demonstrates the use of the classes in udouble.h.

// By Evan Manning (manning@alumni.caltech.edu).


#include <cmath>
#include <iostream>
#include <sstream>
#include "UDoubleTest2.hpp"

// logistic function is chaotic for some inputs
#define LOGISTIC_ITERATIONS 50

template<class T>
T logistic(T r)
{
  unsigned int i;
  T x(0.6);

  for (i = 0; i < LOGISTIC_ITERATIONS; i++)
    x = r * x * (1.0 - x);
  return x;
}

void infix_check()
{
  UDoubleTest a, b(0.0, 0.02, "b (0.0 +/- 0.02)");
  UDoubleTest c(5.0, 0.03, "c (5.00 +/- 0.03)"), t;

  a = UDoubleTest(1.0, 0.1, "a (1.00 +/ 0.1)");

  std::cout << std::endl << "UNARY/INFIX" << std::endl;
  std::cout << "Values for unary & infix tests:" << std::endl;
  std::cout << "a = " << a << std::endl;
  std::cout << "b = " << b << std::endl;
  std::cout << "c = " << c << std::endl;
  std::cout << std::endl;

  std::cout << "These Measurements of correlation show teh size of random "
               "correlations for" << std::endl
            << "the different ensemble sizes.  In general correlations will "
               "be smaller for" << std::endl
            << "larger ensembles, leading to better results."
            << std::endl << std::endl;
  std::cout << "a:b ";
  a.print_correlation(b);
  std::cout << "a:c ";
  a.print_correlation(c);
  std::cout << "b:c ";
  b.print_correlation(c);
  std::cout << std::endl << "UNARY OPERATORS" << std::endl << std::endl;
  std::cout << "Unary operators always give the same result for all classes:"
            << std::endl;
  std::cout << "+a = " << +a << std::endl;
  std::cout << "-a = " << -a << std::endl;
  std::cout << "preincrement to " << ++a << std::endl;
  std::cout << "postincrement leaves at " << a++ << std::endl;
  std::cout << "a ends as " << a << std::endl << std::endl;
  std::cout << "predecrement to " << --a << std::endl;
  std::cout << "postdecrement leaves at " << a-- << std::endl;
  std::cout << "a ends as " << a << std::endl << std::endl;

  std::cout << "INFIX OPERATORS" << std::endl << std::endl;

  std::cout << "Operations where only one operand is uncertain give the same "
               "result for" << std::endl
            << "all models (except where division by a near-zero uncertain number "
               "causes" << std::endl
            << "higher-order effects):" << std::endl;

  std::cout << "2 + a = " << (2.0 + a) << std::endl;
  std::cout << "a + 1 = " << (a + 1.0) << std::endl;
  std::cout << "a += 3 = " << (a += 3.0) << std::endl;

  std::cout << "2 - a = " << (2.0 - a) << std::endl;
  std::cout << "a - 1 = " << (a - 1.0) << std::endl;
  std::cout << "a -= 3 = " << (a -= 3.0) << std::endl;

  std::cout << "2 * a = " << (2.0 * a) << std::endl;
  std::cout << "a * 10 = " << (a * 10.0) << std::endl;
  std::cout << "a *= 3 = " << (a *= 3.0) << std::endl;

  std::cout << "2 / c = " << (2.0 / c) << std::endl;
  std::cout << "a / 10 = " << (a / 10.0) << std::endl;
  std::cout << "a /= 3 = " << (a /= 3.0) << std::endl;

  std::cout << std::endl << "Here are some cases where division by an uncertain "
                            "number gives a distribution" << std::endl
            << "which diverges increasingly from the simple Gaussian model.  "
               "At first all" << std::endl
            << "models give the same results, then, as the number of deviations "
               "from zero" << std::endl
            << "gradually decreases, the \"better\" classes are needed to track "
               "this effect." << std::endl
            << "Soon even this second-order correction is not sufficient.  Finally "
               "even the" << std::endl
            << "two different sizes of ensembles do not give the same answer."
            << std::endl;
  std::cout << "5.0 / (a + 4.0)  [~1.0 / (1.000 +/- 0.020)] = "
            << (5.0 / (a + 4.0)) << std::endl;
  std::cout << "4.0 / (a + 3.0)  [~1.0 / (1.000 +/- 0.025)] =" << std::endl
            << (4.0 / (a + 3.0)) << std::endl;
  std::cout << "3.0 / (a + 2.0)  [~1.0 / (1.000 +/- 0.033)] =" << std::endl
            << (3.0 / (a + 2.0)) << std::endl;
  std::cout << "2.0 / (a + 1.0)  [~1.0 / (1.00 +/- 0.05)] =" << std::endl
            << (2.0 / (a + 1.0)) << std::endl;
  std::cout << "1.0 / a  [~1.0 / (1.00 +/- 0.10)] =" << std::endl
            << (1.0 / a) << std::endl;
  std::cout << "0.8 / (a - 0.2)  [~1.0 / (1.00 +/- 0.13)] =" << std::endl
            << (0.8 / (a - 0.2)) << std::endl;
  std::cout << "0.6 / (a - 0.4)  [~1.0 / (1.00 +/- 0.17)] =" << std::endl
            << (0.6 / (a - 0.4)) << std::endl;
  std::cout << "0.4 / (a - 0.6)  [~1.0 / (1.00 +/- 0.25)] =" << std::endl
            << (0.4 / (a - 0.6)) << std::endl;
  std::cout << "0.2 / (a - 0.8)  [~1.0 / (1.00 +/- 0.50)] =" << std::endl
            << (0.2 / (a - 0.8)) << std::endl;
  std::cout << std::endl;

  std::cout << "When both operands are the same, the correlated answer is right:"
            << std::endl;
  std::cout << "a + a =" << std::endl << (a + a) << std::endl;
  std::cout << "a - a =" << std::endl << (a - a) << std::endl;
  std::cout << "(squaring introduces second-order effects)" << std::endl;
  std::cout << "a * a =" << std::endl << (a * a) << std::endl;
  std::cout << "a / a =" << std::endl << (a / a) << std::endl;
  t = a;
  std::cout << "a += a =" << std::endl << (t += a) << std::endl;
  t = a;
  std::cout << "a -= a =" << std::endl << (t -= a) << std::endl;
  t = a;
  std::cout << "a *= a =" << std::endl << (t *= a) << std::endl;
  t = a;
  std::cout << "a /= a =" << std::endl << (t /= a) << std::endl;

  std::cout << std::endl << "When the operands are independent, the uncorrelated answer "
            << "is right:" << std::endl;
  std::cout << "a + b =" << std::endl << (a + b) << std::endl;
  std::cout << "a - b =" << std::endl << (a - b) << std::endl;
  std::cout << "a * b =" << std::endl << (a * b) << std::endl;
  std::cout << "a * c =" << std::endl << (a * c) << std::endl;
  std::cout << "c / a =" << std::endl << (c / a) << std::endl;

  std::cout << "For these more complicated expressions, neither the uncorrelated"
            << std::endl
            << "nor the correlated uncertainty is right but correlation tracking"
            << std::endl
            << "comes through.  Correlation tracking tracking classes can also give"
            << std::endl
            << "precise information on sources of uncertainty, while the ensemble"
            << std::endl
            << "class can give pretty good estimates."
            << std::endl;
  std::cout << "a + b - a + c =" << std::endl << (t = a + b - a + c);
  t.print_uncertain_sources();
  std::cout << "(a + b) / (a + b) =" << std::endl << (t = (a + b) / (a + b)) << std::endl;
  std::cout << "a + a + b = " << std::endl << (t = a + a + b) << std::endl;
  t.print_uncertain_sources();
  std::cout << "(c - a) / (c + a) =" << std::endl << (t = (c - a) / (c + a)) << std::endl;
  t.print_uncertain_sources();
}

void io_check()
{
  std::string input = "3 +/- 4 1+/-2 2.+/-1 3+/-.3";
  std::stringstream ibuf;
  ibuf.str(input);
  UDoubleTest a, b, c, d;

  std::cout << "EXTRACTOR" << std::endl;
  std::cout << "From \"" << input << "\" extracted: " << std::endl;
  ibuf >> a;
  std::cout << a << "," << std::endl;
  ibuf >> b;
  std::cout << b << "," << std::endl;
  ibuf >> c >> d;
  std::cout << c << "," << std::endl << d << std::endl;
}

int main()
{
  UDoubleInit udi_main;

  io_check();
  UDoubleTest::new_epoch();

  infix_check();
  UDoubleTest::new_epoch();

  UDoubleTest a, b(0.0, 0.02, "b (0.0 +/- 0.02)");
  UDoubleTest c(5.0, 0.03, "c (5.00 +/- 0.03)"), d, t;

  a = UDoubleTest(1.0, 0.1, "a (1.00 +/ 0.1)");
  d = UDoubleTest(1.0, 0.1, "d (1.00 +/ 0.1)");

  std::cout << "Here is a simple demonstration from the article of what happens"
            << std::endl
            << "when two uncertain numbers are added." << std::endl;
  // make "t" be half correlated with "a" and half with "d"
  // renormalized to be 1.00 +/- 0.10
  t = (a + d) / sqrt(2.0) + 1.0 - sqrt(2.0);
  std::cout << "*** Correlated: " << a << "  +  " << a << "  =  " << std::endl
            << (a + a) << std::endl;
  std::cout << "*** Uncorrelated: " << a << "  +  " << d << "  =  " << std::endl
            << (a + d) << std::endl;
  std::cout << "*** Half correlated: " << a << "  +  1.00 +/- 0.10  =  " << std::endl
            << (a + t) << std::endl;

  std::cout << std::endl << "MATH LIBRARY FUNCTIONS" << std::endl << std::endl;

  std::cout << "Math library functions with just one argument have no important "
               "effects" << std::endl
            << "from correlation." << std::endl;

  std::cout << "sqrt(" << c << ") = " << sqrt(c) << std::endl;
  std::cout << "sin(" << b << ") = " << sin(b) << std::endl;
  std::cout << "cos(" << c << ") = " << cos(c) << std::endl;
  std::cout << "ldexp(" << (c + 0.5) << ", " << 5 << ") = " << ldexp(c + 0.5, 5) << std::endl;
  std::cout << "ldexp(" << a << ", " << 3 << ") = " << ldexp(a, 3) << std::endl;

  std::cout << "In some cases the \"better\" classes are needed to model "
               "higher-order" << std::endl
            << "effects, as seen in their agreement with the ensembles." << std::endl;
  std::cout << "sin(" << c << ") = " << std::endl << sin(c) << std::endl;
  std::cout << "sin(" << a << ") = " << std::endl << sin(a) << std::endl;
  std::cout << "cos(" << b << ") = " << std::endl << cos(b) << std::endl;
  std::cout << "cos(" << (0.02 + b) << ") = " << std::endl << cos(0.02 + b) << std::endl;
  std::cout << "cos(" << (a) << ") = " << std::endl << cos(a) << std::endl;

  int int_temp;
  t = frexp(c + 0.1, &int_temp);
  std::cout << "frexp(" << (c + 0.1) << ", i) =" << std::endl << t << std::endl;
  std::cout << std::endl;

  std::cout << "For inherently discontinuous functions, different values should "
               "be expected" << std::endl
            << "from different models (and these functions should rarely be used)."
            << std::endl;
  std::cout << "ceil(" << c << ") =" << std::endl << (t = ceil(c)) << std::endl;
  std::cout << "floor(" << c << ") =" << std::endl << (t = floor(c)) << std::endl;
  std::cout << "fabs(" << b << ") =" << std::endl << (t = fabs(b)) << std::endl;
  std::cout << "fabs(" << (b + 0.02) << ") =" << std::endl << (t = fabs(b + 0.02)) << std::endl;
  std::cout << std::endl;

  std::cout << "Library functions with two arguments raise the same issues as "
               "infix" << std::endl
            << "operators." << std::endl;
  std::cout << "pow(" << c << ", " << 2.0 * a << ") =" << std::endl << pow(c, 2.0 * a)
            << std::endl;
  std::cout << "pow(" << c << ", " << c << ") =" << std::endl << pow(c, c) << std::endl;
  std::cout << "atan2(" << c << ", " << c << ") =" << std::endl << atan2(c, c) << std::endl;
  std::cout << "atan2(" << a << ", " << c << ") =" << std::endl << (t = atan2(a, c));
  t.print_uncertain_sources();
  std::cout << "fmod(" << (c + 0.5) << ", " << 1.0 << ") = "
            << fmod(c + 0.5, 1.0) << std::endl;
  std::cout << "fmod(" << 5.5 << ", " << (0.89 + 0.11 * a) << ") = "
            << fmod(5.5, 0.89 + 0.11 * a) << std::endl;
  double dbl_temp;
  t = modf(c, &dbl_temp);
  std::cout << "modf(" << c << ", x) =" << std::endl
            << t << "(" << dbl_temp << ")" << std::endl;

  std::cout << std::endl;
  std::cout << "Now we test some formulas to see if we get correct values.  These "
               "expressions" << std::endl
            << "frequently trip up the \"better\" classes because they apply "
               "second-order" << std::endl
            << "corrections too soon." << std::endl;
  std::cout << "sin^2 + cos^2 (" << c << ") =" << std::endl
            << (cos(c) * cos(c) + sin(c) * sin(c)) << std::endl;
  std::cout << "sin^2 + cos^2 (" << b << ") =" << std::endl
            << (cos(b) * cos(b) + sin(b) * sin(b)) << std::endl;
  t = a + b;
  std::cout << "sin^2 + cos^2 (" << t << ") =" << std::endl
            << (cos(t) * cos(t) + sin(t) * sin(t)) << std::endl;

  std::cout << "checking trig and exp functions: should be 1.0 +/- 0.0" << std::endl;
  std::cout << "(asin(sin(atan(tan(acos(cos(log(exp(a)))))))) / a) =" << std::endl;
  std::cout << (asin(sin(atan(tan(acos(cos(log(exp(a)))))))) / a) << std::endl;

  std::cout << "checking log10: should be 1.0 +/- 0.0" << std::endl;
  std::cout << "(log10(a) * log(10.0) / log(a)) = " << std::endl;
  print_multi(std::cout, (log10(a + 10.0) * log(10.0) / log(a + 10.0)),
              (log10(b + 2.0) * log(10.0) / log(b + 2.0)),
              (log10(c) * log(10.0) / log(c)));

  std::cout << std::endl;
  std::cout << "checking hyperbolic trig functions: should be 1.0 +/- 0.0" << std::endl;
  std::cout << "(cosh(a)*cosh(a)/sqrt((1.0+sinh(a)*sinh(a)) / "
            << "(1.0 - tanh(a)*tanh(a)))) = " << std::endl;
  print_multi(std::cout,
              (cosh(a) * cosh(a) / sqrt((1.0 + sinh(a) * sinh(a)) / (1.0 - tanh(a) * tanh(a)))),
              (cosh(b) * cosh(b) / sqrt((1.0 + sinh(b) * sinh(b)) / (1.0 - tanh(b) * tanh(b)))),
              (cosh(c) * cosh(c) / sqrt((1.0 + sinh(c) * sinh(c)) / (1.0 - tanh(c) * tanh(c)))));

  std::cout << "checking pow: should be 1.0 +/- 0.0" << std::endl;
  std::cout << (sqrt(c) / pow(c, 0.5)) << std::endl;
  std::cout << (exp(c) / pow(exp(1.0), c)) << std::endl;

  std::cout << "checking atan2: should be 1.0 +/- 0.0" << std::endl;
  std::cout << (atan(a) / atan2(a, 1.0)) << std::endl;

  UDoubleTest::new_epoch();

  std::cout << std::endl;
  std::cout << "Some functions can become chaotic, showing sensitive dependence "
               "on initial" << std::endl
            << "conditions.  Only the ensemble model can really handle this case, "
               "but the" << std::endl
            << "onset of chaos can also be seen by looking at the sudden increase "
               "in the" << std::endl
            << "uncertainty as reported by the correlated and correlation tracking "
               "classes." << std::endl;
  std::cout << std::endl;

  std::cout << LOGISTIC_ITERATIONS << " iterations of the logistic function:"
            << std::endl;
  a = UDoubleTest(2.9, 0.000001);
  t = logistic(a);
  std::cout << "before 1st bifurcation:" << std::endl << t;
  a.print_correlation(t);

  b = UDoubleTest(3.1, 0.000001);
  t = logistic(b);
  std::cout << "after 1st bifurcation:" << std::endl << t;
  b.print_correlation(t);

  c = UDoubleTest(3.7, 0.000001);
  t = logistic(c);
  std::cout << "chaotic regime:" << std::endl << t;
  c.print_correlation(t);

  d = UDoubleTest(3.84, 0.000001);
  t = logistic(d);
  std::cout << "island of order in chaos:" << std::endl << t;
  d.print_correlation(t);

  return 0;
}

