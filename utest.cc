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

// These lint flags are for Gimpel Flexi-Lint
//lint -save
//lint -e762 -e578 -e1917 -e1918 -e1510 -e1712 -e1730 -e655 -e1909 -e578
//lint -e1511 -e1907 -e36 -e1908 -e763 -e534
#include <math.h>
#include <iostream.h>
//lint -restore
#include "udouble.h"


// logistic function is chaotic for some inputs
#define LOGISTIC_ITERATIONS 50

template <class T>
T logistic(T r)
{
   unsigned int i;
   T x(0.6);

   for (i = 0; i < LOGISTIC_ITERATIONS; i++)
      x = r * x * (1.0 - x);
   return x;
}

void infix_check (void) {
   UDoubleTest     a, b(0.0, 0.02, "b (0.0 +/- 0.02)");
   UDoubleTest     c(5.0, 0.03, "c (5.00 +/- 0.03)"), t;

   a = UDoubleTest(1.0, 0.1, "a (1.00 +/ 0.1)");

   cout << endl << "UNARY/INFIX" << endl;
   cout << "Values for unary & infix tests:" << endl;
   cout << "a = " << a << endl;
   cout << "b = " << b << endl;
   cout << "c = " << c << endl;
   cout << endl;

   cout << "These Measurements of correlation show teh size of random "
           "correlations for" << endl
        << "the different ensemble sizes.  In general correlations will "
           "be smaller for" << endl
        << "larger ensembles, leading to better results."
        << endl << endl;
   cout << "a:b "; a.print_correlation(b);
   cout << "a:c "; a.print_correlation(c);
   cout << "b:c "; b.print_correlation(c);
   cout << endl << "UNARY OPERATORS" << endl << endl;
   cout << "Unary operators always give the same result for all classes:"
        << endl;
   cout << "+a = " << +a << endl;
   cout << "-a = " << -a << endl;
   cout << "preincrement to " << ++a << endl;
   cout << "postincrement leaves at " << a++ << endl;
   cout << "a ends as " << a << endl << endl;
   cout << "predecrement to " << --a << endl;
   cout << "postdecrement leaves at " << a-- << endl;
   cout << "a ends as " << a << endl << endl;

   cout << "INFIX OPERATORS" << endl << endl;

   cout << "Operations where only one operand is uncertain give the same "
           "result for" << endl
        << "all models (except where division by a near-zero uncertain number "
           "causes" << endl
        << "higher-order effects):" << endl;

   cout << "2 + a = " << (2.0 + a) << endl;
   cout << "a + 1 = " << (a + 1.0) << endl;
   cout << "a += 3 = " << (a += 3.0) << endl;

   cout << "2 - a = " << (2.0 - a) << endl;
   cout << "a - 1 = " << (a - 1.0) << endl;
   cout << "a -= 3 = " << (a -= 3.0) << endl;

   cout << "2 * a = " << (2.0 * a) << endl;
   cout << "a * 10 = " << (a * 10.0) << endl;
   cout << "a *= 3 = " << (a *= 3.0) << endl;

   cout << "2 / c = " << (2.0 / c) << endl;
   cout << "a / 10 = " << (a / 10.0) << endl;
   cout << "a /= 3 = " << (a /= 3.0) << endl;

   cout << endl << "Here are some cases where division by an uncertain "
           "number gives a distribution" << endl
        << "which diverges increasingly from the simple Gaussian model.  "
           "At first all" << endl
        << "models give the same results, then, as the number of deviations "
           "from zero" << endl
        << "gradually decreases, the \"better\" classes are needed to track "
           "this effect." << endl
        << "Soon even this second-order correction is not sufficient.  Finally "
           "even the" << endl
        << "two different sizes of ensembles do not give the same answer."
        << endl;
   cout << "5.0 / (a + 4.0)  [~1.0 / (1.000 +/- 0.020)] = "
        << (5.0 / (a + 4.0)) << endl;
   cout << "4.0 / (a + 3.0)  [~1.0 / (1.000 +/- 0.025)] =" << endl
        << (4.0 / (a + 3.0)) << endl;
   cout << "3.0 / (a + 2.0)  [~1.0 / (1.000 +/- 0.033)] =" << endl
        << (3.0 / (a + 2.0)) << endl;
   cout << "2.0 / (a + 1.0)  [~1.0 / (1.00 +/- 0.05)] =" << endl
        << (2.0 / (a + 1.0)) << endl;
   cout << "1.0 / a  [~1.0 / (1.00 +/- 0.10)] =" << endl
        << (1.0 / a) << endl;
   cout << "0.8 / (a - 0.2)  [~1.0 / (1.00 +/- 0.13)] =" << endl
        << (0.8 / (a - 0.2)) << endl;
   cout << "0.6 / (a - 0.4)  [~1.0 / (1.00 +/- 0.17)] =" << endl
        << (0.6 / (a - 0.4)) << endl;
   cout << "0.4 / (a - 0.6)  [~1.0 / (1.00 +/- 0.25)] =" << endl
        << (0.4 / (a - 0.6)) << endl;
   cout << "0.2 / (a - 0.8)  [~1.0 / (1.00 +/- 0.50)] =" << endl
        << (0.2 / (a - 0.8)) << endl;
   cout << endl;

   cout << "When both operands are the same, the correlated answer is right:"
        << endl;
   cout << "a + a =" << endl << (a + a) << endl;
   cout << "a - a =" << endl << (a - a) << endl;
   cout << "(squaring introduces second-order effects)" << endl;
   cout << "a * a =" << endl << (a * a) << endl;
   cout << "a / a =" << endl << (a / a) << endl;
   t = a;
   cout << "a += a =" << endl << (t += a) << endl;
   t = a;
   cout << "a -= a =" << endl << (t -= a) << endl;
   t = a;
   cout << "a *= a =" << endl << (t *= a) << endl;
   t = a;
   cout << "a /= a =" << endl << (t /= a) << endl;

   cout << endl << "When the operands are independent, the uncorrelated answer "
        << "is right:" << endl;
   cout << "a + b =" << endl << (a + b) << endl;
   cout << "a - b =" << endl << (a - b) << endl;
   cout << "a * b =" << endl << (a * b) << endl;
   cout << "a * c =" << endl << (a * c) << endl;
   cout << "c / a =" << endl << (c / a) << endl;

   cout << "For these more complicated expressions, neither the uncorrelated"
        << endl
        << "nor the correlated uncertainty is right but correlation tracking"
        << endl
        << "comes through.  Correlation tracking tracking classes can also give"
        << endl
        << "precise information on sources of uncertainty, while the ensemble"
        << endl
        << "class can give pretty good estimates."
        << endl;
   cout << "a + b - a + c =" << endl << (t = a + b - a + c);
   t.print_uncertain_sources();
   cout << "(a + b) / (a + b) =" << endl << (t = (a + b) / (a + b)) << endl;
   cout << "a + a + b = " << endl << (t = a + a + b) << endl;
   t.print_uncertain_sources();
   cout << "(c - a) / (c + a) =" << endl << (t = (c - a) / (c + a)) << endl;
   t.print_uncertain_sources();
}

void io_check(void) {
   const int size = 80;
   char input[size] = "3 +/- 4 1+/-2 2.+/-1 3+/-.3";
   istrstream ibuf(input, size);
   UDoubleTest a, b, c, d;

   cout << "EXTRACTOR" << endl;
   cout << "From \"" << input << "\" extracted: " << endl;
   ibuf >> a;
   cout << a << "," << endl;
   ibuf >> b;
   cout << b << "," << endl;
   ibuf >> c >> d;
   cout << c << "," << endl << d << endl;
}


int main (void) {
   UDoubleInit udi_main;

   io_check();
   UDoubleTest::new_epoch();

   infix_check();
   UDoubleTest::new_epoch();

   UDoubleTest     a, b(0.0, 0.02, "b (0.0 +/- 0.02)");
   UDoubleTest     c(5.0, 0.03, "c (5.00 +/- 0.03)"), d, t;

   a = UDoubleTest(1.0, 0.1, "a (1.00 +/ 0.1)");
   d = UDoubleTest(1.0, 0.1, "d (1.00 +/ 0.1)");

   cout << "Here is a simple demonstration from the article of what happens"
        << endl
        << "when two uncertain numbers are added." << endl;
   // make "t" be half correlated with "a" and half with "d"
   // renormalized to be 1.00 +/- 0.10
   t = (a + d) / sqrt(2.0) + 1.0 - sqrt(2.0);
   cout << "*** Correlated: " << a << "  +  " << a << "  =  " << endl
        << (a + a) << endl;
   cout << "*** Uncorrelated: " << a << "  +  " << d << "  =  " << endl
        << (a + d) << endl;
   cout << "*** Half correlated: " << a << "  +  1.00 +/- 0.10  =  " << endl
        << (a + t) << endl;

   cout << endl << "MATH LIBRARY FUNCTIONS" << endl << endl;

   cout << "Math library functions with just one argument have no important "
           "effects" << endl
        << "from correlation." << endl;

   cout << "sqrt(" << c << ") = " << sqrt(c) << endl;
   cout << "sin(" << b << ") = " << sin(b) << endl;
   cout << "cos(" << c << ") = " << cos(c) << endl;
   cout << "ldexp(" << (c+0.5) << ", " << 5 << ") = " << ldexp(c+0.5, 5) << endl;
   cout << "ldexp(" << a << ", " << 3 << ") = " << ldexp(a, 3) << endl;

   cout << "In some cases the \"better\" classes are needed to model "
           "higher-order" << endl
        << "effects, as seen in their agreement with the ensembles." << endl;
   cout << "sin(" << c << ") = " << endl << sin(c) << endl;
   cout << "sin(" << a << ") = " << endl << sin(a) << endl;
   cout << "cos(" << b << ") = " << endl << cos(b) << endl;
   cout << "cos(" << (0.02+b) << ") = " << endl << cos(0.02+b) << endl;
   cout << "cos(" << (a) << ") = " << endl << cos(a) << endl;

   int int_temp;
   t = frexp(c+0.1, &int_temp);
   cout << "frexp(" << (c+0.1) << ", i) =" << endl << t << endl;
   cout << endl;


   cout << "For inherently discontinuous functions, different values should "
           "be expected" << endl
        << "from different models (and these functions should rarely be used)."
        << endl;
   cout << "ceil(" << c << ") =" << endl << (t = ceil(c)) << endl;
   cout << "floor(" << c << ") =" << endl << (t = floor(c)) << endl;
   cout << "fabs(" << b << ") =" << endl << (t = fabs(b)) << endl;
   cout << "fabs(" << (b + 0.02) << ") =" << endl << (t = fabs(b+0.02)) << endl;
   cout << endl;

   cout << "Library functions with two arguments raise the same issues as "
           "infix" << endl
        << "operators." << endl;
   cout << "pow(" << c << ", " << 2.0*a << ") =" << endl << pow(c, 2.0*a)
        << endl;
   cout << "pow(" << c << ", " << c << ") =" << endl << pow(c, c) << endl;
   cout << "atan2(" << c << ", " << c << ") =" << endl <<  atan2(c, c) << endl;
   cout << "atan2(" << a << ", " << c << ") =" << endl << (t = atan2(a, c));
   t.print_uncertain_sources();
   cout << "fmod(" << (c+0.5) << ", " << 1.0 << ") = "
        << fmod(c+0.5, 1.0) << endl;
   cout << "fmod(" << 5.5 << ", " << (0.89 + 0.11*a) << ") = "
        << fmod(5.5, 0.89 + 0.11*a) << endl;
   double dbl_temp;
   t = modf(c, &dbl_temp);
   cout << "modf(" << c << ", x) =" << endl
        << t << "(" << dbl_temp << ")" << endl;


   cout << endl;
   cout << "Now we test some formulas to see if we get correct values.  These "
           "expressions" << endl
        << "frequently trip up the \"better\" classes because they apply "
           "second-order" << endl
        << "corrections too soon." << endl;
   cout << "sin^2 + cos^2 (" << c << ") =" << endl
        << (cos(c) * cos(c) + sin(c) * sin(c)) << endl;
   cout << "sin^2 + cos^2 (" << b << ") =" << endl
        << (cos(b) * cos(b) + sin(b) * sin(b)) << endl;
   t = a + b;
   cout << "sin^2 + cos^2 (" << t << ") =" << endl
        << (cos(t) * cos(t) + sin(t) * sin(t)) << endl;

   cout << "checking trig and exp functions: should be 1.0 +/- 0.0" << endl;
   cout << "(asin(sin(atan(tan(acos(cos(log(exp(a)))))))) / a) =" << endl;
   cout << (asin(sin(atan(tan(acos(cos(log(exp(a)))))))) / a) << endl;

   cout << "checking log10: should be 1.0 +/- 0.0" << endl;
   cout << "(log10(a) * log(10.0) / log(a)) = " << endl;
   print_multi(cout, (log10(a + 10.0) * log(10.0) / log(a + 10.0)),
                     (log10(b + 2.0) * log(10.0) / log(b + 2.0)),
                     (log10(c) * log(10.0) / log(c)));
 
   cout << endl;
   cout << "checking hyperbolic trig functions: should be 1.0 +/- 0.0" << endl;
   cout << "(cosh(a)*cosh(a)/sqrt((1.0+sinh(a)*sinh(a)) / "
        << "(1.0 - tanh(a)*tanh(a)))) = " << endl;
   print_multi(cout,
        (cosh(a)*cosh(a)/sqrt((1.0+sinh(a)*sinh(a)) / (1.0 - tanh(a)*tanh(a)))),
        (cosh(b)*cosh(b)/sqrt((1.0+sinh(b)*sinh(b)) / (1.0 - tanh(b)*tanh(b)))),
        (cosh(c)*cosh(c)/sqrt((1.0+sinh(c)*sinh(c)) / (1.0 - tanh(c)*tanh(c)))));


   cout << "checking pow: should be 1.0 +/- 0.0" << endl;
   cout << (sqrt(c) / pow(c, 0.5)) << endl;
   cout << (exp(c) / pow(exp(1.0), c)) << endl;

   cout << "checking atan2: should be 1.0 +/- 0.0" << endl;
   cout << (atan(a) / atan2(a, 1.0)) << endl;

   UDoubleTest::new_epoch();

   cout << endl;
   cout << "Some functions can become chaotic, showing sensitive dependence "
           "on initial" << endl
        << "conditions.  Only the ensemble model can really handle this case, "
           "but the" << endl
        << "onset of chaos can also be seen by looking at the sudden increase "
           "in the" << endl
        << "uncertainty as reported by the correlated and correlation tracking "
           "classes." << endl;
   cout << endl;

   cout << LOGISTIC_ITERATIONS << " iterations of the logistic function:"
        << endl;
   a = UDoubleTest(2.9, 0.000001);
   t = logistic(a);
   cout << "before 1st bifurcation:" << endl << t;
   a.print_correlation(t);

   b = UDoubleTest(3.1, 0.000001);
   t = logistic(b);
   cout << "after 1st bifurcation:" << endl << t;
   b.print_correlation(t);

   c = UDoubleTest(3.7, 0.000001);
   t = logistic(c);
   cout << "chaotic regime:" << endl << t;
   c.print_correlation(t);

   d = UDoubleTest(3.84, 0.000001);
   t = logistic(d);
   cout << "island of order in chaos:" << endl << t;
   d.print_correlation(t);

   return 0;
}

