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


#pragma once

#include <uncertain/UDoubleMS.hpp>

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
  uncertain::UDoubleMSUncorr msu;
  uncertain::UDoubleMSCorr msc;

 public:
  UDoubleTest(const double val = 0.0, const double unc = 0.0)
      : msu(val, unc), msc(val, unc) {}

  UDoubleTest(const UDoubleTest& ud)
      : msu(ud.msu), msc(ud.msc) {}

  UDoubleTest operator+() const
  {
    UDoubleTest retval;

    retval.msu = +msu;
    retval.msc = +msc;
    return retval;
  }

  UDoubleTest operator-() const
  {
    UDoubleTest retval;

    retval.msu = -msu;
    retval.msc = -msc;
    return retval;
  }

  friend UDoubleTest operator+(const UDoubleTest& a, const UDoubleTest& b)
  {
    UDoubleTest retval;

    retval.msu = a.msu + b.msu;
    retval.msc = a.msc + b.msc;
    return retval;
  }

  friend UDoubleTest operator+(const double& a, const UDoubleTest& b)
  {
    UDoubleTest retval;
    retval.msu = a + b.msu;
    retval.msc = a + b.msc;
    return retval;
  }

  friend UDoubleTest operator+(const UDoubleTest& a, const double& b)
  {
    UDoubleTest retval;

    retval.msu = a.msu + b;
    retval.msc = a.msc + b;
    return retval;
  }

  friend UDoubleTest operator-(const UDoubleTest& a, const UDoubleTest& b)
  {
    UDoubleTest retval;

    retval.msu = a.msu - b.msu;
    retval.msc = a.msc - b.msc;
    return retval;
  }

  friend UDoubleTest operator-(const double& a, const UDoubleTest& b)
  {
    UDoubleTest retval;

    retval.msu = a - b.msu;
    retval.msc = a - b.msc;
    return retval;
  }

  friend UDoubleTest operator-(const UDoubleTest& a, const double& b)
  {
    UDoubleTest retval;

    retval.msu = a.msu - b;
    retval.msc = a.msc - b;
    return retval;
  }

  friend UDoubleTest operator*(const UDoubleTest& a, const UDoubleTest& b)
  {
    UDoubleTest retval;

    retval.msu = a.msu * b.msu;
    retval.msc = a.msc * b.msc;
    return retval;
  }

  friend UDoubleTest operator*(const double& a, const UDoubleTest& b)
  {
    UDoubleTest retval;

    retval.msu = a * b.msu;
    retval.msc = a * b.msc;
    return retval;
  }

  friend UDoubleTest operator*(const UDoubleTest& a, const double& b)
  {
    UDoubleTest retval;

    retval.msu = a.msu * b;
    retval.msc = a.msc * b;
    return retval;
  }

  friend UDoubleTest operator/(const UDoubleTest& a, const UDoubleTest& b)
  {
    UDoubleTest retval;

    retval.msu = a.msu / b.msu;
    retval.msc = a.msc / b.msc;
    return retval;
  }

  friend UDoubleTest operator/(const double& a, const UDoubleTest& b)
  {
    UDoubleTest retval;

    retval.msu = a / b.msu;
    retval.msc = a / b.msc;

    return retval;
  }

  friend UDoubleTest operator/(const UDoubleTest& a, const double& b)
  {
    UDoubleTest retval;

    retval.msu = a.msu / b;
    retval.msc = a.msc / b;
    return retval;
  }

  // preincrement and predecrement operators return value after changes
  UDoubleTest operator++()
  {
    ++msu;
    ++msc;
    return *this;
  }

  UDoubleTest operator--()
  {
    --msu;
    --msc;
    return *this;
  }

  // postincrement and postdecrement operators return value before changes
  UDoubleTest operator++(int)
  {
    UDoubleTest retval(*this);

    msu++;
    msc++;
    return retval;
  }

  UDoubleTest operator--(int)
  {
    UDoubleTest retval(*this);

    msu--;
    msc--;
    return retval;
  }

  UDoubleTest& operator+=(const UDoubleTest& ud)
  {
    msu += ud.msu;
    msc += ud.msc;
    return *this;
  }

  UDoubleTest& operator-=(const UDoubleTest& ud)
  {
    msu -= ud.msu;
    msc -= ud.msc;
    return *this;
  }

  UDoubleTest& operator*=(const UDoubleTest& ud)
  {
    msu *= ud.msu;
    msc *= ud.msc;
    return *this;
  }

  UDoubleTest& operator/=(const UDoubleTest& ud)
  {
    msu /= ud.msu;
    msc /= ud.msc;
    return *this;
  }

  UDoubleTest& operator+=(const double& d)
  {
    msu += d;
    msc += d;
    return *this;
  }

  UDoubleTest& operator-=(const double& d)
  {
    msu -= d;
    msc -= d;
    return *this;
  }

  UDoubleTest& operator*=(const double& d)
  {
    msu *= d;
    msc *= d;
    return *this;
  }

  UDoubleTest& operator/=(const double& d)
  {
    msu /= d;
    msc /= d;
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& os, const UDoubleTest& ud)
  {
    std::stringstream osu, osc;
    std::string up, cp;
    osu << ud.msu;
    osc << ud.msc;
    if ((up = osu.str()) != (cp = osc.str()))
      os << "Uncorrelated: " << up << "  Correlated: " << cp;
    else
      os << up;
    return os;
  }

#define UDoubleTestfunc1(func) \
      UDoubleTest func(UDoubleTest arg) \
      { \
         std::stringstream os, alt_os; \
         std::string str, slope_str; \
         uncertain::UDoubleMSUncorr alt_msu = PropagateUncertaintiesBySlope(func, arg.msu);\
         arg.msu = func(arg.msu); \
         os << arg.msu; \
         alt_os << alt_msu; \
         if ((slope_str = alt_os.str()) != (str = os.str())) \
            std::cerr << "Warning: different values for " << #func "(): " \
                 << str << " vs. " << slope_str << std::endl; \
         arg.msc = func(arg.msc); \
         return arg; \
      }
#define UDoubleTestfunc2(func) \
      UDoubleTest func(const UDoubleTest& arg1, const UDoubleTest& arg2) \
      { \
         UDoubleTest retval; \
         std::stringstream os, alt_os; \
         std::string str, slope_str; \
         uncertain::UDoubleMSUncorr alt_msu = PropagateUncertaintiesBySlope(func, \
                                                  arg1.msu, arg2.msu); \
         retval.msu = func(arg1.msu, arg2.msu); \
         os << retval.msu; \
         alt_os << alt_msu; \
         if ((slope_str = alt_os.str()) != (str = os.str())) \
            std::cerr << "Warning: different values for " << #func << "(): " \
                 << str << " vs. " << slope_str << std::endl; \
         retval.msc = func(arg1.msc, arg2.msc); \
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
    std::stringstream os, alt_os;
    std::string str, slope_str;
    GlobalInt = intarg;
    uncertain::UDoubleMSUncorr alt_msu = PropagateUncertaintiesBySlope(my_ldexp, arg.msu);
    arg.msu = ldexp(arg.msu, intarg);
    os << arg.msu;
    alt_os << alt_msu;
    if ((slope_str = alt_os.str()) != (str = os.str()))
      std::cerr << "Warning: different values for " << "ldexp(): "
                << str << " vs. " << slope_str << std::endl;
    arg.msc = ldexp(arg.msc, intarg);
    return arg;
  }

  friend UDoubleTest frexp(UDoubleTest arg, int* intarg)
  {
    std::stringstream os, alt_os;
    std::string str, slope_str;
    uncertain::UDoubleMSUncorr alt_msu = PropagateUncertaintiesBySlope(my_frexp, arg.msu);
    arg.msu = frexp(arg.msu, intarg);
    os << arg.msu;
    alt_os << alt_msu;
    if ((slope_str = alt_os.str()) != (str = os.str()))
      std::cerr << "Warning: different values for " << "frexp(): "
                << str << " vs. " << slope_str << std::endl;
    arg.msc = frexp(arg.msc, intarg);
    return arg;
  }

  friend UDoubleTest modf(UDoubleTest arg, double* dblarg)
  {
    std::stringstream os, alt_os;
    std::string str, slope_str;
    uncertain::UDoubleMSUncorr alt_msu = PropagateUncertaintiesBySlope(my_modf, arg.msu);
    arg.msu = modf(arg.msu, dblarg);
    os << arg.msu;
    alt_os << alt_msu;
    if ((slope_str = alt_os.str()) != (str = os.str()))
      std::cerr << "Warning: different values for " << "modf(): "
                << str << " vs. " << slope_str << std::endl;
    arg.msc = modf(arg.msc, dblarg);
    return arg;
  }

};
