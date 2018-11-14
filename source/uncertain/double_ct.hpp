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

#include <uncertain/source_set.hpp>
#include <uncertain/simple_array.hpp>
#include <uncertain/scaled_array.hpp>

namespace uncertain
{

// Correlation tracking class keeps an array of uncertainty
// components from various sources.  (Array implementation
// is specified by template parameter.)
template<class T, size_t size>
class UDoubleCT
{
 private:
  double value;
  T unc_components;
  size_t epoch;
  static SourceSet<size> sources;

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
      size_t new_source_num = sources.get_new_source(source_name);
      unc_components.setelement(new_source_num, unc);
    }
  }

  // copy constructor does not create a new independent uncertainty element
  UDoubleCT(const UDoubleCT& ud) : value(ud.value), epoch(ud.epoch)
  {
    unc_components = ud.unc_components;
  }
  // operator= ?

  ~UDoubleCT()
  {}

  static void new_epoch()
  { sources.new_epoch(); }

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

  UDoubleCT operator+() const
  { return *this; }

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

  UDoubleCT& operator+=(double b)
  {
    value += b;
    return *this;
  }

  friend UDoubleCT operator+(UDoubleCT a, const UDoubleCT& b)
  { return a += b; }

  friend UDoubleCT operator+(UDoubleCT a, double b)
  { return a += b; }

  friend UDoubleCT operator+(double b, UDoubleCT a)
  { return a += b; }

  UDoubleCT& operator-=(const UDoubleCT& b)
  {
    sources.check_epoch(epoch);
    sources.check_epoch(b.epoch);
    unc_components -= b.unc_components;
    value -= b.value;
    return *this;
  }

  UDoubleCT& operator-=(double b)
  {
    value -= b;
    return *this;
  }

  friend UDoubleCT operator-(UDoubleCT a, const UDoubleCT& b)
  { return a -= b; }

  friend UDoubleCT operator-(UDoubleCT a, double b)
  { return a -= b; }

  friend UDoubleCT operator-(double b, UDoubleCT a)
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

  UDoubleCT& operator*=(double b)
  {
    unc_components *= b;
    value *= b;
    return *this;
  }

  UDoubleCT operator++()
  { return (*this += 1.0); }

  UDoubleCT operator--()
  { return (*this -= 1.0); }

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

  friend UDoubleCT operator*(UDoubleCT a, const UDoubleCT& b)
  { return a *= b; }

  friend UDoubleCT operator*(UDoubleCT a, double b)
  { return a *= b; }

  friend UDoubleCT operator*(double b, UDoubleCT a)
  { return a *= b; }

  UDoubleCT& operator/=(const UDoubleCT& b)
  {
    sources.check_epoch(epoch);
    sources.check_epoch(b.epoch);
    unc_components /= b.value;
    unc_components -= b.unc_components * (value / (b.value * b.value));
    value /= b.value;
    return *this;
  }

  UDoubleCT& operator/=(double b)
  {
    unc_components /= b;
    value /= b;
    return *this;
  }

  friend UDoubleCT operator/(UDoubleCT a, const UDoubleCT& b)
  { return a /= b; }

  friend UDoubleCT operator/(UDoubleCT a, double b)
  { return a /= b; }

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
    ud = UDoubleCT<T, size>(mean, sigma, source_name);
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
      UDoubleCT<T, size>::sources.check_epoch(arg1.epoch); \
      UDoubleCT<T, size>::sources.check_epoch(arg2.epoch); \
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

  double mean() const
  {
    return value;
  }

  double deviation() const
  {
    return unc_components.norm();
  }
};

}
