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

#include <uncertain/functions.hpp>
#include <stdlib.h>
#include <cmath>


// \todo add exceptions

namespace uncertain
{

// Specialized array class has only those members needed to be an array
// of uncertainty elements used as the template parameter in UDoubleCT<>.
// This is the simplest possible implementation.
template<size_t size>
class SimpleArray
{
 private:
  double element[size];

 public:
  SimpleArray(double initval = 0.0)
  {
    for (size_t i = 0; i < size; i++)
      element[i] = initval;
  }

  SimpleArray(const SimpleArray& a)
  {
    for (size_t i = 0; i < size; i++)
      element[i] = a.element[i];
  }

  ~SimpleArray() {}

  SimpleArray operator-() const
  {
    SimpleArray retval;

    for (size_t i = 0; i < size; i++)
      retval.element[i] = -element[i];
    return retval;
  }

  SimpleArray& operator+=(const SimpleArray& b)
  {
    for (size_t i = 0; i < size; i++)
      element[i] += b.element[i];
    return *this;
  }

  friend SimpleArray operator+(SimpleArray a, const SimpleArray& b) { return a += b; }

  SimpleArray& operator-=(const SimpleArray& b)
  {
    for (size_t i = 0; i < size; i++)
      element[i] -= b.element[i];
    return *this;
  }

  SimpleArray& operator*=(double b)
  {
    for (size_t i = 0; i < size; i++)
      element[i] *= b;
    return *this;
  }

  friend SimpleArray operator*(SimpleArray a, double b) { return a *= b; }

  SimpleArray& operator/=(double b)
  {
    for (size_t i = 0; i < size; i++)
      element[i] /= b;
    return *this;
  }

  double operator[](size_t subscript)
  {
    if (subscript >= size)
    {
      throw std::runtime_error("Error: oversize subscript: "
                                   + std::to_string(subscript)
                                   + " Greater than " + std::to_string(size - 1));
    }
    return element[subscript];
  }

  void setelement(size_t subscript, double value)
  {
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
    for (size_t i = 0; i < size; i++)
      tot += sqr(element[i]);
    return std::sqrt(tot);
  }
};

}
