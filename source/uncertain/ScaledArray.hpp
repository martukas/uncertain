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



namespace uncertain
{

// Specialized array class has only those members needed to be an array
// of uncertainty elements used as the template parameter in UDoubleCT<>.
// This class is like SimpleArray but adds an undistributed factor for
// greater efficiency in the common case when all members of an array
// are to be multiplied by some factor.
template<size_t size>
class ArrayWithScale
{
 private:
  double element[size];
  double scale;

 public:
  ArrayWithScale(double initval = 0.0)
  {
    for (size_t i = 0; i < size; i++)
      element[i] = initval;
    scale = 1.0;
  }

  ArrayWithScale(const ArrayWithScale& a)
  {
    for (size_t i = 0; i < size; i++)
      element[i] = a.element[i];
    scale = a.scale;
  }

  ~ArrayWithScale() {}

  ArrayWithScale operator-() const
  {
    ArrayWithScale retval;

    for (size_t i = 0; i < size; i++)
      retval.element[i] = element[i];
    retval.scale = -scale;
    return retval;
  }

  ArrayWithScale& operator+=(const ArrayWithScale& b)
  {
    if (scale)
    {
      double scale_factor = b.scale / scale;
      for (size_t i = 0; i < size; i++)
        element[i] += b.element[i] * scale_factor;
    }
    else
    {
      scale = b.scale;
      for (size_t i = 0; i < size; i++)
        element[i] = b.element[i];
    }
    return *this;
  }

  friend ArrayWithScale operator+(ArrayWithScale a, const ArrayWithScale& b) { return a += b; }

  ArrayWithScale& operator-=(const ArrayWithScale& b)
  {
    if (scale)
    {
      double scale_factor = b.scale / scale;
      for (size_t i = 0; i < size; i++)
        element[i] -= b.element[i] * scale_factor;
    }
    else
    {
      scale = -b.scale;
      for (size_t i = 0; i < size; i++)
        element[i] = b.element[i];
    }
    return *this;
  }

  ArrayWithScale& operator*=(double b)
  {
    scale *= b;
    return *this;
  }

  friend ArrayWithScale operator*(ArrayWithScale a, double b) { return a *= b; }

  ArrayWithScale& operator/=(double b)
  {
    scale /= b;
    return *this;
  }

  double operator[](size_t subscript)
  {
    if (subscript >= size)
    {
      throw std::runtime_error("Error: oversize subscript: " + std::to_string(subscript) +
          " Greater than " + std::to_string(size - 1));
    }
    return element[subscript] * scale;
  }

  void setelement(size_t subscript, double value)
  {
    if (subscript >= size)
    {
      throw std::runtime_error("Error: oversize subscript: " + std::to_string(subscript) +
          " Greater than " + std::to_string(size - 1));
    }
    if (scale != 0.0)
      element[subscript] = value / scale;
  }

  double norm() const
  {
    if (scale == 0.0)
      return 0.0;
    double tot = 0.0;
    for (size_t i = 0; i < size; i++)
      tot += sqr(element[i]);
    return std::sqrt(tot * sqr(scale));
  }
};

}
