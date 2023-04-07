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

#pragma once

#include <uncertain/functions.hpp>
#include <vector>

namespace uncertain {

// Specialized array class has only those members needed to be an array
// of uncertainty elements used as the template parameter in UDoubleCT<>.
// This class is like SimpleArray but adds an undistributed factor for
// greater efficiency in the common case when all members of an array
// are to be multiplied by some factor.
class ScaledArray {
 private:
  std::vector<double> elements;
  double scale{1.0};

 public:
  ScaledArray() = default;

  ScaledArray(const ScaledArray &a) = default;

  ~ScaledArray() = default;

  ScaledArray operator-() const {
    ScaledArray retval = *this;
    retval.scale = -retval.scale;
    return retval;
  }

  ScaledArray &operator+=(const ScaledArray &b) {
    if (scale != 0.0) {
      if (elements.size() < b.elements.size()) elements.resize(b.elements.size(), 0.0);
      double scale_factor = b.scale / scale;
      for (size_t i = 0; i < b.elements.size(); i++) elements[i] += b.elements[i] * scale_factor;
    } else {
      scale = b.scale;
      elements = b.elements;
    }
    return *this;
  }

  friend ScaledArray operator+(ScaledArray a, const ScaledArray &b) { return a += b; }

  ScaledArray &operator-=(const ScaledArray &b) {
    if (scale != 0.0) {
      if (elements.size() < b.elements.size()) elements.resize(b.elements.size(), 0.0);
      double scale_factor = b.scale / scale;
      for (size_t i = 0; i < b.elements.size(); i++) elements[i] -= b.elements[i] * scale_factor;
      return *this;
    } else {
      scale = -b.scale;
      elements = b.elements;
    }
    return *this;
  }

  ScaledArray &operator*=(double b) {
    scale *= b;
    return *this;
  }

  friend ScaledArray operator*(ScaledArray a, double b) { return a *= b; }

  ScaledArray &operator/=(double b) {
    scale /= b;
    return *this;
  }

  double operator[](size_t subscript) const { return elements[subscript] * scale; }

  void set_element(size_t idx, double value) {
    if (scale == 0.0) return;
    if (idx >= elements.size()) elements.resize(idx + 1, 0.0);
    elements[idx] = value / scale;
  }

  double norm() const {
    if (scale == 0.0) return 0.0;
    double tot = 0.0;
    for (const auto &e : elements) tot += sqr(e);
    return std::sqrt(tot * sqr(scale));
  }
};

}  // namespace uncertain
