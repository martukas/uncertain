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

#include <uncertain/functions.hpp>

// \todo add exceptions

namespace uncertain
{

// prints uncertainty to 2 digits and value to same precision
void uncertain_print(double mean, double sigma, std::ostream& os)
{
  auto original_precision = os.precision();
  auto original_format = os.flags(std::ios::showpoint);

  // std::cerr << "<" << mean << " +/- " << sigma << "> " << std::endl;
  int precision;
  // special cases for zero, NaN, and Infinities (positive & negative)
  if ((sigma == 0.0) || (sigma != sigma) || (1.0 / sigma == 0.0))
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
  os << std::setprecision(precision)
     << mean << " +/- "
     << std::setprecision(2)
     << sigma
     << std::setprecision(original_precision);
  os.flags(original_format);
}

// reads uncertainty as mean +/- sigma
void uncertain_read(double& mean, double& sigma, std::istream& is)
{
  char plus, slash, minus;
  is >> mean >> plus >> slash >> minus >> sigma;
  if ((plus != '+') || (slash != '/') || (minus != '-'))
  {
    throw std::runtime_error("Error: illegal characters encountered in reading mean +/- sigma");
  }
}

// \todo include skewing of distribution
void gauss_loss(const double& uncertainty, const double& disc_dist,
                const discontinuity_type& disc_type,
                std::string id_string,
                std::string func_str,
                const double& disc_thresh)
{
  double scaled_disc_dist = fabs(disc_dist / uncertainty);
  if ((scaled_disc_dist < disc_thresh) && (disc_type != none))
  {
    int original_precision = std::cerr.precision();
    std::cerr << std::setprecision(2);
    std::cerr << func_str << "is " << scaled_disc_dist << " sigmas"
              << id_string;
    if (disc_type == step)
      std::cerr << " from a step discontinuity" << std::endl;
    else if (disc_type == infinite_wrap)
      std::cerr << " from an infinite wrap discontinuity" << std::endl;
    else if (disc_type == infinite_then_undef)
      std::cerr << " from an infinite "
                << "discontinuity beyond which it is undefined" << std::endl;
    else if (disc_type == slope_only)
      std::cerr << " from a discontinuity in slope" << std::endl;
    else if (disc_type == undefined_beyond)
      std::cerr << " from a point beyond which it is undefined" << std::endl;
    else
      std::cerr << " from unknown discontinuity " << disc_type << std::endl;
    std::cerr << std::setprecision(original_precision);
  }
}

}