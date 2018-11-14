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

// \todo get rid of this include
#include <string.h>

#include <uncertain/UDoubleMS.hpp>
#include <uncertain/SourceSet.hpp>

namespace uncertain
{

// Ensemble uncertainty class.  Represents a distribution by a
// set of n=esize possible values distributed at intervals of
// uniform probability throughout the
// distribution.  The order of the possible values is randomized
// so that all operations can be done element-by-element with
// a minimum of spurious correlations between independent sources
// of uncertainty.  Depending on esize this class can be anywhere
// from very expensive computationally to unusably expensive.
// But for small problems and big esizes it gives "perfect"
// answers.
template<size_t esize>
class UDoubleEnsemble
{
 public:
  static double src_ensemble[MAX_UNC_ELEMENTS][esize];

 private:
  double ensemble[esize];
  static UncertainSourceSet sources;
  unsigned long epoch;

 public:
  // The main constructor initializes a new source of uncertainty
  // (if there is uncertainty).
  UDoubleEnsemble(const double val = 0.0, const double unc = 0.0,
                  std::string name = "")
      : epoch(sources.get_epoch())
  {
    if (unc < 0.0)
    {
      throw std::runtime_error("Error: negative uncertainty: " + std::to_string(unc));
    }
    if (unc != 0.0)
    {
      static double gauss_ensemble[esize];
      static int gauss_ensemble_inited = 0;

      // The base ensemble of n=esize points needs be initialized only
      // once for each ensemble size.  Once it is initialized, each new
      // independent uncertainty element can be made by copying & shuffling
      // this array then scaling it to the appropriate uncertainty and
      // translating it to the appropriate mean.
      if (!gauss_ensemble_inited)
      {
        gauss_ensemble_inited = 1;
        if (esize & 1) // odd ensemble size
        {
          for (size_t i = 0; i < esize / 2; i++)
          {
            double deviate = inverse_gaussian_density(
                (2.0 * (i + 1.0)) / (2.0 * esize));
            gauss_ensemble[2 * i] = deviate;
            gauss_ensemble[2 * i + 1] = -deviate;
          }
          gauss_ensemble[esize - 1] = 0.0;
        }
        else
        {
          for (size_t i = 0; i < esize / 2; i++)
          {
            double deviate = 0.0;
            double k = (2.0 * i + 1.0) / (2.0 * esize);
            for (unsigned j = 0; j < 100; j++)
              deviate += inverse_gaussian_density(k + (j - 49.5)
                  / (100.0 * esize));
            deviate /= 100.0;
            gauss_ensemble[2 * i] = deviate;
            gauss_ensemble[2 * i + 1] = -deviate;
          }
        }
        // Move the points a little to make all the first 5 moments give
        // exact values.
        PerfectEnsemble(gauss_ensemble, esize);
      }
      for (size_t i = 0; i < esize; i++)
        ensemble[i] = val + gauss_ensemble[i] * unc;
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
      unsigned long source_num = sources.get_num_sources();
      this->shuffle();
      if (sources.can_get_new_source())
      {
        memcpy(src_ensemble[source_num], ensemble, sizeof(ensemble));
        source_num = sources.get_new_source(source_name);
      }
    }
    else  // uncertainty is zero
      for (unsigned long i = 0; i < esize; i++)
        ensemble[i] = val;
  }

  // copy constructor does not introduce a new uncertainty element
  UDoubleEnsemble(const UDoubleEnsemble& ud)
      : epoch(ud.epoch)
  {
    for (size_t i = 0; i < esize; i++)
      ensemble[i] = ud.ensemble[i];
  }

  // constructor from an ensemble.
  // \todo add similar function that shuffles its input
  UDoubleEnsemble(const double* const newensemble,
                  std::string name = "")
      : epoch(sources.get_epoch())
  {
    for (size_t i = 0; i < esize; i++)
      ensemble[i] = newensemble[i];
    std::string source_name;
    if (!name.empty())
    {
      source_name = name;
    }
    else
    {
      source_name = "anon from ensemble: " + std::to_string(ensemble[0]);
    }
    unsigned long source_num = sources.get_num_sources();
    memcpy(src_ensemble[source_num], ensemble, sizeof(ensemble));
    source_num = sources.get_new_source(source_name);
  }
// \todo add constructors with other distributions.

  ~UDoubleEnsemble() {}

  UDoubleEnsemble<esize> operator+() const
  {
    return *this;
  }

  UDoubleEnsemble<esize> operator-() const
  {
    UDoubleEnsemble<esize> retval;
    for (size_t i = 0; i < esize; i++)
      retval.ensemble[i] = -ensemble[i];
    retval.epoch = epoch;
    return retval;
  }

  friend UDoubleEnsemble<esize> operator+(UDoubleEnsemble<esize> a,
                                          const UDoubleEnsemble<esize>& b) { return a += b; }

  friend UDoubleEnsemble<esize> operator+(UDoubleEnsemble<esize> a,
                                          double b) { return a += b; }

  friend UDoubleEnsemble<esize> operator+(double b,
                                          UDoubleEnsemble<esize> a) { return a += b; }

  friend UDoubleEnsemble<esize> operator-(UDoubleEnsemble<esize> a,
                                          const UDoubleEnsemble<esize>& b) { return a -= b; }

  friend UDoubleEnsemble<esize> operator-(UDoubleEnsemble<esize> a,
                                          double b) { return a -= b; }

  friend UDoubleEnsemble<esize> operator-(double b,
                                          UDoubleEnsemble<esize> a) { return -(a -= b); }

  UDoubleEnsemble<esize> operator++() { return (*this += 1.0); }

  UDoubleEnsemble<esize> operator--() { return (*this -= 1.0); }

  UDoubleEnsemble<esize> operator++(int)
  {
    UDoubleEnsemble<esize> retval(*this);
    *this += 1.0;
    return retval;
  }

  UDoubleEnsemble<esize> operator--(int)
  {
    UDoubleEnsemble<esize> retval(*this);
    *this -= 1.0;
    return retval;
  }

  friend UDoubleEnsemble<esize> operator*(UDoubleEnsemble<esize> a,
                                          const UDoubleEnsemble<esize>& b) { return a *= b; }

  friend UDoubleEnsemble<esize> operator*(UDoubleEnsemble<esize> a,
                                          double b) { return a *= b; }

  friend UDoubleEnsemble<esize> operator*(double b,
                                          UDoubleEnsemble<esize> a) { return a *= b; }

  friend UDoubleEnsemble<esize> operator/(UDoubleEnsemble<esize> a,
                                          const UDoubleEnsemble<esize>& b) { return a /= b; }

  friend UDoubleEnsemble<esize> operator/(UDoubleEnsemble<esize> a,
                                          double b) { return a /= b; }

  // this one promotes a to UDoubleEnsemble
  friend UDoubleEnsemble<esize> operator/(double a,
                                          const UDoubleEnsemble<esize>& b)
  {
    UDoubleEnsemble<esize> uda(a);
    return uda /= b;
  }

  UDoubleEnsemble<esize>& operator+=(const UDoubleEnsemble<esize>& ud)
  {
    sources.check_epoch(epoch);
    sources.check_epoch(ud.epoch);

    for (size_t i = 0; i < esize; i++)
      ensemble[i] += ud.ensemble[i];
    return *this;
  }

  UDoubleEnsemble<esize>& operator+=(double d)
  {
    for (size_t i = 0; i < esize; i++)
      ensemble[i] += d;
    return *this;
  }

  UDoubleEnsemble<esize>& operator-=(const UDoubleEnsemble<esize>& ud)
  {
    sources.check_epoch(epoch);
    sources.check_epoch(ud.epoch);

    for (size_t i = 0; i < esize; i++)
      ensemble[i] -= ud.ensemble[i];
    return *this;
  }

  UDoubleEnsemble<esize>& operator-=(double d)
  {
    for (size_t i = 0; i < esize; i++)
      ensemble[i] -= d;
    return *this;
  }

  UDoubleEnsemble<esize>& operator*=(const UDoubleEnsemble<esize>& ud)
  {
    sources.check_epoch(epoch);
    sources.check_epoch(ud.epoch);

    for (size_t i = 0; i < esize; i++)
      ensemble[i] *= ud.ensemble[i];
    return *this;
  }

  UDoubleEnsemble<esize>& operator*=(double d)
  {
    for (size_t i = 0; i < esize; i++)
      ensemble[i] *= d;
    return *this;
  }

  UDoubleEnsemble<esize>& operator/=(const UDoubleEnsemble<esize>& ud)
  {
    sources.check_epoch(epoch);
    sources.check_epoch(ud.epoch);

    for (size_t i = 0; i < esize; i++)
      ensemble[i] /= ud.ensemble[i];
    return *this;
  }

  UDoubleEnsemble<esize>& operator/=(double d)
  {
    for (size_t i = 0; i < esize; i++)
      ensemble[i] /= d;
    return *this;
  }

// \todo add procedures to make persistant
  friend std::ostream& operator<<(std::ostream& os, const UDoubleEnsemble<esize>& ud)
  {
    double mean, sigma, skew, kurtosis, m5;
    moments(ud.ensemble, mean, sigma, skew, kurtosis, m5, esize);
    uncertain_print(mean, sigma, os);

    if (sigma != 0.0)
    {
      auto original_precision = os.precision();
      auto original_format = os.flags(std::ios::showpoint);
      os << std::setprecision(2)
         << " [" << skew << " : " << kurtosis << " : " << m5 << "]"
         << std::setprecision(original_precision);
      os.flags(original_format);
    }
    return os;
  }

  friend std::istream& operator>>(std::istream& is, UDoubleEnsemble<esize>& ud)
  {
    double mean, sigma;
    uncertain_read(mean, sigma, is);
    ud = UDoubleEnsemble<esize>(mean, sigma);
    return is;
  }

#define UDoubleEnsemblefunc1(func) \
      UDoubleEnsemble<esize> func(UDoubleEnsemble<esize> arg) \
      { \
         for (size_t i = 0; i < esize; i++) \
            arg.ensemble[i] = func(arg.ensemble[i]); \
         return arg; \
      }
#define UDoubleEnsemblefunc2(func) \
      UDoubleEnsemble<esize> func(const UDoubleEnsemble<esize>& arg1, \
                                  const UDoubleEnsemble<esize>& arg2) \
      { \
         UDoubleEnsemble<esize> retval(arg1); \
         for (size_t i = 0; i < esize; i++) \
            retval.ensemble[i] = func(arg1.ensemble[i], arg2.ensemble[i]); \
         return retval; \
      }

  friend UDoubleEnsemblefunc1(sqrt)

  friend UDoubleEnsemblefunc1(sin)

  friend UDoubleEnsemblefunc1(cos)

  friend UDoubleEnsemblefunc1(tan)

  friend UDoubleEnsemblefunc1(asin)

  friend UDoubleEnsemblefunc1(acos)

  friend UDoubleEnsemblefunc1(atan)

  friend UDoubleEnsemblefunc2(atan2)

  friend UDoubleEnsemblefunc1(ceil)

  friend UDoubleEnsemblefunc1(floor)

  friend UDoubleEnsemblefunc1(fabs)

  friend UDoubleEnsemblefunc2(fmod)

  friend UDoubleEnsemblefunc1(exp)

  friend UDoubleEnsemblefunc1(log)

  friend UDoubleEnsemblefunc1(log10)

  friend UDoubleEnsemblefunc1(sinh)

  friend UDoubleEnsemblefunc1(cosh)

  friend UDoubleEnsemblefunc1(tanh)

  friend UDoubleEnsemblefunc2(pow)

  friend UDoubleEnsemble<esize> ldexp(UDoubleEnsemble<esize> arg,
                                      const int intarg)
  {
    for (size_t i = 0; i < esize; i++)
      arg.ensemble[i] = ldexp(arg.ensemble[i], intarg);
    return arg;
  }

  friend UDoubleEnsemble<esize> frexp(UDoubleEnsemble<esize> arg, int* intarg)
  {
    // use library frexp on mean to get value of return in second arg
    frexp(arg.mean(), intarg);
    for (size_t i = 0; i < esize; i++)
    {
      int tempint;  // ignore return in second arg in loop
      arg.ensemble[i] = frexp(arg.ensemble[i], &tempint);
    }
    return arg;
  }

  friend UDoubleEnsemble<esize> modf(UDoubleEnsemble<esize> arg,
                                     double* dblarg)
  {
    // use library modf on mean to get value of return in second arg
    modf(arg.mean(), dblarg);
    for (size_t i = 0; i < esize; i++)
    {
      double tempdbl;  // ignore return in second arg in loop
      arg.ensemble[i] = modf(arg.ensemble[i], &tempdbl);
    }
    return arg;
  }

  double mean() const
  {
    double sum = 0.0;

    for (size_t i = 0; i < esize; i++)
      sum += ensemble[i];
    return sum / esize;
  }

  double deviation() const
  {
    size_t i;
    double diff, sum_2_diff = 0.0; // watch overflow!
    double value = this->mean();

    for (i = 0; i < esize; i++)
    {
      diff = ensemble[i] - value;
      sum_2_diff += sqr(diff);
    }
    return std::sqrt(sum_2_diff / esize);
  }

  static void new_epoch() { sources.new_epoch(); }

  void print_uncertain_sources(std::ostream& os = std::cout)
  {
    if (deviation() == 0.0)
      os << "No uncertainty";
    else
    {
      double unaccounted_uncertainty = 1.0;
      for (size_t i = 0; i < sources.get_num_sources(); i++)
      {
        double unc_portion = this->correlation(src_ensemble[i]);
        unc_portion *= unc_portion;
        unaccounted_uncertainty -= unc_portion;
        os << sources.get_source_name(i) << ": "
           << int_percent(unc_portion) << "%" << std::endl;
      }
      os << "other: " << int_percent(unaccounted_uncertainty) << "%" << std::endl;
    }
    os << std::endl;
  }

  double correlation(const UDoubleEnsemble<esize>& ud,
                     const size_t offset = 0) const
  {
    size_t i;
    double diff, diff_ud;
    double value = this->mean();
    double ud_value = ud.mean();
    // watch overflow!
    double sum_2_diff = 0.0, sum_2_diff_ud = 0.0, sum_prod_diff = 0.0;

    for (i = 0; i < esize; i++)
    {
      diff = ensemble[i] - value;
      sum_2_diff += diff * diff;
      size_t j = (i + offset) % esize;
      diff_ud = ud.ensemble[j] - ud_value;
      sum_2_diff_ud += diff_ud * diff_ud;
      sum_prod_diff += diff * diff_ud;
    }
    if (!sum_2_diff)
      return 0.0;
    if (!sum_2_diff_ud)
      return 0.0;
    if (!sum_prod_diff)
      return 0.0;
    return sum_prod_diff / std::sqrt(sum_2_diff * sum_2_diff_ud);
  }

  double correlation(const double ens[esize],
                     const size_t offset = 0) const
  {
    size_t i;
    double diff, diff_ud;
    double value = this->mean();
    double ud_value = 0.0;
    // watch overflow!
    double sum_2_diff = 0, sum_2_diff_ud = 0, sum_prod_diff = 0;

    for (i = 0; i < esize; i++)
      ud_value += ens[i];
    ud_value /= esize;
    for (i = 0; i < esize; i++)
    {
      diff = ensemble[i] - value;
      sum_2_diff += diff * diff;
      size_t j = (i + offset) % esize;
      diff_ud = ens[j] - ud_value;
      sum_2_diff_ud += diff_ud * diff_ud;
      sum_prod_diff += diff * diff_ud;
    }
    if (!sum_2_diff)
      return 0.0;
    if (!sum_2_diff_ud)
      return 0.0;
    if (!sum_prod_diff)
      return 0.0;
    return sum_prod_diff / std::sqrt(sum_2_diff * sum_2_diff_ud);
  }

  // \todo add function that gives a description
  // \todo allow superimposing histograms
  void print_histogram(std::ostream& os = std::cout) const
  {
    // centered bins for each 0.5 sigmas from -4 sigmas to +4 sigmas
    // outliers go in the outer bins
    int bin[17] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    size_t i;
    double value = this->mean();
    double sigma = this->deviation();

    if (sigma == 0.0)
    {
      os << "No histogram when no uncertainty" << std::endl;
      return;
    }
    for (i = 0; i < esize; i++)
    {
      double normval = (ensemble[i] - value) / sigma;
      int intval = int(floor(2.0 * normval + 0.5) + 8.0);
      if (intval < 0)
        bin[0]++;
      else if (intval > 16)
        bin[16]++;
      else
        bin[intval]++;
    }
    int binmax = 1;
    for (i = 0; i < 17; i++)
      if (bin[i] >= binmax)
        binmax = bin[i];
    int scale_divisor = 1;
    while (binmax / scale_divisor > 74)
      scale_divisor++;
    os << "Histogram:  (each * represents ";
    if (scale_divisor == 1)
      os << "1 point)";
    else
      os << scale_divisor << " points)";
    os << setiosflags(std::ios::showpos) << std::endl;
    int bin_display[17]; // number of characters displayed in each bin
    for (i = 0; i < 17; i++)
    {
      bin_display[i] = int((bin[i] + 0.5) / scale_divisor);
    }
    size_t first_display_bin = 0, last_display_bin = 16;
    while (bin_display[first_display_bin] == 0)
      first_display_bin++;
    while (bin_display[last_display_bin] == 0)
      last_display_bin--;
    for (i = first_display_bin; i <= last_display_bin; i++)
    {
      if (i & 1)
        os << "    | ";
      else
        os << std::setw(3) << (i / 2 - 4) << " + ";
      for (int j = 0; j < bin_display[i]; j++)
        os << "*";
      os << std::endl;
    }
    os << resetiosflags(std::ios::showpos) << std::endl;
  }

  friend UDoubleEnsemble<esize> Invoke(double (* certainfunc)(double),
                                       const UDoubleEnsemble<esize>& arg)
  {
    UDoubleEnsemble<esize> retval;

    for (unsigned i = 0; i < esize; i++)
      retval.ensemble[i] = certainfunc(arg.ensemble[i]);
    return retval;
  }

  friend UDoubleEnsemble<esize> Invoke(double (* certainfunc)(double, double),
                                       const UDoubleEnsemble<esize>& arg1,
                                       const UDoubleEnsemble<esize>& arg2)
  {
    UDoubleEnsemble<esize> retval;

    for (unsigned i = 0; i < esize; i++)
      retval.ensemble[i] = certainfunc(arg1.ensemble[i], arg2.ensemble[i]);
    return retval;
  }

  void shuffle()
  {
    for (size_t i = 0; i < esize - 1; i++)
    {
      size_t j = i + (size_t) rand() % (esize - i);
      if (j != i)
      {
        double temp = ensemble[i];
        ensemble[i] = ensemble[j];
        ensemble[j] = temp;
      }
    }
  }

  // figure the moments (sigma, skew, kurtosis, & 5th moment) from an
  // ensemble given the mean.
  static void moments_fixed_mean(double const* const ensemble,
                                 double mean, double& sigma, double& skew,
                                 double& kurtosis, double& m5,
                                 const size_t ens_size)
  {
    size_t i;
    double gddiff2 = 0.0, gddiff3 = 0.0, gddiff4 = 0.0, gddiff5 = 0.0;
    double* gddiff = new double[ens_size];

    // Sorting ensemble by absolute value first increases accuracy
    // but at a cost to performance.  Sorting is O(nlog(n)), whereas
    // the rest of what's being done here is O(n) in ens_size.  Operationally
    // this routine should be called mostly only when there is output,
    // so the cost should be okay.
    for (i = 0; i < ens_size; i++)
      gddiff[i] = ensemble[i] - mean;
    qsort(gddiff, ens_size, sizeof(double), abs_double_compare);
    for (i = 0; i < ens_size; i++)
    {
      gddiff2 += gddiff[i] * gddiff[i];
      gddiff3 += gddiff[i] * gddiff[i] * gddiff[i];
      gddiff4 += gddiff[i] * gddiff[i] * gddiff[i] * gddiff[i];
      gddiff5 += gddiff[i] * gddiff[i] * gddiff[i] * gddiff[i] * gddiff[i];
    }
    double var = gddiff2 / ens_size;
    sigma = std::sqrt(var);
    skew = gddiff3 / (var * sigma * ens_size);
    kurtosis = gddiff4 / (var * var * ens_size) - 3;
    m5 = gddiff5 / (var * var * sigma * ens_size);
    // std::cerr << mean << " +/- " << sigma << " [s: " << skew << ", k: ";
    // std::cerr << kurtosis << " m5: " << m5 << "]" << std::endl;
    delete[] gddiff;
  }

  // figure the moments (mean, sigma, skew, kurtosis, & 5th moment) from an
  // ensemble
  static void moments(double const* const ensemble,
                      double& mean, double& sigma, double& skew,
                      double& kurtosis, double& m5, const size_t ens_size)
  {
    size_t i;
    double gdsum = 0.0;
    double* gddiff = new double[ens_size];

    for (i = 0; i < ens_size; i++)
      gdsum += ensemble[i];
    mean = gdsum / ens_size;
    // Sorting ensemble by absolute value first increases accuracy.
    // See note in moments_fixed_mean().
    for (i = 0; i < ens_size; i++)
      gddiff[i] = ensemble[i] - mean;
    qsort(gddiff, ens_size, sizeof(double), abs_double_compare);
    gdsum = 0.0;
    for (i = 0; i < ens_size; i++)
      gdsum += gddiff[i];
    mean += gdsum / ens_size;

    moments_fixed_mean(ensemble, mean, sigma, skew, kurtosis, m5, ens_size);
    delete[] gddiff;
  }

  // This function moves points a little bit so the first 5 moments all
  // get measured at precisely the expected values.
  static void PerfectEnsemble(double* ensemble, const size_t ens_size)
  {
    unsigned long i;
    double value, sigma, skew, kurtosis, m5;

    double* test_ensemble = new double[ens_size];
    for (int j = 0; j < 3; j++)
    {
      moments(ensemble, value, sigma, skew, kurtosis, m5, ens_size);

      for (i = 0; i < ens_size; i++)
        ensemble[i] -= value;
      for (i = 0; i < ens_size; i++)
        ensemble[i] /= sigma;
      // future work: improve kurtosis correction
      double kurtfact = 0.045;
      for (size_t k = 0; k < 5; k++)
      {
        for (i = 0; i < ens_size; i++)
          test_ensemble[i] = ensemble[i]
              - kurtfact * kurtosis * ensemble[i]
                  * ensemble[i] * ensemble[i];
        double test_value, test_sigma, test_skew, test_kurtosis;
        double test_m5;
        moments(test_ensemble, test_value, test_sigma, test_skew,
                test_kurtosis, test_m5, ens_size);
        kurtfact /= 1 - test_kurtosis / kurtosis;
      }
      for (i = 0; i < ens_size; i++)
        ensemble[i] -= kurtosis * kurtfact * ensemble[i]
            * ensemble[i] * ensemble[i];

    }
    delete[] test_ensemble;
    moments(ensemble, value, sigma, skew, kurtosis, m5, ens_size);
    for (i = 0; i < ens_size; i++)
      ensemble[i] -= value;
    for (i = 0; i < ens_size; i++)
      ensemble[i] /= sigma;
  }
};

}
