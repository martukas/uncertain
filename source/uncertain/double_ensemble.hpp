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

#pragma once

#include <functional>
#include <iomanip>
#include <uncertain/source_set.hpp>
#include <vector>

// \todo get rid of this include
#include <string.h>

namespace uncertain {

// Ensemble uncertainty class.  Represents a distribution by a
// set of n=ensemble_size possible values distributed at intervals of
// uniform probability throughout the distribution.  The order of the
// possible values is randomized so that all operations can be done
// element-by-element with a minimum of spurious correlations between
// independent sources of uncertainty.  Depending on ensemble_size this
// class can be anywhere from very expensive computationally to unusably
// expensive. But for small problems and big ensemble_sizes it gives
// "perfect" answers.
template <size_t ensemble_size>
class UDoubleEnsemble {
 public:
  static std::vector<std::vector<double>> src_ensemble;
  static SourceSet sources;
  static std::vector<double> gauss_ensemble;

 private:
  size_t epoch;
  std::vector<double> ensemble;

  //  static SourceSet sources;

  typedef double (*fun1type)(double);
  typedef double (*fun2type)(double, double);

 public:
  // The main constructor initializes a new source of uncertainty
  // (if there is uncertainty).
  UDoubleEnsemble(double val = 0.0, double unc = 0.0, const std::string &name = {})
      : epoch(sources.get_epoch()) {
    if (unc < 0.0) {
      throw std::runtime_error("Error: negative uncertainty: " + std::to_string(unc));
    }

    ensemble.resize(ensemble_size);

    if (unc != 0.0) {
      // The base ensemble of n=ensemble_size points needs be initialized only
      // once for each ensemble size.  Once it is initialized, each new
      // independent uncertainty element can be made by copying & shuffling
      // this array then scaling it to the appropriate uncertainty and
      // translating it to the appropriate mean.
      if (gauss_ensemble.size() != ensemble_size) {
        gauss_ensemble.resize(ensemble_size);
        if (ensemble_size & 1)  // odd ensemble size
        {
          for (size_t i = 0; i < ensemble_size / 2; i++) {
            double deviate = inverse_gaussian_density((2.0 * (i + 1.0)) / (2.0 * ensemble_size));
            gauss_ensemble[2 * i] = deviate;
            gauss_ensemble[2 * i + 1] = -deviate;
          }
          gauss_ensemble[ensemble_size - 1] = 0.0;
        } else {
          for (size_t i = 0; i < ensemble_size / 2; i++) {
            double deviate = 0.0;
            double k = (2.0 * i + 1.0) / (2.0 * ensemble_size);
            for (unsigned j = 0; j < 100; j++)
              deviate += inverse_gaussian_density(k + (j - 49.5) / (100.0 * ensemble_size));
            deviate /= 100.0;
            gauss_ensemble[2 * i] = deviate;
            gauss_ensemble[2 * i + 1] = -deviate;
          }
        }
        // Move the points a little to make all the first 5 moments give
        // exact values.
        PerfectEnsemble(gauss_ensemble);
      }
      for (size_t i = 0; i < ensemble_size; i++) ensemble[i] = val + gauss_ensemble[i] * unc;
      std::string source_name;
      if (!name.empty()) {
        source_name = name;
      } else {
        std::stringstream os;
        os << "anon: ";
        uncertain_print(val, unc, os);
        source_name = os.str();
      }
      this->shuffle();
      auto source_num = sources.get_new_source(source_name);
      if (source_num >= src_ensemble.size()) src_ensemble.resize(source_num + 1);
      src_ensemble[source_num] = ensemble;
    } else  // uncertainty is zero
      for (size_t i = 0; i < ensemble_size; i++) ensemble[i] = val;
  }

  // copy constructor does not introduce a new uncertainty element
  UDoubleEnsemble(const UDoubleEnsemble &ud) {
    epoch = ud.epoch;
    ensemble = ud.ensemble;
  }

  // constructor from an ensemble.
  // \todo add similar function that shuffles its input
  UDoubleEnsemble(const std::vector<double> &newensemble, const std::string &name = {})
      : epoch(sources.get_epoch()) {
    if (newensemble.size() != ensemble_size) {
      throw std::runtime_error("Cannot construct from wrong ensemble size");
    }
    ensemble = newensemble;
    std::string source_name;
    if (!name.empty()) {
      source_name = name;
    } else {
      source_name = "anon from ensemble: " + std::to_string(ensemble[0]);
    }
    size_t source_num = sources.get_new_source(source_name);
    if (source_num >= src_ensemble.size()) src_ensemble.resize(source_num + 1);
    src_ensemble[source_num] = ensemble;
  }

  // \todo add constructors with other distributions.

  ~UDoubleEnsemble() = default;

  double mean() const {
    double sum{0.0};
    for (const auto &e : ensemble) sum += e;
    return sum / ensemble_size;
  }

  double deviation() const {
    double sum_2_diff{0.0};  // watch overflow!
    double value = this->mean();

    for (const auto &e : ensemble) sum_2_diff += sqr(e - value);
    return std::sqrt(sum_2_diff / ensemble_size);
  }

  UDoubleEnsemble<ensemble_size> operator+() const { return *this; }

  UDoubleEnsemble<ensemble_size> operator-() const {
    UDoubleEnsemble<ensemble_size> retval;
    for (size_t i = 0; i < ensemble_size; i++) retval.ensemble[i] = -ensemble[i];
    retval.epoch = epoch;
    return retval;
  }

  friend UDoubleEnsemble<ensemble_size> operator+(UDoubleEnsemble<ensemble_size> a,
                                                  const UDoubleEnsemble<ensemble_size> &b) {
    return a += b;
  }

  friend UDoubleEnsemble<ensemble_size> operator+(UDoubleEnsemble<ensemble_size> a, double b) {
    return a += b;
  }

  friend UDoubleEnsemble<ensemble_size> operator+(double b, UDoubleEnsemble<ensemble_size> a) {
    return a += b;
  }

  friend UDoubleEnsemble<ensemble_size> operator-(UDoubleEnsemble<ensemble_size> a,
                                                  const UDoubleEnsemble<ensemble_size> &b) {
    return a -= b;
  }

  friend UDoubleEnsemble<ensemble_size> operator-(UDoubleEnsemble<ensemble_size> a, double b) {
    return a -= b;
  }

  friend UDoubleEnsemble<ensemble_size> operator-(double b, UDoubleEnsemble<ensemble_size> a) {
    return -(a -= b);
  }

  UDoubleEnsemble<ensemble_size> operator++() { return (*this += 1.0); }

  UDoubleEnsemble<ensemble_size> operator--() { return (*this -= 1.0); }

  UDoubleEnsemble<ensemble_size> operator++(int) {
    UDoubleEnsemble<ensemble_size> retval(*this);
    *this += 1.0;
    return retval;
  }

  UDoubleEnsemble<ensemble_size> operator--(int) {
    UDoubleEnsemble<ensemble_size> retval(*this);
    *this -= 1.0;
    return retval;
  }

  friend UDoubleEnsemble<ensemble_size> operator*(UDoubleEnsemble<ensemble_size> a,
                                                  const UDoubleEnsemble<ensemble_size> &b) {
    return a *= b;
  }

  friend UDoubleEnsemble<ensemble_size> operator*(UDoubleEnsemble<ensemble_size> a, double b) {
    return a *= b;
  }

  friend UDoubleEnsemble<ensemble_size> operator*(double b, UDoubleEnsemble<ensemble_size> a) {
    return a *= b;
  }

  friend UDoubleEnsemble<ensemble_size> operator/(UDoubleEnsemble<ensemble_size> a,
                                                  const UDoubleEnsemble<ensemble_size> &b) {
    return a /= b;
  }

  friend UDoubleEnsemble<ensemble_size> operator/(UDoubleEnsemble<ensemble_size> a, double b) {
    return a /= b;
  }

  // this one promotes a to UDoubleEnsemble
  friend UDoubleEnsemble<ensemble_size> operator/(double a,
                                                  const UDoubleEnsemble<ensemble_size> &b) {
    UDoubleEnsemble<ensemble_size> uda(a);
    return uda /= b;
  }

  UDoubleEnsemble<ensemble_size> &operator+=(const UDoubleEnsemble<ensemble_size> &ud) {
    sources.check_epoch(epoch);
    sources.check_epoch(ud.epoch);

    for (size_t i = 0; i < ensemble_size; i++) ensemble[i] += ud.ensemble[i];
    return *this;
  }

  UDoubleEnsemble<ensemble_size> &operator+=(double d) {
    for (size_t i = 0; i < ensemble_size; i++) ensemble[i] += d;
    return *this;
  }

  UDoubleEnsemble<ensemble_size> &operator-=(const UDoubleEnsemble<ensemble_size> &ud) {
    sources.check_epoch(epoch);
    sources.check_epoch(ud.epoch);

    for (size_t i = 0; i < ensemble_size; i++) ensemble[i] -= ud.ensemble[i];
    return *this;
  }

  UDoubleEnsemble<ensemble_size> &operator-=(double d) {
    for (size_t i = 0; i < ensemble_size; i++) ensemble[i] -= d;
    return *this;
  }

  UDoubleEnsemble<ensemble_size> &operator*=(const UDoubleEnsemble<ensemble_size> &ud) {
    sources.check_epoch(epoch);
    sources.check_epoch(ud.epoch);

    for (size_t i = 0; i < ensemble_size; i++) ensemble[i] *= ud.ensemble[i];
    return *this;
  }

  UDoubleEnsemble<ensemble_size> &operator*=(double d) {
    for (size_t i = 0; i < ensemble_size; i++) ensemble[i] *= d;
    return *this;
  }

  UDoubleEnsemble<ensemble_size> &operator/=(const UDoubleEnsemble<ensemble_size> &ud) {
    sources.check_epoch(epoch);
    sources.check_epoch(ud.epoch);

    for (size_t i = 0; i < ensemble_size; i++) ensemble[i] /= ud.ensemble[i];
    return *this;
  }

  UDoubleEnsemble<ensemble_size> &operator/=(double d) {
    for (size_t i = 0; i < ensemble_size; i++) ensemble[i] /= d;
    return *this;
  }

  // \todo add procedures to make persistent
  friend std::ostream &operator<<(std::ostream &os, const UDoubleEnsemble<ensemble_size> &ud) {
    double mean, sigma, skew, kurtosis, m5;
    moments(ud.ensemble, mean, sigma, skew, kurtosis, m5);
    uncertain_print(mean, sigma, os);

    if (sigma != 0.0) {
      auto original_precision = os.precision();
      auto original_format = os.flags(std::ios::showpoint);
      os << std::setprecision(2) << " [" << skew << " : " << kurtosis << " : " << m5 << "]"
         << std::setprecision(original_precision);
      os.flags(original_format);
    }
    return os;
  }

  friend std::istream &operator>>(std::istream &is, UDoubleEnsemble<ensemble_size> &ud) {
    double mean, sigma;
    uncertain_read(mean, sigma, is);
    ud = UDoubleEnsemble<ensemble_size>(mean, sigma);
    return is;
  }

  static UDoubleEnsemble<ensemble_size> func1(std::function<double(double)> func,
                                              UDoubleEnsemble<ensemble_size> arg) {
    for (size_t i = 0; i < ensemble_size; i++) arg.ensemble[i] = func(arg.ensemble[i]);
    return arg;
  }

  static UDoubleEnsemble<ensemble_size> func2(std::function<double(double, double)> func,
                                              const UDoubleEnsemble<ensemble_size> &arg1,
                                              const UDoubleEnsemble<ensemble_size> &arg2) {
    UDoubleEnsemble<ensemble_size> retval(arg1);
    for (size_t i = 0; i < ensemble_size; i++)
      retval.ensemble[i] = func(arg1.ensemble[i], arg2.ensemble[i]);
    return retval;
  }

  friend UDoubleEnsemble<ensemble_size> sqrt(UDoubleEnsemble<ensemble_size> arg) {
    return func1(static_cast<fun1type>(&std::sqrt), arg);
  }

  friend UDoubleEnsemble<ensemble_size> sin(UDoubleEnsemble<ensemble_size> arg) {
    return func1(static_cast<fun1type>(&std::sin), arg);
  }

  friend UDoubleEnsemble<ensemble_size> cos(UDoubleEnsemble<ensemble_size> arg) {
    return func1(static_cast<fun1type>(&std::cos), arg);
  }

  friend UDoubleEnsemble<ensemble_size> tan(UDoubleEnsemble<ensemble_size> arg) {
    return func1(static_cast<fun1type>(&std::tan), arg);
  }

  friend UDoubleEnsemble<ensemble_size> asin(UDoubleEnsemble<ensemble_size> arg) {
    return func1(static_cast<fun1type>(&std::asin), arg);
  }

  friend UDoubleEnsemble<ensemble_size> acos(UDoubleEnsemble<ensemble_size> arg) {
    return func1(static_cast<fun1type>(&std::acos), arg);
  }

  friend UDoubleEnsemble<ensemble_size> atan(UDoubleEnsemble<ensemble_size> arg) {
    return func1(static_cast<fun1type>(&std::atan), arg);
  }

  friend UDoubleEnsemble<ensemble_size> ceil(UDoubleEnsemble<ensemble_size> arg) {
    return func1(static_cast<fun1type>(&std::ceil), arg);
  }

  friend UDoubleEnsemble<ensemble_size> floor(UDoubleEnsemble<ensemble_size> arg) {
    return func1(static_cast<fun1type>(&std::floor), arg);
  }

  friend UDoubleEnsemble<ensemble_size> fabs(UDoubleEnsemble<ensemble_size> arg) {
    return func1(static_cast<fun1type>(&std::fabs), arg);
  }

  friend UDoubleEnsemble<ensemble_size> exp(UDoubleEnsemble<ensemble_size> arg) {
    return func1(static_cast<fun1type>(&std::exp), arg);
  }

  friend UDoubleEnsemble<ensemble_size> log(UDoubleEnsemble<ensemble_size> arg) {
    return func1(static_cast<fun1type>(&std::log), arg);
  }

  friend UDoubleEnsemble<ensemble_size> log10(UDoubleEnsemble<ensemble_size> arg) {
    return func1(static_cast<fun1type>(&std::log10), arg);
  }

  friend UDoubleEnsemble<ensemble_size> sinh(UDoubleEnsemble<ensemble_size> arg) {
    return func1(static_cast<fun1type>(&std::sinh), arg);
  }

  friend UDoubleEnsemble<ensemble_size> cosh(UDoubleEnsemble<ensemble_size> arg) {
    return func1(static_cast<fun1type>(&std::cosh), arg);
  }

  friend UDoubleEnsemble<ensemble_size> tanh(UDoubleEnsemble<ensemble_size> arg) {
    return func1(static_cast<fun1type>(&std::tanh), arg);
  }

  friend UDoubleEnsemble<ensemble_size> fmod(const UDoubleEnsemble<ensemble_size> &arg1,
                                             const UDoubleEnsemble<ensemble_size> &arg2) {
    return func2(static_cast<fun2type>(&std::fmod), arg1, arg2);
  }

  friend UDoubleEnsemble<ensemble_size> atan2(const UDoubleEnsemble<ensemble_size> &arg1,
                                              const UDoubleEnsemble<ensemble_size> &arg2) {
    return func2(static_cast<fun2type>(&std::atan2), arg1, arg2);
  }

  friend UDoubleEnsemble<ensemble_size> pow(const UDoubleEnsemble<ensemble_size> &arg1,
                                            const UDoubleEnsemble<ensemble_size> &arg2) {
    return func2(static_cast<fun2type>(&std::pow), arg1, arg2);
  }

  friend UDoubleEnsemble<ensemble_size> ldexp(UDoubleEnsemble<ensemble_size> arg,
                                              const int intarg) {
    for (size_t i = 0; i < ensemble_size; i++)
      arg.ensemble[i] = std::ldexp(arg.ensemble[i], intarg);
    return arg;
  }

  friend UDoubleEnsemble<ensemble_size> frexp(UDoubleEnsemble<ensemble_size> arg, int *intarg) {
    // use library frexp on mean to get value of return in second arg
    std::frexp(arg.mean(), intarg);
    for (size_t i = 0; i < ensemble_size; i++) {
      int tempint;  // ignore return in second arg in loop
      arg.ensemble[i] = std::frexp(arg.ensemble[i], &tempint);
    }
    return arg;
  }

  friend UDoubleEnsemble<ensemble_size> modf(UDoubleEnsemble<ensemble_size> arg, double *dblarg) {
    // use library modf on mean to get value of return in second arg
    std::modf(arg.mean(), dblarg);
    for (size_t i = 0; i < ensemble_size; i++) {
      double tempdbl;  // ignore return in second arg in loop
      arg.ensemble[i] = std::modf(arg.ensemble[i], &tempdbl);
    }
    return arg;
  }

  static void new_epoch() {
    sources.new_epoch();
    src_ensemble = {};
    sources = SourceSet("???");
  }

  void print_uncertain_sources(std::ostream &os = std::cout) {
    if (deviation() == 0.0)
      os << "No uncertainty";
    else {
      double unaccounted_uncertainty = 1.0;
      for (size_t i = 0; i < sources.get_num_sources(); i++) {
        double unc_portion = this->correlation(src_ensemble[i]);
        unc_portion *= unc_portion;
        unaccounted_uncertainty -= unc_portion;
        os << sources.get_source_name(i) << ": " << int_percent(unc_portion) << "%" << std::endl;
      }
      os << "other: " << int_percent(unaccounted_uncertainty) << "%" << std::endl;
    }
    os << std::endl;
  }

  double correlation(const UDoubleEnsemble<ensemble_size> &ud, const size_t offset = 0) const {
    size_t i;
    double diff, diff_ud;
    double value = this->mean();
    double ud_value = ud.mean();
    // watch overflow!
    double sum_2_diff = 0.0, sum_2_diff_ud = 0.0, sum_prod_diff = 0.0;

    for (i = 0; i < ensemble_size; i++) {
      diff = ensemble[i] - value;
      sum_2_diff += diff * diff;
      size_t j = (i + offset) % ensemble_size;
      diff_ud = ud.ensemble[j] - ud_value;
      sum_2_diff_ud += diff_ud * diff_ud;
      sum_prod_diff += diff * diff_ud;
    }
    if (!sum_2_diff) return 0.0;
    if (!sum_2_diff_ud) return 0.0;
    if (!sum_prod_diff) return 0.0;
    return sum_prod_diff / std::sqrt(sum_2_diff * sum_2_diff_ud);
  }

  double correlation(const std::vector<double> &ens, const size_t offset = 0) const {
    size_t i;
    double diff, diff_ud;
    double value = this->mean();
    double ud_value = 0.0;
    // watch overflow!
    double sum_2_diff = 0, sum_2_diff_ud = 0, sum_prod_diff = 0;

    for (i = 0; i < ens.size(); i++) ud_value += ens[i];
    ud_value /= ens.size();
    for (i = 0; i < ens.size(); i++) {
      diff = ensemble[i] - value;
      sum_2_diff += diff * diff;
      size_t j = (i + offset) % ens.size();
      diff_ud = ens[j] - ud_value;
      sum_2_diff_ud += diff_ud * diff_ud;
      sum_prod_diff += diff * diff_ud;
    }
    if (!sum_2_diff) return 0.0;
    if (!sum_2_diff_ud) return 0.0;
    if (!sum_prod_diff) return 0.0;
    return sum_prod_diff / std::sqrt(sum_2_diff * sum_2_diff_ud);
  }

  // \todo add function that gives a description
  // \todo allow superimposing histograms
  void print_histogram(std::ostream &os = std::cout) const {
    // centered bins for each 0.5 sigmas from -4 sigmas to +4 sigmas
    // outliers go in the outer bins
    int bin[17] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int i;
    double value = this->mean();
    double sigma = this->deviation();

    if (sigma == 0.0) {
      os << "No histogram when no uncertainty" << std::endl;
      return;
    }
    for (i = 0; i < ensemble_size; i++) {
      double normval = (ensemble[i] - value) / sigma;
      int intval = int(std::floor(2.0 * normval + 0.5) + 8.0);
      if (intval < 0)
        bin[0]++;
      else if (intval > 16)
        bin[16]++;
      else
        bin[intval]++;
    }
    int binmax = 1;
    for (i = 0; i < 17; i++)
      if (bin[i] >= binmax) binmax = bin[i];
    int scale_divisor = 1;
    while (binmax / scale_divisor > 74) scale_divisor++;
    os << "Histogram:  (each * represents ";
    if (scale_divisor == 1)
      os << "1 point)";
    else
      os << scale_divisor << " points)";
    os << setiosflags(std::ios::showpos) << std::endl;
    int bin_display[17];  // number of characters displayed in each bin
    for (i = 0; i < 17; i++) {
      bin_display[i] = int((bin[i] + 0.5) / scale_divisor);
    }
    size_t first_display_bin = 0, last_display_bin = 16;
    while (bin_display[first_display_bin] == 0) first_display_bin++;
    while (bin_display[last_display_bin] == 0) last_display_bin--;
    for (i = first_display_bin; i <= last_display_bin; i++) {
      if (i & 1)
        os << "    | ";
      else
        os << std::setw(3) << (i / 2 - 4) << " + ";
      for (int j = 0; j < bin_display[i]; j++) os << "*";
      os << std::endl;
    }
    os << resetiosflags(std::ios::showpos) << std::endl;
  }

  friend UDoubleEnsemble<ensemble_size> Invoke(double (*certainfunc)(double),
                                               const UDoubleEnsemble<ensemble_size> &arg) {
    UDoubleEnsemble<ensemble_size> retval;

    for (unsigned i = 0; i < ensemble_size; i++) retval.ensemble[i] = certainfunc(arg.ensemble[i]);
    return retval;
  }

  friend UDoubleEnsemble<ensemble_size> Invoke(double (*certainfunc)(double, double),
                                               const UDoubleEnsemble<ensemble_size> &arg1,
                                               const UDoubleEnsemble<ensemble_size> &arg2) {
    UDoubleEnsemble<ensemble_size> retval;

    for (unsigned i = 0; i < ensemble_size; i++)
      retval.ensemble[i] = certainfunc(arg1.ensemble[i], arg2.ensemble[i]);
    return retval;
  }

  void shuffle() {
    for (size_t i = 0; i < ensemble_size - 1; i++) {
      size_t j = i + (size_t)rand() % (ensemble_size - i);
      if (j != i) {
        double temp = ensemble[i];
        ensemble[i] = ensemble[j];
        ensemble[j] = temp;
      }
    }
  }

  // figure the moments (sigma, skew, kurtosis, & 5th moment) from an
  // ensemble given the mean.
  static void moments_fixed_mean(const std::vector<double> &ens, double mean, double &sigma,
                                 double &skew, double &kurtosis, double &m5) {
    size_t i;
    double gddiff2 = 0.0, gddiff3 = 0.0, gddiff4 = 0.0, gddiff5 = 0.0;
    std::vector<double> gddiff;
    gddiff.resize(ens.size());

    // Sorting ensemble by absolute value first increases accuracy
    // but at a cost to performance.  Sorting is O(nlog(n)), whereas
    // the rest of what's being done here is O(n) in ens_size.  Operationally
    // this routine should be called mostly only when there is output,
    // so the cost should be okay.
    for (i = 0; i < ens.size(); i++) gddiff[i] = ens[i] - mean;
    qsort(gddiff.data(), gddiff.size(), sizeof(double), abs_double_compare);
    for (i = 0; i < ens.size(); i++) {
      gddiff2 += gddiff[i] * gddiff[i];
      gddiff3 += gddiff[i] * gddiff[i] * gddiff[i];
      gddiff4 += gddiff[i] * gddiff[i] * gddiff[i] * gddiff[i];
      gddiff5 += gddiff[i] * gddiff[i] * gddiff[i] * gddiff[i] * gddiff[i];
    }
    double var = gddiff2 / ens.size();
    sigma = std::sqrt(var);
    skew = gddiff3 / (var * sigma * ens.size());
    kurtosis = gddiff4 / (var * var * ens.size()) - 3;
    m5 = gddiff5 / (var * var * sigma * ens.size());
    // std::cerr << mean << " +/- " << sigma << " [s: " << skew << ", k: ";
    // std::cerr << kurtosis << " m5: " << m5 << "]" << std::endl;
  }

  // figure the moments (mean, sigma, skew, kurtosis, & 5th moment) from an
  // ensemble
  static void moments(const std::vector<double> &ens, double &mean, double &sigma, double &skew,
                      double &kurtosis, double &m5) {
    size_t i;
    double gdsum = 0.0;
    std::vector<double> gddiff;
    gddiff.resize(ens.size());

    for (i = 0; i < ens.size(); i++) gdsum += ens[i];
    mean = gdsum / ens.size();
    // Sorting ensemble by absolute value first increases accuracy.
    // See note in moments_fixed_mean().
    for (i = 0; i < ens.size(); i++) gddiff[i] = ens[i] - mean;
    qsort(gddiff.data(), gddiff.size(), sizeof(double), abs_double_compare);
    gdsum = 0.0;
    for (i = 0; i < ens.size(); i++) gdsum += gddiff[i];
    mean += gdsum / ens.size();

    moments_fixed_mean(ens, mean, sigma, skew, kurtosis, m5);
  }

  // This function moves points a little bit so the first 5 moments all
  // get measured at precisely the expected values.
  static void PerfectEnsemble(std::vector<double> &ens) {
    size_t i;
    double value, sigma, skew, kurtosis, m5;

    std::vector<double> test_ensemble;
    test_ensemble.resize(ens.size());
    for (int j = 0; j < 3; j++) {
      moments(ens, value, sigma, skew, kurtosis, m5);

      for (i = 0; i < ens.size(); i++) ens[i] -= value;
      for (i = 0; i < ens.size(); i++) ens[i] /= sigma;
      // future work: improve kurtosis correction
      double kurtfact = 0.045;
      for (size_t k = 0; k < 5; k++) {
        for (i = 0; i < ens.size(); i++)
          test_ensemble[i] = ens[i] - kurtfact * kurtosis * ens[i] * ens[i] * ens[i];
        double test_value, test_sigma, test_skew, test_kurtosis;
        double test_m5;
        moments(test_ensemble, test_value, test_sigma, test_skew, test_kurtosis, test_m5);
        kurtfact /= 1 - test_kurtosis / kurtosis;
      }
      for (i = 0; i < ens.size(); i++) ens[i] -= kurtosis * kurtfact * ens[i] * ens[i] * ens[i];
    }
    moments(ens, value, sigma, skew, kurtosis, m5);
    for (i = 0; i < ens.size(); i++) ens[i] -= value;
    for (i = 0; i < ens.size(); i++) ens[i] /= sigma;
  }
};

}  // namespace uncertain
