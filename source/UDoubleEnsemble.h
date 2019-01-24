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

#include "UDouble.h"

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
template <unsigned int esize>
class UDoubleEnsemble
{
private:
   double ensemble[esize];
   static UncertainSourceSet sources;
   static double src_ensemble[MAX_UNC_ELEMENTS][esize];
   unsigned long epoch;

public:
   // The main constructor initializes a new source of uncertainty
   // (if there is uncertainty).
   UDoubleEnsemble(const double val = 0.0, const double unc = 0.0,
                   const char * const name = 0)
              : epoch(sources.get_epoch())
   {
      if (unc < 0.0)
      {
         std::cerr << "Error: negative uncertainty: " << unc << std::endl;
         exit(EXIT_FAILURE);
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
               for (unsigned int i = 0; i < esize / 2; i++)
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
               for (unsigned int i = 0; i < esize / 2; i++)
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
         for (unsigned int i = 0; i < esize; i++)
            ensemble[i] = val + gauss_ensemble[i] * unc;
         char source_name[MAX_SRC_NAME + 1];
         if (name && name[0])
         {
            strncpy(source_name, name, MAX_SRC_NAME);
            source_name[MAX_SRC_NAME] = 0;
         }
         else
         {
            std::ostrstream os;
            os << "anon: ";
            uncertain_print(val, unc, os);
            os << std::ends;
            strncpy(source_name, os.str(), MAX_SRC_NAME);
            os.rdbuf()->freeze(0);
            source_name[MAX_SRC_NAME] = 0;
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
      for (unsigned int i = 0; i < esize; i++)
         ensemble[i] = ud.ensemble[i];
   }
   // constructor from an ensemble.
   // future direction: add similar function that shuffles its input
   UDoubleEnsemble(const double * const newensemble,
                   const char * const name = "")
              : epoch(sources.get_epoch())
   {
      for (unsigned int i = 0; i < esize; i++)
         ensemble[i] = newensemble[i];
      char source_name[MAX_SRC_NAME + 1];
      if (name && name[0])
      {
         strncpy(source_name, name, MAX_SRC_NAME);
         source_name[MAX_SRC_NAME] = 0;
      }
      else
      {
         std::ostrstream os;
         os << "anon from ensemble: " << ensemble[0] << std::ends;
         strncpy(source_name, os.str(), MAX_SRC_NAME);
         os.rdbuf()->freeze(0);
         source_name[MAX_SRC_NAME] = 0;
      }
      unsigned long source_num = sources.get_num_sources();
      memcpy(src_ensemble[source_num], ensemble, sizeof(ensemble));
      source_num = sources.get_new_source(source_name);
   }
// future direction: add constructors with other distributions.

   ~UDoubleEnsemble(void) {}


   UDoubleEnsemble<esize> operator +(void) const
   {
      return *this;
   }
   UDoubleEnsemble<esize> operator -(void) const
   {
      UDoubleEnsemble<esize> retval;
      unsigned int i;

      for (i = 0; i < esize; i++)
         retval.ensemble[i] = -ensemble[i];
      retval.epoch = epoch;
      return retval;
   }
   friend UDoubleEnsemble<esize> operator +(UDoubleEnsemble<esize> a,
                                     const UDoubleEnsemble<esize>& b)
      { return a += b; }
   friend UDoubleEnsemble<esize> operator +(UDoubleEnsemble<esize> a,
                                     const double& b)
      { return a += b; }
   friend UDoubleEnsemble<esize> operator +(const double& b,
                                     UDoubleEnsemble<esize> a)
      { return a += b; }
   friend UDoubleEnsemble<esize> operator -(UDoubleEnsemble<esize> a,
                                     const UDoubleEnsemble<esize>& b)
      { return a -= b; }
   friend UDoubleEnsemble<esize> operator -(UDoubleEnsemble<esize> a,
                                     const double& b)
      { return a -= b; }
   friend UDoubleEnsemble<esize> operator -(const double& b,
                                     UDoubleEnsemble<esize> a)
      { return -(a -= b); }
   UDoubleEnsemble<esize> operator ++(void)
      { return (*this += 1.0); }
   UDoubleEnsemble<esize> operator --(void)
      { return (*this -= 1.0); }
   UDoubleEnsemble<esize> operator ++(int)
   {
      UDoubleEnsemble<esize> retval(*this);
      *this += 1.0;
      return retval;
   }
   UDoubleEnsemble<esize> operator --(int)
   {
      UDoubleEnsemble<esize> retval(*this);
      *this -= 1.0;
      return retval;
   }
   friend UDoubleEnsemble<esize> operator *(UDoubleEnsemble<esize> a,
                                     const UDoubleEnsemble<esize>& b)
      { return a *= b; }
   friend UDoubleEnsemble<esize> operator *(UDoubleEnsemble<esize> a,
                                     const double& b)
      { return a *= b; }
   friend UDoubleEnsemble<esize> operator *(const double& b,
                                     UDoubleEnsemble<esize> a)
      { return a *= b; }
   friend UDoubleEnsemble<esize> operator /(UDoubleEnsemble<esize> a,
                                     const UDoubleEnsemble<esize>& b)
      { return a /= b; }
   friend UDoubleEnsemble<esize> operator /(UDoubleEnsemble<esize> a,
                                     const double& b)
      { return a /= b; }
   // this one promotes a to UDoubleEnsemble
   friend UDoubleEnsemble<esize> operator /(const double& a,
                                     const UDoubleEnsemble<esize>& b)
      { UDoubleEnsemble<esize> uda(a); return uda /= b; }
   UDoubleEnsemble<esize> &operator +=(const UDoubleEnsemble<esize>& ud)
   {
      sources.check_epoch(epoch);
      sources.check_epoch(ud.epoch);

      for (unsigned int i = 0; i < esize; i++)
         ensemble[i] += ud.ensemble[i];
      return *this;
   }
   UDoubleEnsemble<esize> &operator +=(const double& d)
   {
      for (unsigned int i = 0; i < esize; i++)
         ensemble[i] += d;
      return *this;
   }
   UDoubleEnsemble<esize> &operator -=(const UDoubleEnsemble<esize>& ud)
   {
      sources.check_epoch(epoch);
      sources.check_epoch(ud.epoch);

      for (unsigned int i = 0; i < esize; i++)
         ensemble[i] -= ud.ensemble[i];
      return *this;
   }
   UDoubleEnsemble<esize> &operator -=(const double& d)
   {
      for (unsigned int i = 0; i < esize; i++)
         ensemble[i] -= d;
      return *this;
   }
   UDoubleEnsemble<esize> &operator *=(const UDoubleEnsemble<esize>& ud)
   {
      sources.check_epoch(epoch);
      sources.check_epoch(ud.epoch);

      for (unsigned int i = 0; i < esize; i++)
         ensemble[i] *= ud.ensemble[i];
      return *this;
   }
   UDoubleEnsemble<esize> &operator *=(const double& d)
   {
      for (unsigned int i = 0; i < esize; i++)
         ensemble[i] *= d;
      return *this;
   }
   UDoubleEnsemble<esize> &operator /=(const UDoubleEnsemble<esize>& ud)
   {
      sources.check_epoch(epoch);
      sources.check_epoch(ud.epoch);

      for (unsigned int i = 0; i < esize; i++)
         ensemble[i] /= ud.ensemble[i];
      return *this;
   }
   UDoubleEnsemble<esize> &operator /=(const double& d)
   {
      for (unsigned int i = 0; i < esize; i++)
         ensemble[i] /= d;
      return *this;
   }
// future direction: add procedures to make persistant
   friend std::ostream& operator <<(std::ostream &os, const UDoubleEnsemble<esize> &ud)
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
   friend std::istream& operator >>(std::istream &is, UDoubleEnsemble<esize> &ud)
   {
      double mean, sigma;
      uncertain_read(mean, sigma, is);
      ud = UDoubleEnsemble<esize>(mean, sigma);
      return is;
   }
   #define UDoubleEnsemblefunc1(func) \
      UDoubleEnsemble<esize> func(UDoubleEnsemble<esize> arg) \
      { \
         for (unsigned int i = 0; i < esize; i++) \
            arg.ensemble[i] = func(arg.ensemble[i]); \
         return arg; \
      }
   #define UDoubleEnsemblefunc2(func) \
      UDoubleEnsemble<esize> func(const UDoubleEnsemble<esize>& arg1, \
                                  const UDoubleEnsemble<esize>& arg2) \
      { \
         UDoubleEnsemble<esize> retval(arg1); \
         for (unsigned int i = 0; i < esize; i++) \
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
      for (unsigned int i = 0; i < esize; i++)
         arg.ensemble[i] = ldexp(arg.ensemble[i], intarg);
      return arg;
   }
   friend UDoubleEnsemble<esize> frexp(UDoubleEnsemble<esize> arg, int * intarg)
   {
      // use library frexp on mean to get value of return in second arg
      frexp(arg.mean(), intarg);
      for (unsigned int i = 0; i < esize; i++)
      {
         int tempint;  // ignore return in second arg in loop
         arg.ensemble[i] = frexp(arg.ensemble[i], &tempint);
      }
      return arg;
   }
   friend UDoubleEnsemble<esize> modf(UDoubleEnsemble<esize> arg,
                                      double * dblarg)
   {
      // use library modf on mean to get value of return in second arg
      modf(arg.mean(), dblarg);
      for (unsigned int i = 0; i < esize; i++)
      {
         double tempdbl;  // ignore return in second arg in loop
         arg.ensemble[i] = modf(arg.ensemble[i], &tempdbl);
      }
      return arg;
   }


   double mean(void) const
   {
      double sum = 0.0;

      for (unsigned int i = 0; i < esize; i++)
         sum += ensemble[i];
      return sum / esize;
   }
   double deviation(void) const
   {
      unsigned int i;
      double diff, sum_2_diff = 0.0; // watch overflow!
      double value = this->mean();

      for (i = 0; i < esize; i++) {
         diff = ensemble[i] - value;
         sum_2_diff += diff * diff;
      }
      return sqrt(sum_2_diff / esize);
   }
   static void new_epoch(void) { sources.new_epoch(); }
   void print_uncertain_sources(std::ostream &os = std::cout)
   {
      if (deviation() == 0.0)
         os << "No uncertainty";
      else
      {
         double unaccounted_uncertainty = 1.0;
         for (unsigned int i = 0; i < sources.get_num_sources(); i++)
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
                      const unsigned int offset = 0) const
   {
      unsigned int i;
      double diff, diff_ud;
      double value = this->mean();
      double ud_value = ud.mean();
      // watch overflow!
      double sum_2_diff = 0.0, sum_2_diff_ud = 0.0, sum_prod_diff = 0.0;

      for (i = 0; i < esize; i++) {
         diff = ensemble[i] - value;
         sum_2_diff += diff * diff;
         int j = (i + offset) % esize;
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
      return sum_prod_diff / sqrt(sum_2_diff * sum_2_diff_ud);
   }
   double correlation(const double ens[esize],
                      const unsigned int offset = 0) const
   {
      unsigned int i;
      double diff, diff_ud;
      double value = this->mean();
      double ud_value = 0.0;
      // watch overflow!
      double sum_2_diff = 0, sum_2_diff_ud = 0, sum_prod_diff = 0;

      for (i = 0; i < esize; i++)
         ud_value += ens[i];
      ud_value /= esize;
      for (i = 0; i < esize; i++) {
         diff = ensemble[i] - value;
         sum_2_diff += diff * diff;
         int j = (i + offset) % esize;
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
      return sum_prod_diff / sqrt(sum_2_diff * sum_2_diff_ud);
   }
   // future direction: add function that gives a description
   // future direction: allow superimposing histograms
   void print_histogram(std::ostream& os = std::cout) const
   {
      // centered bins for each 0.5 sigmas from -4 sigmas to +4 sigmas
      // outliers go in the outer bins
      int bin[17] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      unsigned int i;
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
      int first_display_bin = 0, last_display_bin = 16;
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
   friend UDoubleEnsemble<esize> Invoke(double (*certainfunc)(double),
                          const UDoubleEnsemble<esize>& arg)
   {
      UDoubleEnsemble<esize> retval;

      for (unsigned i = 0; i < esize; i++)
         retval.ensemble[i] = certainfunc(arg.ensemble[i]);
      return retval;
   }
   friend UDoubleEnsemble<esize> Invoke(double (*certainfunc)(double, double),
                          const UDoubleEnsemble<esize>& arg1,
                          const UDoubleEnsemble<esize>& arg2)
   {
      UDoubleEnsemble<esize> retval;

      for (unsigned i = 0; i < esize; i++)
         retval.ensemble[i] = certainfunc(arg1.ensemble[i], arg2.ensemble[i]);
      return retval;
   }
   void shuffle(void)
   {
      for (unsigned i = 0; i < esize - 1; i++)
      {
         unsigned j = i + (unsigned long)rand() % (esize - i);
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
   static void moments_fixed_mean(double const * const ensemble,
                               const double& mean, double& sigma, double& skew,
                               double& kurtosis, double& m5,
                               const unsigned int ens_size)
   {
      unsigned int i;
      double gddiff2 = 0.0, gddiff3 = 0.0, gddiff4 = 0.0, gddiff5 = 0.0;
      double * gddiff = new double[ens_size];

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
      sigma = sqrt(var);
      skew = gddiff3 / (var * sigma * ens_size);
      kurtosis = gddiff4 / (var * var * ens_size) - 3;
      m5 = gddiff5 / (var * var * sigma * ens_size);
      // std::cerr << mean << " +/- " << sigma << " [s: " << skew << ", k: ";
      // std::cerr << kurtosis << " m5: " << m5 << "]" << std::endl;
      delete [] gddiff;
   }

   // figure the moments (mean, sigma, skew, kurtosis, & 5th moment) from an
   // ensemble
   static void moments(double const * const ensemble,
                    double& mean, double& sigma, double& skew,
                    double& kurtosis, double& m5, const unsigned int ens_size)
   {
      unsigned int i;
      double gdsum = 0.0;
      double * gddiff = new double[ens_size];

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
      delete [] gddiff;
   }


  // This function moves points a little bit so the first 5 moments all
  // get measured at precisely the expected values.
  static void PerfectEnsemble(double * ensemble, const unsigned int ens_size)
  {
     unsigned long i;
     double value, sigma, skew, kurtosis, m5;

     double * test_ensemble = new double[ens_size];
     for (int j = 0; j < 3; j++)
     {
        moments(ensemble, value, sigma, skew, kurtosis, m5, ens_size);

        for (i = 0; i < ens_size; i++)
           ensemble[i] -= value;
        for (i = 0; i < ens_size; i++)
           ensemble[i] /= sigma;
        // future work: improve kurtosis correction
        double kurtfact = 0.045;
        for (unsigned int k = 0; k < 5; k++)
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
     delete [] test_ensemble;
     moments(ensemble, value, sigma, skew, kurtosis, m5, ens_size);
     for (i = 0; i < ens_size; i++)
        ensemble[i] -= value;
     for (i = 0; i < ens_size; i++)
        ensemble[i] /= sigma;
  }
};

// Gimpel lint doesn't like it with consts, so for lint use defines
#define ens_a_size 128u
#define ens_b_size 1024u
// static const unsigned int ens_a_size = 128u;
// static const unsigned int ens_b_size = 1024u;



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


// test class has members of all of above classes
// future direction: add timing of calls to subclasses
// future direction: add method to compare to values gotten with +=,
//                   propogateby...()
class UDoubleTest
{
private:
   UDoubleMSUncorr msu;
   UDoubleMSCorr msc;
   UDoubleMSC<0> mscu;
   UDoubleMSC<1> mscc;
   UDoubleCTSA ctsa;
   UDoubleCTAA ctaa;
   UDoubleEnsemble<ens_a_size> ens_a;
   UDoubleEnsemble<ens_b_size> ens_b;

public:
   UDoubleTest(const double val = 0.0, const double unc = 0.0,
              const char * const name = "")
   : msu(val, unc), msc(val, unc), mscu(val, unc), mscc(val, unc),
     ctsa(val, unc, name), ctaa(val, unc, name),
     ens_a(val, unc, name), ens_b(val, unc, name) {}
   UDoubleTest(const UDoubleTest& ud)
   : msu(ud.msu), msc(ud.msc), mscu(ud.mscu), mscc(ud.mscc),
     ctsa(ud.ctsa), ctaa(ud.ctaa),
     ens_a(ud.ens_a), ens_b(ud.ens_b) {}
   ~UDoubleTest(void) {}

   UDoubleTest operator +(void) const
   {
      UDoubleTest retval;

      retval.msu = +msu;
      retval.msc = +msc;
      retval.mscu = +mscu;
      retval.mscc = +mscc;
      retval.ctsa = +ctsa;
      retval.ctaa = +ctaa;
      retval.ens_a = +ens_a;
      retval.ens_b = +ens_b;
      return retval;
   }
   UDoubleTest operator -(void) const
   {
      UDoubleTest retval;

      retval.msu = -msu;
      retval.msc = -msc;
      retval.mscu = -mscu;
      retval.mscc = -mscc;
      retval.ctsa = -ctsa;
      retval.ctaa = -ctaa;
      retval.ens_a = -ens_a;
      retval.ens_b = -ens_b;
      return retval;
   }
   friend UDoubleTest operator +(const UDoubleTest& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu + b.msu;
      retval.msc = a.msc + b.msc;
      retval.mscu = a.mscu + b.mscu;
      retval.mscc = a.mscc + b.mscc;
      retval.ctsa = a.ctsa + b.ctsa;
      retval.ctaa = a.ctaa + b.ctaa;
      retval.ens_a = a.ens_a + b.ens_a;
      retval.ens_b = a.ens_b + b.ens_b;
      return retval;
   }
   friend UDoubleTest operator +(const double& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a + b.msu;
      retval.msc = a + b.msc;
      retval.mscu = a + b.mscu;
      retval.mscc = a + b.mscc;
      retval.ctsa = a + b.ctsa;
      retval.ctaa = a + b.ctaa;
      retval.ens_a = a + b.ens_a;
      retval.ens_b = a + b.ens_b;
      return retval;
   }
   friend UDoubleTest operator +(const UDoubleTest& a, const double& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu + b;
      retval.msc = a.msc + b;
      retval.mscu = a.mscu + b;
      retval.mscc = a.mscc + b;
      retval.ctsa = a.ctsa + b;
      retval.ctaa = a.ctaa + b;
      retval.ens_a = a.ens_a + b;
      retval.ens_b = a.ens_b + b;
      return retval;
   }
   friend UDoubleTest operator -(const UDoubleTest& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu - b.msu;
      retval.msc = a.msc - b.msc;
      retval.mscu = a.mscu - b.mscu;
      retval.mscc = a.mscc - b.mscc;
      retval.ctsa = a.ctsa - b.ctsa;
      retval.ctaa = a.ctaa - b.ctaa;
      retval.ens_a = a.ens_a - b.ens_a;
      retval.ens_b = a.ens_b - b.ens_b;
      return retval;
   }
   friend UDoubleTest operator -(const double& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a - b.msu;
      retval.msc = a - b.msc;
      retval.mscu = a - b.mscu;
      retval.mscc = a - b.mscc;
      retval.ctsa = a - b.ctsa;
      retval.ctaa = a - b.ctaa;
      retval.ens_a = a - b.ens_a;
      retval.ens_b = a - b.ens_b;
      return retval;
   }
   friend UDoubleTest operator -(const UDoubleTest& a, const double& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu - b;
      retval.msc = a.msc - b;
      retval.mscu = a.mscu - b;
      retval.mscc = a.mscc - b;
      retval.ctsa = a.ctsa - b;
      retval.ctaa = a.ctaa - b;
      retval.ens_a = a.ens_a - b;
      retval.ens_b = a.ens_b - b;
      return retval;
   }
   friend UDoubleTest operator *(const UDoubleTest& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu * b.msu;
      retval.msc = a.msc * b.msc;
      retval.mscu = a.mscu * b.mscu;
      retval.mscc = a.mscc * b.mscc;
      retval.ctsa = a.ctsa * b.ctsa;
      retval.ctaa = a.ctaa * b.ctaa;
      retval.ens_a = a.ens_a * b.ens_a;
      retval.ens_b = a.ens_b * b.ens_b;
      return retval;
   }
   friend UDoubleTest operator *(const double& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a * b.msu;
      retval.msc = a * b.msc;
      retval.mscu = a * b.mscu;
      retval.mscc = a * b.mscc;
      retval.ctsa = a * b.ctsa;
      retval.ctaa = a * b.ctaa;
      retval.ens_a = a * b.ens_a;
      retval.ens_b = a * b.ens_b;
      return retval;
   }
   friend UDoubleTest operator *(const UDoubleTest& a, const double& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu * b;
      retval.msc = a.msc * b;
      retval.mscu = a.mscu * b;
      retval.mscc = a.mscc * b;
      retval.ctsa = a.ctsa * b;
      retval.ctaa = a.ctaa * b;
      retval.ens_a = a.ens_a * b;
      retval.ens_b = a.ens_b * b;
      return retval;
   }
   friend UDoubleTest operator /(const UDoubleTest& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu / b.msu;
      retval.msc = a.msc / b.msc;
      retval.mscu = a.mscu / b.mscu;
      retval.mscc = a.mscc / b.mscc;
      retval.ctsa = a.ctsa / b.ctsa;
      retval.ctaa = a.ctaa / b.ctaa;
      retval.ens_a = a.ens_a / b.ens_a;
      retval.ens_b = a.ens_b / b.ens_b;
      return retval;
   }
   friend UDoubleTest operator /(const double& a, const UDoubleTest& b)
   {
      UDoubleTest retval;

      retval.msu = a / b.msu;
      retval.msc = a / b.msc;
      retval.mscu = a / b.mscu;
      retval.mscc = a / b.mscc;
      retval.ctsa = a / b.ctsa;
      retval.ctaa = a / b.ctaa;
      retval.ens_a = a / b.ens_a;
      retval.ens_b = a / b.ens_b;
      return retval;
   }
   friend UDoubleTest operator /(const UDoubleTest& a, const double& b)
   {
      UDoubleTest retval;

      retval.msu = a.msu / b;
      retval.msc = a.msc / b;
      retval.mscu = a.mscu / b;
      retval.mscc = a.mscc / b;
      retval.ctsa = a.ctsa / b;
      retval.ctaa = a.ctaa / b;
      retval.ens_a = a.ens_a / b;
      retval.ens_b = a.ens_b / b;
      return retval;
   }
   // preincrement and predecrement operators return value after changes
   UDoubleTest operator ++(void)
   {
      ++msu;
      ++msc;
      ++mscu;
      ++mscc;
      ++ctsa;
      ++ctaa;
      ++ens_a;
      ++ens_b;
      return *this;
   }
   UDoubleTest operator --(void)
   {
      --msu;
      --msc;
      --mscu;
      --mscc;
      --ctsa;
      --ctaa;
      --ens_a;
      --ens_b;
      return *this;
   }
   // postincrement and postdecrement operators return value before changes
   UDoubleTest operator ++(int)
   {
      UDoubleTest retval(*this);

      msu++;
      msc++;
      mscu++;
      mscc++;
      ctsa++;
      ctaa++;
      ens_a++;
      ens_b++;
      return retval;
   }
   UDoubleTest operator --(int)
   {
      UDoubleTest retval(*this);

      msu--;
      msc--;
      mscu--;
      mscc--;
      ctsa--;
      ctaa--;
      ens_a--;
      ens_b--;
      return retval;
   }
   UDoubleTest &operator +=(const UDoubleTest& ud)
   {
      msu += ud.msu;
      msc += ud.msc;
      mscu += ud.mscu;
      mscc += ud.mscc;
      ctsa += ud.ctsa;
      ctaa += ud.ctaa;
      ens_a += ud.ens_a;
      ens_b += ud.ens_b;
      return *this;
   }
   UDoubleTest &operator -=(const UDoubleTest& ud)
   {
      msu -= ud.msu;
      msc -= ud.msc;
      mscu -= ud.mscu;
      mscc -= ud.mscc;
      ctsa -= ud.ctsa;
      ctaa -= ud.ctaa;
      ens_a -= ud.ens_a;
      ens_b -= ud.ens_b;
      return *this;
   }
   UDoubleTest &operator *=(const UDoubleTest& ud)
   {
      msu *= ud.msu;
      msc *= ud.msc;
      mscu *= ud.mscu;
      mscc *= ud.mscc;
      ctsa *= ud.ctsa;
      ctaa *= ud.ctaa;
      ens_a *= ud.ens_a;
      ens_b *= ud.ens_b;
      return *this;
   }
   UDoubleTest &operator /=(const UDoubleTest& ud)
   {
      msu /= ud.msu;
      msc /= ud.msc;
      mscu /= ud.mscu;
      mscc /= ud.mscc;
      ctsa /= ud.ctsa;
      ctaa /= ud.ctaa;
      ens_a /= ud.ens_a;
      ens_b /= ud.ens_b;
      return *this;
   }
   UDoubleTest &operator +=(const double& d)
   {
      msu += d;
      msc += d;
      mscu += d;
      mscc += d;
      ctsa += d;
      ctaa += d;
      ens_a += d;
      ens_b += d;
      return *this;
   }
   UDoubleTest &operator -=(const double& d)
   {
      msu -= d;
      msc -= d;
      mscu -= d;
      mscc -= d;
      ctsa -= d;
      ctaa -= d;
      ens_a -= d;
      ens_b -= d;
      return *this;
   }
   UDoubleTest &operator *=(const double& d)
   {
      msu *= d;
      msc *= d;
      mscu *= d;
      mscc *= d;
      ctsa *= d;
      ctaa *= d;
      ens_a *= d;
      ens_b *= d;
      return *this;
   }
   UDoubleTest &operator /=(const double& d)
   {
      msu /= d;
      msc /= d;
      mscu /= d;
      mscc /= d;
      ctsa /= d;
      ctaa /= d;
      ens_a /= d;
      ens_b /= d;
      return *this;
   }
   friend std::istrstream& operator >>(std::istrstream &is, UDoubleTest &ud)
   {
      auto init_pos = is.tellg();        // remember initial position
      std::ios::iostate init_state = is.rdstate(); // remember initial state
      is >> ud.msu;
      auto final_pos = is.tellg();       // remember final position
      std::ios::iostate final_state = is.rdstate();// remember final state

      is.seekg(init_pos);                     // return to initial position
      is.clear(init_state);                   // return to initial state

      is >> ud.msc;
      if (final_pos != is.tellg())
         std::cerr << "msc does not leave the input stream in the same "
                 "position as msu" << std::endl;
      if (final_state != is.rdstate())
         std::cerr << "msc does not leave the input stream in the same "
                 "state as msu" << std::endl;

      is.seekg(init_pos);
      is.clear(init_state);

      is >> ud.mscu;
      if (final_pos != is.tellg())
         std::cerr << "mscu does not leave the input stream in the same "
                 "position as msu" << std::endl;
      if (final_state != is.rdstate())
         std::cerr << "mscu does not leave the input stream in the same "
                 "state as msu" << std::endl;

      is.seekg(init_pos);
      is.clear(init_state);

      is >> ud.mscc;
      if (final_pos != is.tellg())
         std::cerr << "mscc does not leave the input stream in the same "
                 "position as msu" << std::endl;
      if (final_state != is.rdstate())
         std::cerr << "mscc does not leave the input stream in the same "
                 "state as msu" << std::endl;

      is.seekg(init_pos);
      is.clear(init_state);

      is >> ud.ctsa;
      if (final_pos != is.tellg())
         std::cerr << "ctsa does not leave the input stream in the same "
                 "position as msu" << std::endl;
      if (final_state != is.rdstate())
         std::cerr << "ctsa does not leave the input stream in the same "
                 "state as msu" << std::endl;

      is.seekg(init_pos);
      is.clear(init_state);

      is >> ud.ctaa;
      if (final_pos != is.tellg())
         std::cerr << "ctaa does not leave the input stream in the same "
                 "position as msu" << std::endl;
      if (final_state != is.rdstate())
         std::cerr << "ctaa does not leave the input stream in the same "
                 "state as msu" << std::endl;

      is.seekg(init_pos);
      is.clear(init_state);

      is >> ud.ens_a;
      if (final_pos != is.tellg())
         std::cerr << "ens_a does not leave the input stream in the same "
                 "position as msu" << std::endl;
      if (final_state != is.rdstate())
         std::cerr << "ens_a does not leave the input stream in the same "
                 "state as msu" << std::endl;

      is.seekg(init_pos);
      is.clear(init_state);

      is >> ud.ens_b;
      if (final_pos != is.tellg())
         std::cerr << "ens_b does not leave the input stream in the same "
                 "position as msu" << std::endl;
      if (final_state != is.rdstate())
         std::cerr << "ens_b does not leave the input stream in the same "
                 "state as msu" << std::endl;

      return is;
   }
   friend std::istream& operator >>(std::istream &is, UDoubleTest &ud)
   {
      static bool warned = false;
      if (!warned)
      {
         warned = true;
         std::cerr << "warning: UDoubleTest std::istream extractor does not test"
                 " the extractors of the" << std::endl
              << "constituent classes.  Use the std::istrstream extractor for "
                 "those tests." << std::endl;
      }
      double mean, sigma;
      uncertain_read(mean, sigma, is);
      ud = UDoubleTest(mean, sigma);
      return is;
   }
   void print_nonrandom_part(std::ostream &os) const
   {
      os << "Uncorrelated:   " << msu << std::endl;
      os << "Correlated:     " << msc << std::endl;
      os << "Better Uncorr:  " << mscu << std::endl;
      os << "Better Corr:    " << mscc << std::endl;
      os << "Corr Track Simp:" << ctsa << std::endl;
      os << "Corr Track Adv: " << ctaa << std::endl;
   }
   // future direction: add shortprint() only prints what's different
   friend std::ostream& operator <<(std::ostream &os, const UDoubleTest &ud)
   {
      std::ostrstream os_msu, os_msc, os_mscu, os_mscc, os_ctsa, os_ctaa;
      std::ostrstream os_ens_a, os_ens_b;

      os_msu << ud.msu << std::ends; char *str_msu = os_msu.str();
      os_msc << ud.msc << std::ends; char *str_msc = os_msc.str();
      os_mscu << ud.mscu << std::ends; char *str_mscu = os_mscu.str();
      os_mscc << ud.mscc << std::ends; char *str_mscc = os_mscc.str();
      os_ctsa << ud.ctsa << std::ends; char *str_ctsa = os_ctsa.str();
      os_ctaa << ud.ctaa << std::ends; char *str_ctaa = os_ctaa.str();
      os_ens_a << ud.ens_a << std::ends; char *str_ens_a = os_ens_a.str();
      os_ens_b << ud.ens_b << std::ends; char *str_ens_b = os_ens_b.str();

      size_t len = strlen(str_msu);  // will ignore higher moments in ensembles
      if (strncmp(str_msu, str_msc, len)
         || strncmp(str_msu, str_mscu, len) || strncmp(str_msu, str_mscc, len)
         || strncmp(str_msu, str_ctsa, len) || strncmp(str_msu, str_ctaa, len)
         || strncmp(str_msu, str_ens_a, len) || strncmp(str_msu, str_ens_b, len))
      {
         ud.print_nonrandom_part(os);
         os << "Ensemble<" << ens_a_size << ">: " << ud.ens_a << std::endl;
         os << "Ensemble<" << ens_b_size << ">: " << ud.ens_b << std::endl;
      }
      else
      {  // if all the same print just one
         os << ud.msc;
      }
      os_msu.rdbuf()->freeze(0); os_msc.rdbuf()->freeze(0);
      os_mscu.rdbuf()->freeze(0); os_mscc.rdbuf()->freeze(0);
      os_ctsa.rdbuf()->freeze(0); os_ctaa.rdbuf()->freeze(0);
      os_ens_a.rdbuf()->freeze(0); os_ens_b.rdbuf()->freeze(0);

      return os;
   }
   // future direction: add check for equality of ud*.ms*, ud*.cta
   friend void print_multi(std::ostream &os,
                           const UDoubleTest &uda,
                           const UDoubleTest &udb,
                           const UDoubleTest &udc)
   {
      uda.print_nonrandom_part(os);
      os << "Ensemble<" << ens_a_size << ">: " << uda.ens_a << std::endl;
      os << "                " << udb.ens_a << std::endl;
      os << "                " << udc.ens_a << std::endl;
      os << "Ensemble<" << ens_b_size << ">: " << udb.ens_b << std::endl;
      os << "                " << udb.ens_b << std::endl;
      os << "                " << udc.ens_b << std::endl;
   }
   #define UDoubleTestfunc1(func) \
      UDoubleTest func(UDoubleTest arg) \
      { \
         std::ostrstream os, alt_os; \
         std::ostrstream osc, alt_osc; \
         std::ostrstream osa, alt_osa; \
         char *str, *slope_str, *a_str, *alt_a_str; \
         UDoubleMSUncorr alt_msu =PropagateUncertaintiesBySlope(func, arg.msu);\
         arg.msu = func(arg.msu); \
         os << arg.msu << std::ends; \
         alt_os << alt_msu << std::ends; \
         if (strcmp(slope_str = alt_os.str(), str = os.str())) \
            std::cerr << "Warning: different values for " << #func "(): " \
                 << str << " vs. " << slope_str << std::endl; \
         arg.msc = func(arg.msc); \
         os.rdbuf()->freeze(0); alt_os.rdbuf()->freeze(0); \
         \
         UDoubleMSC<0> alt_mscu = PropagateUncertaintiesBySlope(func, arg.mscu);\
         arg.mscu = func(arg.mscu); \
         osc << arg.mscu << std::ends; \
         alt_osc << alt_mscu << std::ends; \
         if (strcmp(slope_str = alt_osc.str(), str = osc.str())) \
            std::cerr << "Warning: different values for curved " << #func "(): " \
                 << str << " vs. " << slope_str << std::endl; \
         arg.mscc = func(arg.mscc); \
         osc.rdbuf()->freeze(0); alt_osc.rdbuf()->freeze(0); \
         \
         arg.ctsa = func(arg.ctsa); \
         arg.ctaa = func(arg.ctaa); \
         UDoubleEnsemble<ens_a_size> alt_ens_a = Invoke(func, arg.ens_a); \
         arg.ens_a = func(arg.ens_a); \
         osa << arg.ens_a << std::ends; \
         alt_osa << alt_ens_a << std::ends; \
         if (strcmp(a_str = alt_osa.str(), alt_a_str = osa.str())) \
            std::cerr << "Warning: different values for ensemble<" << ens_a_size \
                 << "> " << #func "(): " \
                 << a_str << " vs. " << alt_a_str << std::endl; \
         arg.ens_b = func(arg.ens_b); \
         return arg; \
      }
   #define UDoubleTestfunc2(func) \
      UDoubleTest func(const UDoubleTest& arg1, const UDoubleTest& arg2) \
      { \
         UDoubleTest retval; \
         std::ostrstream os, alt_os; \
         std::ostrstream osc, alt_osc; \
         char *str, *slope_str; \
         UDoubleMSUncorr alt_msu =PropagateUncertaintiesBySlope(func, arg1.msu, \
                                                                arg2.msu);\
         retval.msu = func(arg1.msu, arg2.msu); \
         os << retval.msu << std::ends; \
         alt_os << alt_msu << std::ends; \
         if (strcmp(slope_str = alt_os.str(), str = os.str())) \
            std::cerr << "Warning: different values for " << #func "(): " \
                 << str << " vs. " << slope_str << std::endl; \
         retval.msc = func(arg1.msc, arg2.msc); \
         os.rdbuf()->freeze(0); alt_os.rdbuf()->freeze(0); \
         \
         UDoubleMSC<0> alt_mscu =PropagateUncertaintiesBySlope(func, arg1.mscu, \
                                                               arg2.mscu);\
         retval.mscu = func(arg1.mscu, arg2.mscu); \
         osc << retval.mscu << std::ends; \
         alt_osc << alt_mscu << std::ends; \
         if (strcmp(slope_str = alt_osc.str(), str = osc.str())) \
            std::cerr << "Warning: different values for curved " << #func "(): " \
                 << str << " vs. " << slope_str << std::endl; \
         retval.mscc = func(arg1.mscc, arg2.mscc); \
         osc.rdbuf()->freeze(0); alt_osc.rdbuf()->freeze(0); \
         \
         retval.ctsa = func(arg1.ctsa, arg2.ctsa); \
         retval.ctaa = func(arg1.ctaa, arg2.ctaa); \
         retval.ens_a = func(arg1.ens_a, arg2.ens_a); \
         retval.ens_b = func(arg1.ens_b, arg2.ens_b); \
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
      std::ostrstream os, alt_os;
      std::ostrstream osc, alt_osc;
      char *str, *slope_str;

      GlobalInt = intarg;
      UDoubleMSUncorr alt_msu =PropagateUncertaintiesBySlope(my_ldexp, arg.msu);
      arg.msu = ldexp(arg.msu, intarg);
      os << arg.msu << std::ends;
      alt_os << alt_msu << std::ends;
      if (strcmp(slope_str = alt_os.str(), str = os.str()))
         std::cerr << "Warning: different values for ldexp(): "
              << str << " vs. " << slope_str << std::endl;
      arg.msc = ldexp(arg.msc, intarg);
      os.rdbuf()->freeze(0); alt_os.rdbuf()->freeze(0);

      UDoubleMSC<0> alt_mscu =PropagateUncertaintiesBySlope(my_ldexp, arg.mscu);
      arg.mscu = ldexp(arg.mscu, intarg);
      osc << arg.mscu << std::ends;
      alt_osc << alt_mscu << std::ends;
      if (strcmp(slope_str = alt_osc.str(), str = osc.str()))
         std::cerr << "Warning: different values for curved ldexp(): "
              << str << " vs. " << slope_str << std::endl;
      arg.mscc = ldexp(arg.mscc, intarg);
      osc.rdbuf()->freeze(0); alt_osc.rdbuf()->freeze(0);

      arg.ctsa = ldexp(arg.ctsa, intarg);
      arg.ctaa = ldexp(arg.ctaa, intarg);
      arg.ens_a = ldexp(arg.ens_a, intarg);
      arg.ens_b = ldexp(arg.ens_b, intarg);
      return arg;
   }
   friend UDoubleTest frexp(UDoubleTest arg, int * intarg)
   {
      std::ostrstream os, alt_os;
      std::ostrstream osc, alt_osc;
      char *str, *slope_str;

      UDoubleMSUncorr alt_msu =PropagateUncertaintiesBySlope(my_frexp, arg.msu);
      alt_os << alt_msu << " (" << GlobalInt << ")" << std::ends;
      arg.msu = frexp(arg.msu, intarg);
      os << arg.msu << " (" << *intarg << ")" << std::ends;
      if (strcmp(slope_str = alt_os.str(), str = os.str()))
         std::cerr << "Warning: different values for frexp(): "
              << str << " vs. " << slope_str << std::endl;
      arg.msc = frexp(arg.msc, intarg);
      os.rdbuf()->freeze(0); alt_os.rdbuf()->freeze(0);

      UDoubleMSC<0> alt_mscu =PropagateUncertaintiesBySlope(my_frexp, arg.mscu);
      alt_osc << alt_mscu << " (" << GlobalInt << ")" << std::ends;
      arg.mscu = frexp(arg.mscu, intarg);
      osc << arg.mscu << " (" << *intarg << ")" << std::ends;
      if (strcmp(slope_str = alt_osc.str(), str = osc.str()))
         std::cerr << "Warning: different values for curved frexp(): "
              << str << " vs. " << slope_str << std::endl;
      arg.mscc = frexp(arg.mscc, intarg);
      osc.rdbuf()->freeze(0); alt_osc.rdbuf()->freeze(0);

      arg.ctsa = frexp(arg.ctsa, intarg);
      arg.ctaa = frexp(arg.ctaa, intarg);
      arg.ens_a = frexp(arg.ens_a, intarg);
      arg.ens_b = frexp(arg.ens_b, intarg);
      return arg;
   }
   friend UDoubleTest modf(UDoubleTest arg, double * dblarg)
   {
      std::ostrstream os, alt_os;
      std::ostrstream osc, alt_osc;
      char *str, *slope_str;

      UDoubleMSUncorr alt_msu =PropagateUncertaintiesBySlope(my_modf, arg.msu);
      alt_os << alt_msu << " (" << GlobalDouble << ")" << std::ends;
      arg.msu = modf(arg.msu, dblarg);
      os << arg.msu << " (" << *dblarg << ")" << std::ends;
      if (strcmp(slope_str = alt_os.str(), str = os.str()))
         std::cerr << "Warning: different values for modf(): "
              << str << " vs. " << slope_str << std::endl;
      arg.msc = modf(arg.msc, dblarg);
      os.rdbuf()->freeze(0); alt_os.rdbuf()->freeze(0);

      UDoubleMSC<0> alt_mscu =PropagateUncertaintiesBySlope(my_modf, arg.mscu);
      alt_osc << alt_mscu << " (" << GlobalDouble << ")" << std::ends;
      arg.mscu = modf(arg.mscu, dblarg);
      osc << arg.mscu << " (" << *dblarg << ")" << std::ends;
      if (strcmp(slope_str = alt_osc.str(), str = osc.str()))
         std::cerr << "Warning: different values for curved modf(): "
              << str << " vs. " << slope_str << std::endl;
      arg.mscc = modf(arg.mscc, dblarg);
      osc.rdbuf()->freeze(0); alt_osc.rdbuf()->freeze(0);

      arg.ctsa = modf(arg.ctsa, dblarg);
      arg.ctaa = modf(arg.ctaa, dblarg);
      arg.ens_a = modf(arg.ens_a, dblarg);
      arg.ens_b = modf(arg.ens_b, dblarg);
      return arg;
   }

   static void new_epoch(void)
   {
      UDoubleCTSA::new_epoch();
      UDoubleCTAA::new_epoch();
      UDoubleEnsemble<ens_a_size>::new_epoch();
      UDoubleEnsemble<ens_b_size>::new_epoch();
   }
   void print_uncertain_sources(std::ostream &os = std::cout)
   {
      os << "Sources of uncertainty:" << std::endl;
      os << "Corr Tracking Simple:  " << std::endl;
      ctsa.print_uncertain_sources(os);
      os << "Corr Tracking Advanced:  " << std::endl;
      ctaa.print_uncertain_sources(os);
      os << "Ensemble<" << ens_a_size << ">: " << std::endl;
      ens_a.print_uncertain_sources(os);
      os << "Ensemble<" << ens_b_size << ">: " << std::endl;
      ens_b.print_uncertain_sources(os);
      os << std::endl;
   }
   void print_correlation(const UDoubleTest ud, std::ostream &os = std::cout) const
   {
      os << "Measurements of correlation:" << std::endl;
      os << "Ensemble<" << ens_a_size << ">: ";
      os << ens_a.correlation(ud.ens_a) << std::endl;
      os << "Ensemble<" << ens_b_size << ">: ";
      os << ens_b.correlation(ud.ens_b) << std::endl;
      os << std::endl;
   }
};

template <>
UncertainSourceSet UDoubleEnsemble<ens_a_size>::sources("Small Ensemble");

template <>
double UDoubleEnsemble<ens_a_size>::src_ensemble[MAX_UNC_ELEMENTS][ens_a_size];

template <>
UncertainSourceSet UDoubleEnsemble<ens_b_size>::sources("Large Ensemble");

template <>
double UDoubleEnsemble<ens_b_size>::src_ensemble[MAX_UNC_ELEMENTS][ens_b_size];

class UDoubleInit {
private:
   static unsigned short count; // # of UDoubleInit objects existing
public:
   UDoubleInit(void)
   {
      if (count++ == 0) {
         // if time() fails, it returns -1.  In that case we just go
         // ahead and seed the random # generator with (unsigned int)-1
         // for lack of anything better.
         srand((unsigned int)time(0));
      }
   }
   ~UDoubleInit(void)
   {
      if (--count == 0) {
         // nothing to be done here now
      }
   }
};

static UDoubleInit udi; // forces construction in every including file
unsigned short UDoubleInit::count;

