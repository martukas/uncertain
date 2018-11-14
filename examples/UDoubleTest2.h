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

#include <uncertain/UDoubleMS.h>
#include <uncertain/UDoubleMSC.h>
#include <uncertain/UDoubleCT.h>
#include <uncertain/UDoubleEnsemble.h>

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
// \todo add timing of calls to subclasses
// \todo add method to compare to values gotten with +=,
//                   propogateby...()
class UDoubleTest
{
 private:
  UDoubleMSUncorr msu;
  UDoubleMSCorr msc;
  UDoubleMSC<false> mscu;
  UDoubleMSC<true> mscc;
  UDoubleCTSA ctsa;
  UDoubleCTAA ctaa;
  UDoubleEnsemble<ens_a_size> ens_a;
  UDoubleEnsemble<ens_b_size> ens_b;

 public:
  UDoubleTest(const double val = 0.0, const double unc = 0.0,
              std::string name = "")
      : msu(val, unc)
        , msc(val, unc)
        , mscu(val, unc)
        , mscc(val, unc)
        , ctsa(val, unc, name)
        , ctaa(val, unc, name)
        , ens_a(val, unc, name)
        , ens_b(val, unc, name) {}

  UDoubleTest(const UDoubleTest& ud)
      : msu(ud.msu)
        , msc(ud.msc)
        , mscu(ud.mscu)
        , mscc(ud.mscc)
        , ctsa(ud.ctsa)
        , ctaa(ud.ctaa)
        , ens_a(ud.ens_a)
        , ens_b(ud.ens_b) {}

  ~UDoubleTest() {}

  UDoubleTest operator+() const
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

  UDoubleTest operator-() const
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

  friend UDoubleTest operator+(const UDoubleTest& a, const UDoubleTest& b)
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

  friend UDoubleTest operator+(const double& a, const UDoubleTest& b)
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

  friend UDoubleTest operator+(const UDoubleTest& a, const double& b)
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

  friend UDoubleTest operator-(const UDoubleTest& a, const UDoubleTest& b)
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

  friend UDoubleTest operator-(const double& a, const UDoubleTest& b)
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

  friend UDoubleTest operator-(const UDoubleTest& a, const double& b)
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

  friend UDoubleTest operator*(const UDoubleTest& a, const UDoubleTest& b)
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

  friend UDoubleTest operator*(const double& a, const UDoubleTest& b)
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

  friend UDoubleTest operator*(const UDoubleTest& a, const double& b)
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

  friend UDoubleTest operator/(const UDoubleTest& a, const UDoubleTest& b)
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

  friend UDoubleTest operator/(const double& a, const UDoubleTest& b)
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

  friend UDoubleTest operator/(const UDoubleTest& a, const double& b)
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
  UDoubleTest operator++()
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

  UDoubleTest operator--()
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
  UDoubleTest operator++(int)
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

  UDoubleTest operator--(int)
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

  UDoubleTest& operator+=(const UDoubleTest& ud)
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

  UDoubleTest& operator-=(const UDoubleTest& ud)
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

  UDoubleTest& operator*=(const UDoubleTest& ud)
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

  UDoubleTest& operator/=(const UDoubleTest& ud)
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

  UDoubleTest& operator+=(const double& d)
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

  UDoubleTest& operator-=(const double& d)
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

  UDoubleTest& operator*=(const double& d)
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

  UDoubleTest& operator/=(const double& d)
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

  friend std::istream& operator>>(std::istream& is, UDoubleTest& ud)
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

  friend std::stringstream& operator>>(std::stringstream& is, UDoubleTest& ud)
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

  void print_nonrandom_part(std::ostream& os) const
  {
    os << "Uncorrelated:   " << msu << std::endl;
    os << "Correlated:     " << msc << std::endl;
    os << "Better Uncorr:  " << mscu << std::endl;
    os << "Better Corr:    " << mscc << std::endl;
    os << "Corr Track Simp:" << ctsa << std::endl;
    os << "Corr Track Adv: " << ctaa << std::endl;
  }

  // \todo add shortprint() only prints what's different
  friend std::ostream& operator<<(std::ostream& os, const UDoubleTest& ud)
  {
    std::stringstream os_msu, os_msc, os_mscu, os_mscc, os_ctsa, os_ctaa;
    std::stringstream os_ens_a, os_ens_b;

    os_msu << ud.msu;
    std::string str_msu = os_msu.str();
    os_msc << ud.msc;
    std::string str_msc = os_msc.str();
    os_mscu << ud.mscu;
    std::string str_mscu = os_mscu.str();
    os_mscc << ud.mscc;
    std::string str_mscc = os_mscc.str();
    os_ctsa << ud.ctsa;
    std::string str_ctsa = os_ctsa.str();
    os_ctaa << ud.ctaa;
    std::string str_ctaa = os_ctaa.str();
    os_ens_a << ud.ens_a;
    std::string str_ens_a = os_ens_a.str();
    os_ens_b << ud.ens_b;
    std::string str_ens_b = os_ens_b.str();

    size_t len = str_msu.length();  // will ignore higher moments in ensembles
    if (str_msu.compare(0, len, str_msc, 0, len)
        || str_msu.compare(0, len, str_mscu, 0, len)
        || str_msu.compare(0, len, str_mscc, 0, len)
        || str_msu.compare(0, len, str_ctsa, 0, len)
        || str_msu.compare(0, len, str_ctaa, 0, len)
        || str_msu.compare(0, len, str_ens_a, 0, len)
        || str_msu.compare(0, len, str_ens_b, 0, len))
    {
      ud.print_nonrandom_part(os);
      os << "Ensemble<" << ens_a_size << ">: " << ud.ens_a << std::endl;
      os << "Ensemble<" << ens_b_size << ">: " << ud.ens_b << std::endl;
    }
    else
    {  // if all the same print just one
      os << ud.msc;
    }

    return os;
  }

  // \todo add check for equality of ud*.ms*, ud*.cta
  friend void print_multi(std::ostream& os,
                          const UDoubleTest& uda,
                          const UDoubleTest& udb,
                          const UDoubleTest& udc)
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
         std::stringstream os, alt_os; \
         std::stringstream osc, alt_osc; \
         std::stringstream osa, alt_osa; \
         std::string str, slope_str, a_str, alt_a_str; \
         UDoubleMSUncorr alt_msu =PropagateUncertaintiesBySlope(func, arg.msu);\
         arg.msu = func(arg.msu); \
         os << arg.msu; \
         alt_os << alt_msu; \
         if ((slope_str = alt_os.str()) != (str = os.str())) \
            std::cerr << "Warning: different values for " << #func "(): " \
                 << str << " vs. " << slope_str << std::endl; \
         arg.msc = func(arg.msc); \
         \
         UDoubleMSC<false> alt_mscu = PropagateUncertaintiesBySlope(func, arg.mscu);\
         arg.mscu = func(arg.mscu); \
         osc << arg.mscu; \
         alt_osc << alt_mscu; \
         if ((slope_str = alt_osc.str()) != (str = osc.str())) \
            std::cerr << "Warning: different values for curved " << #func "(): " \
                 << str << " vs. " << slope_str << std::endl; \
         arg.mscc = func(arg.mscc); \
         \
         arg.ctsa = func(arg.ctsa); \
         arg.ctaa = func(arg.ctaa); \
         UDoubleEnsemble<ens_a_size> alt_ens_a = Invoke(func, arg.ens_a); \
         arg.ens_a = func(arg.ens_a); \
         osa << arg.ens_a; \
         alt_osa << alt_ens_a; \
         if ((a_str = alt_osa.str()) != (alt_a_str = osa.str())) \
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
         std::stringstream os, alt_os; \
         std::stringstream osc, alt_osc; \
         std::string str, slope_str; \
         UDoubleMSUncorr alt_msu =PropagateUncertaintiesBySlope(func, arg1.msu, \
                                                                arg2.msu);\
         retval.msu = func(arg1.msu, arg2.msu); \
         os << retval.msu; \
         alt_os << alt_msu; \
         if ((slope_str = alt_os.str()) != (str = os.str())) \
            std::cerr << "Warning: different values for " << #func "(): " \
                 << str << " vs. " << slope_str << std::endl; \
         retval.msc = func(arg1.msc, arg2.msc); \
         \
         UDoubleMSC<false> alt_mscu =PropagateUncertaintiesBySlope(func, arg1.mscu, \
                                                               arg2.mscu);\
         retval.mscu = func(arg1.mscu, arg2.mscu); \
         osc << retval.mscu; \
         alt_osc << alt_mscu; \
         if ((slope_str = alt_osc.str()) != (str = osc.str())) \
            std::cerr << "Warning: different values for curved " << #func "(): " \
                 << str << " vs. " << slope_str << std::endl; \
         retval.mscc = func(arg1.mscc, arg2.mscc); \
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
    std::stringstream os, alt_os;
    std::stringstream osc, alt_osc;
    std::string str, slope_str;

    GlobalInt = intarg;
    UDoubleMSUncorr alt_msu = PropagateUncertaintiesBySlope(my_ldexp, arg.msu);
    arg.msu = ldexp(arg.msu, intarg);
    os << arg.msu;
    alt_os << alt_msu;
    if ((slope_str = alt_os.str()) != (str = os.str()))
      std::cerr << "Warning: different values for ldexp(): "
                << str << " vs. " << slope_str << std::endl;
    arg.msc = ldexp(arg.msc, intarg);

    UDoubleMSC<false> alt_mscu = PropagateUncertaintiesBySlope(my_ldexp, arg.mscu);
    arg.mscu = ldexp(arg.mscu, intarg);
    osc << arg.mscu;
    alt_osc << alt_mscu;
    if ((slope_str = alt_osc.str()) != (str = osc.str()))
      std::cerr << "Warning: different values for curved ldexp(): "
                << str << " vs. " << slope_str << std::endl;
    arg.mscc = ldexp(arg.mscc, intarg);

    arg.ctsa = ldexp(arg.ctsa, intarg);
    arg.ctaa = ldexp(arg.ctaa, intarg);
    arg.ens_a = ldexp(arg.ens_a, intarg);
    arg.ens_b = ldexp(arg.ens_b, intarg);
    return arg;
  }

  friend UDoubleTest frexp(UDoubleTest arg, int* intarg)
  {
    std::stringstream os, alt_os;
    std::stringstream osc, alt_osc;
    std::string str, slope_str;

    UDoubleMSUncorr alt_msu = PropagateUncertaintiesBySlope(my_frexp, arg.msu);
    alt_os << alt_msu << " (" << GlobalInt << ")";
    arg.msu = frexp(arg.msu, intarg);
    os << arg.msu << " (" << *intarg << ")";
    if ((slope_str = alt_os.str()) != (str = os.str()))
      std::cerr << "Warning: different values for frexp(): "
                << str << " vs. " << slope_str << std::endl;
    arg.msc = frexp(arg.msc, intarg);

    UDoubleMSC<false> alt_mscu = PropagateUncertaintiesBySlope(my_frexp, arg.mscu);
    alt_osc << alt_mscu << " (" << GlobalInt << ")";
    arg.mscu = frexp(arg.mscu, intarg);
    osc << arg.mscu << " (" << *intarg << ")";
    if ((slope_str = alt_osc.str()) != (str = osc.str()))
      std::cerr << "Warning: different values for curved frexp(): "
                << str << " vs. " << slope_str << std::endl;
    arg.mscc = frexp(arg.mscc, intarg);

    arg.ctsa = frexp(arg.ctsa, intarg);
    arg.ctaa = frexp(arg.ctaa, intarg);
    arg.ens_a = frexp(arg.ens_a, intarg);
    arg.ens_b = frexp(arg.ens_b, intarg);
    return arg;
  }

  friend UDoubleTest modf(UDoubleTest arg, double* dblarg)
  {
    std::stringstream os, alt_os;
    std::stringstream osc, alt_osc;
    std::string str, slope_str;

    UDoubleMSUncorr alt_msu = PropagateUncertaintiesBySlope(my_modf, arg.msu);
    alt_os << alt_msu << " (" << GlobalDouble << ")";
    arg.msu = modf(arg.msu, dblarg);
    os << arg.msu << " (" << *dblarg << ")";
    if ((slope_str = alt_os.str()) != (str = os.str()))
      std::cerr << "Warning: different values for modf(): "
                << str << " vs. " << slope_str << std::endl;
    arg.msc = modf(arg.msc, dblarg);

    UDoubleMSC<false> alt_mscu = PropagateUncertaintiesBySlope(my_modf, arg.mscu);
    alt_osc << alt_mscu << " (" << GlobalDouble << ")";
    arg.mscu = modf(arg.mscu, dblarg);
    osc << arg.mscu << " (" << *dblarg << ")";
    if ((slope_str = alt_osc.str()) != (str = osc.str()))
      std::cerr << "Warning: different values for curved modf(): "
                << str << " vs. " << slope_str << std::endl;
    arg.mscc = modf(arg.mscc, dblarg);

    arg.ctsa = modf(arg.ctsa, dblarg);
    arg.ctaa = modf(arg.ctaa, dblarg);
    arg.ens_a = modf(arg.ens_a, dblarg);
    arg.ens_b = modf(arg.ens_b, dblarg);
    return arg;
  }

  static void new_epoch()
  {
    UDoubleCTSA::new_epoch();
    UDoubleCTAA::new_epoch();
    UDoubleEnsemble<ens_a_size>::new_epoch();
    UDoubleEnsemble<ens_b_size>::new_epoch();
  }

  void print_uncertain_sources(std::ostream& os = std::cout)
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

  void print_correlation(const UDoubleTest ud, std::ostream& os = std::cout) const
  {
    os << "Measurements of correlation:" << std::endl;
    os << "Ensemble<" << ens_a_size << ">: ";
    os << ens_a.correlation(ud.ens_a) << std::endl;
    os << "Ensemble<" << ens_b_size << ">: ";
    os << ens_b.correlation(ud.ens_b) << std::endl;
    os << std::endl;
  }
};

template<>
UncertainSourceSet UDoubleEnsemble<ens_a_size>::sources("Small Ensemble");

template<>
double UDoubleEnsemble<ens_a_size>::src_ensemble[MAX_UNC_ELEMENTS][ens_a_size];

template<>
UncertainSourceSet UDoubleEnsemble<ens_b_size>::sources("Large Ensemble");

template<>
double UDoubleEnsemble<ens_b_size>::src_ensemble[MAX_UNC_ELEMENTS][ens_b_size];

class UDoubleInit
{
 private:
  static unsigned short count; // # of UDoubleInit objects existing
 public:
  UDoubleInit()
  {
    if (count++ == 0)
    {
      // if time() fails, it returns -1.  In that case we just go
      // ahead and seed the random # generator with (unsigned int)-1
      // for lack of anything better.
      srand((unsigned int) time(0));
    }
  }

  ~UDoubleInit()
  {
    if (--count == 0)
    {
      // nothing to be done here now
    }
  }
};

static UDoubleInit udi; // forces construction in every including file
unsigned short UDoubleInit::count;

