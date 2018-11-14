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
#include <iostream>
#include <sstream>
#include <iomanip>
#include <time.h>

// \todo add exceptions

// \todo make template variables or class consts
#define MAX_UNC_ELEMENTS 5

namespace uncertain
{

// A class for sources of uncertainties.
class UncertainSourceSet
{
 private:
  unsigned long num_sources;
  unsigned long source_epoch;
  std::string source_name[MAX_UNC_ELEMENTS];
  std::string class_name;
 public:
  UncertainSourceSet(std::string cname = "") : num_sources(0), source_epoch(0)
  {
    for (unsigned int i = 0; i < MAX_UNC_ELEMENTS; i++)
      source_name[i][0] = 0;
    class_name = cname;
  }

  unsigned long get_epoch() const { return source_epoch; }

  void check_epoch(const unsigned long& epoch) const
  {
    if (epoch != source_epoch)
    {
      throw std::runtime_error("Bad epoch: "
                                   + std::to_string(epoch) + " expected: " + std::to_string(source_epoch)
                                   + " in class " + class_name);
    }
  }

  void new_epoch()
  {
    for (unsigned int i = 0; i < num_sources; i++)
      source_name[i][0] = 0;
    source_epoch++;
    num_sources = 0;
  }

  int can_get_new_source() const
  {
    if (num_sources >= MAX_UNC_ELEMENTS)
    {
      std::stringstream ss;
      ss << "Already at maximum number of permissible uncertainty elements: "
         << MAX_UNC_ELEMENTS << "(" << num_sources << ")" << "in class "
         << class_name << ".Change the value of MAX_UNC_ELEMENTS and recompile"
         << " or use new_epoch().";
      throw std::runtime_error(ss.str());
    }
    return 1;
  }

  unsigned long get_new_source(std::string name)
  {
    // std::cerr << "trying to get new source (" << num_sources << ") for "
    //      << name << std::endl;
    (void) can_get_new_source();
    source_name[num_sources] = name;
    return num_sources++;
  }

  unsigned long get_num_sources() const
  {
    return num_sources;
  }

  std::string get_source_name(const unsigned long i) const
  {
    if (i >= num_sources)
    {
      throw std::runtime_error("get_source_name called with illegal source number: "
                                   + std::to_string(i));
    }
    return source_name[i];
  }
};

}