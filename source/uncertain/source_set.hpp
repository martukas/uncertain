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

#include <sstream>
#include <uncertain/functions.hpp>
#include <vector>

namespace uncertain {

// A class for sources of uncertainties.
class SourceSet {
 private:
  size_t source_epoch{0};
  std::vector<std::string> source_names;
  std::string class_name;

 public:
  SourceSet(const std::string &cname = {}) : class_name(cname) {}

  // \todo reserve

  size_t get_epoch() const { return source_epoch; }

  void check_epoch(size_t epoch) const {
    if (epoch != source_epoch) {
      throw std::runtime_error("Wrong epoch: " + std::to_string(epoch) + " expected: " +
                               std::to_string(source_epoch) + " in class " + class_name);
    }
  }

  void new_epoch() {
    source_names.clear();
    source_epoch++;
  }

  size_t get_new_source(std::string name) {
    source_names.push_back(name);
    return source_names.size() - 1;
  }

  size_t get_num_sources() const { return source_names.size(); }

  std::string get_source_name(size_t i) const {
    if (i >= source_names.size()) {
      throw std::runtime_error("get_source_name called with illegal source number: " +
                               std::to_string(i));
    }
    return source_names[i];
  }
};

}  // namespace uncertain
