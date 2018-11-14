//
// Authors: Martin Shetty <martin.shetty@esss.se>
//          Morten Jagd Christensen <Morten.JagdChristensen@esss.se>
//

#pragma  once

#include <gtest/gtest.h>

namespace testing
{
namespace internal
{
enum GTestColor {
  COLOR_DEFAULT,
  COLOR_RED,
  COLOR_GREEN,
  COLOR_YELLOW
};

extern void ColoredPrintf(GTestColor color, const char* fmt, ...);
}
}
#define PRINTF(...)  do { testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, "[          ] "); testing::internal::ColoredPrintf(testing::internal::COLOR_YELLOW, __VA_ARGS__); } while(0)

// C++ stream interface
class TestCout : public std::stringstream
{
 public:
  ~TestCout()
  {
    PRINTF("%s\n",str().c_str());
  }
};

#define TEST_COUT  TestCout()