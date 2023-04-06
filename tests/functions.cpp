#include <uncertain/functions.hpp>
#include "test_lib/gtest_print.hpp"

TEST(Functions, HalfPi)
{
  EXPECT_EQ(uncertain::kHalfPi, uncertain::kPi / 2);
}

TEST(Functions, hypot3)
{
  EXPECT_EQ(uncertain::hypot(1, 2, 2), 3);
  EXPECT_EQ(uncertain::hypot(2, 3, 6), 7);
}

TEST(Functions, sqr)
{
  EXPECT_EQ(uncertain::sqr(2), 4);
}
