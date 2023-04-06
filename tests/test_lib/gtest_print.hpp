#pragma once
#include <gtest/gtest.h>
#include "test_lib/color_bash.hpp"

class TestBase : public ::testing::Test
{
protected:
  class Message : public std::stringstream
  {
  public:
    ~Message()
    {
      std::cout << col(BashColor::GREEN) <<  "[          ] ";
      std::cout << col(BashColor::YELLOW) <<  str();
    }
  };
#define MESSAGE Message
};
