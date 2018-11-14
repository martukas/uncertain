#include <uncertain/double_msc.hpp>
#include "gtest_print.hpp"

class UDoubleMSCTest : public TestBase
{
};

TEST_F(UDoubleMSCTest, ConstructCorrelatedPositive)
{
  uncertain::UDoubleMSCCorr ud(2.0, 1.0);
  EXPECT_FLOAT_EQ(ud.mean(), 2.0);
  EXPECT_FLOAT_EQ(ud.deviation(), 1.0);
}

TEST_F(UDoubleMSCTest, ConstructCorrelatedNegative)
{
  uncertain::UDoubleMSCCorr ud(2.0, -1.0);
  EXPECT_FLOAT_EQ(ud.mean(), 2.0);
  EXPECT_FLOAT_EQ(ud.deviation(), 1.0);
//  EXPECT_FLOAT_EQ(ud.deviation(), -1.0);
}

TEST_F(UDoubleMSCTest, ConstructUncorrelatedPositive)
{
  uncertain::UDoubleMSCUncorr ud(2.0, 1.0);
  EXPECT_FLOAT_EQ(ud.mean(), 2.0);
  EXPECT_FLOAT_EQ(ud.deviation(), 1.0);
}

TEST_F(UDoubleMSCTest, UncorrelatedNegativeThrows)
{
  EXPECT_ANY_THROW(uncertain::UDoubleMSCUncorr(2.0, -1.0));
}

TEST_F(UDoubleMSCTest, Copy)
{
  uncertain::UDoubleMSCCorr ud(2.0, -1.0);
  uncertain::UDoubleMSCCorr ud2(ud);
  EXPECT_FLOAT_EQ(ud2.mean(), 2.0);
  EXPECT_FLOAT_EQ(ud2.deviation(), 1.0);
//  EXPECT_FLOAT_EQ(ud2.deviation(), -1.0);
}

TEST_F(UDoubleMSCTest, UnaryPlus)
{
  uncertain::UDoubleMSCCorr ud(2.0, -1.0);
  uncertain::UDoubleMSCCorr ud2 = +ud;
  EXPECT_FLOAT_EQ(ud2.mean(), 2.0);
  EXPECT_FLOAT_EQ(ud2.deviation(), 1.0);
//  EXPECT_FLOAT_EQ(ud2.deviation(), -1.0);
}

TEST_F(UDoubleMSCTest, UnaryNegateCorrelated)
{
  uncertain::UDoubleMSCCorr ud(2.0, -1.0);
  uncertain::UDoubleMSCCorr ud2 = -ud;
  EXPECT_FLOAT_EQ(ud2.mean(), -2.0);
  EXPECT_FLOAT_EQ(ud2.deviation(), 1.0);
}

TEST_F(UDoubleMSCTest, UnaryNegateUncorrelated)
{
  uncertain::UDoubleMSCUncorr ud(2.0, 1.0);
  uncertain::UDoubleMSCUncorr ud2 = -ud;
  EXPECT_FLOAT_EQ(ud2.mean(), -2.0);
  EXPECT_FLOAT_EQ(ud2.deviation(), 1.0);
}


TEST_F(UDoubleMSCTest, PlusEquals)
{
  uncertain::UDoubleMSCCorr ud(2.0, 1.0);
  uncertain::UDoubleMSCCorr ud2(3.0, 0.5);

  ud += ud2;
  EXPECT_FLOAT_EQ(ud.mean(), 5.0);
  EXPECT_FLOAT_EQ(ud.deviation(), 1.5);
}

TEST_F(UDoubleMSCTest, PlusEqualsUncorr)
{
  uncertain::UDoubleMSCUncorr ud(2.0, 3.0);
  uncertain::UDoubleMSCUncorr ud2(3.0, 4.0);

  ud += ud2;
  EXPECT_FLOAT_EQ(ud.mean(), 5.0);
  EXPECT_FLOAT_EQ(ud.deviation(), 5);
}

TEST_F(UDoubleMSCTest, MinusEquals)
{
  uncertain::UDoubleMSCCorr ud(3.0, 1.0);
  uncertain::UDoubleMSCCorr ud2(1.0, 0.5);

  ud -= ud2;
  EXPECT_FLOAT_EQ(ud.mean(), 2.0);
  EXPECT_FLOAT_EQ(ud.deviation(), 0.5);
}

TEST_F(UDoubleMSCTest, MinusEqualsUncorr)
{
  uncertain::UDoubleMSCUncorr ud(3.0, 3.0);
  uncertain::UDoubleMSCUncorr ud2(2.0, 4.0);

  ud -= ud2;
  EXPECT_FLOAT_EQ(ud.mean(), 1.0);
  EXPECT_FLOAT_EQ(ud.deviation(), 5);
}

TEST_F(UDoubleMSCTest, DivEquals)
{
  uncertain::UDoubleMSCCorr ud(4.0, 2.0);

  ud /= uncertain::UDoubleMSCCorr(2.0, 1.0);
  EXPECT_FLOAT_EQ(ud.mean(), 2.0);
  EXPECT_FLOAT_EQ(ud.deviation(), 0);
}

TEST_F(UDoubleMSCTest, DivEqualsReciprocal)
{
  auto ud = uncertain::UDoubleMSCCorr(1.0, 0.0);
  ud /= uncertain::UDoubleMSCCorr(2.0, 1.0);
  EXPECT_FLOAT_EQ(ud.mean(), 0.625);
//  EXPECT_FLOAT_EQ(ud.mean(), 0.5);
  EXPECT_FLOAT_EQ(ud.deviation(), 0.30618623);
//  EXPECT_FLOAT_EQ(ud.deviation(), -0.25);
}

TEST_F(UDoubleMSCTest, DivEqualsUncorr)
{
  uncertain::UDoubleMSCUncorr ud(8.0, 6.0);

  ud /= uncertain::UDoubleMSCUncorr(2.0, 2.0);
  EXPECT_FLOAT_EQ(ud.mean(), 8.0);
//  EXPECT_FLOAT_EQ(ud.mean(), 4.0);
  EXPECT_FLOAT_EQ(ud.deviation(), 9.16515);
//  EXPECT_FLOAT_EQ(ud.deviation(), 5.0);
}

TEST_F(UDoubleMSCTest, DivEqualsUncorrReciprocal)
{
  uncertain::UDoubleMSCUncorr ud(1.0, 0.0);

  ud /= uncertain::UDoubleMSCUncorr(2.0, 2.0);
  EXPECT_FLOAT_EQ(ud.mean(), 1.0);
//  EXPECT_FLOAT_EQ(ud.mean(), 0.5);
  EXPECT_FLOAT_EQ(ud.deviation(), 0.86602539);
//  EXPECT_FLOAT_EQ(ud.deviation(), 0.5);
}

TEST_F(UDoubleMSCTest, TimesEquals)
{
  uncertain::UDoubleMSCCorr ud(1.0, 0.0);

  ud /= uncertain::UDoubleMSCCorr(2.0, 1.0);
  EXPECT_FLOAT_EQ(ud.mean(), 0.625);
//  EXPECT_FLOAT_EQ(ud.mean(), 0.5);
  EXPECT_FLOAT_EQ(ud.deviation(), 0.30618623);
//  EXPECT_FLOAT_EQ(ud.deviation(), -0.25);

  auto ud2 = uncertain::UDoubleMSCCorr(2.0, 0.0);
  ud2 /= ud;
  EXPECT_FLOAT_EQ(ud2.mean(), 3.9679999);
//  EXPECT_FLOAT_EQ(ud2.mean(), 4.0);
  EXPECT_FLOAT_EQ(ud2.deviation(), 1.90715710);
//  EXPECT_FLOAT_EQ(ud2.deviation(), 2.0);
}

TEST_F(UDoubleMSCTest, TimesEqualsUncorr)
{
  uncertain::UDoubleMSCUncorr ud(1.0, 0.0);

  ud /= uncertain::UDoubleMSCUncorr(2.0, 2.0);
  EXPECT_FLOAT_EQ(ud.mean(), 1.0);
//  EXPECT_FLOAT_EQ(ud.mean(), 0.5);
  EXPECT_FLOAT_EQ(ud.deviation(), 0.86602539);
//  EXPECT_FLOAT_EQ(ud.deviation(), 0.5);

  auto ud2 = uncertain::UDoubleMSCUncorr(4.0, 5.0);
  ud2 /= ud;
  EXPECT_FLOAT_EQ(ud2.mean(), 7.0);
//  EXPECT_FLOAT_EQ(ud2.mean(), 8.0);
  EXPECT_FLOAT_EQ(ud2.deviation(), 10.093314);
//  EXPECT_FLOAT_EQ(ud2.deviation(), 6.0);
}

TEST_F(UDoubleMSCTest, Ceiling)
{
  uncertain::UDoubleMSCUncorr ud(2.5, 1.0);
  auto ud2 = ceil(ud);
  EXPECT_FLOAT_EQ(ud2.mean(), 3.0);
  EXPECT_FLOAT_EQ(ud2.deviation(), 0.0);
}
