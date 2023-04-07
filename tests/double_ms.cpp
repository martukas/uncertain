#include <uncertain/double_ms.hpp>

#include "test_lib/gtest_print.hpp"

class UDoubleMSTest : public TestBase {};

TEST_F(UDoubleMSTest, ConstructCorrelatedPositive) {
  uncertain::UDoubleMSCorr ud(2.0, 1.0);
  EXPECT_FLOAT_EQ(ud.mean(), 2.0);
  EXPECT_FLOAT_EQ(ud.deviation(), 1.0);
}

TEST_F(UDoubleMSTest, ConstructCorrelatedNegative) {
  uncertain::UDoubleMSCorr ud(2.0, -1.0);
  EXPECT_FLOAT_EQ(ud.mean(), 2.0);
  EXPECT_FLOAT_EQ(ud.deviation(), 1.0);
  //  EXPECT_FLOAT_EQ(ud.deviation(), -1.0);
}

TEST_F(UDoubleMSTest, ConstructUncorrelatedPositive) {
  uncertain::UDoubleMSUncorr ud(2.0, 1.0);
  EXPECT_FLOAT_EQ(ud.mean(), 2.0);
  EXPECT_FLOAT_EQ(ud.deviation(), 1.0);
}

TEST_F(UDoubleMSTest, UncorrelatedNegativeThrows) {
  EXPECT_ANY_THROW(uncertain::UDoubleMSUncorr(2.0, -1.0));
}

TEST_F(UDoubleMSTest, Copy) {
  uncertain::UDoubleMSCorr ud(2.0, -1.0);
  uncertain::UDoubleMSCorr ud2(ud);
  EXPECT_FLOAT_EQ(ud2.mean(), 2.0);
  EXPECT_FLOAT_EQ(ud2.deviation(), 1.0);
  //  EXPECT_FLOAT_EQ(ud2.deviation(), -1.0);
}

TEST_F(UDoubleMSTest, UnaryPlus) {
  uncertain::UDoubleMSCorr ud(2.0, -1.0);
  uncertain::UDoubleMSCorr ud2 = +ud;
  EXPECT_FLOAT_EQ(ud2.mean(), 2.0);
  EXPECT_FLOAT_EQ(ud2.deviation(), 1.0);
  //  EXPECT_FLOAT_EQ(ud2.deviation(), -1.0);
}

TEST_F(UDoubleMSTest, UnaryNegateCorrelated) {
  uncertain::UDoubleMSCorr ud(2.0, -1.0);
  uncertain::UDoubleMSCorr ud2 = -ud;
  EXPECT_FLOAT_EQ(ud2.mean(), -2.0);
  EXPECT_FLOAT_EQ(ud2.deviation(), 1.0);
}

TEST_F(UDoubleMSTest, UnaryNegateUncorrelated) {
  uncertain::UDoubleMSUncorr ud(2.0, 1.0);
  uncertain::UDoubleMSUncorr ud2 = -ud;
  EXPECT_FLOAT_EQ(ud2.mean(), -2.0);
  EXPECT_FLOAT_EQ(ud2.deviation(), 1.0);
}

TEST_F(UDoubleMSTest, PlusEquals) {
  uncertain::UDoubleMSCorr ud(2.0, 1.0);
  uncertain::UDoubleMSCorr ud2(3.0, 0.5);

  ud += ud2;
  EXPECT_FLOAT_EQ(ud.mean(), 5.0);
  EXPECT_FLOAT_EQ(ud.deviation(), 1.5);
}

TEST_F(UDoubleMSTest, PlusEqualsUncorr) {
  uncertain::UDoubleMSUncorr ud(2.0, 3.0);
  uncertain::UDoubleMSUncorr ud2(3.0, 4.0);

  ud += ud2;
  EXPECT_FLOAT_EQ(ud.mean(), 5.0);
  EXPECT_FLOAT_EQ(ud.deviation(), 5);
}

TEST_F(UDoubleMSTest, MinusEquals) {
  uncertain::UDoubleMSCorr ud(3.0, 1.0);
  uncertain::UDoubleMSCorr ud2(1.0, 0.5);

  ud -= ud2;
  EXPECT_FLOAT_EQ(ud.mean(), 2.0);
  EXPECT_FLOAT_EQ(ud.deviation(), 0.5);
}

TEST_F(UDoubleMSTest, MinusEqualsUncorr) {
  uncertain::UDoubleMSUncorr ud(3.0, 3.0);
  uncertain::UDoubleMSUncorr ud2(2.0, 4.0);

  ud -= ud2;
  EXPECT_FLOAT_EQ(ud.mean(), 1.0);
  EXPECT_FLOAT_EQ(ud.deviation(), 5);
}

TEST_F(UDoubleMSTest, DivEquals) {
  uncertain::UDoubleMSCorr ud(4.0, 2.0);

  ud /= uncertain::UDoubleMSCorr(2.0, 1.0);
  EXPECT_FLOAT_EQ(ud.mean(), 2.0);
  EXPECT_FLOAT_EQ(ud.deviation(), 0);
}

TEST_F(UDoubleMSTest, DivEqualsReciprocal) {
  auto ud = uncertain::UDoubleMSCorr(1.0, 0.0);
  ud /= uncertain::UDoubleMSCorr(2.0, 1.0);
  EXPECT_FLOAT_EQ(ud.mean(), 0.5);
  EXPECT_FLOAT_EQ(ud.deviation(), 0.25);
  //  EXPECT_FLOAT_EQ(ud.deviation(), -0.25);
}

TEST_F(UDoubleMSTest, DivEqualsUncorr) {
  uncertain::UDoubleMSUncorr ud(8.0, 6.0);

  ud /= uncertain::UDoubleMSUncorr(2.0, 2.0);
  EXPECT_FLOAT_EQ(ud.mean(), 4.0);
  EXPECT_FLOAT_EQ(ud.deviation(), 5.0);
}

TEST_F(UDoubleMSTest, DivEqualsUncorrReciprocal) {
  uncertain::UDoubleMSUncorr ud(1.0, 0.0);

  ud /= uncertain::UDoubleMSUncorr(2.0, 2.0);
  EXPECT_FLOAT_EQ(ud.mean(), 0.5);
  EXPECT_FLOAT_EQ(ud.deviation(), 0.5);
}

TEST_F(UDoubleMSTest, TimesEquals) {
  uncertain::UDoubleMSCorr ud(1.0, 0.0);

  ud /= uncertain::UDoubleMSCorr(2.0, 1.0);
  EXPECT_FLOAT_EQ(ud.mean(), 0.5);
  EXPECT_FLOAT_EQ(ud.deviation(), 0.25);
  //  EXPECT_FLOAT_EQ(ud.deviation(), -0.25);

  auto ud2 = uncertain::UDoubleMSCorr(2.0, 0.0);
  ud2 /= ud;
  EXPECT_FLOAT_EQ(ud2.mean(), 4.0);
  EXPECT_FLOAT_EQ(ud2.deviation(), 2.0);
}

TEST_F(UDoubleMSTest, TimesEqualsUncorr) {
  uncertain::UDoubleMSUncorr ud(1.0, 0.0);

  ud /= uncertain::UDoubleMSUncorr(2.0, 2.0);
  EXPECT_FLOAT_EQ(ud.mean(), 0.5);
  EXPECT_FLOAT_EQ(ud.deviation(), 0.5);

  auto ud2 = uncertain::UDoubleMSUncorr(4.0, 5.0);
  ud2 /= ud;
  EXPECT_FLOAT_EQ(ud2.mean(), 8.0);
  //  EXPECT_FLOAT_EQ(ud2.deviation(), 6.0);
}

TEST_F(UDoubleMSTest, Ceiling) {
  uncertain::UDoubleMSUncorr ud(2.5, 1.0);
  auto ud2 = ceil(ud);
  EXPECT_FLOAT_EQ(ud2.mean(), 3.0);
  EXPECT_FLOAT_EQ(ud2.deviation(), 0.0);
}

TEST_F(UDoubleMSTest, SqrtCorr) {
  uncertain::UDoubleMSCorr ud(4.0, 2.0);

  auto ud2 = sqrt(ud);
  EXPECT_FLOAT_EQ(ud2.mean(), 2.0);
  EXPECT_FLOAT_EQ(ud2.deviation(), 0.5);
}

TEST_F(UDoubleMSTest, SqrtUncorr) {
  uncertain::UDoubleMSUncorr ud(4.0, 2.0);

  auto ud2 = sqrt(ud);
  EXPECT_FLOAT_EQ(ud2.mean(), 2.0);
  EXPECT_FLOAT_EQ(ud2.deviation(), 0.5);
}

TEST_F(UDoubleMSTest, PowCorr) {
  uncertain::UDoubleMSCorr ud(4.0, 2.0);
  uncertain::UDoubleMSCorr ud2(2.0, 0.1);

  auto ud3 = pow(ud, ud2);
  EXPECT_FLOAT_EQ(ud3.mean(), 16.0);
  EXPECT_FLOAT_EQ(ud3.deviation(), 18.218071);
}

TEST_F(UDoubleMSTest, PowUncorr) {
  uncertain::UDoubleMSUncorr ud(4.0, 2.0);
  uncertain::UDoubleMSUncorr ud2(2.0, 0.1);

  auto ud3 = pow(ud, ud2);
  EXPECT_FLOAT_EQ(ud3.mean(), 16.0);
  EXPECT_FLOAT_EQ(ud3.deviation(), 16.153013);
}
