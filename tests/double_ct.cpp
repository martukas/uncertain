#include <uncertain/double_ct.hpp>
#include <uncertain/scaled_array.hpp>
#include <uncertain/simple_array.hpp>

#include "test_lib/gtest_print.hpp"

namespace uncertain {

using UDoubleCTSA = UDoubleCT<SimpleArray>;
using UDoubleCTAA = UDoubleCT<ScaledArray>;

template <>
SourceSet UDoubleCTSA::sources("Simple Array");

template <>
SourceSet UDoubleCTAA::sources("Scaled Array");

}  // namespace uncertain

class UDoubleCTTest : public TestBase {
  virtual void SetUp() {
    uncertain::UDoubleCTSA::new_epoch();
    uncertain::UDoubleCTAA::new_epoch();
  }

  virtual void TearDown() {}
};

TEST_F(UDoubleCTTest, ConstructSimple) {
  uncertain::UDoubleCTSA ud(2.0, 1.0);
  EXPECT_DOUBLE_EQ(ud.mean(), 2.0);
  EXPECT_DOUBLE_EQ(ud.deviation(), 1.0);
}

TEST_F(UDoubleCTTest, SimpleNegativeThrows) { EXPECT_ANY_THROW(uncertain::UDoubleCTSA(2.0, -1.0)); }

TEST_F(UDoubleCTTest, ConstructScaled) {
  uncertain::UDoubleCTAA ud(2.0, 1.0);
  EXPECT_DOUBLE_EQ(ud.mean(), 2.0);
  EXPECT_DOUBLE_EQ(ud.deviation(), 1.0);
}

TEST_F(UDoubleCTTest, ScaledNegativeThrows) { EXPECT_ANY_THROW(uncertain::UDoubleCTAA(2.0, -1.0)); }

TEST_F(UDoubleCTTest, Copy) {
  uncertain::UDoubleCTSA ud(2.0, 1.0);
  uncertain::UDoubleCTSA ud2(ud);
  EXPECT_DOUBLE_EQ(ud2.mean(), 2.0);
  EXPECT_DOUBLE_EQ(ud2.deviation(), 1.0);
}

TEST_F(UDoubleCTTest, UnaryPlus) {
  uncertain::UDoubleCTSA ud(2.0, 1.0);
  uncertain::UDoubleCTSA ud2 = +ud;
  EXPECT_DOUBLE_EQ(ud2.mean(), 2.0);
  EXPECT_DOUBLE_EQ(ud2.deviation(), 1.0);
}

TEST_F(UDoubleCTTest, UnaryNegateSimple) {
  uncertain::UDoubleCTSA ud(2.0, 1.0);
  uncertain::UDoubleCTSA ud2 = -ud;
  EXPECT_DOUBLE_EQ(ud2.mean(), -2.0);
  EXPECT_DOUBLE_EQ(ud2.deviation(), 1.0);
}

TEST_F(UDoubleCTTest, UnaryNegateScaled) {
  uncertain::UDoubleCTAA ud(2.0, 1.0);
  uncertain::UDoubleCTAA ud2 = -ud;
  EXPECT_DOUBLE_EQ(ud2.mean(), -2.0);
  EXPECT_DOUBLE_EQ(ud2.deviation(), 1.0);
}

TEST_F(UDoubleCTTest, PlusEquals) {
  uncertain::UDoubleCTSA ud(2.0, 1.0);
  uncertain::UDoubleCTSA ud2(3.0, 0.5);

  ud += ud2;
  EXPECT_DOUBLE_EQ(ud.mean(), 5.0);
  EXPECT_DOUBLE_EQ(ud.deviation(), 1.118034);
}

TEST_F(UDoubleCTTest, PlusEqualsScaled) {
  uncertain::UDoubleCTAA ud(2.0, 3.0);
  uncertain::UDoubleCTAA ud2(3.0, 4.0);

  ud += ud2;
  EXPECT_DOUBLE_EQ(ud.mean(), 5.0);
  EXPECT_DOUBLE_EQ(ud.deviation(), 5);
}

TEST_F(UDoubleCTTest, MinusEquals) {
  uncertain::UDoubleCTSA ud(3.0, 1.0);
  uncertain::UDoubleCTSA ud2(1.0, 0.5);

  ud -= ud2;
  EXPECT_DOUBLE_EQ(ud.mean(), 2.0);
  EXPECT_DOUBLE_EQ(ud.deviation(), 1.118034);
}

TEST_F(UDoubleCTTest, MinusEqualsScaled) {
  uncertain::UDoubleCTAA ud(3.0, 3.0);
  uncertain::UDoubleCTAA ud2(2.0, 4.0);

  ud -= ud2;
  EXPECT_DOUBLE_EQ(ud.mean(), 1.0);
  EXPECT_DOUBLE_EQ(ud.deviation(), 5);
}

TEST_F(UDoubleCTTest, DivEquals) {
  uncertain::UDoubleCTSA ud(4.0, 2.0);

  ud /= uncertain::UDoubleCTSA(2.0, 1.0);
  EXPECT_DOUBLE_EQ(ud.mean(), 2.0);
  EXPECT_DOUBLE_EQ(ud.deviation(), 1.4142135);
}

TEST_F(UDoubleCTTest, DivEqualsReciprocal) {
  auto ud = uncertain::UDoubleCTSA(1.0, 0.0);
  ud /= uncertain::UDoubleCTSA(2.0, 1.0);
  EXPECT_DOUBLE_EQ(ud.mean(), 0.5);
  EXPECT_DOUBLE_EQ(ud.deviation(), 0.25);
}

TEST_F(UDoubleCTTest, DivEqualsScaled) {
  uncertain::UDoubleCTAA ud(8.0, 6.0);

  ud /= uncertain::UDoubleCTAA(2.0, 2.0);
  EXPECT_DOUBLE_EQ(ud.mean(), 4.0);
  EXPECT_DOUBLE_EQ(ud.deviation(), 5.0);
}

TEST_F(UDoubleCTTest, DivEqualsScaledReciprocal) {
  uncertain::UDoubleCTAA ud(1.0, 0.0);

  ud /= uncertain::UDoubleCTAA(2.0, 2.0);
  EXPECT_DOUBLE_EQ(ud.mean(), 0.5);
  EXPECT_DOUBLE_EQ(ud.deviation(), 0.5);
}

TEST_F(UDoubleCTTest, TimesEqualsSimple) {
  uncertain::UDoubleCTSA ud(1.0, 0.0);

  ud /= uncertain::UDoubleCTSA(2.0, 1.0);
  EXPECT_DOUBLE_EQ(ud.mean(), 0.5);
  EXPECT_DOUBLE_EQ(ud.deviation(), 0.25);

  auto ud2 = uncertain::UDoubleCTSA(4.0, 1.0);
  ud2 /= ud;
  EXPECT_DOUBLE_EQ(ud2.mean(), 8.0);
  EXPECT_DOUBLE_EQ(ud2.deviation(), 4.472136);
}

TEST_F(UDoubleCTTest, TimesEqualsScaled) {
  uncertain::UDoubleCTAA ud(1.0, 0.0);

  ud /= uncertain::UDoubleCTAA(2.0, 1.0);
  EXPECT_DOUBLE_EQ(ud.mean(), 0.5);
  EXPECT_DOUBLE_EQ(ud.deviation(), 0.25);

  auto ud2 = uncertain::UDoubleCTAA(4.0, 1.0);
  ud2 /= ud;
  EXPECT_DOUBLE_EQ(ud2.mean(), 8.0);
  EXPECT_DOUBLE_EQ(ud2.deviation(), 4.472136);
}

TEST_F(UDoubleCTTest, Ceiling) {
  uncertain::UDoubleCTAA ud(2.5, 1.0);
  auto ud2 = ceil(ud);
  EXPECT_DOUBLE_EQ(ud2.mean(), 3.0);
  EXPECT_DOUBLE_EQ(ud2.deviation(), 0.0);
}

TEST_F(UDoubleCTTest, SqrtSimple) {
  uncertain::UDoubleCTSA ud(4.0, 2.0);

  auto ud2 = sqrt(ud);
  EXPECT_DOUBLE_EQ(ud2.mean(), 2.0);
  EXPECT_DOUBLE_EQ(ud2.deviation(), 0.5);
}

TEST_F(UDoubleCTTest, SqrtScaled) {
  uncertain::UDoubleCTAA ud(4.0, 2.0);

  auto ud2 = sqrt(ud);
  EXPECT_DOUBLE_EQ(ud2.mean(), 2.0);
  EXPECT_DOUBLE_EQ(ud2.deviation(), 0.5);
}

TEST_F(UDoubleCTTest, PowSimple) {
  uncertain::UDoubleCTSA ud(4.0, 2.0);
  uncertain::UDoubleCTSA ud2(2.0, 0.1);

  auto ud3 = pow(ud, ud2);
  EXPECT_DOUBLE_EQ(ud3.mean(), 16);
  EXPECT_DOUBLE_EQ(ud3.deviation(), 16.153013);
}

TEST_F(UDoubleCTTest, PowScaled) {
  uncertain::UDoubleCTAA ud(4.0, 2.0);
  uncertain::UDoubleCTAA ud2(2.0, 0.1);

  auto ud3 = pow(ud, ud2);
  EXPECT_DOUBLE_EQ(ud3.mean(), 16);
  EXPECT_DOUBLE_EQ(ud3.deviation(), 16.153013);
}
