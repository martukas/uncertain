#include <uncertain/double_ensemble.hpp>

#include "test_lib/gtest_print.hpp"

static constexpr size_t ens_a_size = 128u;
static constexpr size_t ens_b_size = 1024u;

namespace uncertain
{

using EnsembleSmall = UDoubleEnsemble<ens_a_size>;
using EnsembleLarge = UDoubleEnsemble<ens_b_size>;

template<>
SourceSet EnsembleSmall::sources("Small Ensemble");

template<>
SourceSet EnsembleLarge::sources("Large Ensemble");

template<>
std::vector<std::vector<double>> EnsembleSmall::src_ensemble = {};

template<>
std::vector<double> EnsembleSmall::gauss_ensemble = {};

template<>
std::vector<std::vector<double>> EnsembleLarge::src_ensemble = {};

template<>
std::vector<double> EnsembleLarge::gauss_ensemble = {};

}

class UDoubleEnsembleTest : public TestBase
{
  virtual void SetUp()
  {
    uncertain::EnsembleSmall::new_epoch();
    uncertain::EnsembleLarge::new_epoch();
  }

  virtual void TearDown()
  {

  }
};

TEST_F(UDoubleEnsembleTest, ConstructSmall)
{
  uncertain::EnsembleSmall ud(2.0, 1.0);
  EXPECT_DOUBLE_EQ(ud.mean(), 2.0);
  EXPECT_DOUBLE_EQ(ud.deviation(), 1.0);
}

TEST_F(UDoubleEnsembleTest, SmallNegativeThrows)
{
  EXPECT_ANY_THROW(uncertain::EnsembleSmall(2.0, -1.0));
}

TEST_F(UDoubleEnsembleTest, ConstructLarge)
{
  uncertain::EnsembleLarge ud(2.0, 1.0);
  EXPECT_DOUBLE_EQ(ud.mean(), 2.0);
  EXPECT_DOUBLE_EQ(ud.deviation(), 1.0);
}

TEST_F(UDoubleEnsembleTest, LargeNegativeThrows)
{
  EXPECT_ANY_THROW(uncertain::EnsembleLarge(2.0, -1.0));
}

TEST_F(UDoubleEnsembleTest, Copy)
{
  uncertain::EnsembleSmall ud(2.0, 1.0);
  EXPECT_DOUBLE_EQ(ud.mean(), 2.0);
  EXPECT_DOUBLE_EQ(ud.deviation(), 1.0);

  uncertain::EnsembleSmall ud2(ud);
  EXPECT_DOUBLE_EQ(ud2.mean(), 2.0);
  EXPECT_DOUBLE_EQ(ud2.deviation(), 1.0);
}

TEST_F(UDoubleEnsembleTest, UnaryPlus)
{
  uncertain::EnsembleSmall ud(2.0, 1.0);
  uncertain::EnsembleSmall ud2 = +ud;
  EXPECT_DOUBLE_EQ(ud2.mean(), 2.0);
  EXPECT_DOUBLE_EQ(ud2.deviation(), 1.0);
}

TEST_F(UDoubleEnsembleTest, UnaryNegateSmall)
{
  uncertain::EnsembleSmall ud(2.0, 1.0);
  uncertain::EnsembleSmall ud2 = -ud;
  EXPECT_DOUBLE_EQ(ud2.mean(), -2.0);
  EXPECT_DOUBLE_EQ(ud2.deviation(), 1.0);
}

TEST_F(UDoubleEnsembleTest, UnaryNegateLarge)
{
  uncertain::EnsembleLarge ud(2.0, 1.0);
  uncertain::EnsembleLarge ud2 = -ud;
  EXPECT_DOUBLE_EQ(ud2.mean(), -2.0);
  EXPECT_DOUBLE_EQ(ud2.deviation(), 1.0);
}

TEST_F(UDoubleEnsembleTest, PlusEquals)
{
  uncertain::EnsembleSmall ud(2.0, 1.0);
  uncertain::EnsembleSmall ud2(3.0, 0.5);

  ud += ud2;
  EXPECT_DOUBLE_EQ(ud.mean(), 5.0);
  EXPECT_DOUBLE_EQ(ud.deviation(), 1.1530961);
}

TEST_F(UDoubleEnsembleTest, PlusEqualsLarge)
{
  uncertain::EnsembleLarge ud(2.0, 3.0);
  uncertain::EnsembleLarge ud2(3.0, 4.0);

  ud += ud2;
  EXPECT_DOUBLE_EQ(ud.mean(), 5.0);
  EXPECT_DOUBLE_EQ(ud.deviation(), 5.1418414);
}

TEST_F(UDoubleEnsembleTest, MinusEquals)
{
  uncertain::EnsembleSmall ud(3.0, 1.0);
  uncertain::EnsembleSmall ud2(1.0, 0.5);

  ud -= ud2;
  EXPECT_DOUBLE_EQ(ud.mean(), 2.0);
  EXPECT_DOUBLE_EQ(ud.deviation(), 1.1260865);
}

TEST_F(UDoubleEnsembleTest, MinusEqualsLarge)
{
  uncertain::EnsembleLarge ud(3.0, 3.0);
  uncertain::EnsembleLarge ud2(2.0, 4.0);

  ud -= ud2;
  EXPECT_DOUBLE_EQ(ud.mean(), 1.0);
  EXPECT_DOUBLE_EQ(ud.deviation(), 5.1139555);
}

TEST_F(UDoubleEnsembleTest, DivEquals)
{
  uncertain::EnsembleSmall ud(4.0, 2.0);

  ud /= uncertain::EnsembleSmall(2.0, 1.0);
  EXPECT_DOUBLE_EQ(ud.mean(), 2.3744812);
  EXPECT_DOUBLE_EQ(ud.deviation(), 7.9007583);
}

TEST_F(UDoubleEnsembleTest, DivEqualsReciprocal)
{
  auto ud = uncertain::EnsembleSmall(1.0, 0.0);
  ud /= uncertain::EnsembleSmall(2.0, 1.0);
  EXPECT_DOUBLE_EQ(ud.mean(), 0.63049471);
  EXPECT_DOUBLE_EQ(ud.deviation(), 1.8518502);
}

TEST_F(UDoubleEnsembleTest, DivEqualsLarge)
{
  uncertain::EnsembleLarge ud(8.0, 6.0);

  ud /= uncertain::EnsembleLarge(2.0, 2.0);
  EXPECT_DOUBLE_EQ(ud.mean(), 12.625203);
  EXPECT_DOUBLE_EQ(ud.deviation(), 272.90494);
}

TEST_F(UDoubleEnsembleTest, DivEqualsLargeReciprocal)
{
  uncertain::EnsembleLarge ud(1.0, 0.0);

  ud /= uncertain::EnsembleLarge(2.0, 2.0);
  EXPECT_DOUBLE_EQ(ud.mean(), 0.81952477);
  EXPECT_DOUBLE_EQ(ud.deviation(), 18.990055);
}

TEST_F(UDoubleEnsembleTest, TimesEquals)
{
  uncertain::EnsembleSmall ud(1.0, 0.0);

  ud /= uncertain::EnsembleSmall(2.0, 1.0);
  EXPECT_DOUBLE_EQ(ud.mean(), 0.63049471);
  EXPECT_DOUBLE_EQ(ud.deviation(), 1.8518502);

  auto ud2 = uncertain::EnsembleSmall(2.0, 0.0);
  ud2 /= ud;
  EXPECT_DOUBLE_EQ(ud2.mean(), 4.0);
  EXPECT_DOUBLE_EQ(ud2.deviation(), 2.0);
}

TEST_F(UDoubleEnsembleTest, TimesEqualsLarge)
{
  uncertain::EnsembleLarge ud(1.0, 0.0);

  ud /= uncertain::EnsembleLarge(2.0, 2.0);
  EXPECT_DOUBLE_EQ(ud.mean(), 0.81952477);
  EXPECT_DOUBLE_EQ(ud.deviation(), 18.990055);

  auto ud2 = uncertain::EnsembleLarge(4.0, 5.0);
  ud2 /= ud;
  EXPECT_DOUBLE_EQ(ud2.mean(), 7.330265);
  EXPECT_DOUBLE_EQ(ud2.deviation(), 15.943852);
//  EXPECT_DOUBLE_EQ(ud2.deviation(), 6.0);
}

TEST_F(UDoubleEnsembleTest, Ceiling)
{
  uncertain::EnsembleLarge ud(2.5, 1.0);
  auto ud2 = ceil(ud);
  EXPECT_DOUBLE_EQ(ud2.mean(), 3.0);
  EXPECT_DOUBLE_EQ(ud2.deviation(), 1.0364453);
}

TEST_F(UDoubleEnsembleTest, SqrtSmall)
{
  uncertain::EnsembleSmall ud(64.0, 1.0);

  auto ud2 = sqrt(ud);
  EXPECT_DOUBLE_EQ(ud2.mean(), 7.9997559);
  EXPECT_DOUBLE_EQ(ud2.deviation(), 0.062506683);
}

TEST_F(UDoubleEnsembleTest, SqrtLarge)
{
  uncertain::EnsembleLarge ud(64.0, 2.0);

  auto ud2 = sqrt(ud);
  EXPECT_DOUBLE_EQ(ud2.mean(), 7.9990225);
  EXPECT_DOUBLE_EQ(ud2.deviation(), 0.12505354);
}

TEST_F(UDoubleEnsembleTest, PowSmall)
{
  uncertain::EnsembleSmall ud(8.0, 1.0);
  uncertain::EnsembleSmall ud2(2.0, 0.1);

  auto ud3 = pow(ud, ud2);
  EXPECT_DOUBLE_EQ(ud3.mean(), 66.276054);
  EXPECT_DOUBLE_EQ(ud3.deviation(), 21.70063);
}

TEST_F(UDoubleEnsembleTest, PowLarge)
{
  uncertain::EnsembleLarge ud(8.0, 1.0);
  uncertain::EnsembleLarge ud2(2.0, 0.1);

  auto ud3 = pow(ud, ud2);
  EXPECT_DOUBLE_EQ(ud3.mean(), 66.119255);
  EXPECT_DOUBLE_EQ(ud3.deviation(), 21.194429);
}
