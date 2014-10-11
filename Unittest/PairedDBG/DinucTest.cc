#include "PairedDBG/Dinuc.h"

#include <gtest/gtest.h>

using namespace std;

const uint8_t A(0), C(1), G(2), T(3);

bool operator==(const Dinuc& a, const Dinuc& b)
{
	return a.toInt() == b.toInt();
}

TEST(Dinuc, Dinuc)
{
	uint8_t CG = C | G << 2;
	Dinuc cg1(C, G);
	Dinuc cg2(CG);

	EXPECT_EQ(cg1, cg2);
	EXPECT_EQ(cg1.toInt(), CG);
	EXPECT_EQ(cg1.a(), C);
	EXPECT_EQ(cg1.b(), G);

	uint8_t GG = G | G << 2;
	++cg2;
	Dinuc gg(GG);
	EXPECT_EQ(cg2, gg);
	EXPECT_TRUE(cg1 < gg);

	Dinuc gc(G, C);
	EXPECT_EQ(cg1, cg1.reverseComplement());
}

TEST(Dinuc, complementNuc)
{
	EXPECT_EQ(Dinuc::complementNuc(A), T);
	EXPECT_EQ(Dinuc::complementNuc(T), A);
	EXPECT_EQ(Dinuc::complementNuc(C), G);
	EXPECT_EQ(Dinuc::complementNuc(G), C);
}

/* Tests mask() and operator==() as well complement(). */
TEST(DinucSet, general)
{
	Dinuc AT(A, T);
	Dinuc CG(C, G);
	Dinuc GT(G, T);
	Dinuc CC(C, C);

	DinucSet ds;
	ASSERT_FALSE(ds.hasExtension());

	ds.setBase(AT);
	ASSERT_TRUE(ds.hasExtension());
	ds.setBase(CG);
	ds.setBase(GT);

	ASSERT_EQ(ds.outDegree(), 3);
	ASSERT_TRUE(ds.checkBase(AT));
	ASSERT_TRUE(ds.checkBase(CG));
	ASSERT_TRUE(ds.checkBase(GT));
	ASSERT_FALSE(ds.checkBase(CC));

	uint16_t x = 1 << AT.toInt() | 1 << CG.toInt() | 1 << GT.toInt();

	Dinuc AC(A, C);	
	uint16_t y = 1 << AT.toInt() | 1 << CG.toInt() | 1 << AC.toInt();

	DinucSet ds_new = DinucSet::mask(x);
	ASSERT_EQ(ds_new.outDegree(), 3);
	DinucSet ds_new_rc = DinucSet::mask(y);
	EXPECT_EQ(ds, ds_new);
	EXPECT_EQ(ds.complement(), ds_new_rc);

	ds.clear();
	ASSERT_FALSE(ds.hasExtension());

	ds.setBase(AT);
	ds_new.clear(ds);
	ASSERT_EQ(ds_new.outDegree(), 2);
}
