#include <gtest/gtest.h>
#include "../src/imageMat.h"

TEST(Basic_Test, Equality_Test) {
    ASSERT_EQ(42, int(41.0)) << "Test Didn't Pass";
}