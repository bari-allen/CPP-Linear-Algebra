#include "../lib/imageVec.h"
#include <gtest/gtest.h>
#include <math.h>

TEST(SQRT, 2) {
    float expected = 1.4142135623730;
    float actual = fast_sqrt(2);

    ASSERT_TRUE(std::abs(expected - actual) < 0.0001);
}

TEST(SQRT, 125) {
    float expected = 11.180339887498;
    float actual = fast_sqrt(125);

    ASSERT_TRUE(std::abs(expected - actual) < 0.0001);
}

TEST(SQRT, 5678123) {
    float expected = sqrt(5678123);
    float actual = fast_sqrt(5678123);

    ASSERT_TRUE(std::abs(expected - actual) < 0.0001);
}

TEST(SQRT, 8) {
    float expected = sqrt(2147395599);
    float actual = fast_sqrt(2147395599);

    ASSERT_TRUE(std::abs(expected - actual) < 0.0001);
}