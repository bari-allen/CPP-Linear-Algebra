#include "../include/imageVec.h"
#include <gtest/gtest.h>
#include <math.h>

TEST(SQRT, 2) {
    double expected = 1.4142135623730;
    double actual = fast_sqrt(2);

    ASSERT_TRUE(std::abs(expected - actual) < 0.001);
}

TEST(SQRT, 125) {
    double expected = 11.180339887498;
    double actual = fast_sqrt(125);

    ASSERT_TRUE(std::abs(expected - actual) < 0.001);
}

TEST(SQRT, 5678123) {
    double expected = sqrt(5678123);
    double actual = fast_sqrt(5678123);

    ASSERT_TRUE(std::abs(expected - actual) < 0.001);
}

TEST(SQRT, 2147395599) {
    double expected = sqrt(2147395599);
    double actual = fast_sqrt(2147395599);

    ASSERT_TRUE(std::abs(expected - actual) < 0.001);
}

TEST(SQRT, 2147395600) {
    double expected = sqrt(2147395600);
    double actual = fast_sqrt(2147395600);

    ASSERT_TRUE(std::abs(expected - actual) < 0.001);
}