#include "../include/imageVec.h"
#include <cmath>
#include <gtest/gtest.h>
#include <math.h>

constexpr float THRESHOLD = 0.000001;

TEST(SQRT, 2) {
    double expected = 1.4142135623730;
    double actual = fast_sqrt(2);

    ASSERT_TRUE(std::abs(expected - actual) < THRESHOLD);
}

TEST(SQRT, 125) {
    double expected = 11.180339887498;
    double actual = fast_sqrt(125);

    ASSERT_TRUE(std::abs(expected - actual) < THRESHOLD);
}

TEST(SQRT, 5678123) {
    double expected = 2382.88123917;
    double actual = fast_sqrt(5678123);

    ASSERT_TRUE(std::abs(expected - actual) < THRESHOLD);
}

TEST(SQRT, 2147395599) {
    double expected = 46339.9999892;
    double actual = fast_sqrt(2147395599);

    ASSERT_TRUE(std::abs(expected - actual) < THRESHOLD);
}

TEST(SQRT, 2147395600) {
    double expected = 46340;
    double actual = fast_sqrt(2147395600);

    ASSERT_TRUE(std::abs(expected - actual) < THRESHOLD);
}

TEST(SQRT, DENORMALIZED) {
    constexpr double expected = 1e-155;
    double actual = fast_sqrt(1e-310);

    ASSERT_TRUE(std::abs(expected - actual) < THRESHOLD);
}

TEST(SQRT, 0) {
    constexpr double expected = 0.0;
    double actual = fast_sqrt(0);

    ASSERT_TRUE(std::abs(expected - actual) < THRESHOLD);
}