#include "../lib/imageMat.h"
#include <gtest/gtest.h>
#include <memory>

TEST(DOT_PRODUCT, SIMPLE) {
    int lhs_data[3] = {1, 2, 5};
    auto unique_lhs_data = std::make_unique<int[]>(3);
    std::memcpy(unique_lhs_data.get(), lhs_data, sizeof(int) * 3);
    I_Matrix<int> lhs(1, 3, unique_lhs_data);

    int rhs_data[3] = {10, 20, 1};
    auto unique_rhs_data = std::make_unique<int[]>(3);
    std::memcpy(unique_rhs_data.get(), rhs_data, sizeof(int) * 3);
    I_Matrix<int> rhs(3, 1, unique_rhs_data);

    int expected_data[1] = {55};
    auto unique_expected_data = std::make_unique<int[]>(1);
    std::memcpy(unique_expected_data.get(), expected_data, sizeof(int) * 1);
    I_Matrix<int> expected(1, 1, unique_expected_data);

    auto actual = lhs * rhs;
    ASSERT_EQ(actual, expected);
}

TEST(DOT_PRODUCT, INTERMEDIATE) {
    int lhs_data[6] = {1, 2, 5, 10, 0, 20};
    auto unique_lhs_data = std::make_unique<int[]>(6);
    std::memcpy(unique_lhs_data.get(), lhs_data, sizeof(int) * 6);
    I_Matrix<int> lhs(2, 3, unique_lhs_data);

    int rhs_data[6] = {10, 20, 20, 0, 1, 2};
    auto unique_rhs_data = std::make_unique<int[]>(6);
    std::memcpy(unique_rhs_data.get(), rhs_data, sizeof(int) * 6);
    I_Matrix<int> rhs(3, 2, unique_rhs_data);

    int expected_data[4] = {55, 30, 120, 240};
    auto unique_expected_data = std::make_unique<int[]>(4);
    std::memcpy(unique_expected_data.get(), expected_data, sizeof(int) * 4);
    I_Matrix<int> expected(2, 2, unique_expected_data);

    auto actual = lhs * rhs;
    ASSERT_EQ(actual, expected);
}

TEST(DOT_PRODUCT, HARD) {
    int lhs_data[9] = {9, 2, 7, 5, 8, -1, 0, 0, -3};
    auto unique_lhs_data = std::make_unique<int[]>(9);
    std::memcpy(unique_lhs_data.get(), lhs_data, sizeof(int) * 9);
    I_Matrix<int> lhs(3, 3, unique_lhs_data);

    int rhs_data[9] = {-1, 2, 13, -3, 0, 7, -4, -1, 12};
    auto unique_rhs_data = std::make_unique<int[]>(9);
    std::memcpy(unique_rhs_data.get(), rhs_data, sizeof(int) * 9);
    I_Matrix<int> rhs(3, 3, unique_rhs_data);

    int expected_data[9] = {-43, 11, 215, -25, 11, 109, 12, 3, -36};
    auto unique_expected_data = std::make_unique<int[]>(9);
    std::memcpy(unique_expected_data.get(), expected_data, sizeof(int) * 9);
    I_Matrix<int> expected(3, 3, unique_expected_data);

    auto actual = lhs * rhs;
    ASSERT_EQ(actual, expected);
}