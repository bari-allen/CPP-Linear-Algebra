#include "../lib/imageMat.h"
#include <gtest/gtest.h>
#include <memory>

TEST(EQ_TEST, TRUE) {
    int data[5] = {1, 2, 3, 4, 5};

    auto lhs_data = std::make_unique<int[]>(5);
    std::memcpy(lhs_data.get(), data, sizeof(int) * 5);
    I_Matrix<int> lhs(1, 5, lhs_data);

    auto rhs_data = std::make_unique<int[]>(5);
    std::memcpy(rhs_data.get(), data, sizeof(int) * 5);
    I_Matrix<int> rhs(1, 5, rhs_data);

    ASSERT_TRUE(lhs == rhs);
}

TEST(EQ_TEST, FALSE) {
    int data[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

    auto lhs_data = std::make_unique<int[]>(10);
    std::memcpy(lhs_data.get(), data, sizeof(int) * 10);
    I_Matrix<int> lhs(2, 5, lhs_data);

    I_Matrix<int> rhs;

    ASSERT_FALSE(lhs == rhs);
}

TEST(ADDITION, MATRIX) {
    int lhs_data[6] = {0, 1, 2, 3, 4, 5};
    auto unique_lhs_data = std::make_unique<int[]>(6);
    std::memcpy(unique_lhs_data.get(), lhs_data, sizeof(int) * 6);
    I_Matrix<int> lhs(3, 2, unique_lhs_data);

    int rhs_data[6] = {1, 2, 3, 4, 5, 6};
    auto unique_rhs_data = std::make_unique<int[]>(6);
    std::memcpy(unique_rhs_data.get(), rhs_data, sizeof(int) * 6);
    I_Matrix<int> rhs(3, 2, unique_rhs_data);

    int expected_data[6] = {1, 3, 5, 7, 9, 11};
    auto unique_expected = std::make_unique<int[]>(6);
    std::memcpy(unique_expected.get(), expected_data, sizeof(int) * 6);
    I_Matrix<int> expected(3, 2, unique_expected);

    auto actual = lhs + rhs;
    ASSERT_EQ(actual, expected);
}