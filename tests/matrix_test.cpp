#include "../lib/imageMat.cpp"
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

