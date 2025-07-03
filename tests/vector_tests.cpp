#include "../include/imageMat.h"
#include "../include/vector_exception.h"
#include <cmath>
#include <cstdint>
#include <cstring>
#include <gtest/gtest.h>
#include <memory>
#include <sys/types.h>

using fast_math::fast_sqrt;

constexpr float THRESHOLD = 0.000001;

TEST(SQRT, 2) {
    constexpr double expected = 1.4142135623730;
    double actual = fast_sqrt(2);

    ASSERT_TRUE(std::abs(expected - actual) < THRESHOLD);
}

TEST(SQRT, 125) {
    constexpr double expected = 11.180339887498;
    double actual = fast_sqrt(125);

    ASSERT_TRUE(std::abs(expected - actual) < THRESHOLD);
}

TEST(SQRT, 5678123) {
    constexpr double expected = 2382.88123917;
    double actual = fast_sqrt(5678123);

    ASSERT_TRUE(std::abs(expected - actual) < THRESHOLD);
}

TEST(SQRT, 2147395599) {
    constexpr double expected = 46339.9999892;
    double actual = fast_sqrt(2147395599);

    ASSERT_TRUE(std::abs(expected - actual) < THRESHOLD);
}

TEST(SQRT, 2147395600) {
    constexpr double expected = 46340;
    double actual = fast_sqrt(2147395600);

    ASSERT_TRUE(std::abs(expected - actual) < THRESHOLD);
}

TEST(SQRT, DENORMALIZED) {
    constexpr double expected = 0.0;
    double actual = fast_sqrt(1e-310);

    ASSERT_TRUE(std::abs(expected - actual) < THRESHOLD);
}

TEST(SQRT, 0) {
    constexpr double expected = 0.0;
    double actual = fast_sqrt(0);

    ASSERT_TRUE(std::abs(expected - actual) < THRESHOLD);
}

TEST(ADDITION, EXCEPTION) {
    I_Vector<int> lhs(3);
    I_Vector<int> rhs(4);

    ASSERT_THROW(lhs + rhs, vector_exception);
}

TEST(ADDITION, 1D) {
    auto lhs_data = std::make_unique<uint32_t[]>(1);
    auto rhs_data = std::make_unique<uint32_t[]>(1);

    lhs_data[0] = 1;
    rhs_data[0] = 4;

    I_Vector<uint32_t> lhs(1, std::move(lhs_data));
    I_Vector<uint32_t> rhs(1, std::move(rhs_data));
    
    auto result = lhs + rhs;

    ASSERT_EQ(result.get_element(0), 5);
}

TEST(ADDITION, 5D) {
    auto lhs_data = 
        std::unique_ptr<uint32_t[]>(new uint32_t[5]{1, 2, 3, 4, 5});
    auto rhs_data 
        = std::unique_ptr<uint32_t[]>(new uint32_t[5]{2, 3, 4, 5, 6});
    
    I_Vector<uint32_t> lhs(5, std::move(lhs_data));
    I_Vector<uint32_t> rhs(5, std::move(rhs_data));

    auto expected_data = 
        std::unique_ptr<uint32_t[]>(new uint32_t[5]{3, 5, 7, 9, 11});
    I_Vector<uint32_t> expected(5, std::move(expected_data));

    ASSERT_EQ(lhs + rhs, expected);
}

TEST(TRANSPOSE, 3D) {
    uint32_t data[3] = {1, 2, 3};
    auto vec_data = std::make_unique<uint32_t[]>(3);
    std::memcpy(vec_data.get(), data, 3 * sizeof(uint32_t));

    auto mat_data = std::make_unique<uint32_t[]>(3);
    std::memcpy(mat_data.get(), data, 3 * sizeof(uint32_t));

    I_Vector<uint32_t> vec(3, std::move(vec_data));
    I_Matrix<uint32_t> mat(1, 3, std::move(mat_data));

    I_Matrix<uint32_t> transpose = vec.transpose();
    bool test = transpose == mat;
    ASSERT_EQ(transpose, mat);
}

TEST(MULT, 3x) {
    uint32_t data[3] = {1, 5, 20};
    auto vec_data = std::make_unique<uint32_t[]>(3);
    std::memcpy(vec_data.get(), data, 3 * sizeof(uint32_t));
    I_Vector<uint32_t> lhs(3, std::move(vec_data));

    uint32_t expected_data[3] = {3, 15, 60};
    auto unique_expected_data = std::make_unique<uint32_t[]>(3);
    std::memcpy(unique_expected_data.get(), expected_data, 3* sizeof(uint32_t));
    I_Vector<uint32_t> expected(3, std::move(unique_expected_data));

    auto result = lhs * 3;
    ASSERT_EQ(result, expected);
}

TEST(INNER_PRODUCT, 3_DIM) {
    uint32_t lhs_data[3] = {10, 21, 4};
    uint32_t rhs_data[3] = {1, 2, 3};
    uint32_t expected = 64;

    auto unique_lhs_data = std::make_unique<uint32_t[]>(3);
    auto unique_rhs_data = std::make_unique<uint32_t[]>(3);

    std::memcpy(unique_lhs_data.get(), lhs_data, 3 * sizeof(uint32_t));
    std::memcpy(unique_rhs_data.get(), rhs_data, 3 * sizeof(uint32_t));

    I_Vector<uint32_t> lhs(3, std::move(unique_lhs_data));
    I_Vector<uint32_t> rhs(3, std::move(unique_rhs_data));

    uint32_t actual = I_Vector<uint32_t>::dot(lhs, rhs);

    ASSERT_EQ(expected, actual);
}

TEST(INNER_PRODUCT, 0_DIM) {
    I_Vector<uint32_t> lhs{};
    I_Vector<uint32_t> rhs{};

    uint32_t expected = 0;
    uint32_t actual = I_Vector<uint32_t>::dot(lhs, rhs);

    ASSERT_EQ(expected, actual);
}

TEST(INNER_PRODUCT, EXCEPTION) {
    I_Vector<uint32_t> lhs{};
    I_Vector<uint32_t> rhs(4);

    ASSERT_THROW(I_Vector<uint32_t>::dot(lhs, rhs), vector_exception);
}

TEST(OUTER_PRODUCT, 3_DIM) {
    uint32_t lhs_data[3] = {10, 21, 4};
    uint32_t rhs_data[3] = {1, 2, 3};
    uint32_t expected_data[9] = {10, 20, 30, 
                                21, 42, 63, 
                                4, 8, 12};
    
    auto unique_lhs = std::make_unique<uint32_t[]>(3);
    auto unique_rhs = std::make_unique<uint32_t[]>(3);
    auto unique_expected = std::make_unique<uint32_t[]>(9);

    std::memcpy(unique_lhs.get(), lhs_data, 3 * sizeof(uint32_t));
    std::memcpy(unique_rhs.get(), rhs_data, 3 * sizeof(uint32_t));
    std::memcpy(unique_expected.get(), expected_data, 9 * sizeof(uint32_t));

    I_Vector<uint32_t> lhs(3, std::move(unique_lhs));
    I_Matrix<uint32_t> rhs(1, 3, std::move(unique_rhs));
    I_Matrix<uint32_t> expected(3, 3, std::move(unique_expected));

    auto actual = I_Vector<uint32_t>::dot(lhs, rhs);

    ASSERT_EQ(expected, actual);
}

TEST(NORM, 3_DIM) {
    uint32_t data[3] = {10, 21, 4};
    auto unique_data = std::make_unique<uint32_t[]>(3);
    std::memcpy(unique_data.get(), data, 3 * sizeof(uint32_t));
    I_Vector<uint32_t> vec(3, std::move(unique_data));

    double expected = 23.600847;
    double actual = I_Vector<uint32_t>::norm(vec);

    ASSERT_TRUE(std::abs(actual - expected) < 0.0001);
}

TEST(NORM, 4_DIM) {
    uint32_t data[4] = {10, 21, 4, 40};
    auto unique_data = std::make_unique<uint32_t[]>(4);
    std::memcpy(unique_data.get(), data, 4 * sizeof(uint32_t));
    I_Vector<uint32_t> vec(4, std::move(unique_data));

    double expected = 46.443514;
    double actual = I_Vector<uint32_t>::norm(vec);

    ASSERT_TRUE(std::abs(actual - expected) < 0.0001);
}

TEST(DOUBLE_ABS, NEGATIVE) {
    double expected = 0.5;
    double actual = fast_math::abs(-0.5);

    ASSERT_TRUE(std::abs(expected - actual) < 0.001);
}

TEST(DOUBLE_ABS, POSITIVE) {
    double expected = 0.9999;
    double actual = fast_math::abs(0.9999);

    ASSERT_TRUE(std::abs(expected - actual) < 0.001);
}

TEST(INT_ABS, NEGATIVE) {
    uint32_t expected = 10000;
    uint32_t actual = fast_math::abs(-10000);

    ASSERT_TRUE(expected == actual);
}

TEST(INT_ABS, POSITIVE) {
    uint32_t expected = 646;
    uint32_t actual = fast_math::abs(646);

    ASSERT_EQ(expected, actual);
}