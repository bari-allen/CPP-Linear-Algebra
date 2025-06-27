#include "../include/imageVec.h"
#include "../include/vector_exception.h"
#include <cmath>
#include <cstdint>
#include <cstring>
#include <gtest/gtest.h>
#include <memory>
#include <sys/types.h>

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