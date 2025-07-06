#include "../include/imageMat.h"
#include <cstdint>
#include <gtest/gtest.h>
#include <memory>
#include <chrono>
#include <stdexcept>
#include <sys/types.h>

TEST(DOT_PRODUCT, SIMPLE) {
    int lhs_data[3] = {1, 2, 5};
    auto unique_lhs_data = std::make_unique<int[]>(3);
    std::memcpy(unique_lhs_data.get(), lhs_data, sizeof(int) * 3);
    I_Matrix<int> lhs(1, 3, std::move(unique_lhs_data));

    int rhs_data[3] = {10, 20, 1};
    auto unique_rhs_data = std::make_unique<int[]>(3);
    std::memcpy(unique_rhs_data.get(), rhs_data, sizeof(int) * 3);
    I_Matrix<int> rhs(3, 1, std::move(unique_rhs_data));

    int expected_data[1] = {55};
    auto unique_expected_data = std::make_unique<int[]>(1);
    std::memcpy(unique_expected_data.get(), expected_data, sizeof(int) * 1);
    I_Matrix<int> expected(1, 1, std::move(unique_expected_data));

    auto actual = lhs * rhs;
    
    ASSERT_EQ(actual, expected);
}

TEST(DOT_PRODUCT, INTERMEDIATE) {
    int lhs_data[6] = {1, 2, 5, 10, 0, 20};
    auto unique_lhs_data = std::make_unique<int[]>(6);
    std::memcpy(unique_lhs_data.get(), lhs_data, sizeof(int) * 6);
    I_Matrix<int> lhs(2, 3, std::move(unique_lhs_data));

    int rhs_data[6] = {10, 20, 20, 0, 1, 2};
    auto unique_rhs_data = std::make_unique<int[]>(6);
    std::memcpy(unique_rhs_data.get(), rhs_data, sizeof(int) * 6);
    I_Matrix<int> rhs(3, 2, std::move(unique_rhs_data));

    int expected_data[4] = {55, 30, 120, 240};
    auto unique_expected_data = std::make_unique<int[]>(4);
    std::memcpy(unique_expected_data.get(), expected_data, sizeof(int) * 4);
    I_Matrix<int> expected(2, 2, std::move(unique_expected_data));

    auto actual = lhs * rhs;
    ASSERT_EQ(actual, expected);
}

TEST(DOT_PRODUCT, HARD) {
    int lhs_data[9] = {9, 2, 7, 5, 8, -1, 0, 0, -3};
    auto unique_lhs_data = std::make_unique<int[]>(9);
    std::memcpy(unique_lhs_data.get(), lhs_data, sizeof(int) * 9);
    I_Matrix<int> lhs(3, 3, std::move(unique_lhs_data));

    int rhs_data[9] = {-1, 2, 13, -3, 0, 7, -4, -1, 12};
    auto unique_rhs_data = std::make_unique<int[]>(9);
    std::memcpy(unique_rhs_data.get(), rhs_data, sizeof(int) * 9);
    I_Matrix<int> rhs(3, 3, std::move(unique_rhs_data));

    int expected_data[9] = {-43, 11, 215, -25, 11, 109, 12, 3, -36};
    auto unique_expected_data = std::make_unique<int[]>(9);
    std::memcpy(unique_expected_data.get(), expected_data, sizeof(int) * 9);
    I_Matrix<int> expected(3, 3, std::move(unique_expected_data));

    auto actual = lhs * rhs;
    ASSERT_EQ(actual, expected);
}

TEST(DOT_PRODUCT, VERY_HARD) {
    int lhs_data[12] = {-1, 1, -1, 5, 2, -5, 6, -5, 1, -5, 6, 0};
    auto unique_lhs_data = std::make_unique<int[]>(12);
    std::memcpy(unique_lhs_data.get(), lhs_data, sizeof(int) * 12);
    I_Matrix<int> lhs(4, 3, std::move(unique_lhs_data));

    int rhs_data[6] = {6, 5, 5, -6, 6, 0};
    auto unique_rhs_data = std::make_unique<int[]>(6);
    std::memcpy(unique_rhs_data.get(), rhs_data, sizeof(int) * 6);
    I_Matrix<int> rhs(3, 2, std::move(unique_rhs_data));

    int expected_data[8] = {-7, -11, 10, 13, 17, 60, 0, -61};
    auto unique_expected_data = std::make_unique<int[]>(8);
    std::memcpy(unique_expected_data.get(), expected_data, sizeof(int) * 8);
    I_Matrix<int> expected(4, 2, std::move(unique_expected_data));

    auto start = std::chrono::high_resolution_clock::now();
    
    auto actual = lhs * rhs;

    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    ASSERT_EQ(actual, expected);
    std::cout << "Time Taken: " << duration.count() << " microseconds" <<std::endl;
}

TEST(DOT_PRODUCT, EXTREME) {
    auto start = std::chrono::high_resolution_clock::now();

    int lhs_data[81] = {1, 2, 3, 4, 5, 6, 7, 8, 9,
                        10, 11, 12, 13, 14, 15, 16, 17, 18,
                        19, 20, 21, 22, 23, 24, 25, 26, 27,
                        28, 29, 30, 31, 32, 33, 34, 35, 36,
                        37, 38, 39, 40, 41, 42, 43, 44, 45,
                        46, 47, 48, 49, 50, 51, 52, 53, 54,
                        55, 56, 57, 58, 59, 60, 61, 62, 63,
                        64, 65, 66, 67, 68, 69, 70, 71, 72,
                        73, 74, 75, 76, 77, 78, 79, 80, 81};
    auto unique_lhs_data = std::make_unique<int[]>(81);
    std::memcpy(unique_lhs_data.get(), lhs_data, sizeof(int) * 81);
    I_Matrix<int> lhs(9, 9, std::move(unique_lhs_data));

    int rhs_data[81] = {81, 80, 79, 78, 77, 76, 75, 74, 73,
                        72, 71, 70, 69, 68, 67, 66, 65, 64,
                        63, 62, 61, 60, 59, 58, 57, 56, 55,
                        54, 53, 52, 51, 50, 49, 48, 47, 46,
                        45, 44, 43, 42, 41, 40, 39, 38, 37,
                        36, 35, 34, 33, 32, 31, 30, 29, 28,
                        27, 26, 25, 24, 23, 22, 21, 20, 19,
                        18, 17, 16, 15, 14, 13, 12, 11, 10,
                        9, 8, 7, 6, 5, 4, 3, 2, 1};
    auto unique_rhs_data = std::make_unique<int[]>(81);
    std::memcpy(unique_rhs_data.get(), rhs_data, sizeof(int) * 81);
    I_Matrix<int> rhs(9, 9, std::move(unique_rhs_data));

    int expected_data[81] = {1485, 1440, 1395, 1350, 1305, 1260, 1215, 1170, 1125,
                            5130, 5004, 4878, 4752, 4626, 4500, 4374, 4248, 4122,
                            8775, 8568, 8361, 8154, 7947, 7740, 7533, 7326, 7119,
                            12420, 12132, 11844, 11556, 11268, 10980, 10692, 10404, 10116,
                            16065, 15696, 15327, 14958, 14589, 14220, 13851, 13482, 13113,
                            19710, 19260, 18810, 18360, 17910, 17460, 17010, 16560, 16110,
                            23355, 22824, 22293, 21762, 21231, 20700, 20169, 19638, 19107,
                            27000, 26388, 25776, 25164, 24552, 23940, 23328, 22716, 22104,
                            30645, 29952, 29259, 28566, 27873, 27180, 26487, 25794, 25101};
    auto unique_expected_data = std::make_unique<int[]>(81);
    std::memcpy(unique_expected_data.get(), expected_data, sizeof(int) * 81);
    I_Matrix<int> expected(9, 9, std::move(unique_expected_data));

    auto actual = lhs * rhs;

    ASSERT_EQ(actual, expected);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Time Taken: " << duration.count() << " microseconds" <<std::endl;
}

TEST(DOT_PRODUCT, 14x14) {
    auto start = std::chrono::high_resolution_clock::now();

    int lhs_data[196] = {5, 0, 3, 3, 7, 9, 3, 5, 2, 4, 7, 6, 8, 8,
                        1, 6, 7, 7, 8, 1, 5, 9, 8, 9, 4, 3, 0, 3,
                        5, 0, 2, 3, 8, 1, 3, 3, 3, 7, 0, 1, 9, 9,
                        0, 4, 7, 3, 2, 7, 2, 0, 0, 4, 5, 5, 6, 8,
                        4, 1, 4, 9, 8, 1, 1, 7, 9, 9, 3, 6, 7, 2,
                        0, 3, 5, 9, 4, 4, 6, 4, 4, 3, 4, 4, 8, 4,
                        3, 7, 5, 5, 0, 1, 5, 9, 3, 0, 5, 0, 1, 2,
                        4, 2, 0, 3, 2, 0, 7, 5, 9, 0, 2, 7, 2, 9,
                        2, 3, 3, 2, 3, 4, 1, 2, 9, 1, 4, 6, 8, 2,
                        3, 0, 0, 6, 0, 6, 3, 3, 8, 8, 8, 2, 3, 2,
                        0, 8, 8, 3, 8, 2, 8, 4, 3, 0, 4, 3, 6, 9,
                        8, 0, 8, 5, 9, 0, 9, 6, 5, 3, 1, 8, 0, 4,
                        9, 6, 5, 7, 8, 8, 9, 2, 8, 6, 6, 9, 1, 6,
                        8, 8, 3, 2, 3, 6, 3, 6, 5, 7, 0, 8, 4, 6};
    auto unique_lhs_data = std::make_unique<int[]>(196);
    std::memcpy(unique_lhs_data.get(), lhs_data, sizeof(int) * 196);
    I_Matrix<int> lhs(14, 14, std::move(unique_lhs_data));

    int rhs_data[196] = {5, 8, 2, 3, 9, 7, 5, 3, 4, 5, 3, 3, 7, 9,
                        9, 9, 7, 3, 2, 3, 9, 7, 7, 5, 1, 2, 2, 8,
                        1, 5, 8, 4, 0, 2, 5, 5, 0, 8, 1, 1, 0, 3,
                        8, 8, 4, 4, 0, 9, 3, 7, 3, 2, 1, 1, 2, 1,
                        4, 2, 5, 5, 5, 2, 5, 7, 7, 6, 1, 6, 7, 2,
                        3, 1, 9, 5, 9, 9, 2, 0, 9, 1, 9, 0, 6, 0,
                        4, 8, 4, 3, 3, 8, 8, 7, 0, 3, 8, 7, 7, 1,
                        8, 4, 7, 0, 4, 9, 0, 6, 4, 2, 4, 6, 3, 3,
                        7, 8, 5, 0, 8, 5, 4, 7, 4, 1, 3, 3, 9, 2,
                        5, 2, 3, 5, 7, 2, 7, 1, 6, 5, 0, 0, 3, 1,
                        9, 9, 6, 6, 7, 8, 8, 7, 0, 8, 6, 8, 9, 8,
                        3, 6, 1, 7, 4, 9, 2, 0, 8, 2, 7, 8, 4, 4,
                        1, 7, 6, 9, 4, 1, 5, 9, 7, 1, 3, 5, 7, 3,
                        6, 6, 7, 9, 1, 9, 6, 0, 3, 8, 4, 1, 4, 5};
    auto unique_rhs_data = std::make_unique<int[]>(196);
    std::memcpy(unique_rhs_data.get(), rhs_data, sizeof(int) * 196);
    I_Matrix<int> rhs(14, 14, std::move(unique_rhs_data));

    int expected_data[196] = {330, 373, 383, 376, 347, 440, 318, 290, 339, 287, 299, 272, 385, 241,
                            413, 400, 375, 254, 291, 392, 348, 367, 287, 302, 197, 249, 316, 217,
                            244, 288, 274, 293, 237, 270, 272, 250, 258, 232, 152, 178, 284, 171,
                            238, 295, 316, 314, 202, 307, 276, 205, 244, 247, 210, 162, 235, 184,
                            372, 400, 337, 301, 328, 382, 308, 363, 336, 254, 195, 260, 354, 209,
                            311, 378, 338, 298, 229, 366, 296, 334, 262, 226, 225, 230, 299, 177,
                            297, 317, 272, 147, 166, 301, 229, 272, 146, 194, 161, 178, 200, 194,
                            296, 354, 246, 221, 222, 373, 243, 235, 206, 192, 225, 228, 293, 194,
                            242, 314, 266, 234, 251, 274, 228, 260, 252, 164, 194, 203, 290, 175,
                            306, 311, 263, 219, 300, 339, 252, 237, 219, 183, 208, 175, 304, 160,
                            332, 408, 391, 322, 211, 362, 363, 356, 260, 311, 233, 259, 311, 240,
                            315, 385, 299, 261, 272, 408, 301, 306, 248, 283, 235, 281, 321, 219,
                            472, 538, 441, 392, 434, 568, 454, 387, 400, 369, 355, 329, 473, 321,
                            355, 393, 347, 303, 333, 411, 325, 265, 365, 261, 249, 225, 322, 265};
    auto unique_expected_data = std::make_unique<int[]>(196);
    std::memcpy(unique_expected_data.get(), expected_data, sizeof(int) * 196);
    I_Matrix<int> expected(14, 14, std::move(unique_expected_data));

    auto actual = lhs * rhs;
    ASSERT_EQ(actual, expected);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "14 x 14 Test: " << duration.count() << " microseconds" <<std::endl;
}

TEST(DETERMINANT, 3x3) {
    int data[9] = {7, -4, 2, 3, 1, -5, 2, 2, -5};
    auto unique_data = std::make_unique<int[]>(9);
    std::memcpy(unique_data.get(), data, sizeof(int) * 9);
    I_Matrix<int> matrix(3, 3, std::move(unique_data));

    double expected = 23.0;
    double actual = det(matrix);

    ASSERT_EQ(expected, actual);
}

TEST(DETERMINANT, 4x4) {
    int data[16] = {2, 1, 3, 4, 0, -1, 2, 1, 3, 2, 0, 5, -1, 3, 2, 1};
    auto unique_data = std::make_unique<int[]>(16);
    std::memcpy(unique_data.get(), data, sizeof(int) * 16);
    I_Matrix<int> matrix(4, 4, std::move(unique_data));

    double expected = 35.0;
    double actual = det(matrix);

    ASSERT_EQ(actual, expected);
}

TEST(DETERMINANT, ANOTHER_4x4) {
    auto start = std::chrono::high_resolution_clock::now();
    int data[16] = {0, -5, -2, -2, 2, 4, -2, 0, -3, -1, 2, 1, 3, 3, 5, -4};
    auto unique_data = std::make_unique<int[]>(16);
    std::memcpy(unique_data.get(), data, sizeof(int) * 16);
    I_Matrix<int> matrix(4, 4, std::move(unique_data));

    double expected = 176.0;
    double actual = det(matrix);

    ASSERT_EQ(actual, expected);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "4x4 Determinant Test: " << duration.count() << " microseconds" <<std::endl;
}

TEST(DETERMINANT, 5x5) {
    auto start = std::chrono::high_resolution_clock::now();
    int data[25] = {0, 3, 4, 0, -5,
                    -5, -4, 2, 1, 4,
                    -3, -1, 0, -3, -1,
                    5, -3, -1, 2, 2,
                    4, -4, 2, -5, 1};
    auto unique_data = std::make_unique<int[]>(25);
    std::memcpy(unique_data.get(), data, sizeof(int) * 25);
    I_Matrix<int> matrix(5, 5, std::move(unique_data));

    double expected = 3186.0;
    double actual = det(matrix);

    ASSERT_EQ(actual, expected);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "5x5 Determinant Test: " << duration.count() << " microseconds" <<std::endl;
}

TEST(DETERMINANT, 10x10) {
    auto start = std::chrono::high_resolution_clock::now();
    int data[64] = {8, 8, 4, 6, 4, 5, 1, 4,
                    5, 8, 5, 4, 1, 1, 2, 3,
                    6, 5, 6, 6, 4, 6, 3, 4,
                    6, 6, 6, 3, 2, 5, 3, 7,
                    1, 7, 4, 7, 3, 1, 5, 2,
                    1, 9, 6, 3, 9, 7, 5, 3,
                    2, 3, 5, 4, 7, 7, 2, 2,
                    6, 5, 8, 4, 2, 7, 8, 6};
    auto unique_data = std::make_unique<int[]>(64);
    std::memcpy(unique_data.get(), data, sizeof(int) * 64);
    I_Matrix<int> matrix(8, 8, std::move(unique_data));

    double expected = 8968.0;
    double actual = det(matrix);

    ASSERT_EQ(actual, expected);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "10x10 Determinant Test: " << duration.count() << " microseconds" <<std::endl;
}

TEST(DETERMINANT, ERROR) {
    auto unique_data = std::make_unique<int[]>(1);
    I_Matrix<int> matrix(5, 2, std::move(unique_data));

    ASSERT_THROW(det(matrix), std::logic_error);
}

TEST(INVERSE, 3x3) {
    int data[9] = {2, 4, -6, 7, 3, 5, 1, -2, 4};
    auto unique_data = std::make_unique<int[]>(9);
    std::memcpy(unique_data.get(), data, sizeof(int) * 9);
    I_Matrix<int> matrix(3, 3, std::move(unique_data));

    double expected_data[9] = {11/27.0, -2/27.0, 19/27.0, -23/54.0, 
                        7/27.0, -26/27.0, -17/54.0, 4/27.0, -11/27.0};
    auto unique_expected_data = std::make_unique<double[]>(9);
    std::memcpy(unique_expected_data.get(), expected_data, sizeof(double) * 9);
    I_Matrix<double> expected(3, 3, std::move(unique_expected_data));

    auto actual = inv(matrix);
    ASSERT_EQ(expected, actual);
}

TEST(INVERSE, 4x4) {
    int data[16] = {1, 0, 0, 1, 0, 2, 1, 2, 2, 1, 0, 1, 2, 0, 1, 4};
    auto unique_data = std::make_unique<int[]>(16);
    std::memcpy(unique_data.get(), data, sizeof(int) * 16);
    I_Matrix<int> matrix(4, 4, std::move(unique_data));

    double expected_data[16] = {-2.0, -0.5, 1.0, 0.5, 
                                1.0, 0.5, 0, -0.5, 
                                -8.0, -1.0, 2.0, 2.0, 
                                3.0, 0.5, -1.0, -0.5};
    auto unique_expected_data = std::make_unique<double[]>(16);
    std::memcpy(unique_expected_data.get(), expected_data, sizeof(double) * 16);
    I_Matrix<double> expected(4, 4, std::move(unique_expected_data));

    auto actual = inv(matrix);
    ASSERT_EQ(expected, actual);
}

TEST(INVERSE, SINGULAR) {
    int data[9] = {1, 2, 2, 1, 2, 2, 3, 2, -1};
    auto unique_data = std::make_unique<int[]>(9);
    std::memcpy(unique_data.get(), data, sizeof(int) * 9);
    I_Matrix<int> singular_mat(3, 3, std::move(unique_data));

    ASSERT_THROW(inv(singular_mat), std::logic_error);
}

TEST(INVERSE, NEAR_SINGULAR) {
    double data[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0001, 8.0001, 9.0001};
    auto unique_data = std::make_unique<double[]>(9);
    std::memcpy(unique_data.get(), data, sizeof(double) * 9);
    I_Matrix<double> near_singular(3, 3, std::move(unique_data));

    ASSERT_THROW(inv(near_singular), std::logic_error);
}

TEST(TRANSPOSE, SQUARE) {
    int data[4] = {1, 2, 3, 4};
    auto unique_data = std::make_unique<int[]>(4);
    std::memcpy(unique_data.get(), data, sizeof(int) * 4);
    I_Matrix<int> square_mat(2, 2, std::move(unique_data));

    int expected_data[4] = {1, 3, 2, 4};
    auto unique_expected_data = std::make_unique<int[]>(4);
    std::memcpy(unique_expected_data.get(), expected_data, sizeof(int) * 4);
    I_Matrix<int> expected(2, 2, std::move(unique_expected_data));

    ASSERT_EQ(square_mat.transpose(), expected);
}

TEST(TRANSPOSE, NON_SQUARE) {
    int data[6] = {1, 2, 3, 4, 5, 6};
    auto unique_data = std::make_unique<int[]>(6);
    std::memcpy(unique_data.get(), data, sizeof(int) * 6);
    I_Matrix<int> non_square(3, 2, std::move(unique_data));

    int expected_data[6] = {1, 3, 5, 2, 4, 6};
    auto unique_expected_data = std::make_unique<int[]>(6);
    std::memcpy(unique_expected_data.get(), expected_data, sizeof(int) * 6);
    I_Matrix<int> expected(2, 3, std::move(unique_expected_data));

    ASSERT_EQ(non_square.transpose(), expected);
}

TEST(GET_COLUMN, 3x3) {
    uint32_t mat_data[9] = {1, 2, 10, 3, 4, 1, 5, 6, 0};
    auto unique_matrix = std::make_unique<uint32_t[]>(9);
    std::memcpy(unique_matrix.get(), mat_data, 9 * sizeof(uint32_t));
    I_Matrix<uint32_t> mat(3, 3, std::move(unique_matrix));

    uint32_t vec_data[3] = {1, 3, 5};
    auto unique_vec = std::make_unique<uint32_t[]>(3);
    std::memcpy(unique_vec.get(), vec_data, 3 * sizeof(uint32_t));
    I_Vector<uint32_t> vec(3, std::move(unique_vec));

    auto result = mat.get_column(0);

    ASSERT_EQ(vec, result);
}

TEST(GET_ROW, 3x3) {
    uint32_t mat_data[9] = {1, 2, 10, 3, 4, 1, 5, 6, 0};
    auto unique_matrix = std::make_unique<uint32_t[]>(9);
    std::memcpy(unique_matrix.get(), mat_data, 9 * sizeof(uint32_t));
    I_Matrix<uint32_t> mat(3, 3, std::move(unique_matrix));

    uint32_t vec_data[3] {5, 6, 0};
    auto unique_vec = std::make_unique<uint32_t[]>(3);
    std::memcpy(unique_vec.get(), vec_data, 3 * sizeof(uint32_t));
    I_Vector<uint32_t> vec(3, std::move(unique_vec));

    auto result = mat.get_row(2);

    ASSERT_EQ(vec, result);
}

TEST(SET_COLUMN, 4x4) {
    uint32_t mat_data[16] = {1, 7, 0, 5, 
                            2, 10, 0, 72, 
                            3, 42, 1, 95, 
                            5, 1, 0, 11};
    auto unique_matrix = std::make_unique<uint32_t[]>(16);
    std::memcpy(unique_matrix.get(), mat_data, 16 * sizeof(uint32_t));
    I_Matrix<uint32_t> mat(4, 4, std::move(unique_matrix));

    uint32_t vec_data[4] = {71, 17, 82, 28};
    auto unique_vec = std::make_unique<uint32_t[]>(4);
    std::memcpy(unique_vec.get(), vec_data, 4 * sizeof(uint32_t));
    I_Vector<uint32_t> vec(4, std::move(unique_vec));

    uint32_t expected_data[16] = {1, 7, 71, 5, 
                                2, 10, 17, 72, 
                                3, 42, 82, 95, 
                                5, 1, 28, 11};
    auto uniue_expected = std::make_unique<uint32_t[]>(16);
    std::memcpy(uniue_expected.get(), expected_data, 16 * sizeof(uint32_t));
    I_Matrix<uint32_t> expected(4, 4, std::move(uniue_expected));

    mat.set_col(2, vec);

    ASSERT_EQ(mat, expected);
}

TEST(SET_ROW, 4x4) {
    uint32_t mat_data[16] = {1, 7, 0, 5, 
                            2, 10, 0, 72, 
                            3, 42, 1, 95, 
                            5, 1, 0, 11};
    auto unique_matrix = std::make_unique<uint32_t[]>(16);
    std::memcpy(unique_matrix.get(), mat_data, 16 * sizeof(uint32_t));
    I_Matrix<uint32_t> mat(4, 4, std::move(unique_matrix));

    uint32_t vec_data[4] = {71, 17, 82, 28};
    auto unique_vec = std::make_unique<uint32_t[]>(4);
    std::memcpy(unique_vec.get(), vec_data, 4 * sizeof(uint32_t));
    I_Vector<uint32_t> vec(4, std::move(unique_vec));

    uint32_t expected_data[16] = {1, 7, 0, 5, 
                                2, 10, 0, 72, 
                                3, 42, 1, 95, 
                                71, 17, 82, 28};
    auto unique_expected = std::make_unique<uint32_t[]>(16);
    std::memcpy(unique_expected.get(), expected_data, 16 * sizeof(uint32_t));
    I_Matrix<uint32_t> expected(4, 4, std::move(unique_expected));

    mat.set_row(3, vec);

    ASSERT_EQ(mat, expected);
}

TEST(QR_TEST, 2x2) {
    double a_data[4] = {1, 2, 3, 4};
    auto unique_a_data = std::make_unique<double[]>(4);
    std::memcpy(unique_a_data.get(), a_data, 4 * sizeof(double));
    I_Matrix<double> A(2, 2, std::move(unique_a_data));

    std::tuple<I_Matrix<double>, I_Matrix<double>> actualQR = I_Matrix<double>::QR(A);

    auto identity = std::get<0>(actualQR) * std::get<1>(actualQR);

    ASSERT_TRUE(identity == A);
}

TEST(QR_TEST, 3x3) {
    double a_data[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    auto unique_a_data = std::make_unique<double[]>(9);
    std::memcpy(unique_a_data.get(), a_data, 9 * sizeof(double));
    I_Matrix<double> A(3, 3, std::move(unique_a_data));

    auto actualQR = I_Matrix<double>::QR(A);
    auto actual_A = std::get<0>(actualQR) * std::get<1>(actualQR);

    ASSERT_TRUE(actual_A == A);
}


TEST(QR_TEST, 4x4) {
    double a_data[16] = {1, 2, 3, 4, 
                        5, 6, 7, 8, 
                        9, 10, 11, 12, 
                        13, 14, 15, 16};
    auto unique_a_data = std::make_unique<double[]>(16);
    std::memcpy(unique_a_data.get(), a_data, 16 * sizeof(double));
    I_Matrix<double> A(4, 4, std::move(unique_a_data));

    auto actualQR = I_Matrix<double>::QR(A);
    auto actual_A = std::get<0>(actualQR) * std::get<1>(actualQR);

    ASSERT_TRUE(actual_A == A);
}

TEST(QR_TEST, 5x5) {
    double a_data[25] = {3, 7, 2, 6, 9,
                        1, 5, 8, 4, 0,
                        2, 3, 6, 7, 1,
                        9, 4, 3, 5, 8,
                        0, 1, 7, 2, 6};
    auto unique_a_data = std::make_unique<double[]>(25);
    std::memcpy(unique_a_data.get(), a_data, 25 * sizeof(double));
    I_Matrix<double> A(5, 5, std::move(unique_a_data));

    auto actualQR = I_Matrix<double>::QR(A);
    auto actual_A = std::get<0>(actualQR) * std::get<1>(actualQR);

    ASSERT_TRUE(actual_A == A);
}

TEST(QR_TEST, 10x10) {
    constexpr double a_data[100] = {3.14, 2.71, 1.62, 4.66, 5.55, 9.81, 7.77, 8.88, 6.66, 0.01,
                                    1.23, 3.45, 6.78, 9.10, 2.34, 5.67, 8.90, 1.11, 2.22, 3.33,
                                    4.44, 5.55, 6.66, 7.77, 8.88, 9.99, 0.12, 0.23, 0.34, 0.45,
                                    0.56, 0.67, 0.78, 0.89, 1.90, 2.01, 3.12, 4.23, 5.34, 6.45,
                                    7.56, 8.67, 9.78, 0.89, 1.98, 2.87, 3.76, 4.65, 5.54, 6.43,
                                    7.32, 8.21, 9.10, 0.09, 1.18, 2.27, 3.36, 4.45, 5.54, 6.63,
                                    7.72, 8.81, 9.90, 0.01, 1.02, 2.03, 3.04, 4.05, 5.06, 6.07,
                                    7.08, 8.09, 9.10, 1.11, 2.22, 3.33, 4.44, 5.55, 6.66, 7.77,
                                    8.88, 9.99, 0.12, 1.23, 2.34, 3.45, 4.56, 5.67, 6.78, 7.89,
                                    8.90, 9.01, 0.12, 1.23, 2.34, 3.45, 4.56, 5.67, 6.78, 7.89};
    auto unique_a_data = std::make_unique<double[]>(100);
    std::memcpy(unique_a_data.get(), a_data, 100 * sizeof(double));
    I_Matrix<double> A(10, 10, std::move(unique_a_data));

    auto actualQR = I_Matrix<double>::QR(A);
    auto actual_A = std::get<0>(actualQR) * std::get<1>(actualQR);

    ASSERT_EQ(actual_A, A);
}

TEST(QR_TEST, 14x14) {
    constexpr double a_data[196] = {1.11, 2.22, 3.33, 4.44, 5.55, 6.66, 7.77, 8.88, 9.99, 0.11, 1.23, 2.34, 3.45, 4.56,
                                    5.67, 6.78, 7.89, 8.90, 9.01, 0.12, 1.34, 2.46, 3.58, 4.60, 5.72, 6.84, 7.96, 8.08,
                                    9.20, 0.32, 1.44, 2.56, 3.68, 4.80, 5.92, 6.04, 7.16, 8.28, 9.40, 0.52, 1.64, 2.76,
                                    3.88, 5.00, 6.12, 7.24, 8.36, 9.48, 0.60, 1.72, 2.84, 3.96, 5.08, 6.20, 7.32, 8.44,
                                    9.56, 0.68, 1.80, 2.92, 4.04, 5.16, 6.28, 7.40, 8.52, 9.64, 0.76, 1.88, 3.00, 4.12,
                                    5.24, 6.36, 7.48, 8.60, 9.72, 0.84, 1.96, 3.08, 4.20, 5.32, 6.44, 7.56, 8.68, 9.80,
                                    0.92, 2.04, 3.16, 4.28, 5.40, 6.52, 7.64, 8.76, 9.88, 1.00, 2.12, 3.24, 4.36, 5.48,
                                    6.60, 7.72, 8.84, 9.96, 1.08, 2.20, 3.32, 4.44, 5.56, 6.68, 7.80, 8.92, 0.04, 1.16,
                                    2.28, 3.40, 4.52, 5.64, 6.76, 7.88, 9.00, 0.12, 1.24, 2.36, 3.48, 4.60, 5.72, 6.84,
                                    7.96, 9.08, 0.20, 1.32, 2.44, 3.56, 4.68, 5.80, 6.92, 8.04, 9.16, 0.28, 1.40, 2.52,
                                    3.64, 4.76, 5.88, 7.00, 8.12, 9.24, 0.36, 1.48, 2.60, 3.72, 4.84, 5.96, 7.08, 8.20,
                                    9.32, 0.44, 1.56, 2.68, 3.80, 4.92, 6.04, 7.16, 8.28, 9.40, 0.52, 1.64, 2.76, 3.88,
                                    5.00, 6.12, 7.24, 8.36, 9.48, 0.60, 1.72, 2.84, 3.96, 5.08, 6.20, 7.32, 8.44, 9.56};
    auto unique_a_data = std::make_unique<double[]>(196);
    std::memcpy(unique_a_data.get(), a_data, 196 * sizeof(double));
    I_Matrix<double> A(14, 14, std::move(unique_a_data));

    auto QR = I_Matrix<double>::QR(A);
    auto actual_A = std::get<0>(QR) * std::get<1>(QR);

    ASSERT_EQ(actual_A, A);
}

TEST(EIG, 3x3) {
    double a_data[9] = {41, 12, 87, 62, 63, 98, 72, 21, 61};
    auto unique_a_data = std::make_unique<double[]>(9);
    std::memcpy(unique_a_data.get(), a_data, 9 * sizeof(double));
    I_Matrix<double> A(3, 3, std::move(unique_a_data));

    double expected_data[3] = {159.21038910441368, 35.59260509899343, -29.802994203407092};
    auto unique_expected_data = std::make_unique<double[]>(3);
    std::memcpy(unique_expected_data.get(), expected_data, 3 * sizeof(double));
    I_Vector<double> expected(3, std::move(unique_expected_data));

    I_Vector<double> actual = I_Matrix<double>::eig(A);
    ASSERT_TRUE(fast_math::abs(actual.get_element(2) - expected.get_element(2)) < 0.000001);
}

TEST(EIG, ANOTHER_3x3) {
    double a_data[9] = {0, 1, 1, 1, 0, 1, 1, 1, 0};
    auto unique_a_data = std::make_unique<double[]>(9);
    std::memcpy(unique_a_data.get(), a_data, 9 * sizeof(double));
    I_Matrix<double> A(3, 3, std::move(unique_a_data));

    double expected_last = -1;

    I_Vector<double> actual = I_Matrix<double>::eig(A);
    ASSERT_TRUE(fast_math::abs(actual.get_element(2) - expected_last) < 0.000001);
}

TEST(EIG, YET_ANOTHER_3x3) {
    double a_data[9] = {-1, 18, 0, 1, 2, 0, 5, -3, -1};
    auto unique_a_data = std::make_unique<double[]>(9);
    std::memcpy(unique_a_data.get(), a_data, 9 * sizeof(double));
    I_Matrix<double> A(3, 3, std::move(unique_a_data));

    double expected_last = -1;

    I_Vector<double> actual = I_Matrix<double>::eig(A);
    ASSERT_TRUE(fast_math::abs(actual.get_element(2) - expected_last) < 0.000001);
}

TEST(EIG, EXCEPTION) {
    I_Matrix<double> exception(4, 3);

    ASSERT_THROW(I_Matrix<double>::eig(exception), std::invalid_argument);
}