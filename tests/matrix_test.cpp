#include "../lib/imageMat.h"
#include <gtest/gtest.h>
#include <memory>
#include <chrono>

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

TEST(DOT_PRODUCT, VERY_HARD) {
    int lhs_data[12] = {-1, 1, -1, 5, 2, -5, 6, -5, 1, -5, 6, 0};
    auto unique_lhs_data = std::make_unique<int[]>(12);
    std::memcpy(unique_lhs_data.get(), lhs_data, sizeof(int) * 12);
    I_Matrix<int> lhs(4, 3, unique_lhs_data);

    int rhs_data[6] = {6, 5, 5, -6, 6, 0};
    auto unique_rhs_data = std::make_unique<int[]>(6);
    std::memcpy(unique_rhs_data.get(), rhs_data, sizeof(int) * 6);
    I_Matrix<int> rhs(3, 2, unique_rhs_data);

    int expected_data[8] = {-7, -11, 10, 13, 17, 60, 0, -61};
    auto unique_expected_data = std::make_unique<int[]>(8);
    std::memcpy(unique_expected_data.get(), expected_data, sizeof(int) * 8);
    I_Matrix<int> expected(4, 2, unique_expected_data);

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
    I_Matrix<int> lhs(9, 9, unique_lhs_data);

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
    I_Matrix<int> rhs(9, 9, unique_rhs_data);

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
    I_Matrix<int> expected(9, 9, unique_expected_data);

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
    I_Matrix<int> lhs(14, 14, unique_lhs_data);

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
    I_Matrix<int> rhs(14, 14, unique_rhs_data);

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
    I_Matrix<int> expected(14, 14, unique_expected_data);

    auto actual = lhs * rhs;
    ASSERT_EQ(actual, expected);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Time Taken: " << duration.count() << " microseconds" <<std::endl;
}
