#include <gtest/gtest.h>
#include "../lib/imageMat.h"
using std::size_t;

TEST(Matrix_Addition, Simple) {
    auto rows = 3;
    auto cols = 2;
    const double a[6] = {3.0, 2.0, 1.5, 5.3, 6.0, 9.2};
    ImageMat a_mat(rows, cols, a);

    const double b[6] = {1.5, 2.0, 1.0, 6.3, 7.0, 1.1};
    ImageMat b_mat(rows, cols, b);

    const double e[6] = {4.5, 4.0, 2.5, 11.6, 13.0, 10.3};
    ImageMat e_mat(rows, cols, e);

    auto result = a_mat + b_mat;
    EXPECT_EQ(e_mat, result);
}

TEST(Matrix_Mult, Intermediate) {
    constexpr size_t a_rows = 2;
    constexpr size_t a_cols = 3;
    constexpr double a[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    ImageMat a_mat(a_rows, a_cols, a);

    constexpr size_t b_rows = 3;
    constexpr size_t b_cols = 2;
    constexpr double b[6] = {7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
    ImageMat b_mat(b_rows, b_cols, b);

    constexpr double e[4] = {58.0, 64.0, 139.0, 154.0};
    ImageMat e_mat(2, 2, e);

    auto result = a_mat * b_mat;
    EXPECT_EQ(result, e_mat);
}

TEST(Matrix_Swap, Same_Row) {
    constexpr double a[4] = {1.0, 2.0, 3.0, 4.0};
    ImageMat a_mat(2, 2, a);

    ImageMat e_mat(a_mat);

    a_mat.swap_row(0, 0);
    EXPECT_EQ(a_mat, e_mat);
}

TEST(Matrix_Swap, Simple) {
    constexpr size_t a_rows = 2;
    constexpr size_t a_cols = 2;
    constexpr double a[4] = {1.0, 2.0, 3.0, 4.0};
    ImageMat a_mat(a_rows, a_cols, a);

    constexpr size_t e_rows = 2;
    constexpr size_t e_cols = 2;
    constexpr double e[4] = {3.0, 4.0, 1.0, 2.0};
    ImageMat e_mat(e_rows, e_cols, e);

    a_mat.swap_row(0, 1);
    EXPECT_EQ(a_mat, e_mat);
}

TEST(Matrix_Swap, Intermediate) {
    constexpr double a[8] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    ImageMat a_mat(4, 2, a);

    constexpr double e[8] = {1.0, 2.0, 7.0, 8.0, 5.0, 6.0, 3.0, 4.0};
    ImageMat e_mat(4, 2, e);

    a_mat.swap_row(3, 1);
    EXPECT_EQ(a_mat, e_mat);
}