#include <gtest/gtest.h>
#include "../src/imageMat.h"

TEST(Matrix_Multiplication, Simple) {
    int a_rows = 1;
    int a_cols = 3;
    unsigned int a[3] = {1, 2, 5};
    imageMat a_mat(a_rows, a_cols, a);

    int b_rows = 3;
    int b_cols = 1;
    unsigned int b[3] = {10, 20, 1};
    imageMat b_mat(b_rows, b_cols, b);

    int e_rows = 1;
    int e_cols = 1;
    unsigned int e[1] = {55};
    imageMat e_mat(e_rows, e_cols, e);

    imageMat actual_mat = a_mat * b_mat;

    EXPECT_TRUE(e_mat == actual_mat);
}

TEST(Matrix_Multiplication, Intermediate) {
    int a_rows = 2;
    int a_cols = 3;
    unsigned int a[6] = {1, 2, 5, 10, 0, 20};
    imageMat a_mat(a_rows, a_cols, a);

    int b_rows = 3;
    int b_cols = 2;
    unsigned int b[6] = {10, 20, 20, 0, 1, 2};
    imageMat b_mat(b_rows, b_cols, b);

    int e_rows = 2;
    int e_cols = 2;
    unsigned int e[4] = {55, 30, 120, 240};
    imageMat e_mat(e_rows, e_cols, e);

    imageMat actual_mat = a_mat * b_mat;

    EXPECT_TRUE(e_mat == actual_mat);
}

TEST(Matrix_Multiplication, Advanced) {
    int a_rows = 6;
    int a_cols = 3;
    unsigned int a[18] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
    imageMat a_mat(a_rows, a_cols, a);

    int b_rows = 3;
    int b_cols = 6;
    unsigned int b[18] = {1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0};
    imageMat b_mat(b_rows, b_cols, b);

    int e_rows = 6;
    int e_cols = 6;
    unsigned int e[36] = {1, 2, 3, 3, 2, 1, 4, 5, 6, 6, 5, 4, 7, 8, 9, 9, 8, 7, 10, 11, 12, 12, 11, 10, 
        13, 14, 15, 15, 14, 13, 16, 17, 18, 18, 17, 16};
    imageMat e_mat(e_rows, e_cols, e);

    imageMat actual_mat = a_mat * b_mat;

    EXPECT_TRUE(e_mat == actual_mat);
}

TEST(Separate_Matrix, Error) {
    int a_rows = 2;
    int a_cols = 3;
    unsigned int a[6] = {0, 0, 0, 0, 0, 0};
    imageMat a_mat(a_rows, a_cols, a);

    int out_of_bounds_col = 4;
    imageMat* test_1 = new imageMat();
    imageMat* test_2 = new imageMat();

    EXPECT_THROW(a_mat.separate(test_1, test_2, out_of_bounds_col), std::invalid_argument);

    out_of_bounds_col = a_cols;
    EXPECT_THROW(a_mat.separate(test_1, test_2, out_of_bounds_col), std::invalid_argument);

    out_of_bounds_col = 0;
    EXPECT_THROW(a_mat.separate(test_1, test_2, out_of_bounds_col), std::invalid_argument);

    out_of_bounds_col = -10;
    EXPECT_THROW(a_mat.separate(test_1, test_2, out_of_bounds_col), std::invalid_argument);
}

TEST(Separate_Matrix, Simple) {
    unsigned int a[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    imageMat a_mat(3, 4, a);

    imageMat* b_mat = new imageMat();
    imageMat* c_mat = new imageMat();

    unsigned int e_1[6] = {0, 1, 4, 5, 8, 9};
    imageMat e_1_mat(3, 2, e_1);

    unsigned int e_2[6] = {2, 3, 6, 7, 10, 11};
    imageMat e_2_mat(3, 2, e_2);

    a_mat.separate(b_mat, c_mat, 2);

    EXPECT_TRUE(*b_mat == e_1_mat);
    EXPECT_TRUE(*c_mat == e_2_mat);
}

TEST(Separate_Matrix, Advanced) {
    unsigned int a[25] = {1, 1, 1, 1, 1, 0, 2, 4, 10, 3, 20, 10, 15, 10,
        12, 0, 1, 1, 0, 1, 100, 10, 2, 3, 7};
    imageMat a_mat(5, 5, a);

    imageMat* b_mat = new imageMat();
    imageMat* c_mat = new imageMat();

    unsigned int e_1[15] = {1, 1, 1, 0, 2, 4, 20, 10, 15, 0, 1, 1, 100, 10, 2};
    imageMat e_1_mat(5, 3, e_1);

    unsigned int e_2[10] = {1, 1, 10, 3, 10, 12, 0, 1, 3, 7};
    imageMat e_2_mat(5, 2, e_2);

    a_mat.separate(b_mat, c_mat, 3);

    EXPECT_TRUE(*b_mat == e_1_mat);
    EXPECT_TRUE(*c_mat == e_2_mat);
}

TEST(Swap_Rows, Error) {
    unsigned int a[6] = {0, 0, 0, 0, 0, 0};
    imageMat a_mat(3, 2, a);
    int out_of_bounds_row = -1;

    EXPECT_THROW(a_mat.swap_row(out_of_bounds_row, 2), std::invalid_argument);

    out_of_bounds_row = 4;
    EXPECT_THROW(a_mat.swap_row(1, out_of_bounds_row), std::invalid_argument);
}

TEST(Swap_Rows, Simple) {
    unsigned int a[6] = {10, 2, 11, 4, 3, 20};
    imageMat a_mat(3, 2, a);

    unsigned int e[6] = {3, 20, 11, 4, 10, 2};
    imageMat e_mat(3, 2, e);

    a_mat.swap_row(0, 2);
    EXPECT_TRUE(a_mat == e_mat);
}

TEST(Swap_Rows, Intermediate) {
    unsigned int a[12] = {10, 2, 31, 22, 5, 0, 11, 100, 250, 0, 1, 0};
    imageMat a_mat(4, 3, a);

    unsigned int e[12] = {10, 2, 31, 0, 1, 0, 11, 100, 250, 22, 5, 0};
    imageMat e_mat(4, 3, e);

    a_mat.swap_row(1, 3);
    EXPECT_EQ(a_mat, e_mat);
}
