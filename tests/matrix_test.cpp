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
