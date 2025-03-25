#include <gtest/gtest.h>
#include "../src/imageMat.h"

TEST(Matrix, Multiplication_Simple) {
    int a_rows = 1;
    int a_cols = 3;
    unsigned char a[9] = {1,2,1,2,2,1,1,0,1};
    imageMat a_mat(a_rows, a_cols, a);
    delete[] a;

    int b_rows = 3;
    int b_cols = 1;
    unsigned char b[9] = {2,2,2,1,1,0,0,0,1};
    imageMat b_mat(b_rows, b_cols, b);
    delete[] b;

    int e_rows = 1;
    int e_cols = 1;
    unsigned char e[3] = {4, 6, 3};
    imageMat e_mat(e_rows, e_cols, e);
    delete[] e;

    imageMat actual = a_mat * b_mat;

    EXPECT_TRUE(e_rows == actual.get_rows() && e_cols == actual.get_cols());
    EXPECT_TRUE(e_mat == actual);
};

TEST(Matrix, Multiplication_Intermediate) {
    int a_rows = 1;
    int a_cols = 3;
    unsigned char a[9] = {1, 2, 1, 2, 2, 1, 1, 0, 1};
    imageMat a_mat(a_rows, a_cols, a);

    int b_rows = 3;
    int b_cols = 2;
    unsigned char b[18] = {2, 2, 2, 10, 2, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 10, 5};
    imageMat b_mat(b_rows, b_cols, b);

    int e_rows = 1;
    int e_cols = 2;
    unsigned char e[6] = {4, 6, 3, 11, 4, 6};
    imageMat e_mat(e_rows, e_cols, e);

    imageMat actual = a_mat * b_mat;

    EXPECT_TRUE(e_rows == actual.get_rows() && e_cols == actual.get_cols());
    EXPECT_TRUE(e_mat == actual);
}