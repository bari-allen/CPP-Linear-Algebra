#include <gtest/gtest.h>
#include "../lib/imageMat.h"

TEST(Matrix_Addition, Simple) {
    auto rows = 3;
    auto cols = 2;
    const double a[6] = {3.0, 2.0, 1.5, 5.3, 6.0, 9.2};
    imageMat a_mat(rows, cols, a);

    const double b[6] = {1.5, 2.0, 1.0, 6.3, 7.0, 1.1};
    imageMat b_mat(rows, cols, b);

    const double e[6] = {4.5, 4.0, 2.5, 11.6, 13.0, 10.3};
    imageMat e_mat(rows, cols, e);

    auto result = a_mat + b_mat;
    EXPECT_EQ(e_mat, result);
}