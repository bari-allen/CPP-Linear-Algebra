#ifndef IMAGEMAT_H
#define IMAGEMAT_H

#include <cstring>
#include <stdexcept>
#include <omp.h>
#include <math.h>
#include <cmath>
#include <vector>
#include <memory>

/**
 * A matrix class that will include functions to perform simple linear algebra operations
 * The purpose of this class is to allow for PCA transformations of images
 * @author Karl Haidinyak, bari-allen (GitHub)
 */
class ImageMat {
    public:
        //Below are some various constructors that may or may not be useful
        ImageMat();
        ImageMat(int nRows, int nCols);
        //Use an array of unsigned chars because each pixel in BMP is represented
        //and 3 BGR hex values
        ImageMat(int nRows, int nCols, const double* input_data);
        ImageMat(const ImageMat& inputMat);
        ImageMat(int nRows, int nCols, const std::vector<double> *input_data);

        //The destructor
        ~ImageMat();

        //Configuration functions
        void resize(int nRows, int nCols);
        void set_identity();
        void clamp(const double low, const double high);

        //Element access functions
        double get(int row, int col) const;
        bool set(int row, int col, double rgb_data);
        int get_rows() const;
        int get_cols() const;
        std::unique_ptr<double[]> get_data() const;
        std::unique_ptr<unsigned char[]> get_pixel_data() const noexcept(true);


        //Invert the current matrix if it is invertible
        bool inverse();

        //Overload the equals operator
        bool operator== (const ImageMat& rhs) const;
        bool within_tolerance(const ImageMat& comparator, double tolerance) const noexcept;
        bool close_enough(double val1, double val2) const;

        //Overload the +, -, and * operators
        //First takes two matrices as inputs
        //Second takes a scalar and a matrix as inputs
        //Third takes a matrix and a scalar as inputs
        friend ImageMat operator+ (const ImageMat& lhs, const ImageMat& rhs);
        friend ImageMat operator+ (const ImageMat& lhs, const double& rhs);

        friend ImageMat operator- (const ImageMat& lhs, const ImageMat& rhs);
        friend ImageMat operator- (const ImageMat& lhs, const unsigned char& rhs);

        friend ImageMat operator* (const ImageMat& lhs, const ImageMat& rhs);
        friend ImageMat operator* (const ImageMat& lhs, const unsigned char& rhs);

        void separate(ImageMat* const matrix_1, ImageMat* const matrix_2, int split_col);

    public: //TODO: Make these functions private after testing
        int get_linear_index(int row, int col) const noexcept;
        double get_linear(int linear_index) const;
        bool is_square() const noexcept;
        void swap_row(int row1, int row2);
        void mult_add(int addend, int multiplicant, int multiplication_factor);
        void mult_row(int row, int multiplication_factor);
        bool join(const ImageMat& matrix2);
        int row_with_max(int col_num, int starting_row) const;

    private:
        double* m_data;
        int n_rows, n_cols, n_elements;
        

};

#endif