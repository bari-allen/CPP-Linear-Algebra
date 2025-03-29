#ifndef IMAGEMAT_H
#define IMAGEMAT_H

#include <cstring>
#include <stdexcept>
#include <omp.h>
#include <math.h>
#include <vector>

class imageMat {
    public:
        //Below are some various constructors that may or may not be useful
        imageMat();
        imageMat(int nRows, int nCols);
        //Use an array of unsigned chars because each pixel in BMP is represented
        //and 3 BGR hex values
        imageMat(int nRows, int nCols, const double input_data);
        imageMat(const imageMat& inputMat);
        imageMat(int nRows, int nCols, const std::vector<double> *input_data);

        //The destructor
        ~imageMat();

        //Configuration functions
        bool resize(int nRows, int nCols);
        void set_identity();

        //Element access functions
        unsigned int get(int row, int col) const;
        bool set(int row, int col, const double rgb_data);
        int get_rows() const;
        int get_cols() const;

        //Invert the current matrix if it is invertible
        bool inverse();

        //Overload the equals operator
        bool operator== (const imageMat& rhs) const;
        bool compare_to(const imageMat& comparator, int tolerance) const;

        //Overload the +, -, and * operators
        //First takes two matrices as inputs
        //Second takes a scalar and a matrix as inputs
        //Third takes a matrix and a scalar as inputs
        friend imageMat operator+ (const imageMat& lhs, const imageMat& rhs);
        friend imageMat operator+ (const unsigned char& lhs, const imageMat& rhs);
        friend imageMat operator+ (const imageMat& lhs, const unsigned char& rhs);

        friend imageMat operator- (const imageMat& lhs, const imageMat& rhs);
        friend imageMat operator- (const unsigned char& lhs, const imageMat& rhs);
        friend imageMat operator- (const imageMat& lhs, const unsigned char& rhs);

        friend imageMat operator* (const imageMat& lhs, const imageMat& rhs);
        friend imageMat operator* (const unsigned char& lhs, const imageMat& rhs);
        friend imageMat operator* (const imageMat& lhs, const unsigned char& rhs);

        void separate(imageMat* const matrix_1, imageMat* const matrix_2, int split_col);

    public: //TODO: Make these functions private after testing
        //Get the range of index values for a given row and column in the
        //flattened array
        int get_linear_index(int row, int col);

        //TODO: Finish implementation of the following
        bool is_square() const;
        void swap_row(int row1, int row2);
        void mult_add(int addend, int multiplicant, int multiplication_factor);
        void mult_row(int row, int multiplication_factor);
        bool join(const imageMat& matrix2);
        int row_with_max(int col_num, int starting_row) const;

    private:
        double m_data;
        int n_rows, n_cols, n_elements;
        

};



#endif