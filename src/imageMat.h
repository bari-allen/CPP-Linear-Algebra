#ifndef IMAGEMAT_H
#define IMAGEMAT_H

#include <cstring>
#include <stdexcept>
#include <omp.h>
#include <math.h>
#include <vector>
#include <memory>

/**
 * A matrix class that will include functions to perform simple linear algebra operations
 * The purpose of this class is to allow for PCA transformations of images
 * @author Karl Haidinyak, bari-allen (GitHub)
 */
class imageMat {
    public:
        //Below are some various constructors that may or may not be useful
        imageMat();
        imageMat(int nRows, int nCols);
        //Use an array of unsigned chars because each pixel in BMP is represented
        //and 3 BGR hex values
        imageMat(int nRows, int nCols, const double* input_data);
        imageMat(const imageMat& inputMat);
        imageMat(int nRows, int nCols, const std::vector<double> *input_data);

        //The destructor
        ~imageMat();

        //Configuration functions
        void resize(int nRows, int nCols);
        void set_identity();

        //Element access functions
        unsigned int get(int row, int col) const;
        bool set(int row, int col, const double rgb_data);
        int get_rows() const;
        int get_cols() const;
        std::unique_ptr<double> get_data() const;

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
        int get_linear_index(int row, int col);
        bool is_square() const;
        void swap_row(int row1, int row2);
        void mult_add(int addend, int multiplicant, int multiplication_factor);
        void mult_row(int row, int multiplication_factor);
        bool join(const imageMat& matrix2);
        int row_with_max(int col_num, int starting_row) const;

    private:
        double* m_data;
        int n_rows, n_cols, n_elements;
        

};

/**
 * Default constructor with no paramters. Should never typically be used
 */
imageMat::imageMat() {
    n_rows = 0;
    n_cols = 0;
    n_elements = 0;

    m_data = new double[n_elements]();
}

/**
 * Constructor that makes an empty C-Style array with the inputted shape
 * @param rows the number of rows
 * @param cols the number of columns
 */
imageMat::imageMat(int rows, int cols) {
    n_rows = rows;
    n_cols = cols;
    n_elements = n_rows * n_cols;

    m_data = new double[n_elements]();
}

/**
 * Constructor that makes a C-Style array with the inputted shape and populates
 * it with the given data in the input_data parameter
 * @param rows the number of rows
 * @param cols the number of columns
 * @param input_data a pointer to an array of doubles that will be copied into the matrix
 */
imageMat::imageMat(int rows, int cols, const double* input_data) {
    n_rows = rows;
    n_cols = cols;
    n_elements = n_rows * n_cols;

    std::memcpy(m_data, input_data, n_elements * sizeof(double));
}

/**
 * Constructor that makes a copy of another instance of an imageMat
 * @param input_mat the matrix that will be copied
 */
imageMat::imageMat(const imageMat& input_mat) {
    n_rows = input_mat.get_rows();
    n_cols = input_mat.get_cols();
    n_elements = n_rows * n_cols;

    std::memcpy(m_data, input_mat.get_data().get(), n_elements * sizeof(double));
}

/**
 * Constructor that constructs a C-Style array with the given shape and input data
 * @param rows the number of rows
 * @param cols the number of columns
 * @param input_data a pointer to a vector of doubles that will be copied into the matrix
 */
imageMat::imageMat(int rows, int cols, const std::vector<double>* input_data) {
    n_rows = rows;
    n_cols = cols;
    n_elements = n_rows * n_cols;

    if (n_elements != input_data->size()) {
        throw std::invalid_argument("Number of elements must match the size of the vector!");
    }

    std::copy(input_data->begin(), input_data->end(), m_data);
}

imageMat::~imageMat() {
    if (m_data != nullptr) {
        delete[] m_data;
    }
}

/**
 * Resizes the matrix to the inputted shape
 * @param rows the new number of rows for the matrix
 * @param cols the new number of columns for the matrix
 * @throw invalid_argument when the inputted rows or columns is less than 0
 */
void imageMat::resize(int rows, int cols) {
    if (rows <= 0 || cols <= 0) {
        throw std::invalid_argument("The number of rows or columns cannot be less than 1!");
    }

    n_rows = rows;
    n_cols = cols;
    n_elements = n_rows * n_cols;

    m_data = new double[n_elements]();
}

/**
 * Converts the matrix into the identity matrix if the matrix is square
 * @throw std::logic_error when the matrix is not square
 */
void imageMat::set_identity() {
    if (!is_square()) {
        throw std::logic_error("Cannot set a non-square matrix to the identity matrix!");
    }

    std::memset(m_data, 0, n_elements * sizeof(double));
    for (int i = 0; i < n_rows; ++i) {
        int diag_index = get_linear_index(i, i);
        m_data[diag_index] = 1;
    }
}


#endif