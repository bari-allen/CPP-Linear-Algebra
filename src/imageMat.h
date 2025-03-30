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
        void clamp(const double low, const double high);

        //Element access functions
        double get(int row, int col) const;
        bool set(int row, int col, double rgb_data);
        int get_rows() const;
        int get_cols() const;
        std::unique_ptr<double[]> get_data() const;
        std::unique_ptr<unsigned char[]> get_pixel_data() const;


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
        int get_linear_index(int row, int col) const;
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
        m_data = nullptr;
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

/**
 * Clamps the matrix to a certain range
 * @param low the minimum value an item in the matrix can have
 * @param high the maximum value an item in the matrix can have
 * @throws invalid_argument when low is greater than high
 */
void imageMat::clamp(const double low, const double high) {
    if (low > high) {
        throw std::invalid_argument("The low value cannot be greater than the high value");
    }

    for (int i = 0; i < n_elements; ++i) {
        double data = m_data[i];
        m_data[i] = std::max(low, std::max(high, data));
    }
}

/**
 * Returns the item at the given row and column
 * @param row the row of the desired item
 * @param col the column of the desired item
 * @throws invalid_argument when the row and column aren't in the matrix
 * @returns the item at the given row and column
 */
double imageMat::get(int row, int col) const {
    int linear_index = get_linear_index(row, col);

    if (linear_index == -1) {
        throw std::invalid_argument("Invalid row or column inputted!");
    }

    return m_data[linear_index];
}

/**
 * Sets the item at the given row and column to the given RGB value
 * @param row the row to be set
 * @param col the column to be set
 * @param rgb_data the data to be set
 * @returns true when the item is properly set and false when the index doesn't exist
 */
bool imageMat::set(int row, int col, double rgb_data) {
    int linear_index = get_linear_index(row, col);

    if (linear_index == -1) {
        return false;
    }

    m_data[linear_index] = rgb_data;
    return true;
}

/**
 * Returns the number of rows in the matrix
 * @returns the number of rows in the matrix
 */
int imageMat::get_rows() const {
    return n_rows;
}

/**
 * Returns the number of columns in the matrix
 * @returns the number of columns in the matrix
 */
int imageMat::get_cols() const {
    return n_cols;
}

/**
 * Returns a unique pointer with a COPY of the matrix data
 * @returns a unique pointer with a COPY of the matrix data
 */
std::unique_ptr<double[]> imageMat::get_data() const {
    auto copy_data = std::make_unique<double[]>(n_elements);
    std::memcpy(copy_data.get(), m_data, n_elements * sizeof(double));
    
    return copy_data;
}

/**
 * Returns a unique pointer with a COPY of the matrix data clamped and cast into  unsigned chars
 * @returns a unique pointer with a COPY of the matrix data clamped and cast into  unsigned chars
 */
std::unique_ptr<unsigned char[]> imageMat::get_pixel_data() const {
    std::unique_ptr<unsigned char[]> copy_data = std::make_unique<unsigned char[]>(n_elements);
    
    for (int i = 0; i < n_elements; ++i) {
        double clamped_val = std::max(0.0, std::min(255.0, m_data[i]));
        unsigned char char_val = static_cast<unsigned char>(clamped_val);
        copy_data[i] = char_val;
    }

    return copy_data;
}

//TODO: Finish implementation
bool imageMat::inverse() {}

#endif