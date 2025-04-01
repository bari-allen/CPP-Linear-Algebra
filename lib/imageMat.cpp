#include "imageMat.h"

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
    m_data = new double[n_elements];

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
    m_data = new double[n_elements];

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
    m_data = new double[n_elements];

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
std::unique_ptr<unsigned char[]> imageMat::get_pixel_data() const noexcept(true){
    auto copy_data = std::make_unique<unsigned char[]>(n_elements);
    
    for (int i = 0; i < n_elements; ++i) {
        double clamped_val = std::max(0.0, std::min(255.0, m_data[i]));
        auto char_val = static_cast<unsigned char>(clamped_val);
        copy_data[i] = char_val;
    }

    return copy_data;
}

//TODO: Finish implementation
bool imageMat::inverse() {return false;}

/**
 * Determines if two matrices are equal by seeing if each of their corresponding elements
 * are close enough in value (difference is less than 1e-9) to each other
 * @param rhs the other matrix
 * @returns whether the two matrices are equal or close enough to equal (tolerance of 1e-9 difference)
 */
bool imageMat::operator==(const imageMat& rhs) const{
    if (this->n_rows != rhs.get_rows() || this->n_cols != rhs.get_cols()) {
        return false;
    }

    for (std::size_t i = 0; i < n_elements; ++i) {
        auto elem_1 = m_data[i];
        auto elem_2 = rhs.get_linear(i);

        if (!close_enough(elem_1, elem_2)) {
            return false;
        }
    }

    return true;
}

/**
 * Compares two matrices by finding the standard deviation of the differences between their values
 * and checking whether it is less than the inputted tolerance
 * If both the matrices are 1x1, then we use the distance between their only values and check if that
 * is less than the tolerance
 * @param rhs the other matrix being compared to
 * @param tolerance the exclusive maximum that the standard deviation can be
 * @returns whether the two matrices are similar
 */
bool imageMat::within_tolerance(const imageMat& rhs, double tolerance) const noexcept{
    auto rhs_rows = rhs.get_rows();
    auto rhs_cols = rhs.get_cols();

    if(rhs_rows != n_rows || rhs_cols != n_cols) {
        return false;
    }

    if (n_elements == 1) {
        auto elem_1 = m_data[0];
        auto elem_2 = rhs.get_linear(0);

        return fabs(elem_1 - elem_2) < tolerance;
    }

    auto cum_sum = 0.0;
    for (std::size_t i = 0; i < n_elements; ++i) {
        auto elem_1 = m_data[i];
        auto elem_2 = rhs.get_linear(i);
        auto difference = elem_1 - elem_2;
        
        cum_sum += difference * difference;
    }

    auto final_val = sqrt(cum_sum / ((n_elements) - 1));

    return final_val < tolerance;
}

/**
 * Determine whether the distance between the two inputs is less than 1e-9
 * @param val the first value
 * @param other_val the second value
 * @returns whether the distance between the two inputs is less than 1e-9
 */
bool imageMat::close_enough(double val, double other_val) const {
    return fabs(val - other_val) < 1e-9;
}

imageMat operator+(const imageMat& lhs, const imageMat& rhs) {
    if (lhs.get_rows() != rhs.get_rows() || lhs.get_cols() != rhs.get_cols()) {
        throw std::logic_error("Both matrices must be the same shape");
    }

    int n_rows = lhs.get_rows();
    int n_cols = lhs.get_cols();
    int n_elements = n_rows * n_cols;
    double* data = new double[n_elements];

    for (std::size_t i = 0; i < n_elements; ++i) {
        auto lhs_val = lhs.get_linear(i);
        auto rhs_val = rhs.get_linear(i);

        auto new_val = lhs_val + rhs_val;
        data[i] = new_val;
    }

    imageMat return_mat(n_rows, n_cols, data);
    delete[] data;
    data = nullptr;

    return return_mat;
}

imageMat operator+(const imageMat& lhs, const double rhs) {
    int n_rows = lhs.get_rows();
    int n_cols = lhs.get_cols();
    int n_elements = n_rows * n_cols;
    double* data = new double[n_elements];

    for (std::size_t i = 0; i < n_elements; ++i) {
        auto lhs_val = lhs.get_linear(i);
        auto new_val = lhs_val + rhs;

        data[i] = new_val;
    }

    imageMat return_mat(n_rows, n_cols, data);
    delete[] data;
    data = nullptr;

    return return_mat;
}

imageMat operator-(const imageMat& lhs, const imageMat& rhs) {
    if (lhs.get_rows() != rhs.get_rows() || lhs.get_cols() != rhs.get_cols()) {
        throw std::logic_error("Both matrices must be the same shape");
    }

    int n_rows = lhs.get_rows();
    int n_cols = lhs.get_cols();
    int n_elements = n_rows * n_cols;
    double* data = new double[n_elements];

    for (std::size_t i = 0; i < n_elements; ++i) {
        auto lhs_val = lhs.get_linear(i);
        auto rhs_val = rhs.get_linear(i);
        auto new_val = lhs_val - rhs_val;

        data[i] = new_val;
    }

    imageMat return_mat(n_rows, n_cols, data);
    delete[] data;
    data = nullptr;

    return return_mat;
}

imageMat operator-(const imageMat& lhs, const double rhs) {
    int n_rows = lhs.get_rows();
    int n_cols = lhs.get_cols();
    int n_elements = n_rows * n_cols;
    double* data = new double[n_elements];

    for (std::size_t i = 0; i < n_elements; ++i) {
        auto lhs_val = lhs.get_linear(i);
        auto new_val = lhs_val - rhs;

        data[i] = new_val;
    }

    imageMat return_mat(n_rows, n_cols, data);
    delete[] data;
    data = nullptr;

    return return_mat;
}

/**
 * Returns the linear index given the row and column numbers
 * @param row the row of the element
 * @param col the column of the element
 * @returns the linear index to one of the elements
 */
int imageMat::get_linear_index(int row, int col) const noexcept {
    if (row < 0 || col < 0 || row >= n_rows || row >= n_cols) {
        return -1;
    }

    return (row * n_cols) + col;
}

/**
 * Returns the element at the inputted linear index
 * @param linear_index the linear index of the element
 * @returns the element at the inputted linear index
 * @throws invalid_argument when the linear index is greater than the number of elements or less than 0
 */
double imageMat::get_linear(int linear_index) const {
    if (linear_index >= n_elements || linear_index < 0) {
        throw std::invalid_argument("Index out of bounds!");
    }
    return m_data[linear_index];
}

/**
 * Determines whether the matrix is a square matrix
 * @returns whether the matrix is a square matrix
 */
bool imageMat::is_square() const noexcept {
    return n_rows == n_cols;
}