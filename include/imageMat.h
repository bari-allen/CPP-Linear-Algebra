#ifndef IMAGEMAT_H
#define IMAGEMAT_H

#include <memory>
#include <cstring>
#include <stdexcept>
#include <type_traits>
#include <thread>
#include <vector>

/**
 * A matrix class that will include functions to perform simple linear algebra operations
 * The purpose of this class is to allow for PCA transformations of images
 * @author Karl Haidinyak, bari-allen (GitHub)
 */
template <class T>
class I_Matrix {
    //Public Functions
    public:
        //Define the various constructors
        /**The default constructor which should never typically be used */
        I_Matrix();
        /**A constructor which initializes the number of rows and columns with
        the inputted integers */
        I_Matrix(int n_rows, int n_cols);
        /**This constructor initializes the number of rows and columns with the
        inputted integers and also initializes the data to the inputted pointer */
        I_Matrix(int n_rows, int n_cols, std::unique_ptr<T[]> input_data);
        /**The copy constructor */
        I_Matrix(const I_Matrix<T>& input_matrix);

        /**
         * Resizes the matrix to the given number of rows and columns and initializes
         * each element using the default value (calls default constructor)
         * 
         * If a number less than or equal to 0 is inputted then this function automatically
         * returns false
         * 
         * @param n_rows the new number of rows
         * @param n_cols the new number of columns
         * @returns whether the matrix could be resized
         */
        bool resize(int n_rows, int n_cols);

        //Element access functions

        /**
         * Returns the element at the given row and column
         * 
         * Throws an error if the inputted row or column are out of bounds
         * 
         * @param row the row of the desired element
         * @param col the column of the desired element
         * @returns the element at the given row and column
         * @throws invalid_argument is thrown if the inputted row or column are
         *         out of bounds for the matrix
         */
        T get_element(int row, int col) const;

        /**
         * Sets the element at the given row and column to the inputted element
         * 
         * @param row the row of the element to be set
         * @param col the column of the element to be set
         * @param element the new element to be set at the given row and column
         * @returns Returns true if the element was properly set and false if the
         *          inputted row or column are out of bounds
         */
        bool set_element(int row, int col, T element);

        /**
         * Returns the number of rows in the matrix
         * 
         * @returns Returns the number of rows in the matrix
         */
        int rows() const;

        /**
         * Returns the number of columns in the matrix
         * 
         * @returns Returns the number of columns in the matrix
         */
        int cols() const;

        /**
         * Returns a copy of the matrix data
         * 
         * @returns Returns a copy of the matrix data
         */
        std::unique_ptr<T[]> get_elements() const;

        /**
         * Determines if the inputted matrix is equal to this matrix by checking for
         * equality element-wise
         * 
         * The inputted matrix must match the shape of this matrix
         * 
         * @param rhs the inputted matrix
         * @returns Returns true if the inputted matrix is equal to this matrix and 
         *          false otherwise or if the shapes mismatch
         */
        bool operator== (const I_Matrix<T>& rhs) const  ;

        /**
         * Adds two matrices element-wise
         * 
         * The shape of the two matrices must match in order to add them
         * 
         * @param lhs the lefthand matrix of the addition
         * @param rhs the righthand matrix of the addition
         * @return the element-wise addition of the two inputted matrices
         * @throws invalid_argument is thrown if the shape of the two matrices mismatch
         */
        template <class U> friend I_Matrix<U> operator+ (const I_Matrix<U>& lhs, const I_Matrix<U>& rhs);

        /**
         * Adds a constant to each element of the matrix
         * 
         * @param lhs the constant being added to each element of the matrix
         * @param rhs the matrix
         * @return the matrix with the constant added to it
         */
        template <class U> friend I_Matrix<U> operator+ (const U& lhs, const I_Matrix<U>& rhs);

        /**
         * Adds a constant to each element of the matrix
         * 
         * @param lhs the matrix
         * @param rhs the constant being added to each element of the matrix
         * @return the matrix with the constant added to it
         */
        template <class U> friend I_Matrix<U> operator+ (const I_Matrix<U>& lhs, const U& rhs);

        /**
         * Substracts two matrices element-wise
         * 
         * The shape of the two matrices must match in order to subtract them
         * 
         * @param lhs the lefthand matrix of the subtraction
         * @param rhs the righthand matrix of the subtraction
         * @return the element-wise subtraction of the two inputted matrices
         * @throws invalid_argument is thrown if the the shape of the two matrices mismatch
         */
        template <class U> friend I_Matrix<U> operator- (const I_Matrix<U>& lhs, const I_Matrix<U>& rhs);

        /**
         * Subtracts a constant from each element of the matrix
         * 
         * @param lhs the constant being subtracted from each element of the matrix
         * @param rhs the matrix
         * @return the matrix with the constant subtracted from it
         */
        template <class U> friend I_Matrix<U> operator- (const U& lhs, const I_Matrix<U>& rhs);

        /**
         * Subtracts a constant from each element of the matrix
         * 
         * @param lhs the matrix
         * @param rhs the constant being subtracted from each element of the matrix
         * @return the matrix with the constant subtracted from it
         */
        template <class U> friend I_Matrix<U> operator- (const I_Matrix<U>& lhs, const U& rhs);

        /**
        * Computes the dot product between the two inputted matrices 
        * 
        * The number of rows of the lefthand side must match the number of columns in the righthand side
        * 
        * @param lhs the lefthand side of the multiplication
        * @param rhs the righthand side of the multiplication
        * @return the matrix multiplication between the inputted matrices
        * @throws invalid_argument is thrown if the number of rows of the lefthand side do not match
        *         the number of columns of the righthand side
        */
        template <class U> friend I_Matrix<U> operator* (const I_Matrix<U>& lhs, const I_Matrix<U>& rhs);

        /**
         * Multiplies a constant to each element of the matrix
         * 
         * @param lhs the constant being multiplied to each element of the matrix
         * @param rhs the matrix
         * @return the matrix with the constant multiplied to each element
         */
        template <class U> friend I_Matrix<U> operator* (const U& lhs, const I_Matrix<U>& rhs);

        /**
         * Multiplies a constant to each element of the matrix
         * 
         * @param lhs the matrix
         * @param rhs the constant being multiplied to each element of the matrix
         * @return the matrix with the constant multiplied to each element
         */
        template <class U> friend I_Matrix<U> operator* (const I_Matrix<U>& lhs, const U& rhs);

        /**
         * Finds the determinant of the inputted matrix
         * For the documentation about the algorithm, see the det_helper documentation
         * 
         * @param mat the inputted matrix
         * @return the determinant of the matrix
         * @throws throws a logic_error if the inputted matrix is not square
         */
        template <class U> friend double det(const I_Matrix<U>& mat);

        /**
         * Finds the inverse of the inputted matrix
         * For the documentation about the algorithm, see the det_helper documentation
         * 
         * Note that all matrices, no matter the inputted type, will be converted to a matrix of doubles
         * 
         * @param mat the inputted matrix
         * @return the inverse of the matrix
         * @throws throws logic_error if the inputted matrix is not square or is not invertible
         */
        template <class U> friend I_Matrix<double> inv(const I_Matrix<U>& mat);

        /**
         * Returns the transpose of the given matrix
         * @returns the transpose of the given matrix
         */
        I_Matrix<T> transpose(void) const;

    //Private functions
    private:
        /**
         * Returns the linear index given the inputted row and column
         * 
         * @param row the row of the desired index
         * @param col the column of the desired index
         * @returns the linear index given the inputted row and column or -1 if the 
         *          row and/or column are out of bounds
         */
        int linear_index(int row, int col) const;

        /**
         * A recursive helper function that finds the determinant of a square matrix
         * 
         * This function uses a recursive cofactor expansion to shrink larger matrices into
         * 2x2 matrices to find their determinants
         * 
         * This function does not use the I_Matrix object since it is easier to work directly
         * with a pointer
         * 
         * A link to cofactor expansion:
         * https://textbooks.math.gatech.edu/ila/determinants-cofactors.html
         * 
         * @param size the number of rows and columns the square matrix has
         * @param data a pointer to the matrix data
         * @returns the determinant of the sub-matrix
         */
        static double det_helper(const int size, const std::unique_ptr<T[]>& data, const int depth = 0);


        /**
         * Returns the cofactors of the given point in the matrix data excluding the data points 
         * from the excluded row and column
         * 
         * @param size the number of both rows and columns of the matrix
         * @param excluded_col the column whos data is not included in the cofactor data
         * @param data the input data
         * @param excluded_row a optional parameter to exclude data from the inputted row, default in 0
         * @returns the list of cofactors of the given data point
         */
        static std::unique_ptr<T[]> get_cofactors(const int size, const int excluded_col,  
                                                const std::unique_ptr<T[]>& data, const int excluded_row = 0);

        /**
         * A function that assists with inverting the given matrix
         * 
         * Uses the adjugate matrix to calculate the inverse of the matrix then multiplies the matrix by the
         * inverse of the determinant
         * 
         * If the matrix is singular or near-singluar then an exception is thrown
         * 
         * A link to the adjugate algorithm:
         * https://www.mathsisfun.com/algebra/matrix-inverse-minors-cofactors-adjugate.html
         * 
         * @param size the number of rows and columns of the matrix
         * @param data the matrix data
         * @returns the inverse of the matrix if one exists
         * @throws throws logic_error if the matrix is singular or near-singular
         */
        static I_Matrix<double> inv_helper(const int size, const std::unique_ptr<T[]>& data);

    //Private variables
    private:
        std::unique_ptr<T[]> matrix_data;
        int m_rows, m_cols, m_elements;
};

template <class T>
I_Matrix<T>::I_Matrix() {
    m_rows = 0;
    m_cols = 0;
    m_elements = 0;
    matrix_data = std::make_unique<T[]>(0);
}

template <class T>
I_Matrix<T>::I_Matrix(int n_rows, int n_cols) {
    m_rows = n_rows;
    m_cols = n_cols;
    m_elements = m_rows * m_cols;
    matrix_data = std::make_unique<T[]>(m_elements);
}

template <class T>
I_Matrix<T>::I_Matrix(int n_rows, int n_cols, std::unique_ptr<T[]> input_data) {
    m_rows = n_rows;
    m_cols = n_cols;
    m_elements = m_rows * m_cols;
    matrix_data = std::move(input_data);
}

template <class T>
I_Matrix<T>::I_Matrix(const I_Matrix<T>& input_matrix) {
    m_rows = input_matrix.rows();
    m_cols = input_matrix.cols();
    m_elements = m_rows * m_cols;
    std::copy(input_matrix.get_elements().get(), input_matrix.get_elements().get() + m_elements, matrix_data.get());
}

template <class T> 
bool I_Matrix<T>::resize(int n_rows, int n_cols) {
    if (n_rows <= 0 || n_cols <= 0) {
        return false;
    }

    m_rows = n_rows;
    m_cols = n_cols;
    m_elements = m_cols * m_rows;
    matrix_data = std::make_unique<T[]>(m_elements);

    return true;
}

template <class T>
T I_Matrix<T>::get_element(int row, int col) const {
    int index = linear_index(row, col);
    if (index == -1) {
        throw std::invalid_argument("The number of rows or columns is invalid!");
    }

    return matrix_data[index];
}

template <class T>
bool I_Matrix<T>::set_element(int row, int col, T element) {
    int index = linear_index(row, col);
    if (index == -1) return false;

    matrix_data[index] = element;
    return true;
}

template <class T>
int I_Matrix<T>::rows() const {
    return m_rows;
}

template <class T>
int I_Matrix<T>::cols() const {
    return m_cols;
}

template <class T>
std::unique_ptr<T[]> I_Matrix<T>::get_elements() const {
    auto copy_data = std::make_unique<T[]>(m_elements);
    std::copy(matrix_data.get(), matrix_data.get() + m_elements, copy_data.get());

    return copy_data;
}

template <class T>
bool I_Matrix<T>::operator==(const I_Matrix<T>& rhs) const {
    if (m_rows != rhs.m_rows || m_cols != rhs.m_cols) {
        return false;
    }

    std::unique_ptr<T[]> rhs_data = rhs.get_elements();
    for (int i = 0; i < m_elements; ++i) {
        if (!std::is_floating_point<T>::value && matrix_data[i] != rhs_data[i]) {
            return false;
        } else if (matrix_data[i] - rhs_data[i] > 0.0000000001 
                || rhs_data[i] - matrix_data[i] > 0.0000000001) {
                return false;
        }
    }

    return true;
}

template <class T>
I_Matrix<T> operator+(const I_Matrix<T>& lhs, const I_Matrix<T>& rhs) {
    if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
        throw std::invalid_argument("Size mismatch between arguments");
    }

    auto data = std::make_unique<T[]>(lhs.rows() * lhs.cols());

    for (int i = 0; i < (lhs.rows() * lhs.cols()); ++i) {
        data[i] = lhs.matrix_data[i] + rhs.matrix_data[i];
    }

    I_Matrix<T> matrix(lhs.rows(), lhs.cols(), data);
    return matrix;
}

template <class T>
I_Matrix<T> operator+(const T& lhs, const I_Matrix<T>& rhs) {
    int rows = rhs.rows();
    int cols = rhs.cols();
    int n_elements = rows * cols;
    auto data = std::make_unique<T[]>(n_elements);

    for (int i = 0; i < n_elements; ++i) {
        data[i] = lhs + rhs.matrix_data[i];
    }

    I_Matrix<T> matrix(rows, cols, data);
    return matrix;
}

template <class T>
I_Matrix<T> operator+(const I_Matrix<T>& lhs, const T& rhs) {
    return rhs + lhs;
}

template <class T>
I_Matrix<T> operator-(const I_Matrix<T>& lhs, const I_Matrix<T>& rhs) {
    if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
        throw std::invalid_argument("Size mismatch between arguments!");
    }

    auto data = std::make_unique<T[]>(lhs.rows() * lhs.cols());
    auto lhs_data = lhs.matrix_data;
    auto rhs_data = rhs.matrix_data;

    for (int i = 0; i < lhs.rows() * lhs.cols(); ++i) {
        data[i] = lhs_data[i] - rhs_data[i];
    }

    I_Matrix<T> matrix(lhs.rows(), lhs.cols(), data);
    return matrix;
}

template <class T>
I_Matrix<T> operator-(const I_Matrix<T>& lhs, const T& rhs) {
    int rows = lhs.rows();
    int cols = lhs.cols();
    int n_elements = rows * cols;
    auto data = std::make_unique<T[]>(n_elements);

    for (int i = 0; i < n_elements; ++i) {
        data[i] = lhs.matrix_data[i] - rhs;
    }

    I_Matrix<T> matrix(rows, cols, data);
    return matrix;
}

template <class T>
I_Matrix<T> operator-(const T& lhs, const I_Matrix<T>& rhs) {
    int rows = rhs.rows();
    int cols = rhs.cols();
    int n_elements = rows * cols;
    auto data = std::make_unique<T[]>(n_elements);

    for (int i = 0; i < n_elements; ++i) {
        data[i] = lhs - rhs.matrix_data[i];
    }

    I_Matrix<T> matrix(rows, cols, data);
    return matrix;
}

template <class T>
I_Matrix<T> operator*(const I_Matrix<T>& lhs, const I_Matrix<T>& rhs) {
    if (lhs.cols() != rhs.rows()) {
        throw std::invalid_argument("Left-hand side's rows do not match right-hand side's columns");
    }

    int m_rows = lhs.rows();
    int m_cols = rhs.cols();
    int n_elements = m_rows * m_cols;
    auto data = std::make_unique<T[]>(n_elements);

    //Iterates over every row of the lefthand matrix
    for (int lhs_row = 0; lhs_row < lhs.rows(); ++lhs_row) {
        //Iterates over every column of the righthand matrix
        for (int rhs_col = 0; rhs_col < rhs.cols(); ++rhs_col) {

            T cum_sum{}; //The sum of each lefthand row multiplication

            //Iterates over every row of the righthand matrix
            for (int rhs_row = 0; rhs_row < rhs.rows(); ++rhs_row) {
                int lhs_index = rhs_row + (lhs_row * lhs.cols());
                int rhs_index = rhs_col + (rhs_row * rhs.cols());

                T sum = rhs.matrix_data[rhs_index] * lhs.matrix_data[lhs_index];
                cum_sum += sum;
            }

            int data_index = rhs_col + (lhs_row * m_cols);
            data[data_index] = cum_sum;
        }
    }

    I_Matrix<T> matrix(m_rows, m_cols, std::move(data));
    return matrix;
}

template <class T>
I_Matrix<T> operator*(const T& lhs, const I_Matrix<T>& rhs) {
    int rows = rhs.rows();
    int cols = rhs.cols();
    int n_elements = rows * cols;
    auto data = std::make_unique<T[]>(n_elements);

    for (int i = 0; i < n_elements; ++i) {
        data[i] = rhs.matrix_data[i] * lhs;
    }

    I_Matrix<T> matrix(rows, cols, std::move(data));
    return matrix;
}

template <class T>
I_Matrix<T> operator*(const I_Matrix<T>& lhs, const T& rhs) {
    int rows = lhs.rows();
    int cols = lhs.cols();
    int n_elements = rows * cols;
    auto data = std::make_unique<T[]>(n_elements);

    for (int i = 0; i < n_elements; ++i) {
        data[i] = lhs.matrix_data[i] * rhs;
    }

    I_Matrix<T> matrix(rows, cols, data);
    return matrix;
}

template <class T>
I_Matrix<T> I_Matrix<T>::transpose() const {
    auto transpose_data = std::make_unique<T[]>(m_rows * m_cols);

    for (int row = 0; row < m_rows; ++row) {
        for (int col = 0; col < m_cols; ++col) {
            int index = col + (row * m_cols);
            int trans_index = row + (col * m_rows);

            transpose_data[trans_index] = matrix_data[index];
        }
    }

    I_Matrix<T> transpose_mat(m_cols, m_rows, std::move(transpose_data));

    return transpose_mat;
}

template <class T>
std::unique_ptr<T[]> I_Matrix<T>::get_cofactors(const int size, const int excluded_col, 
                                    const std::unique_ptr<T[]>& data, 
                                    const int excluded_row) {
    int i = 0;
    auto cofactors = std::make_unique<T[]>((size - 1) * (size - 1));

    for (int row = 0; row < size; ++row) {
        if (row == excluded_row) {continue;}
        for (int col = 0; col < size; ++col) {
            if (col == excluded_col) {continue;}

            int index = col + (row * size);
            cofactors[i++] = data[index];
        }
    }

    return cofactors;
}

constexpr int MAX_THREAD_DEPTH = 1;

template <class T>
double I_Matrix<T>::det_helper(const int size, const std::unique_ptr<T[]>& data, const int depth) {
    //If the matrix is a 1x1 matrix we return the only value in the matrix
    if (size == 1) {
        return static_cast<double>(data[0]);
    }

    //If the matrix is a 2x2 matrix we compute the determinant like normal
    if (size == 2) {
        auto value = (data[0] * data[3]) - (data[1] * data[2]);
        return static_cast<double>(value);
    }

    std::vector<std::thread> threads;
    std::vector<double> results(size, 0.0);
    double det = 0.0;

    for (int excluded = 0; excluded < size; ++excluded) {
        auto cofactors = I_Matrix<T>::get_cofactors(size, excluded, data);

        int multiplier = (excluded % 2 == 0) ? 1 : -1;
        auto value = data[excluded];
        
        if (depth >= MAX_THREAD_DEPTH || size <= 5) {
            results[excluded] = static_cast<double>((multiplier * value) 
                * det_helper(size - 1, cofactors, depth + 1));
            continue;
        }

        threads.emplace_back([=, cofactors = std::move(cofactors), &results]() {
            results[excluded] = static_cast<double>((multiplier * value) 
                * det_helper(size - 1, cofactors, depth + 1));
        });
    }

    for (auto& thread : threads) {
        if (thread.joinable())
            thread.join();
    }

    for (double partial : results) {
        det += partial;
    }

    return det;
}


template <class T>
double det(const I_Matrix<T>& mat) {
    static_assert(std::is_arithmetic<T>::value, "Type must be an arithmetic type");

    if (mat.rows() != mat.cols()) {
        throw std::logic_error("Cannot compute the determinant of a non-square matrix!");
    }

    if (mat.rows() == 0 || mat.cols() == 0) {
        throw std::logic_error("Cannot compute the determinant of a matrix with size 0");
    }

    return I_Matrix<T>::det_helper(mat.rows(), mat.matrix_data);
}

template <class T>
I_Matrix<double> I_Matrix<T>::inv_helper(const int size, const std::unique_ptr<T[]>& data) {
    auto inv_data = std::make_unique<double[]>(size * size);
    double det = det_helper(size, data);
        
    if (det == 0.0 || det < 0.0000001) {
        throw std::logic_error("Matrix is singular or near-singular!");
    }


    if (size == 1) {
        inv_data[0] = 1 / data[0];
        I_Matrix<double> mat(size, size, std::move(inv_data));
    }

    if (size == 2) {
        inv_data[0] = static_cast<double>(data[3]);
        inv_data[1] = static_cast<double>(-data[1]);
        inv_data[2] = static_cast<double>(-data[2]);
        inv_data[3] = static_cast<double>(data[0]);

        I_Matrix<double> mat(size, size, std::move(inv_data));
        return (1 / det) * mat;
    }

    for (int excluded_row = 0; excluded_row < size; ++excluded_row) {
        for (int excluded_col = 0; excluded_col < size; ++excluded_col) {

            auto cofactors = get_cofactors(size, excluded_col, data, excluded_row);
            double cofactor_det = det_helper(size - 1, cofactors);
            int multiplier = (excluded_col + excluded_row) % 2 == 0 ? 1 : -1;

            int index = excluded_col + (excluded_row * size);
            inv_data[index] = multiplier * cofactor_det;
        }
    }

    I_Matrix<double> mat(size, size, std::move(inv_data));
    auto adjugate_mat = mat.transpose();
    return (1.0 / det) * adjugate_mat;
}

template <class T>
I_Matrix<double> inv(const I_Matrix<T>& mat) {
    static_assert(std::is_arithmetic<T>::value, "Type must be an arithmetic type");

    if (mat.rows() != mat.cols()) {
        throw std::logic_error("Cannot compute the inverse of a non-square matrix!");
    }

    if (mat.rows() == 0 || mat.cols() == 0) {
        throw std::logic_error("Cannot compute the inverse of a matrix with size 0!");
    }

    return I_Matrix<T>::inv_helper(mat.rows(), mat.matrix_data);
}

template <class T>
int I_Matrix<T>::linear_index(int row, int col) const{
    if (row > m_rows || col > m_cols || row < 0 || col < 0) {
        return -1;
    }

    return col + (row * m_cols);
}

#endif