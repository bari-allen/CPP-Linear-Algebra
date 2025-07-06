#ifndef IMAGEMAT_H
#define IMAGEMAT_H

#include <memory>
#include <cstring>
#include <stdexcept>
#include <type_traits>
#include <thread>
#include <vector>
#include "vector_exception.h"
#include "../src/fast_math.cpp"
#include <tuple>

template <class T>
class I_Vector;

/**
 * A matrix class that will include functions to perform simple linear algebra operations
 * @author Karl Haidinyak, bari-allen (GitHub)
 */

/*******************************************************************************
*I_Matrix class definitions
*******************************************************************************/
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
        /**Constructor to cast the data from a matrix to type T and copy it */
        template <class U> I_Matrix(const I_Matrix<U>& input_matrix);
        /**Move copying */
        I_Matrix<T>& operator=(I_Matrix<T>&& other) noexcept = default;

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
         * Sets the column at the given column number to the inputted vector
         * 
         * @param col the column number to change
         * @param col_data the vector whos data will be set at the given column
         * @throws  throws invalid_argument if the given column is out of bounds
         */
        void set_col(const uint32_t col, const I_Vector<T>& col_data); 

        /**
         * Sets the row at the given row number to the inputted vector
         * 
         * @param row the row number to change
         * @param row_data the vector whos data will be set at the given row
         * @throws  throws invalid_argument if the given row is out of bounds
         */
        void set_row(const uint32_t row, const I_Vector<T>& row_data);

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
         * Returns a vector of the elements in the inputted column
         * 
         * Throws an error if the inputted column is out of bounds
         * 
         * @param col the column of the desired element
         * @returns a vector of the column data
         * @throws invalid_argument is thrown if the inputted column is
         *         out of bounds for the matrix
         */
        I_Vector<T> get_column(const uint32_t col) const;

        /**
         * Returns a vector of the elements in the inputted row
         * 
         * Throws an error if the inputted row is out of bounds
         * 
         * @param row the row of the desired element
         * @returns a vector of the column data
         * @throws invalid_argument is thrown if the inputted row is
         *         out of bounds for the matrix
         */
        I_Vector<T> get_row(const uint32_t row) const;

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
         * Returns a square, identity matrix with the given size 
         * 
         * @param size the number of rows and columns in the matrix
         * @returns a square, identity matrix*/
        static I_Matrix<T> eye(const uint32_t size);

        /**
         * Returns the transpose of the given matrix
         * @returns the transpose of the given matrix
         */
        I_Matrix<T> transpose(void) const;

        /**
         * Performs QR factorization on the given matrix given that the matrix is square
         * 
         * This function will throw an exception if the inputted matrix is not square
         * 
         * This function performs QR factorization using Givens Rotations due to its
         * numberical stability and speed
         *  
         * Algorithm used is adapted from this video:
         * https://youtu.be/kyG8YMIfNA0
         * 
         * @param  A the matrix to be factored
         * @returns a tuple containing both the Q and R matrices
         * @throws throws logic_error if the inputted matrix is non-square
         */
        static std::tuple<I_Matrix<double>, I_Matrix<double>> QR(const I_Matrix<T>& A);

        /**
         * @brief Computes the eigenvalues of the given matrix
         * 
         * This function computes the eigenvalues using the QR algorithm with 
         * Wilkison-style shifts where the last element on the diagonal is used
         * for the shift matrix
         * 
         * A link to the Wikipedia page for the QR algorithm:
         * https://en.wikipedia.org/wiki/QR_algorithm
         * 
         * @param A the input matrix
         * @param tolerance the percision of the last eigenvalue in the matrix
         * @return the eigenvalues of the matrix
         * @throws throws invalid_argument exception if the A matrix is non-square
         */
        static I_Vector<double> eig(const I_Matrix<T>& A, const double tolerance = 1e-32);

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
        
        /**
         * Returns the sine and cosine of the given sides of the triangle
         * 
         * @param a the opposite side of the given angle
         * @param b the adjacent side of the given angle
         * @returns the sine and cosine from the given sides */
        static std::unique_ptr<double[]> get_angles(const T& a, const  T& b);
        
        /**
         * @brief Computes the RQ matrix of the given matrix
         * 
         * If A = QR then  this function computes RQ which should the converge
         * to a matrix with the eigenvalues along the diagonal
         * 
         * @param A The input matrix
         * @return The RQ matrix of the given matrix
         */
        static I_Matrix<double> make_similar(const I_Matrix<T>& A); 

    //Private variables
    private:
        std::unique_ptr<T[]> matrix_data;
        int m_rows, m_cols, m_elements;
};

/*******************************************************************************
*I_Vector class definitions
*******************************************************************************/

template <class T>
class I_Vector {
    public:
        /** The constructor with default values of 0 and nullptr */
        I_Vector(const uint32_t dimensions = 0, std::unique_ptr<T[]> input_data = nullptr) noexcept;

        /**
         * Returns the number of elements in the vector
         * 
         * @returns the number of elements in the vector*/
        size_t get_dims(void) const noexcept;

        /**
         * Returns the element at the given index
         * 
         * Throws an exception if the index is out of bounds
         * 
         * @param index the desired element's index
         * @throws throws invalid_argument if the given index is out of bounds*/
        T get_element(const uint32_t index) const;

        /**
         * @brief Adds the vector element-wise
         * 
         * @param rhs the right-hand vector
         * @return the element-wise addition of this and the rhs vector
         * @throws throws vector-exception if the two vectors have different dimensions
         */
        I_Vector<T> operator+(const I_Vector<T>& rhs) const;

        /**
         * @brief Subtracts the vector element-wise
         * 
         * @param rhs the right-hand vector
         * @return the element-wise subtraction of this and the rhs vector
         * @throws throws vector-exception if the two vectors have different dimensions
         */
        I_Vector<T> operator-(const I_Vector<T>& rhs) const;

        /**
         * @brief Multiplies this vector with a scaler
         * 
         * @param rhs the right-hand scaler 
         * @return this vector multiplied by the scaler rhs value
         */
        I_Vector<T> operator*(const T& rhs) const noexcept;

        /**
         * @brief Multiplies the rhs vector by the lhs scaler
         * 
         * @tparam U 
         * @param lhs the scaler value
         * @param rhs the vector
         * @return the rhs vector multiplied by the scaler lhs value
         */
        template <class U> friend I_Vector<U> operator*(const U& lhs, const I_Vector<U>& rhs) noexcept;

        /**
         * @brief Checks whether the two matrices are equal
         *
         * Checks for equality by checking for equality element-wise or checks if the absolute value
         * difference is less than a tolerance if the U parameter is a floating point
         * 
         * @tparam U 
         * @param lhs the left-hand side vector
         * @param rhs the right-hand side vector
         * @return true if the two vectors are equal
         * @return false if the two vectors aren't equal
         */
        template <class U> friend bool operator==(const I_Vector<U>& lhs, const I_Vector<U>& rhs) noexcept;

        /**
         * @brief Computes the inner-product of the two matrices
         * 
         * @param lhs the left-hand side vector
         * @param rhs the right-hand side vector
         * @return the inner-product of the vectors
         * @throws throws vector_exception if the two vectors have different dimensions
         */
        static T dot(const I_Vector<T>& lhs, const I_Vector<T>& rhs);

        /**
         * @brief Computes the outer-product of the two matrices
         * 
         * The rhs is a matrix instead of a vector since it is a row-vector instead
         * of a column vector
         * 
         * @param lhs the left-hand side column vector
         * @param rhs the right-hand side row vector
         * @return the outer-product of the two vectors
         * @throws throws vector_exception if the two vectors have different dimensions
         */
        static I_Matrix<T> dot(const I_Vector<T>& lhs, const I_Matrix<T>& rhs);

        /**
         * @brief Computes the cross product of the two vectors
         * 
         * @param lhs the left-hand side vector
         * @param rhs the right-hand side vector
         * @return the cross product of the two vectors
         * @throws throws vector_exception if the lhs and rhs either have a dimension
         *         mismatch or don't have a dimension of 3
         */
        static I_Vector<T> cross(const I_Vector<T>& lhs, const I_Vector<T>& rhs);

        /**
         * @brief Computes the 2-norm of the given vector
         * 
         * @param vec the input vector
         * @return the 2-norm of the input vector
         */
        static double norm(const I_Vector<T>& vec) noexcept;
        
        /**
         * @brief Computes the transpose of the given vector
         * 
         * Since the transpose of a column vector is a row vector an I_Matrix is
         * returned 
         * 
         * @return the transpose of the given vector
         */
        I_Matrix<T> transpose(void) noexcept;

    
    private:
        std::unique_ptr<T[]> m_data;
        int m_dims;
};

/*******************************************************************************
*I_Matrix function implementations
*******************************************************************************/

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
    matrix_data = std::make_unique<T[]>(m_elements);

    for (uint32_t i = 0; i < m_elements; ++i) {
        matrix_data[i] = static_cast<T>(input_matrix.get_elements()[i]);
    }
}

template <class T>
template <class U>
I_Matrix<T>::I_Matrix(const I_Matrix<U>& input_matrix) {
    static_assert(std::is_convertible<U, T>::value, "Types must be convertible!");

    m_rows = input_matrix.rows();
    m_cols = input_matrix.cols();
    m_elements = m_rows * m_cols;
    matrix_data = std::make_unique<T[]>(m_elements);

    for (uint32_t i = 0; i < m_elements; ++i) {
        matrix_data[i] = static_cast<T>(input_matrix.get_elements()[i]);
    }
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
void I_Matrix<T>::set_col(const uint32_t col, const I_Vector<T>& col_vec) {
    if (col >= m_cols) {
        throw std::invalid_argument("Inputted column out of bounds");
    }

    uint32_t data_index{};
    for (uint32_t row = 0; row < m_rows; ++row) {
        uint32_t index = linear_index(row, col);
        matrix_data[index] = col_vec.get_element(data_index++);
    }
}

template <class T>
void I_Matrix<T>::set_row(const uint32_t row, const I_Vector<T>& row_vec) {
    if (row >= m_rows) {
        throw std::invalid_argument("Inputted row out of bounds");
    }

    uint32_t data_index{};
    for (uint32_t col = 0; col < m_cols; ++col) {
        uint32_t index = linear_index(row, col);
        matrix_data[index] = row_vec.get_element(data_index++);
    }
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
I_Vector<T> I_Matrix<T>::get_column(const uint32_t col) const {
    if (col >= m_cols) {
        throw std::invalid_argument("Inputted column out of bounds");
    }

    auto vec_data = std::make_unique<T[]>(m_rows);
    uint32_t data_index{};

    for (uint32_t row = 0; row < m_rows; ++row) {
        uint32_t index = linear_index(row, col);
        vec_data[data_index++] = matrix_data[index];
    }

    I_Vector<T> vec(m_rows, std::move(vec_data));
    return vec;
}

template <class T>
I_Vector<T> I_Matrix<T>::get_row(const uint32_t row) const {
    if (row >= m_rows) {
        throw std::invalid_argument("Inputted row out of bounds");
    }

    auto vec_data = std::make_unique<T[]>(m_cols);
    uint32_t data_index{};

    for (uint32_t col = 0; col < m_cols; ++col) {
        uint32_t index = linear_index(row, col);
        vec_data[data_index++] = matrix_data[index];
    }

    I_Vector<T> vec(m_cols, std::move(vec_data));
    return vec;
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
        }

        double difference = static_cast<double>(matrix_data[i] - rhs_data[i]);
        
        if (std::is_floating_point<T>::value && fast_math::abs(difference) > 0.0000001) {
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

    I_Matrix<T> matrix(lhs.rows(), lhs.cols(), std::move(data));
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

    I_Matrix<T> matrix(rows, cols, std::move(data));
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

    for (int i = 0; i < lhs.rows() * lhs.cols(); ++i) {
        data[i] = lhs.matrix_data[i] - rhs.matrix_data[i];
    }

    I_Matrix<T> matrix(lhs.rows(), lhs.cols(), std::move(data));
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

    I_Matrix<T> matrix(rows, cols, std::move(data));
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

    I_Matrix<T> matrix(rows, cols, std::move(data));
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

    I_Matrix<T> matrix(rows, cols, std::move(data));
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
I_Matrix<T> I_Matrix<T>::eye(const uint32_t size) {
    auto identity_data = std::make_unique<T[]>(size * size);

    for (uint32_t i = 0; i < size; ++i) {
        uint32_t index = i + (i * size);
        identity_data[index] = 1;
    }

    I_Matrix<T> identity_mat(size, size, std::move(identity_data));
    return identity_mat;
}

template <class T>
std::unique_ptr<double[]> I_Matrix<T>::get_angles(const T& a, const T& b) {
    static_assert(std::is_arithmetic<T>::value, "Type must be a number");

    double d_a = static_cast<double>(a);
    double d_b = static_cast<double>(b);

    double hyp = fast_math::fast_sqrt((d_a * d_a) + (d_b * d_b));
    double cos = d_a / hyp;
    double sin = -d_b / hyp;

    auto angles = std::make_unique<double[]>(2);
    angles[0] = cos;
    angles[1] = sin;

    return angles;
}

template <class T>
std::tuple<I_Matrix<double>, I_Matrix<double>> I_Matrix<T>::QR(const I_Matrix<T>& A) {
    if (A.cols() != A.rows()) {
        throw std::logic_error("Cannot return a QR matrix of a non-square matrix");
    }

    uint32_t shape = A.cols();
    I_Matrix<double> R(A);
    I_Matrix<double> Q = I_Matrix<T>::eye(shape);
    std::tuple<I_Matrix<double>, I_Matrix<double>> return_tup;

    for (uint32_t i = 0; i < shape - 1; ++i) {
        for (uint32_t j = i + 1; j < shape; ++j) {
            auto angles = get_angles(R.get_element(i, i), R.get_element(j, i));
            double cos = angles[0];
            double sin = angles[1];

            I_Vector<double> ith_row = R.get_row(i);
            I_Vector<double> jth_row = R.get_row(j);

            R.set_row(i, (cos * ith_row) + ((-sin) * jth_row));
            R.set_row(j, (sin * ith_row) + (cos * jth_row));

            I_Vector<double> ith_col = Q.get_column(i);
            I_Vector<double> jth_col = Q.get_column(j);

            Q.set_col(i, (cos * ith_col) + ((-sin) * jth_col));
            Q.set_col(j, (sin * ith_col) + (cos * jth_col));    
        }
    }

    return_tup = std::make_tuple(std::move(Q), std::move(R));

    return return_tup;
}

template <class T>
I_Matrix<double> I_Matrix<T>::make_similar(const I_Matrix<T>& A) {
    std::tuple<I_Matrix<double>, I_Matrix<double>> QR = I_Matrix<T>::QR(A);
    I_Matrix<double> Q = std::get<0>(QR);
    I_Matrix<double> R = std::get<1>(QR);

    I_Matrix<double> B = R * Q;
    return B;
}

template <class T>
I_Vector<double> I_Matrix<T>::eig(const I_Matrix<T>& A, const double tolerance) {
    static_assert(std::is_arithmetic<T>::value, "Type must be an arithmetic type");

    if (A.rows() != A.cols()) {
        throw std::invalid_argument("Cannot get the eigenvalues of a non-square matrix");
    }

    I_Matrix<double> B = make_similar(A);
    uint32_t iters = 0;
    double last_eig = B.get_element(B.rows() - 1, B.cols() - 1);
    double diff = 1;
    const uint32_t last_row = B.rows() - 1;
    const uint32_t last_col = B.cols() - 1;

    while (diff > tolerance && iters < 100) {
        I_Matrix<double> shift = eye(B.rows()) * B.get_element(last_row, last_col);
        I_Matrix<double> C = B - shift;
        B = make_similar(C) + shift;
        ++iters;
        diff = fast_math::abs(last_eig - B.get_element(last_row, last_col));
        last_eig = B.get_element(last_row, last_col);
    }

    auto eig_data = std::make_unique<double[]>(B.rows());

    for (uint32_t i = 0; i < B.rows(); ++i) {
        eig_data[i] = B.get_element(i, i);
    }

    I_Vector<double> eig(B.rows(), std::move(eig_data));
    return eig;
}


template <class T>
int I_Matrix<T>::linear_index(int row, int col) const{
    if (row > m_rows || col > m_cols || row < 0 || col < 0) {
        return -1;
    }

    return col + (row * m_cols);
}

/*******************************************************************************
*I_Vector function implementations
*******************************************************************************/

template <class T>
I_Vector<T>::I_Vector(const uint32_t dimensions, std::unique_ptr<T[]> input_data) noexcept
    : m_dims(dimensions)  {
        if (input_data != nullptr) {
            m_data = std::move(input_data);
        } else {
            m_data = std::make_unique<T[]>(m_dims);
        }

}

template <class T>
size_t I_Vector<T>::get_dims() const noexcept{
    return m_dims;
}

template <class T>
T I_Vector<T>::get_element(const uint32_t index) const{
    if (index >= m_dims)
        throw std::invalid_argument("Index out of bounds");

    return m_data[index];
}



template <class T>
I_Vector<T> I_Vector<T>::operator+(const I_Vector<T>& rhs) const{
    if (m_dims != rhs.get_dims()) {
        throw vector_exception("Size mismatch!");
    }

    auto new_data = std::make_unique<T[]>(m_dims);

    for (uint32_t i = 0; i < m_dims; ++i) {
        new_data[i] = m_data[i] + rhs.get_element(i);
    }

    I_Vector<T> new_vector(m_dims, std::move(new_data));
    return new_vector;
}

template <class T>
I_Vector<T> I_Vector<T>::operator-(const I_Vector<T>& rhs) const {
    if (m_dims != rhs.get_dims()) {
        throw vector_exception("Size mismatch!");
    }

    auto new_data = std::make_unique<T[]>(m_dims);

    for (uint32_t i = 0; i < m_dims; ++i) {
        new_data[i] = m_data[i] - rhs.get_element(i);
    }

    I_Vector<T> new_vector(m_dims, std::move(new_data));
    return new_vector;
}

template <class T>
I_Vector<T> I_Vector<T>::operator*(const T& rhs) const noexcept{
    auto new_data = std::make_unique<T[]>(m_dims);

    for (uint32_t i = 0; i < m_dims; ++i) {
        new_data[i] = m_data[i] * rhs;
    }

    I_Vector<T> new_vector(m_dims, std::move(new_data));
    return new_vector;
}

template <class T>
I_Vector<T> operator*(const T& lhs, const I_Vector<T>& rhs) noexcept {
    return rhs * lhs;
}

template <class T>
bool operator==(const I_Vector<T>& lhs, const I_Vector<T>& rhs) noexcept {
    if (lhs.get_dims() != rhs.get_dims()) {return false; }

    for (size_t i = 0; i < lhs.get_dims(); ++i) {
        if (lhs.get_element(i) != rhs.get_element(i)) {return false; }
    }

    return true;
}

template <class T>
T I_Vector<T>::dot(const I_Vector<T>& lhs, const I_Vector<T>& rhs) {
    if (lhs.get_dims() != rhs.get_dims()) {
        throw vector_exception("Size mismatch!");
    }

    T sum{};

    for (uint32_t i = 0; i < lhs.get_dims(); ++i) {
        sum += lhs.get_element(i) * rhs.get_element(i);
    }

    return sum;
}

template <class T>
I_Matrix<T> I_Vector<T>::dot(const I_Vector<T>& lhs, const I_Matrix<T>& rhs) {
    if (rhs.rows() != 1) {
        throw vector_exception("Inputted matrix must be a row vector!");
    }

    auto output_data = std::make_unique<T[]>(lhs.get_dims() * rhs.cols());

    for (uint32_t i = 0; i < lhs.get_dims(); ++i) {
        for (uint32_t j = 0; j < rhs.cols(); ++j) {
            uint32_t index = j + (i * rhs.cols());
            output_data[index] = lhs.get_element(i) * rhs.get_element(0, j);
        }
    }

    I_Matrix<T> outer_product(lhs.get_dims(), rhs.cols(), std::move(output_data));
    return outer_product;
}

template <class T>
I_Vector<T> I_Vector<T>::cross(const I_Vector<T>& lhs, const I_Vector<T>& rhs) {
    if (lhs.get_dims() != 3 || rhs.get_dims() != 3) {
        throw vector_exception("Cross product is only defined for 3-dimensional vectors!");
    }

    if (lhs.get_dims() != rhs.get_dims()) {
        throw vector_exception("Dimension mismatch");
    }

    auto cross_data = std::make_unique<T[]>(3);

    cross_data[0] = lhs.get_element(1) * rhs.get_element(2) - lhs.get_element(2) * rhs.get_element(1);
    cross_data[1] = lhs.get_element(2) * rhs.get_element(0) - lhs.get_element(0) * rhs.get_element(2);
    cross_data[2] = lhs.get_element(0) * rhs.get_element(1) - lhs.get_element(1) * rhs.get_element(0);

    I_Vector<T> cross(3, std::move(cross_data));
    return cross;
}

template <class T>
double I_Vector<T>::norm(const I_Vector<T>& vec) noexcept{
    static_assert(std::is_arithmetic<T>::value, "Type must be a number!");
    double square_sum{};

    for (uint32_t i = 0; i < vec.get_dims(); ++i) {
        double val = static_cast<double>(vec.get_element(i));
        square_sum += val * val;
    }

    return fast_math::fast_sqrt(square_sum);
}

template <class T>
I_Matrix<T> I_Vector<T>::transpose(void) noexcept{
    auto trans_data = std::make_unique<T[]>(m_dims);
    std::copy(m_data.get(), m_data.get() + m_dims, trans_data.get());

    I_Matrix<T> transpose(1, m_dims, std::move(trans_data));
    return transpose;
}

#endif