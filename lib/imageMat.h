#ifndef IMAGEMAT_H
#define IMAGEMAT_H

#include <memory>
#include <cstring>
#include <stdexcept>

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
        I_Matrix(int n_rows, int n_cols, std::unique_ptr<T[]>& input_data);
        /**The copy constructor */
        I_Matrix(const I_Matrix<T>& input_matrix);

        //Configuration functions
        bool resize(int n_rows, int n_cols);

        //Element access functions
        T get_element(int row, int col) const;
        bool set_element(int row, int col, T element);
        int rows() const;
        int cols() const;
        std::unique_ptr<T[]> get_elements() const;

        //Overload == operator
        bool operator== (const I_Matrix<T>& rhs) const  ;

        //Overload +, -, and * operators (friends)

        /**Adds the two matrices element-wise
         * Both matrices must have the same shape
         */
        template <class U> friend I_Matrix<U> operator+ (const I_Matrix<U>& lhs, const I_Matrix<U>& rhs);
        template <class U> friend I_Matrix<U> operator+ (const U& lhs, const I_Matrix<U>& rhs);
        template <class U> friend I_Matrix<U> operator+ (const I_Matrix<U>& lhs, const U& rhs);

        template <class U> friend I_Matrix<U> operator- (const I_Matrix<U>& lhs, const I_Matrix<U>& rhs);
        template <class U> friend I_Matrix<U> operator- (const U& lhs, const I_Matrix<U>& rhs);
        template <class U> friend I_Matrix<U> operator- (const I_Matrix<U>& lhs, const U& rhs);

        template <class U> friend I_Matrix<U> operator* (const I_Matrix<U>& lhs, const I_Matrix<U>& rhs);
        template <class U> friend I_Matrix<U> operator* (const U& lhs, const I_Matrix<U>& rhs);
        template <class U> friend I_Matrix<U> operator* (const I_Matrix<U>& lhs, const U& rhs);

    //Private functions
    private:
        int linear_index(int row, int col);

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
I_Matrix<T>::I_Matrix(int n_rows, int n_cols, std::unique_ptr<T[]>& input_data) {
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
        if (matrix_data[i] != rhs_data[i]) {
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

    for (int lhs_row = 0; lhs_row < lhs.rows(); ++lhs_row) {
        for (int rhs_col = 0; rhs_col < rhs.cols(); ++rhs_col) {
            T cum_sum{};
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

    I_Matrix<T> matrix(m_rows, m_cols, data);
    return matrix;
}

template <class T>
int I_Matrix<T>::linear_index(int row, int col) {
    if (row > m_rows || col > m_cols || row < 0 || col < 0) {
        return -1;
    }

    return col + (row * m_cols);
}

#endif