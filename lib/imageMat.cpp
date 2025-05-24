#include "imageMat.h"

template <class T>
I_Matrix<T>::I_Matrix() {
    m_rows = 0;
    m_cols = 0;
    m_elements = 0;
    matrix_data = nullptr;
}

template <class T>
I_Matrix<T>::I_Matrix(int n_rows, int n_cols) {
    m_rows = n_rows;
    m_cols = n_cols;
    m_elements = m_rows * m_cols;
    matrix_data = new T[m_elements];
}

template <class T>
I_Matrix<T>::I_Matrix(int n_rows, int n_cols, std::unique_ptr<T> input_data) {
    m_rows = n_rows;
    m_cols = n_cols;
    m_elements = m_rows * m_cols;
    matrix_data = input_data.get();
}

template <class T>
I_Matrix<T>::I_Matrix(const I_Matrix<T>& input_matrix) {
    m_rows = input_matrix->m_rows;
    m_cols = input_matrix->m_cols;
    m_elements = m_rows * m_cols;
    matrix_data = std::memcpy(matrix_data, input_matrix->matrix_data, sizeof(T) * m_elements);
}

template <class T>
I_Matrix<T>::~I_Matrix() {
    delete[] matrix_data;
}

template <class T> 
bool I_Matrix<T>::resize(int n_rows, int n_cols) {
    m_rows = n_rows;
    m_cols = n_cols;
    m_elements = m_cols * m_rows;
    delete[] matrix_data;

    matrix_data = new(std::nothrow) T[m_elements];
    return matrix_data == nullptr;
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
int I_Matrix<T>::get_num_rows() const {
    return m_rows;
}

template <class T>
int I_Matrix<T>::get_num_cols() const {
    return m_cols;
}

template <class T>
std::unique_ptr<T> I_Matrix<T>::get_elements() const {
    T* copy = new T[m_elements];
    std::memcpy(copy, matrix_data, sizeof(T) * m_elements);
    return std::make_unique(copy);
}

template <class T>
int I_Matrix<T>::linear_index(int row, int col) {
    if (row > m_rows || col > m_cols || row < 0 || col < 0) {
        return -1;
    }

    return row + (col * m_cols);
}