#ifndef IMAGEVEC_H
#define IMAGEVEC_H

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <cstdint>
#include <type_traits>
#include "vector_exception.h"
#include "imageMat.h"
#include "../src/fast_math.cpp"

template <class T>
class I_Vector {
    public:
        I_Vector(const uint32_t dimensions = 0, std::unique_ptr<T[]> input_data = nullptr) noexcept;

        size_t get_dims(void) const noexcept;

        T get_element(const uint32_t index) const;

        I_Vector<T> operator+(const I_Vector<T>& rhs) const;
        I_Vector<T> operator-(const I_Vector<T>& rhs) const;
        I_Vector<T> operator*(const T& rhs) const noexcept;

        template <class U> friend I_Vector<U> operator*(const U& lhs, const I_Vector<U>& rhs) noexcept;
        template <class U> friend bool operator==(const I_Vector<U>& lhs, const I_Vector<U>& rhs) noexcept;

        static T dot(const I_Vector<T>& lhs, const I_Vector<T>& rhs);
        static I_Matrix<T> dot(const I_Vector<T>& lhs, const I_Matrix<T>& rhs);
        static I_Vector<T> cross(const I_Vector<T>& lhs, const I_Vector<T>& rhs);
        static double norm(const I_Vector<T>& vec);
        
        I_Matrix<T> transpose(void);

    
    private:
        std::unique_ptr<T[]> m_data;
        int m_dims;
};

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

    auto cross_data = std::make_unique<T[]>(3);

    cross_data[0] = lhs.get_element(1) * rhs.get_element(2) - lhs.get_element(2) * rhs.get_element(1);
    cross_data[1] = lhs.get_element(2) * rhs.get_element(0) - lhs.get_element(0) * rhs.get_element(2);
    cross_data[2] = lhs.get_element(0) * rhs.get_element(1) - lhs.get_element(1) * rhs.get_element(0);

    I_Vector<T> cross(3, std::move(cross_data));
    return cross;
}

template <class T>
double I_Vector<T>::norm(const I_Vector<T>& vec) {
    static_assert(std::is_arithmetic<T>::value, "Type must be a number!");
    double square_sum{};

    for (uint32_t i = 0; i < vec.get_dims(); ++i) {
        double val = static_cast<double>(vec.get_element(i));
        square_sum += val * val;
    }

    return fast_math::fast_sqrt(square_sum);
}

template <class T>
I_Matrix<T> I_Vector<T>::transpose(void) {
    auto trans_data = std::make_unique<T[]>(m_dims);
    std::copy(m_data.get(), m_data.get() + m_dims, trans_data.get());

    I_Matrix<T> transpose(1, m_dims, std::move(trans_data));
    return transpose;
}

#endif