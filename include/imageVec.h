#ifndef IMAGEVEC_H
#define IMAGEVEC_H

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <cstdint>
#include <bit>
#include <limits>
#include "vector_exception.h"
#include <imageMat.h>

template <class T>
class I_Vector {
    public:
        I_Vector(const uint32_t dimensions = 0, std::unique_ptr<T[]> input_data = nullptr) noexcept;

        size_t get_dims(void) const noexcept;

        T get_element(const uint32_t index) const;

        I_Vector<T> operator+(const I_Vector<T>& rhs) const;
        I_Vector<T> operator-(const I_Vector<T>& rhs) const;
        I_Vector<T> operator*(const T& rhs) const;

        template <class U> friend I_Vector<U> operator*(const U& lhs, const I_Vector<U>& rhs);
        template <class U> friend bool operator==(const I_Vector<U>& lhs, const I_Vector<U>& rhs);

        static T dot(const I_Vector<T>& lhs, const I_Vector<T>& rhs);
        static I_Matrix<T> dot(const I_Vector<T>& lhs, const I_Matrix<T>& rhs);
        static I_Vector<T> cross(const I_Vector<T>& lhs, const I_Vector<T>& rhs);
        
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
bool operator==(const I_Vector<T>& lhs, const I_Vector<T>& rhs) {
    if (lhs.get_dims() != rhs.get_dims()) {return false; }

    for (size_t i = 0; i < lhs.get_dims(); ++i) {
        if (lhs.get_element(i) != rhs.get_element(i)) {return false; }
    }

    return true;
}

template <class T>
I_Matrix<T> I_Vector<T>::transpose(void) {
    auto trans_data = std::make_unique<T[]>(m_dims);
    std::copy(m_data.get(), m_data.get() + m_dims, trans_data.get());

    I_Matrix<T> transpose(1, m_dims, std::move(trans_data));
    return transpose;
}



inline double fast_sqrt(const double number) {
    static_assert(std::numeric_limits<float>::is_iec559, "IEEE 754 must be used!");

    if (number < 0) {return std::numeric_limits<double>::quiet_NaN();}

    constexpr double smallest_normal = std::numeric_limits<double>::min();
    if (number < smallest_normal) { return 0;}


    double d_number = number;
    double half = d_number * 0.5;
    //magic = 3/2 * 2^52 * (1023 - 0.043)
    constexpr uint64_t magic = 0x5FE339A0219FF02D;
    constexpr double three_halves = 1.5;

    //Evil bit hack
    uint64_t reinterpret = std::bit_cast<uint64_t>(d_number);
    reinterpret = magic - (reinterpret >> 1);
    d_number = std::bit_cast<double>(reinterpret);

    //Newton iterations to get a closer approximation
    d_number = d_number * (three_halves - (half * d_number * d_number));
    d_number = d_number * (three_halves - (half * d_number * d_number));
    d_number = d_number * (three_halves - (half * d_number * d_number));
    d_number = d_number * (three_halves - (half * d_number * d_number));

    //S * S ^ -1/2 = S ^ 1/2
    return d_number * number;
}

inline double fast_sqrt(const int number) {
    static_assert(std::numeric_limits<float>::is_iec559, "IEEE 754 must be used!");

    if (number < 0) {return std::numeric_limits<double>::quiet_NaN();}


    double d_number = static_cast<double>(number);
    double half = d_number * 0.5;
    //magic = 3/2 * 2^52 * (1023 - 0.043)
    constexpr uint64_t magic = 0x5FE339A0219FF02D;
    constexpr double three_halves = 1.5;

    //Evil bit hack
    uint64_t reinterpret = std::bit_cast<uint64_t>(d_number);
    reinterpret = magic - (reinterpret >> 1);
    d_number = std::bit_cast<double>(reinterpret);

    //Newton iterations to get a closer approximation
    d_number = d_number * (three_halves - (half * d_number * d_number));
    d_number = d_number * (three_halves - (half * d_number * d_number));
    d_number = d_number * (three_halves - (half * d_number * d_number));
    d_number = d_number * (three_halves - (half * d_number * d_number));

    //S * S ^ -1/2 = S ^ 1/2
    return d_number * static_cast<double>(number);
}

#endif