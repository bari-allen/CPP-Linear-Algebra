#ifndef IMAGEVEC_H
#define IMAGEVEC_H

#include <memory>
#include <stdexcept>
#include <type_traits>
#include <cstdint>
#include <bit>
#include <cmath>
#include <limits>

template <class T>
class I_Vector {
    public:
        I_Vector(const uint32_t dimensions = 0, const std::unique_ptr<T[]>& input_data = nullptr) noexcept;

        uint32_t get_dims(void) const;

        T get_element(const uint32_t index) const;

        I_Vector<T> operator+(const I_Vector<T>& rhs) const;
        I_Vector<T> operator-(const I_Vector<T>& rhs) const;
        I_Vector<T> operator*(const T& rhs) const;

        template <class U> friend I_Vector<U> operator*(const U& lhs, const I_Vector<U>& rhs);
        template <class U> friend bool operator==(const I_Vector<U>& lhs, const I_Vector<U>& rhs);

        static T dot(const I_Vector<T>& lhs, const I_Vector<T>& rhs);
        static I_Vector<T> cross(const I_Vector<T>& lhs, const I_Vector<T>& rhs);        
    
    private:
        std::unique_ptr<T[]> m_data;
        int m_dims;
};

template <class T>
I_Vector<T>::I_Vector(const uint32_t dimensions, const std::unique_ptr<T[]>& input_data) noexcept
    : m_dims(dimensions)  {
        if (input_data != nullptr) {
            m_data = std::move(input_data);
        } else {
            m_data = std::make_unique<T[]>(m_dims);
        }

}

template <class T>
uint32_t I_Vector<T>::get_dims() const{
    return m_dims;
}

template <class T>
T I_Vector<T>::get_element(const uint32_t index) const{
    if (index >= m_dims)
        throw std::invalid_argument("Index out of bounds");

    return m_data[index];
}



template <class T>
double fast_sqrt(const T number) {
    static_assert(std::is_arithmetic<T>::value, "Inputted type must be a number!");
    static_assert(std::numeric_limits<float>::is_iec559, "IEEE 754 must be used!");

    if (number < 0) [[unlikely]]
        throw std::invalid_argument("The input cannot be negative!");

    if (std::is_floating_point_v<T> && std::fpclassify(number) != FP_NORMAL) {
        return std::sqrt(number);
    }

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