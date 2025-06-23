#ifndef IMAGEVEC_H
#define IMAGEVEC_H

#include <memory>
#include <stdexcept>
#include <type_traits>
#include <cstdint>
#include <bit>

template <class T>
class I_Vector {
    public:
        I_Vector(void);

        I_Vector(const uint32_t dimensions, const std::unique_ptr<T[]>& input_data);

        uint32_t get_dims(void) const;

        T get_element(const uint32_t index) const;

        I_Vector<T> operator+(const I_Vector<T>& rhs) const;
        I_Vector<T> operator-(const I_Vector<T>& rhs) const;
        I_Vector<T> operator*(const T& rhs) const;

        template <class U> friend I_Vector<U> operator*(const U& lhs, const I_Vector<U>& rhs);
        template <class U> friend bool operator==(const I_Vector<U>& lhs, const I_Vector<U>& rhs);

        static T dot(const I_Vector<T>& lhs, const I_Vector<T>& rhs);
        static I_Vector<T> cross(const I_Vector<T>& lhs, const I_Vector<T>& rhs);
        static double fast_sqrt(const T number);

    private:
        std::unique_ptr<T[]> vec_data;
        int dims;
};

template <class T>
I_Vector<T>::I_Vector() {
    dims = 0;
    vec_data = std::make_unique<T[]>(dims);
}

template <class T>
I_Vector<T>::I_Vector(const uint32_t dimensions, const std::unique_ptr<T[]>& input_data) {
    dims = dimensions;
    vec_data = std::move(input_data);
}

template <class T>
uint32_t I_Vector<T>::get_dims() const{
    return dims;
}

template <class T>
T I_Vector<T>::get_element(const uint32_t index) const{
    if (index >= dims) 
        throw std::invalid_argument("Index out of bounds");

    return vec_data[index];
}



template <class T>
float fast_sqrt(const T number) {
    static_assert(std::is_arithmetic<T>::value, "Inputted type must be a number!");

    if (number < 0)
        throw std::invalid_argument("The input cannot be negative!");

    float f_number = static_cast<float>(number);
    float saved = static_cast<float>(number);
    uint32_t reinterpret;
    
    reinterpret = std::bit_cast<uint32_t>(f_number);
    reinterpret = 0x1FBD3F7D + (reinterpret >> 1);
    f_number = std::bit_cast<float>(reinterpret);

    f_number = 0.5 * (f_number + (saved / f_number));
    f_number = 0.5 * (f_number + (saved / f_number));
    f_number = 0.5 * (f_number + (saved / f_number));

    return f_number;
}

#endif