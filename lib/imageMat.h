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
        I_Matrix(int n_rows, int n_cols, std::unique_ptr<T> input_data);
        /**The copy constructor */
        I_Matrix(const I_Matrix<T>& input_matrix);

        //The destructor
        ~I_Matrix();

        //Configuration functions
        bool resize(int n_rows, int n_cols);

        //Element access functions
        T get_element(int row, int col) const;
        bool set_element(int row, int col, T element);
        int get_num_rows() const;
        int get_num_cols() const;
        std::unique_ptr<T> get_elements() const;

        //Overload == operator
        bool operator== (const I_Matrix<T>& rhs);

        //Overload +, -, and * operators (friends)
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
        T* matrix_data;
        int m_rows, m_cols, m_elements;
};

#endif