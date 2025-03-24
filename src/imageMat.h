#ifndef IMAGEMAT_H
#define IMAGEMAT_H

#include <cstring>
#include <stdexcept>
#include <array>
#include <omp.h>
#include <gtest/gtest.h>

class imageMat {
    public:
        //Below are some various constructors that may or may not be useful
        imageMat();
        imageMat(int nRows, int nCols);
        //Use an array of unsigned chars because each pixel in BMP is represented
        //and 3 BGR hex values
        imageMat(int nRows, int nCols, const unsigned char* input_data);
        imageMat(const imageMat& inputMat);

        //The destructor
        ~imageMat();

        //Configuration functions
        bool resize(int nRows, int nCols);

        //Element access functions
        std::array<unsigned char, 3> get(int row, int col);
        bool set(int row, int col, unsigned char* rgb_data);
        int get_rows();
        int get_cols();

        //Overload the equals operator
        bool operator== (const imageMat& rhs);

        //Overload the +, -, and * operators
        //First takes two matrices as inputs
        //Second takes a scalar and a matrix as inputs
        //Third takes a matrix and a scalar as inputs
        friend imageMat operator+ (const imageMat& lhs, const imageMat& rhs);
        friend imageMat operator+ (const unsigned char& lhs, const imageMat& rhs);
        friend imageMat operator+ (const imageMat& lhs, const unsigned char& rhs);

        friend imageMat operator- (const imageMat& lhs, const imageMat& rhs);
        friend imageMat operator- (const unsigned char& lhs, const imageMat& rhs);
        friend imageMat operator- (const imageMat& lhs, const unsigned char& rhs);

        friend imageMat operator* (const imageMat& lhs, const imageMat& rhs);
        friend imageMat operator* (const unsigned char& lhs, const imageMat& rhs);
        friend imageMat operator* (const imageMat& lhs, const unsigned char& rhs);

    private:
        //Get the range of index values for a given row and column in the
        //flattened array
        int get_range_start(int row, int col);
        
    private:
        unsigned char* m_data;
        int n_rows, n_cols, n_elements;
        

};

//The default constructor
imageMat::imageMat() {
    n_rows = 0;
    n_cols = 0;
    n_elements = 0;
    m_data = new unsigned char[n_elements];
}

//Constructor that creates an a flattened matrix containing only white pixels
imageMat::imageMat(int input_rows, int input_cols) {
    n_rows = input_rows;
    n_cols = input_cols;
    n_elements = n_rows * n_cols * 3;
    m_data = new unsigned char[n_elements];
    unsigned char white_pixel[3] = {0xFF, 0xFF, 0xFF};

    for (int i = 0; i < n_elements - 3; i += 3) {
        std::memcpy(m_data + i, white_pixel, 3 * sizeof(unsigned char));
    }
}

//Constructor that creates a flattened matrix containing the input_data
imageMat::imageMat(int input_rows, int input_cols, const unsigned char* input_data) {
    if (input_data == nullptr) {
        throw std::invalid_argument("The input_data cannot be null!");
    }

    n_rows = input_rows;
    n_cols = input_cols;
    n_elements = n_rows * n_cols * 3;

    m_data = new unsigned char[n_elements];

    std::memcpy(m_data, input_data, n_elements * sizeof(unsigned char));
}

//Copy constructor
imageMat::imageMat(const imageMat& inputMat) {
    n_rows = inputMat.n_rows;
    n_cols = inputMat.n_cols;
    n_elements = inputMat.n_elements;
    m_data = new unsigned char[n_elements];

    for (int i = 0; i < n_elements; ++i) {
        m_data[i] = inputMat.m_data[i];
    }

}

//Destructor
imageMat::~imageMat() {
    if (m_data != nullptr) {
        delete[] m_data;
    }
}

//Deletes all the elements in the matrix, resizes, and the makes each pixel white
bool imageMat::resize(int r_rows, int r_cols) {
    n_rows = r_rows;
    n_cols = r_cols;
    n_elements = n_rows * n_cols * 3;

    if (m_data != nullptr) {
        delete[] m_data;
    }

    unsigned char white_pixel[3] = {0xFF, 0xFF, 0xFF};

    m_data = new unsigned char[n_elements];
    if (m_data != nullptr) {
        for (int i = 0; i < n_elements - 3; i += 3) {
            std::memcpy(m_data + i, white_pixel, 3 * sizeof(unsigned char));
        }

        return true;
    } else {
        return false;
    }
}

std::array<unsigned char, 3> imageMat::get(int row, int col) {
    int linear_index = get_range_start(row, col);

    std::array<unsigned char, 3> pixel = {0x00, 0x00, 0x00};

    if (linear_index >= 0) {
        std::copy(m_data + linear_index, m_data + linear_index + 3, pixel.begin());
    }

    return pixel;
}

bool imageMat::set(int row, int col, unsigned char* rgb_data) {
    int linear_index = get_range_start(row, col);

    if (linear_index >= 0) {
        std::memcpy(m_data + linear_index, rgb_data, 3 * sizeof(unsigned char));
        return true;
    }

    return false;
}

int imageMat::get_rows() {
    return n_rows;
}

int imageMat::get_cols() {
    return n_cols;
}

imageMat operator+(const imageMat& lhs, const imageMat& rhs) {
    if (lhs.n_cols != rhs.n_cols || lhs.n_rows != rhs.n_rows) {
        throw std::invalid_argument("Both sides must have the same size!");
    }

    int n_rows = lhs.n_rows;
    int n_cols = lhs.n_cols;
    int n_elements = n_rows * n_cols * 3;
    unsigned char* temp_results = new unsigned char[n_elements];

    for (int i = 0; i < n_elements; ++i) {
        temp_results[i] = lhs.m_data[i] + rhs.m_data[i];
    }

    imageMat result = imageMat(n_rows, n_cols, temp_results);
    delete[] temp_results;
    return result;
}

imageMat operator+(const unsigned char& lhs, const imageMat& rhs) {
    int n_rows = rhs.n_rows;
    int n_cols = rhs.n_cols;
    int n_elements = n_rows * n_cols * 3;
    unsigned char* temp_results = new unsigned char[n_elements];

    for (int i = 0; i < n_elements; ++i) {
        temp_results[i] = lhs + rhs.m_data[i];
    }

    imageMat result(n_rows, n_cols, temp_results);
    delete[] temp_results;
    return result;
}

imageMat operator+(const imageMat& lhs, const unsigned char& rhs) {
    int n_rows = lhs.n_rows;
    int n_cols = lhs.n_cols;
    int n_elements = n_rows * n_cols * 3;
    unsigned char* temp_results = new unsigned char[n_elements];

    for (int i = 0; i < n_elements; ++i) {
        temp_results [i] = lhs.m_data[i] + rhs;
    }

    imageMat result(n_rows, n_cols, temp_results);
    delete[] temp_results;
    return result;
}

imageMat operator*(const imageMat& lhs, const imageMat& rhs) {
    int r_rows = rhs.n_rows;
    int r_cols = rhs.n_cols;
    int l_rows = lhs.n_rows;
    int l_cols = lhs.n_cols;

    if (l_cols != r_rows) {
        throw std::invalid_argument("Dimension mismatch!");
    }

    unsigned char* temp_result = new unsigned char[lhs.n_rows * rhs.n_cols * 3];

    #pragma omp parallel for
    for (int lhs_row = 0; lhs_row < l_rows; lhs_row += 3) {
        for (int rhs_col = 0; rhs_col < r_cols; rhs_col += 3) {

            unsigned char blue_result = 0x00;
            unsigned char green_result = 0x00;
            unsigned char red_result = 0x00;

            for (int lhs_col = 0; lhs_col < l_cols; lhs_col += 3) {
                int lhs_linear_index = 3 * ((lhs_row * l_cols) + lhs_col);
                int rhs_linear_index = 3 * ((lhs_col * r_cols) + rhs_col);

                unsigned char blue_val = (lhs.m_data[lhs_linear_index] * rhs.m_data[rhs_linear_index]);
                unsigned char green_val = (lhs.m_data[lhs_linear_index + 1] * rhs.m_data[rhs_linear_index + 1]);
                unsigned char red_val = (lhs.m_data[lhs_linear_index + 2] * rhs.m_data[rhs_linear_index + 2]);

                blue_result = std::min(255, int(blue_result) + int(blue_val));
                green_result = std::min(255, int(green_result) + int(green_val));
                red_result = std::min(255, int(red_result) + int(red_val));
            }

            int result_linear_index = 3 * ((lhs_row * r_cols) + rhs_col);
                temp_result[result_linear_index] = blue_result;
                temp_result[result_linear_index + 1] = green_result;
                temp_result[result_linear_index + 2] = red_result;
        }
    }

    imageMat result(l_rows, r_cols, temp_result);
    delete[] temp_result;
    return result;
}



//Returns the flattened index of the 3D matrix
//Similar to the formula index = (row * n_cols) + col, but multiplies by 3
//because each element in the 2D matrix contains 3 elements
int imageMat::get_range_start(int row, int col) {
    if ((row < n_rows) && (row >= 0) && (col < n_cols) && (col >= 0)) {
        return 3 * ((row * n_cols) + col);
    } else {
        return -1;
    }
}


#endif