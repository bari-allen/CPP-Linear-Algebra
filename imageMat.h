#ifndef IMAGEMAT_H
#define IMAGEMAT_H

#include <cstring>
#include <stdexcept>
#include <array>

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
    n_elements = n_rows * n_cols;
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
    n_elements = n_rows * n_cols;

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
    n_elements = n_rows * n_cols;

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
    int n_elements = n_rows * n_cols;
    unsigned char* temp_results = new unsigned char[n_elements];

    for (int i = 0; i < n_elements; ++i) {
        temp_results[i] = lhs.m_data[i] + rhs.m_data[i];
    }

    imageMat result = imageMat(n_rows, n_cols, temp_results);
    delete[] temp_results;
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