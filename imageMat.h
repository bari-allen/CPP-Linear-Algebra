#ifndef IMAGEMAT_H
#define IMAGEMAT_H

#include <cstring>
#include <stdexcept>

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
        unsigned char* get(int row, int col);
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
        int get_range(int row, int col);
        
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



#endif