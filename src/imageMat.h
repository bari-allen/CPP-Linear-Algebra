#ifndef IMAGEMAT_H
#define IMAGEMAT_H

#include <cstring>
#include <stdexcept>
#include <omp.h>

class imageMat {
    public:
        //Below are some various constructors that may or may not be useful
        imageMat();
        imageMat(int nRows, int nCols);
        //Use an array of unsigned chars because each pixel in BMP is represented
        //and 3 BGR hex values
        imageMat(int nRows, int nCols, const unsigned int* input_data);
        imageMat(const imageMat& inputMat);

        //The destructor
        ~imageMat();

        //Configuration functions
        bool resize(int nRows, int nCols);

        //Element access functions
        unsigned int get(int row, int col);
        bool set(int row, int col, unsigned int rgb_data);
        int get_rows();
        int get_cols();

        //Overload the equals operator
        bool operator== (const imageMat& rhs) const;

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
        unsigned int* m_data;
        int n_rows, n_cols, n_elements;
        

};

//The default constructor
imageMat::imageMat() {
    n_rows = 0;
    n_cols = 0;
    n_elements = 0;
    m_data = new unsigned int[n_elements];
}

//Constructor that creates an a flattened matrix containing only white pixels
imageMat::imageMat(int input_rows, int input_cols) {
    n_rows = input_rows;
    n_cols = input_cols;
    n_elements = n_rows * n_cols;
    m_data = new unsigned int[n_elements];
    unsigned int white_pixel = 0xFFFFFF;

    for (int i = 0; i < n_elements; ++i) {
        m_data[i] = white_pixel;
    }
}

//Constructor that creates a flattened matrix containing the input_data
imageMat::imageMat(int input_rows, int input_cols, const unsigned int* input_data) {
    if (input_data == nullptr) {
        throw std::invalid_argument("The input_data cannot be null!");
    }

    n_rows = input_rows;
    n_cols = input_cols;
    n_elements = n_rows * n_cols;

    m_data = new unsigned int[n_elements];

    std::memcpy(m_data, input_data, n_elements * sizeof(unsigned int));
}

//Copy constructor
imageMat::imageMat(const imageMat& inputMat) {
    n_rows = inputMat.n_rows;
    n_cols = inputMat.n_cols;
    n_elements = inputMat.n_elements;
    m_data = new unsigned int[n_elements];

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

    unsigned int white_pixel = 0xFFFFFF;

    m_data = new unsigned int[n_elements];
    if (m_data != nullptr) {
        for (int i = 0; i < n_elements; ++i) {
            m_data[i] = white_pixel;
        }

        return true;
    } else {
        return false;
    }
}

unsigned int imageMat::get(int row, int col) {
    int linear_index = get_range_start(row, col);

    unsigned int pixel = 0x000000;

    if (linear_index >= 0) {
        pixel = m_data[linear_index];
    }

    return pixel;
}

bool imageMat::set(int row, int col, unsigned int rgb_data) {
    int linear_index = get_range_start(row, col);

    if (linear_index >= 0) {
        m_data[linear_index] = rgb_data;
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

bool imageMat::operator== (const imageMat& rhs) const {
    if (this->n_rows != rhs.n_rows || this->n_cols != rhs.n_cols) {
        return false;
    }

    for (int i = 0; i < this->n_elements; ++i) {
        if (this->m_data[i] != rhs.m_data[i]) {
            return false;
        }
    }

    return true;
}

imageMat operator+(const imageMat& lhs, const imageMat& rhs) {
    if (lhs.n_cols != rhs.n_cols || lhs.n_rows != rhs.n_rows) {
        throw std::invalid_argument("Both sides must have the same size!");
    }

    int n_rows = lhs.n_rows;
    int n_cols = lhs.n_cols;
    int n_elements = n_rows * n_cols;
    unsigned int* temp_results = new unsigned int[n_elements];

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
    int n_elements = n_rows * n_cols;
    unsigned int* temp_results = new unsigned int[n_elements];

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
    int n_elements = n_rows * n_cols;
    unsigned int* temp_results = new unsigned int[n_elements];

    for (int i = 0; i < n_elements; ++i) {
        temp_results [i] = lhs.m_data[i] + rhs;
    }

    imageMat result(n_rows, n_cols, temp_results);
    delete[] temp_results;
    return result;
}

//Takes the dot product of two imageMat objects
imageMat operator*(const imageMat& lhs, const imageMat& rhs) {
    int r_rows = rhs.n_rows;
    int r_cols = rhs.n_cols;
    int l_rows = lhs.n_rows;
    int l_cols = lhs.n_cols;

    if (l_cols != r_rows) {
        throw std::invalid_argument("Dimension mismatch!");
    }

    unsigned int* temp_result = new unsigned int[lhs.n_rows * rhs.n_cols];

    #pragma omp parallel for
    for (int lhs_row = 0; lhs_row < l_rows; ++lhs_row) {
        for (int rhs_col = 0; rhs_col < r_cols; ++rhs_col) {

            unsigned int color_result = 0x000000;

            for (int lhs_col = 0; lhs_col < l_cols; ++lhs_col) {

                //Same liear_indexing algorithm used below
                int lhs_linear_index = ((lhs_row * l_cols) + lhs_col);

                //Different linear indexing algorithm is used since the rhs has
                //the same number of rows as the lhs has columns
                int rhs_linear_index = ((lhs_col * r_cols) + rhs_col);

                unsigned int color_val = (lhs.m_data[lhs_linear_index] * rhs.m_data[rhs_linear_index]);
                
                //Clamps the result to the max value of each pixel (white pixel)
                unsigned int clamp_size = 0xFFFFFF;
                color_result = std::min(clamp_size, color_val + color_result);
            }

            int result_linear_index = ((lhs_row * r_cols) + rhs_col);
                temp_result[result_linear_index] = color_result;
        }
    }

    imageMat result(l_rows, r_cols, temp_result);
    delete[] temp_result;
    return result;
}



//Returns the flattened index of the 3D matrix
int imageMat::get_range_start(int row, int col) {
    if ((row < n_rows) && (row >= 0) && (col < n_cols) && (col >= 0)) {
        return ((row * n_cols) + col);
    } else {
        return -1;
    }
}


#endif