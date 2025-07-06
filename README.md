# CPP-Linear-Algebra<a name="CPP-Linear-Algebra"></a>

A simple linear algebra library that contains basic linear algebra functions to manipulate data

## Building

```bash
mkdir build
cmake -S . -B build
cmake --build build
```

## Usage

```c++
auto matrix_data = std::make_unique<int[]>(4);
//create a 2x2 matrix with data initialized to 0
I_Matrix<int> matrix(2, 2, std::move(matrix_data));

auto vector_data = std::make_unique<int[]>(3);
//create a 3-dimensional vector with data initialized to 0
I_Vector<int> vector(3, std::move(vector_data));
```

## Features 
This project includes various linear algebra functions but some highlights are:
### Eigenvalues
```c++
//I_Matrix<double> matrix;
I_Vector<double> eigenvalues = I_Matrix<double>::eig(matrix);
```
### Matrix Inverses
```c++
//I_Matrix<double> matrix
I_Matrix<double> inverse_matrix = inv(matrix);
```
### Matrix Determinants
```c++
//I_Matrix<double> matrix
double determinant = det(matrix);
```
### Matrix QR Factorization
```c++
//I_Matrix<double> matrix
auto QR_tuple = I_Matrix<double>::QR(matrix);
```

## Testing
After building the project run this to run the matrix tests
```bash
cd build && ./matrix_test && cd../
```
After building the project run this to run the vector tests
```bash
cd build && ./vector_tests && cd ../
```
