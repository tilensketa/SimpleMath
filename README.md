# SimpleMath

A simple, lightweight and header-only math library.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Features](#features)
- [Authors](#authors)

## Installation

Add ```SimpleMath.h``` file in your project and you are good to go.

## Usage

Example of how to use library:

```c++
#include "SimpleMath.h"

// Creating a vector (type float, size 3)
sm::Vector<float, 3> vector;

// Creating a matrix (type int, size 5)
sm::Matrix<int, 5> matrix;

// Manipulating matrix
auto result1 = matrix * vector;
auto result2 = matrix / vector;
matrix[0][1] = 741;

// Manipulating vector
auto result3 = vector * 3;
auto result4 = vector / 5;
vector[2] = 741;
```

## Features
### Vector
```c++
// Create vector
sm::Vector<int, 3> vec1; // elements are initialized with 0
sm::Vector<int, 3> vec2(3); // element are initialized with 3 (value)
sm::Vector<int, 3> vec3({1,2,3}); // elements are initialized with {1,2,3} (initializer list)

// Manipulating vector with other vector, matrix, scalar
sm::Matrix<int, 3> matrix(2);
auto res1 = vec1 + vec2;
auto res2 = vec1 - vec2;
auto res3 = vec1 * vec2;
auto res4 = vec1 / vec2;
auto res5 = matrix * vec1;
auto res6 = matrix / vec1;
auto res7 = vec1 * 2;
auto res8 = vec1 / 2;

// Manipulating elements
vec1[0] = 10;
```
### Matrix (support only square matrix!)
```c++
// Create matrix
sm::Matrix<int, 3> mat1; // elements are initialized with 0
sm::Matrix<int, 3> mat2(8); // element are initialized with 8 (value)
sm::Matrix<int, 3> mat3({1,2,3}); // elements are initialized with {1,2,3,4,5,6,7,8,9} (initializer list)

// Manipulating matrix with other matrix, scalar
auto res1 = mat1 + mat2;
auto res2 = mat1 - mat2;
auto res3 = mat1 * mat2;
auto res4 = mat1 / mat2;
auto res5 = mat1 * 2;
auto res6 = mat1 / 2;

// Manipulating elements
mat1[0][2] = 10;
```

## Authors
me :)