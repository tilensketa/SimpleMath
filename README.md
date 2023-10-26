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
sm::Vector<float> vector(3);

// Creating a matrix (type int, size 5)
sm::Matrix<int> matrix(5);

// Manipulating matrix
auto result1 = matrix * vector;
auto result2 = matrix * 3;
auto result3 = matrix / 5;
matrix[0][1] = 741;

// Manipulating vector
auto result4 = vector * 3;
auto result5 = vector / 5;
vector[2] = 741;
```

## Features
### Vector
```c++
// Create vector
sm::Vector<int> vec1(3); // elements are initialized with 0
sm::Vector<int> vec2(3, 3); // element are initialized with 3 (value)
sm::Vector<int> vec3({1,2,3}); // elements are initialized with {1,2,3}
sm::Vector<int> vec4 = {1,2,3}; // elements are initialized with {1,2,3}

// Manipulating vector with other vector, matrix, scalar
sm::Matrix<int> matrix(3, 2);
auto res1 = vec1 + vec2;
auto res2 = vec1 - vec2;
auto res5 = matrix * vec1;
auto res6 = vec1 * 2;
auto res7 = vec1 / 2;

// Manipulating elements
vec1[0] = 10;
```
### Matrix (support only square matrix!)
```c++
// Create matrix
sm::Matrix<int> mat1(3); // matrix size 3, elements are initialized with 0
sm::Matrix<int> mat2(3, 8); // matrix size 3, element are initialized with 8 (value)
sm::Matrix<int> mat3({1,2,3,4,5,6,7,8,9}); // elements are initialized with {1,2,3,4,5,6,7,8,9}, matrix size is calculated
sm::Matrix<int> mat4 = {1,2,3,4,5,6,7,8,9}; // elements are initialized with {1,2,3,4,5,6,7,8,9}, matrix size is calculated

// Manipulating matrix with other matrix, scalar
auto res1 = mat1 + mat2;
auto res2 = mat1 - mat2;
auto res3 = mat1 * mat2;
auto res4 = mat1 * 2;
auto res5 = mat1 / 2;

// Manipulating elements
mat1[0][2] = 10;
```

## Authors
me :)