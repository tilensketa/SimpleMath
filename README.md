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

// Creating a matrix (type int, 5x3 matrix)
sm::Matrix<int> matrix(5, 3);

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
### Trigonometry
```c++
float angleDeg = 30.0f;
float angleRad = sm::tg::degToRad(angleDeg);
float sine = sm::tg::sin(angleRad);
float cosine = sm::tg::cos(angleRad);
float tangent = sm::tg::tan(angleRad);
```
### Vector
```c++
// Create vector
sm::Vector<int> vec1(3); // elements are initialized with 0
sm::Vector<int> vec2(3, 3); // element are initialized with 3 (value)
sm::Vector<int> vec3({1,2,3}); // elements are initialized with {1,2,3}
sm::Vector<int> vec4 = {1,2,3}; // elements are initialized with {1,2,3}

// Manipulating vector with other vector, matrix, scalar
sm::Matrix<int> matrix(2, 3, 4);
auto res1 = vec1 + vec2;
auto res2 = vec1 - vec2;
auto res5 = matrix * vec1;
auto res6 = vec1 * 2;
auto res7 = vec1 / 2;

auto homogenousVector = vec4.Homogeneous();
auto normalizedVector = vec4.Normalize();
double length = vec4.Length();
bool orthoVector = vec4.IsOrthogonal();

// Manipulating elements
vec1[0] = 10;
```
### Matrix
```c++
// Create matrix
sm::Matrix<int> mat1(3, 3); // 3x3 matrix, elements are initialized with 0
sm::Matrix<int> mat2(3, 8, 8); // 3x8 matrix, element are initialized with 8 (value)
sm::Matrix<int> mat3(3, 3, {1,2,3,4,5,6,7,8,9}); // 3x3 matrix, elements are initialized with {1,2,3,4,5,6,7,8,9}
sm::Matrix<int> mat4(3, sm::Identity) // 3x3 identity matrix (diagonal ones)
sm::Matrix<int> mat5(sm::RotationZ, 30.0f) // rotationZ matrix (rotate z: 30 deg)
sm::Matrix<int> mat6(sm::Translation2D, 2.0f, 1.0f) // translation2D matrix (translate x: 2.0f, y: 1.0f)

// Manipulating matrix with other matrix, scalar
auto res1 = mat1 + mat4;
auto res2 = mat1 - mat4;
auto res3 = mat1 * mat2;
auto res4 = mat1 * 2;
auto res5 = mat1 / 2;

auto homogenousMatrix = mat1.Homogenous();
int determinant = mat1.Determinant();
auto transposedMatrix = mat2.Transpose();
auto subMatrix = mat3.SubMatrix(0, 1); // row and column to ignore
auto inversedMatrix = mat3.Inverse();
bool orthoMatrix = mat3.IsOrthogonal();

// Manipulating elements
mat1[0][2] = 10;
```

### Algorithms
- Thomas algorithm (Tridiagonal matrix algorithm)
- Factorial
- Absolute
- Power
- Taylor series
- Newton method
- Sqrt
- Cross product
- Dot product

## Authors
me :)