# SimpleMath

A simple, lightweight and header-only math library.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Features](#features)
      - [Trigonometry](#trigonometry)
      - [Vector](#vector)
      - [Matrix](#matrix)
- [Authors](#authors)

## Installation

Add ```SimpleMath.h``` file in your project and you are good to go.

## Usage

Example of how to use library:

```c++
#include "SimpleMath.h"

// Create vector with size 3 (elements initialized with 0)
sm::Vector myVector(3);

// Change elements
myVector[0] = 4;
myVector[2] = 1;

// Print vector
std::cout << myVector << std::endl;

// Normalize vector and get vector length
sm::Vector myNormalizedVector = myVector.Normalize();
double myVectorLength = myVector.Length();

// Create matrix 3x3 (elements initialized with initializer list)
sm::Matrix myMatrix(3, 3, { 2,5,8,8,7,6,2,4,6 });

// Change elements
myMatrix[0][0] = 4;
myMatrix[2][1] = 9;

// Print matrix
std::cout << myMatrix << std::endl;

// Multiply matrix and vector
sm::Vector myResult = myMatrix * myVector;
```

## Features
### Trigonometry
```c++
double angle1Deg = 30.0;
double angle2Rad = PI / 3;

// Angle conversion
double angle1Rad = sm::tg::degToRad(angle1Deg);
double angle2Deg = sm::tg::radToDeg(angle2Rad);

// Trigonometric functions
double sine = sm::tg::sin(angle1Rad);
double cosine = sm::tg::cos(angle1Rad);
double tangent = sm::tg::tan(angle1Rad);
double cotangent = sm::tg::cot(angle1Rad);
```
### Vector
```c++
// Constructors
sm::Vector vector1(3);
sm::Vector vector2(3, 5);
sm::Vector vector3({1,2,3});
sm::Vector copyVector = vector1;

// Manipulation
vector1[0] = 4;
sm::Vector vectorResult1 = vector1 + vector2;
sm::Vector vectorResult2 = vector1 - vector2;
sm::Vector vectorResult3 = vector2 * 2.0;
sm::Vector vectorResult4 = vector2 / 2.0;

// Other
sm::Vector homogenousVector = vector3.Homogeneous();
double vectorLength = vector3.Length();
bool vectorOrthogonal = vector3.IsOrthogonal();
sm::Vector normalizedVector = vector3.Normalize();
```
### Matrix
```c++
sm::Matrix matrix1(3, 4);
sm::Matrix matrix2(2, 2, 5);
sm::Matrix matrix3(2, 2, { 1,2,3,4 });
sm::Matrix matrix4(3, 3, { 2,5,8,5,3,5,8,9,7 });
sm::Matrix identityMatrix(3);
sm::Matrix rotationMatrix(sm::RotationZ, 30.0);
sm::Matrix transformationMatrix(sm::Transformation2D, 30.0, 2.0, 3.0);
sm::Matrix copyMatrix = matrix3;

// Manipulation
matrix1[0][1] = 2;
sm::Matrix matrixResult1 = matrix2 + matrix3;
sm::Matrix matrixResult2 = matrix2 - matrix3;
sm::Matrix matrixResult3 = matrix2 * matrix3;
sm::Vector vec({ 1,2 });
sm::Vector vecResult = matrix3 * vec;

// Other
bool matrixOrthogonal = matrix4.IsOrthogonal();
sm::Matrix homogenousMatrix = matrix4.Homogeneous();
sm::Matrix subMatrix = matrix4.SubMatrix(0, 0);
sm::Matrix transposedMatrix = matrix4.Transpose();
double matrixDeterminant = matrix4.Determinant();
sm::Matrix inverseMatrix = matrix4.Inverse();
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