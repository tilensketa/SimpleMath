#pragma once

#include <iostream>
#include <initializer_list>
#include <cassert>
#include <limits>

typedef double Number;
#define MAX_NUMBER std::numeric_limits<Number>::infinity();
#define MIN_NUMBER -std::numeric_limits<Number>::infinity();

#define PI 3.1415926536

#define PRECISION 0.000001
#define COS_SIN_TAYLOR_PRECISION 31
#define NEWTON_PRECISION 0.0001
#define LOG_TAYLOR_PRECISION 1000

namespace sm {
	// Function declarations
	namespace algo {
		Number TaylorCosSin(int precision, bool prime, Number value);
		Number Absolute(Number value);
		Number Power(Number value, int exponent);
		Number Sqrt(Number value);
	}

	namespace tg {

		// Angle conversion
		Number degToRad(Number angleDeg) {
			while(angleDeg >= 360)
				angleDeg -= 360;
			Number result = angleDeg * (PI / 180);
			return result;
		}
		Number radToDeg(Number angleRad) {
			Number result = angleRad * (180 / PI);
			return result;
		}

		// Trigonometric functions
		Number sin(Number angleRad) {
			if (angleRad >= 0 && angleRad < PI / 2)
				return algo::TaylorCosSin(COS_SIN_TAYLOR_PRECISION, false, angleRad);
			else if (angleRad >= PI / 2 && angleRad < PI)
				return algo::TaylorCosSin(COS_SIN_TAYLOR_PRECISION, false, PI - angleRad);
			else if (angleRad >= PI && angleRad < 1.5 * PI)
				return -algo::TaylorCosSin(COS_SIN_TAYLOR_PRECISION, false, angleRad - PI);
			else
				return -algo::TaylorCosSin(COS_SIN_TAYLOR_PRECISION, false, 2 * PI - angleRad);
		}
		Number cos(Number angleRad) {
			if (angleRad >= 0 && angleRad < PI / 2)
				return algo::TaylorCosSin(COS_SIN_TAYLOR_PRECISION, true, angleRad);
			else if (angleRad >= PI / 2 && angleRad < PI)
				return -algo::TaylorCosSin(COS_SIN_TAYLOR_PRECISION, true, PI - angleRad);
			else if (angleRad >= PI && angleRad < 1.5 * PI)
				return -algo::TaylorCosSin(COS_SIN_TAYLOR_PRECISION, true, angleRad - PI);
			else
				return algo::TaylorCosSin(COS_SIN_TAYLOR_PRECISION, true, 2 * PI - angleRad);
		}
		Number tan(Number angleRad) {
			Number result = sin(angleRad) / cos(angleRad);
			return result;
		}
		Number cot(Number angleRad) {
			Number result = cos(angleRad) / sin(angleRad);
			return result;
		}
	}

	struct Vector {
		int m_Size;
		Number* m_Data;

		Vector(int size) : m_Size(size) {
			m_Data = new Number[m_Size];
			for (int i = 0; i < m_Size; i++) {
				m_Data[i] = (Number)0;
			}
		}
		Vector(int size, Number value) : m_Size(size) {
			m_Data = new Number[m_Size];
			for (int i = 0; i < m_Size; i++) {
				m_Data[i] = value;
			}
		}
		Vector(std::initializer_list<Number> values) : m_Size(values.size()) {
			assert(values.size() == m_Size);
			m_Data = new Number[m_Size];
			int i = 0;
			for (const Number& value : values) {
				m_Data[i++] = value;
			}
		}
		Vector(const Vector& other) : m_Size(other.m_Size) {
			m_Data = new Number[m_Size];
			for (int i = 0; i < m_Size; i++) {
				m_Data[i] = other.m_Data[i];
			}
		}

		~Vector() {
			delete[] m_Data;
		}

		Vector& operator=(const Vector& other) {
			if (this != &other) {
				delete[] m_Data;
				m_Size = other.m_Size;
				m_Data = new Number[m_Size];
				for (int i = 0; i < m_Size; i++) {
					m_Data[i] = other.m_Data[i];
				}
			}
			return *this;
		}

		Number& operator[](int index){
			assert(index >= 0 && index < m_Size);
			return m_Data[index];
		}
		const Number& operator[] (int index) const {
			assert(index >= 0 && index < m_Size);
			return m_Data[index];
		}
		
		Vector operator+(const Vector& other) {
			assert(m_Size == other.m_Size);

			Vector result(other.m_Size);
			for (int i = 0; i < m_Size; i++) {
				result[i] = this->m_Data[i] + other[i];
			}
			return result;
		}
		Vector operator-(const Vector& other) {
			assert(m_Size == other.m_Size);

			Vector result(other.m_Size);
			for (int i = 0; i < m_Size; i++) {
				result[i] = this->m_Data[i] - other[i];
			}
			return result;
		}
		bool operator<(const Vector& other) {
			assert(m_Size == other.m_Size);

			for (int i = 0; i < m_Size; i++) {
				if (m_Data[i] > other.m_Data[i])
					return false;
			}
			return true;
		}
		bool operator>(const Vector& other) {
			assert(m_Size == other.m_Size);

			for (int i = 0; i < m_Size; i++) {
				if (m_Data[i] < other.m_Data[i])
					return false;
			}
			return true;
		}

		friend std::ostream& operator<<(std::ostream& os, const Vector& vec) {
			os << "{ ";
			for (int i = 0; i < vec.m_Size; i++) {
				os << vec[i] << " ";
			}
			os << "}";
			return os;
		}

		Vector operator*(const Number& value) {
			Vector result(m_Size);
			for (int i = 0; i < m_Size; i++) {
				result[i] = m_Data[i] * value;
			}
			return result;
		}
		Vector operator/(const Number& value) {
			Vector result(m_Size);
			for (int i = 0; i < m_Size; i++) {
				result[i] = m_Data[i] / value;
			}
			return result;
		}

		Vector Homogeneous() {
			Vector result(m_Size + 1);
			for (int i = 0; i < m_Size; i++) {
				result.m_Data[i] = m_Data[i];
			}
			result.m_Data[m_Size] = (Number)1;
			return result;
		}
		Number Length() {
			Number variable = 0;
			for (int i = 0; i < m_Size; i++) {
				Number value = m_Data[i];
				variable += algo::Power(value, 2);
			}
			return algo::Sqrt(variable);
		}
		bool IsOrthogonal() {
			if (algo::Absolute(Length() - 1.0) < PRECISION)
				return true;
			return false;
		}
		Vector Normalize() {
			Number length = Length();
			Vector normalizedVector(m_Size);
			for (int i = 0; i < m_Size; i++) {
				normalizedVector.m_Data[i] = m_Data[i] / length;
			}
			return normalizedVector;
		}
	};

	enum Type
	{
		RotationX,
		RotationY,
		RotationZ,
		Transformation2D,
	};

	struct Matrix {
		int m_Rows;
		int m_Cols;
		Number* m_Data;

		Matrix(int rows, int cols) : m_Rows(rows), m_Cols(cols) {
			m_Data = new Number[m_Rows * m_Cols];
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					int i = x * m_Cols + y;
					m_Data[i] = (Number)0;
				}
			}
		}
		Matrix(int rows, int cols, Number value) : m_Rows(rows), m_Cols(cols) {
			m_Data = new Number[m_Rows * m_Cols];
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					int i = x * m_Cols + y;
					m_Data[i] = value;
				}
			}
		}
		Matrix(int rows, int cols, std::initializer_list<Number> values) : m_Rows(rows), m_Cols(cols) {
			assert(m_Rows * m_Cols == values.size());

			m_Data = new Number[m_Rows * m_Cols];
			int i = 0;
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					Number value = *(values.begin() + i);
					m_Data[x * m_Cols + y] = value;
					i++;
				}
			}
		}
		Matrix(int size) : m_Rows(size), m_Cols(size) {
			m_Data = new Number[m_Rows * m_Cols];
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					int i = x * m_Cols + y;
					if (x == y)
						m_Data[i] = (Number)1;
					else
						m_Data[i] = (Number)0;
				}
			}
		}
		Matrix(const Type& type, Number first = 0, Number second = 0, Number third = 0) {
			m_Rows = 3;
			m_Cols = 3;
			m_Data = new Number[m_Rows * m_Cols];
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					int i = x * m_Cols + y;
					m_Data[i] = (Number)0;
				}
			}

			Number angleRad = tg::degToRad(first);
			switch (type)
			{
			case RotationX:
				m_Data[0] = (Number)1;
				m_Data[4] = tg::cos(angleRad);
				m_Data[5] = -tg::sin(angleRad);
				m_Data[7] = tg::sin(angleRad);
				m_Data[8] = tg::cos(angleRad);
				break;
			case RotationY:
				m_Data[0] = tg::cos(angleRad);
				m_Data[2] = tg::sin(angleRad);
				m_Data[4] = (Number)1;
				m_Data[6] = -tg::sin(angleRad);
				m_Data[7] = tg::cos(angleRad);
				break;
			case RotationZ:
				m_Data[0] = tg::cos(angleRad);
				m_Data[1] = -tg::sin(angleRad);
				m_Data[3] = tg::sin(angleRad);
				m_Data[4] = tg::cos(angleRad);
				m_Data[8] = (Number)1;
				break;
			case Transformation2D:
				// Rotation
				m_Data[0] = tg::cos(angleRad);
				m_Data[1] = -tg::sin(angleRad);
				m_Data[3] = tg::sin(angleRad);
				m_Data[4] = tg::cos(angleRad);
				m_Data[8] = (Number)1;
				// Translation
				m_Data[2] = second;
				m_Data[5] = third;
				break;
			default:
				assert(0);
				break;
			}
		}
		Matrix(const Matrix& other) : m_Rows(other.m_Rows), m_Cols(other.m_Cols) {
			m_Data = new Number[m_Rows * m_Cols];
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					int i = x * m_Cols + y;
					m_Data[i] = other.m_Data[i];
				}
			}
		}

		~Matrix() {
			delete[] m_Data;
		}

		Matrix& operator=(const Matrix& other) {
			if (this != &other) {
				delete[] m_Data;
				m_Rows = other.m_Rows;
				m_Cols = other.m_Cols;
				m_Data = new Number[m_Rows * m_Cols];
				for (int x = 0; x < m_Rows; x++) {
					for (int y = 0; y < m_Cols; y++) {
						int i = x * m_Cols + y;
						m_Data[i] = other.m_Data[i];
					}
				}
			}
			return *this;
		}

		struct RowProxy {
			int m_Size;
			Number* row;
			RowProxy(int size, Number* rowData) : m_Size(size), row(rowData) {}

			Number& operator[](int y) {
				assert(y >= 0 && y < m_Size);
				return row[y];
			}
			const Number& operator[](int y) const {
				assert(y >= 0 && y < m_Size);
				return row[y];
			}
		};
		RowProxy operator[](int x) {
			assert(x >= 0 && x < m_Rows);
			return RowProxy(m_Cols, &m_Data[x * m_Cols]);
		}
		const RowProxy operator[](int x) const {
			assert(x >= 0 && x < m_Rows);
			return RowProxy(m_Cols, &m_Data[x * m_Cols]);
		}
		/*RowProxy& operator=(const RowProxy& other) {
			for (int col = 0; col < m_Rows; col++) {
				this->row[col] = other.row[col];
			}
			return *this;
		}*/

		Matrix operator+(const Matrix& other) {
			assert(m_Rows == other.m_Rows && m_Cols == other.m_Cols);

			Matrix result(m_Rows, m_Cols);
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					int i = x * m_Cols + y;
					result.m_Data[i] = m_Data[i] + other.m_Data[i];
				}
			}
			return result;
		}
		Matrix operator-(const Matrix& other) {
			assert(m_Rows == other.m_Rows && m_Cols == other.m_Cols);

			Matrix result(m_Rows, m_Cols);
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					int i = x * m_Cols + y;
					result.m_Data[i] = m_Data[i] - other.m_Data[i];
				}
			}
			return result;
		}
		Matrix operator*(const Matrix& other) {
			assert(m_Cols == other.m_Rows);

			int resultRows = m_Rows;
			int resultCols = other.m_Cols;
			Matrix result(resultRows, resultCols);

			for (int rx = 0; rx < resultRows; rx++) {
				for (int ry = 0; ry < resultCols; ry++) {
					int ri = rx * resultCols + ry;
					Number midResult = 0;

					for (int j = 0; j < m_Cols; j++) {
						Number value = m_Data[rx * m_Cols + j];
						Number otherValue = other.m_Data[j * other.m_Cols + ry];
						midResult += value * otherValue;
					}
					result.m_Data[ri] = midResult;
				}
			}
			return result;
		}

		friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
			os << "{ ";
			for (int x = 0; x < matrix.m_Rows; x++) {
				if (x != 0)
					os << "  ";
				for (int y = 0; y < matrix.m_Cols; y++) {
					int i = x * matrix.m_Cols + y;
					os << matrix.m_Data[i] << " ";
				}
				if(x != matrix.m_Rows-1)
					os << std::endl;
			}
			os << "}";
			return os;
		}

		Vector operator*(const Vector& vector) {
			assert(m_Cols == vector.m_Size);

			Vector result(m_Rows);
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					int i = x * m_Cols + y;
					result[x] += m_Data[i] * vector[y];
				}
			}
			return result;
		}

		Matrix operator*(const Number& value) {
			Matrix result(m_Rows, m_Cols);
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					int i = x * m_Cols + y;
					result.m_Data[i] = m_Data[i] * value;
				}
			}
			return result;
		}
		Matrix operator/(const Number& value) {
			Matrix result(m_Rows, m_Cols);
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					int i = x * m_Cols + y;
					result.m_Data[i] = m_Data[i] / value;
				}
			}
			return result;
		}

		bool IsOrthogonal() {
			if (m_Rows != m_Cols)
				return false;

			// Check for horizontal vectors
			for (int x = 0; x < m_Rows; x++) {
				Vector vector(m_Cols);
				for (int y = 0; y < m_Cols; y++) {
					int i = x * m_Cols + y;
					vector.m_Data[y] = m_Data[i];
				}
				if (!vector.IsOrthogonal())
					return false;
			}

			// Check for vertical vectors
			for (int y = 0; y < m_Cols; y++) {
				Vector vector(m_Rows);
				for (int x = 0; x < m_Rows; x++) {
					int i = x * m_Cols + y;
					vector.m_Data[x] = m_Data[i];
				}
				if (!vector.IsOrthogonal())
					return false;
			}
			return true;
		}
		Matrix Homogeneous() {
			assert(m_Rows == m_Cols);

			Matrix result(m_Rows + 1, m_Cols + 1);
			int i = 0;
			for (int x = 0; x < result.m_Rows; x++) {
				for (int y = 0; y < result.m_Cols; y++) {
					int j = x * result.m_Cols + y;
					if (y == m_Cols || x == m_Rows) {
						result.m_Data[j] = (Number)0;
					}
					else {
						result.m_Data[j] = m_Data[i];
						i++;
					}
				}
			}
			result.m_Data[result.m_Rows * result.m_Cols - 1] = (Number)1;
			return result;
		}
		Matrix SubMatrix(int row, int col) {
			assert(m_Rows == m_Cols);

			Matrix result(m_Rows - 1, m_Cols - 1);
			int i = 0;
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					if (x != row && y != col) {
						result.m_Data[i] = m_Data[x * m_Cols + y];
						i++;
					}
				}
			}
			return result;
		}
		Matrix Transpose() {
			Matrix result(m_Cols, m_Rows);
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					int i = x * m_Cols + y;
					int newI = y * m_Rows + x;
					result.m_Data[newI] = m_Data[i];
				}
			}
			return result;
		}
		Number Determinant() {
			assert(m_Rows == m_Cols);

			Number determinant = 0;
			if (m_Rows == 2) {
				determinant = m_Data[0] * m_Data[3] - m_Data[1] * m_Data[2];
			}
			else {
				for (int i = 0; i < m_Cols; i++) {
					Matrix mat = this->SubMatrix(0, i);
					if (i % 2 == 0) {
						determinant += m_Data[i] * mat.Determinant();
					}
					else {
						determinant -= m_Data[i] * mat.Determinant();
					}
				}
			}
			return determinant;
		}
		Matrix Inverse() {
			if (IsOrthogonal())
				return Transpose();

			assert(m_Rows == m_Cols);
			Number determiant = this->Determinant();
			assert(determiant != 0);

			Matrix transposedMatrix = this->Transpose();
			Matrix cofactorMatrix(transposedMatrix.m_Rows, transposedMatrix.m_Cols);

			for (int x = 0; x < cofactorMatrix.m_Rows; x++) {
				for (int y = 0; y < cofactorMatrix.m_Cols; y++) {
					Matrix subMatrix = transposedMatrix.SubMatrix(x, y);
					Number subMatrixDeterminant = subMatrix.Determinant();
					int i = x * cofactorMatrix.m_Cols + y;
					if (i % 2 == 0)
						cofactorMatrix.m_Data[i] = subMatrixDeterminant;
					else
						cofactorMatrix.m_Data[i] = -subMatrixDeterminant;
				}
			}
			Matrix inversedMatrix = cofactorMatrix * (1 / determiant);
			return inversedMatrix;
		}
	};

	namespace algo {

		Vector Thomas(const Vector& a, const Vector& b, const Vector& c, const Vector& z) {
			int m = b.m_Size;
			sm::Vector U(m);
			sm::Vector L(m-1);
			U[0] = b[0];

			for (int i = 0; i < m - 1; i++) {
				L[i] = a[i] / U[i];
				U[i + 1] = b[i + 1] - L[i] * c[i];
			}

			sm::Vector y = z;
			for (int i = 1; i < m; i++) {
				y[i] = z[i] - L[i - 1] * y[i - 1];
			}

			sm::Vector x = y;
			x[m - 1] = x[m - 1] / U[m - 1];
			for (int i = m-2; i > -1; i--) {
				x[i] = (y[i] - c[i] * x[i + 1]) / U[i];
			}
			return x;
		}
		int Factorial(int value) {
			int result = 1;
			if (value == 0)
				return 1;
			else if (value < 0)
				return -Factorial(-value);

			for (int i = 1; i <= value; i++) {
				result *= i;
			}
			return result;
		}
		Number Absolute(Number value) {
			if (value < 0)
				return -value;
			else
				return value;
		}
		Vector Absolute(Vector vector) {
			Vector result = vector;
			for (int i = 0; i < vector.m_Size; i++) {
				if (result[i] < 0)
					result[i] = -result[i];
			}
			return result;
		}
		Number Power(Number value, int exponent) {
			Number result = 1;
			for (int i = 1; i <= Absolute(exponent); i++) {
				result *= value;
			}
			if (exponent < 0)
				return 1 / result;
			else
				return result;
		}
		Number TaylorCosSin(int precision, bool prime, Number value) {
			Number result = (Number)0;
			int j = 0;
			for (int i = 0; i <= precision; i++) {
				if (!prime && i % 2 == 0)
					continue;
				if (prime && i % 2 == 1)
					continue;

				Number upperPart = Power(value, i);
				Number lowerPart = Factorial(i);
				Number midResult = upperPart / lowerPart;

				if (j % 2 == 0)
					result += midResult;
				else
					result -= midResult;
				j++;
			}
			return result;
		}
		Number NewtonMethod(Number value) {
			assert(value >= 0);

			if (value == 0)
				return 0;
			Number initialGuess = value / 2.0;
			Number newGuess;
			while (true) {
				newGuess = 0.5f * (initialGuess + (value / initialGuess));
				Number delta = newGuess - initialGuess;
				if (Absolute(delta) < NEWTON_PRECISION)
					break;
				initialGuess = newGuess;
			}
			return newGuess;
		}
		Number Sqrt(Number value) {
			return NewtonMethod(value);
		}
		Vector CrossProduct(const Vector& a, const Vector& b) {
			assert(a.m_Size == 3 && b.m_Size == 3);

			// Make matrix
			Matrix matrix(3, 3);
			for (int x = 0; x < matrix.m_Rows; x++) {
				for (int y = 0; y < matrix.m_Cols; y++) {
					if (x == 0) {
						if (y % 2 == 0)
							matrix[x][y] = (Number)1;
						else
							matrix[x][y] = -(Number)1;
					}
					else if (x == 1) {
						matrix[x][y] = a[y];
					}
					else if (x == 2) {
						matrix[x][y] = b[y];
					}
				}
			}
			
			Vector crossProduct(a.m_Size);
			for (int i = 0; i < crossProduct.m_Size; i++) {
				Matrix subMatrix = matrix.SubMatrix(0, i);
				if (i % 2 == 0) {
					crossProduct[i] = subMatrix.Determinant();
				}
				else {
					crossProduct[i] = -subMatrix.Determinant();
				}
			}
			return crossProduct;
		}
		Number DotProduct(const Vector& a, const Vector& b) {
			assert(a.m_Size == b.m_Size);

			Number dotProduct = 0;
			for (int i = 0; i < a.m_Size; i++) {
				dotProduct += a[i] * b[i];
			}
			return dotProduct;
		}
		Number Lerp(const Number& a, const Number& b, const Number& t) {
			Number result = a + t * (b - a);
			return result;
		}
		Vector Lerp(const Vector& a, const Vector& b, const Number& t) {
			assert(a.m_Size == b.m_Size);

			Vector result(a.m_Size);
			for (int i = 0; i < a.m_Size; i++) {
				result[i] = a[i] + t * (b[i] - a[i]);
			}
			return result;
		}
		Number InverseLerp(const Number& a, const Number& b, const Number& c) {
			if (a == b)
				return 0;
			return (c - a) / (b - a);
		}
		Vector InverseLerp(const Vector& a, const Vector& b, const Vector& c) {
			Vector result(a.m_Size);
			for (int i = 0; i < a.m_Size; i++) {
				result[i] = (c[i] - a[i]) / (b[i] - a[i]);
			}
			return result;
		}
		Number Max(const Number& a, const Number& b) {
			if (a > b)
				return a;
			else
				return b;
		}
		Number Max(const Vector& a) {
			Number max = MIN_NUMBER;
			for (int i = 0; i < a.m_Size; i++) {
				Number value = a[i];
				if (value > max) {
					max = value;
				}
			}
			return max;
		}
		Number Min(const Number& a, const Number& b) {
			if (a < b)
				return a;
			else
				return b;
		}
		Number Min(const Vector& a) {
			Number min = MAX_NUMBER;
			for (int i = 0; i < a.m_Size; i++) {
				Number value = a[i];
				if (value < min) {
					min = value;
				}
			}
			return min;
		}
		int Sign(const Number& value) {
			if (value >= 0)
				return 1;
			else
				return -1;
		}
		Number TaylorLog(Number x, int precision) {
			assert(x > 0 && precision > 0);

			if (x == 1)
				return 0;
			Number result = 0.0;
			Number term = (x - 1) / (x + 1);
			Number term_squared = term * term;
			Number numerator = term;
			Number denominator = 1;

			for (int i = 1; i <= precision; i++) {
				result += numerator / denominator;
				numerator *= term_squared;
				denominator += 2;
			}
			return 2 * result;
		}
		Number Log(Number x, Number base) {
			assert(x > 0 && base > 0 && base != 1);

			Number natural_log_x = TaylorLog(x, LOG_TAYLOR_PRECISION);
			Number natural_log_base = TaylorLog(base, LOG_TAYLOR_PRECISION);

			return natural_log_x / natural_log_base;
		}
		int Floor(const Number& value) {
			if (value == (int)value)
				return (int)value;
			else if (value > 0)
				return (int)value;
			else
				return (int)value - 1;
		}
		int Ceil(const Number& value) {
			if (value == (int)value)
				return (int)value;
			else if (value > 0)
				return (int)value + 1;
			else
				return (int)value;
		}
		int Round(const Number& value) {
			int floor = Floor(value);
			int ceil = Ceil(value);
			Number deltaUp = Absolute(ceil - value);
			Number deltaDown = Absolute(floor - value);
			if (deltaUp > deltaDown)
				return floor;
			else
				return ceil;
		}
		Number Clamp(const Number& value, const Number& min, const Number& max) {
			return Max(min, Min(value, max));
		}
		Number Clamp01(const Number& value) {
			return Max(0, Min(value, 1));
		}
		Vector CramersRule(Matrix& coefficientMatrix, const Vector& constants) {
			assert(coefficientMatrix.m_Rows == coefficientMatrix.m_Cols && coefficientMatrix.m_Cols == constants.m_Size);

			Vector solutions(coefficientMatrix.m_Rows);
			Number coeffMatDeterminant = coefficientMatrix.Determinant();
			for (int i = 0; i < coefficientMatrix.m_Cols; i++) {
				Matrix mat = coefficientMatrix;
				for (int x = 0; x < mat.m_Rows; x++) {
					mat[x][i] = constants[x];
				}
				Number matDeterminant = mat.Determinant();
				Number solution = matDeterminant / coeffMatDeterminant;
				solutions[i] = solution;
			}
			return solutions;
		}
		Vector GaussJordanElimination(Matrix& coefficientMatrix, const Vector& constants) {
			assert(coefficientMatrix.m_Rows == coefficientMatrix.m_Cols && coefficientMatrix.m_Cols == constants.m_Size);

			Matrix mat(coefficientMatrix.m_Rows, coefficientMatrix.m_Cols + 1);
			// Create matrix
			for (int x = 0; x < mat.m_Rows; x++) {
				for (int y = 0; y < mat.m_Cols; y++) {
					if (y == mat.m_Cols - 1)
						mat[x][y] = constants[x];
					else
						mat[x][y] = coefficientMatrix[x][y];
				}
			}

			// Solve bottom
			for (int x = 0; x < mat.m_Rows; x++) {
				for (int y = 0; y < x+1; y++) {
					Number value = mat[x][y];
					if (x == y) {
						Number divisor = value;
						// 1 -> divide row
						for (int yi = 0; yi < mat.m_Cols; yi++) {
							mat[x][yi] /= divisor;
						}
					}
					else {
						// 0
						Number multiplier = value;
						for (int yi = 0; yi < mat.m_Cols; yi++) {
							mat[x][yi] -= multiplier * mat[y][yi];
						}
					}
				}
			}
			// Solve top
			for (int y = 1; y < mat.m_Cols-1; y++) {
				for (int x = 0; x < y; x++) {
					Number value = mat[x][y];
					// 0
					Number multiplier = value;
					for (int yi = y; yi < mat.m_Cols; yi++) {
						mat[x][yi] -= multiplier * mat[y][yi];
					}
				}
			}
			
			Vector solutions(mat.m_Rows);
			for (int i = 0; i < mat.m_Rows; i++) {
				solutions[i] = mat[i][mat.m_Cols - 1];
			}
			return solutions;
		}
		Vector GaussSeidel(Matrix& coefficientMatrix, const Vector& constants) {
			assert(coefficientMatrix.m_Rows == coefficientMatrix.m_Cols && coefficientMatrix.m_Cols == constants.m_Size);

			Vector prevGuess(coefficientMatrix.m_Rows);
			Vector newGuess = prevGuess;
			int k = 0;
			while (k < 1000) {
				for (int x = 0; x < coefficientMatrix.m_Rows; x++) {
					double guess = 0;
					double divisor = coefficientMatrix[x][x];
					for (int y = 0; y < coefficientMatrix.m_Cols; y++) {
						if (y != x)
							guess -= Sign(divisor) * coefficientMatrix[x][y] * prevGuess[y];
					}
					guess += Sign(divisor) * constants[x];

					guess /= Absolute(divisor);
					newGuess[x] = guess;
				}
				if (Absolute(prevGuess - newGuess) < Vector(3, PRECISION)) {
					return newGuess;
				}
				prevGuess = newGuess;
				k++;
			}
			// Did not converge
			assert(0);
			return newGuess;
		}
	}
}

#define RED "\033[41m"
#define BLUE "\033[46m"
#define RESET "\033[0m"

void PlotData(const sm::Vector& data1, const sm::Vector& data2) {
	std::cout << "data size: " << data1.m_Size << std::endl;
	Number maxValue = sm::algo::Max(data1);
	Number minValue = sm::algo::Min(data2);
	std::cout << "min: " << minValue << std::endl;
	std::cout << "max: " << maxValue << std::endl;

	float delta = 1.0;
	int plotLength = 170;
	int plotHeight = 50;
	float deltaX = data1.m_Size / plotLength;
	float deltaY = (maxValue - minValue) / plotHeight;
	std::cout << "deltaX: " << deltaX << std::endl;
	std::cout << "deltaY: " << deltaY << std::endl;
	
	for (int y = maxValue; y > minValue; y-=deltaY) {
		std::cout << y << " ";
		for (int x = 0; x < data1.m_Size; x+=deltaX) {
			if (sm::algo::Absolute(data1[x] - y) < delta) {
				std::cout << RED << " " << RESET;
			}
			else if (sm::algo::Absolute(data2[x] - y) < delta) {
				std::cout << BLUE << " " << RESET;
			}
			else {
				std::cout << " ";
			}
		}
		std::cout << std::endl;
	}
}