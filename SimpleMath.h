#pragma once

#include <iostream>
#include <initializer_list>
#include <cassert>

#define PI 3.1415926536

#define PRECISION 0.001
#define TRIGONOMETRIC_PRECISION 11
#define NEWTON_PRECISION 0.0001

namespace sm {
	// Function declarations
	namespace algo {
		template <typename T>
		T TaylorSeries(int precision, bool prime, T value);
		template <typename T>
		T Absolute(T value);
		template <typename T>
		T Power(T value, int exponent);
		template <typename T>
		T Sqrt(T value);
	}

	namespace tg {

		// Angle conversion
		template <typename T>
		T degToRad(T angleDeg) {
			T result = angleDeg * (PI / 180);
			return result;
		}
		template <typename T>
		T radToDeg(T angleRad) {
			T result = angleRad * (180 / PI);
			return result;
		}

		// Trigonometric functions
		template <typename T>
		T sin(T angleRad) {
			return algo::TaylorSeries(TRIGONOMETRIC_PRECISION, false, angleRad);
		}
		template <typename T>
		T cos(T angleRad) {
			return algo::TaylorSeries(TRIGONOMETRIC_PRECISION, true, angleRad);
		}
		template <typename T>
		T tan(T angleRad) {
			T result = sin(angleRad) / cos(angleRad);
			return result;
		}
		template <typename T>
		T cot(T angleRad) {
			T result = cos(angleRad) / sin(angleRad);
			return result;
		}
	}

	template <typename T>
	struct Vector {
		int m_Size;
		T* m_Data;

		Vector(int size) : m_Size(size) {
			m_Data = new T[m_Size];
			for (int i = 0; i < m_Size; i++) {
				m_Data[i] = (T)0;
			}
		}
		Vector(int size, T value) : m_Size(size) {
			m_Data = new T[m_Size];
			for (int i = 0; i < m_Size; i++) {
				m_Data[i] = value;
			}
		}
		Vector(std::initializer_list<T> values) : m_Size(values.size()) {
			assert(values.size() == m_Size);
			m_Data = new T[m_Size];
			size_t i = 0;
			for (const T& value : values) {
				m_Data[i++] = value;
			}
		}
		Vector(const Vector& other) : m_Size(other.m_Size) {
			m_Data = new T[m_Size];
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
				m_Data = new T[m_Size];
				for (int i = 0; i < m_Size; i++) {
					m_Data[i] = other.m_Data[i];
				}
			}
			return *this;
		}

		T& operator[](int index){
			assert(index >= 0 && index < m_Size);
			return m_Data[index];
		}
		const T& operator[] (int index) const {
			assert(index >= 0 && index < m_Size);
			return m_Data[index];
		}
		
		Vector<T> operator+(const Vector<T>& other) {
			assert(m_Size == other.m_Size);

			Vector<T> result(other.m_Size);
			for (int i = 0; i < m_Size; i++) {
				result[i] = this->m_Data[i] + other[i];
			}
			return result;
		}
		Vector<T> operator-(const Vector<T>& other) {
			assert(m_Size == other.m_Size);

			Vector<T> result(other.m_Size);
			for (int i = 0; i < m_Size; i++) {
				result[i] = this->m_Data[i] - other[i];
			}
			return result;
		}

		friend std::ostream& operator<<(std::ostream& os, const Vector<T>& vec) {
			os << "{ ";
			for (int i = 0; i < vec.m_Size; i++) {
				os << vec[i] << " ";
			}
			os << "}";
			return os;
		}

		Vector<T> operator*(const T& value) {
			Vector<T> result(m_Size);
			for (int i = 0; i < m_Size; i++) {
				result[i] = m_Data[i] * value;
			}
			return result;
		}
		Vector<T> operator/(const T& value) {
			Vector<T> result(m_Size);
			for (int i = 0; i < m_Size; i++) {
				result[i] = m_Data[i] / value;
			}
			return result;
		}

		Vector<T> Homogeneous() {
			// Homogenous vector is vector that has one additional row of 1
			Vector<T> result(m_Size + 1);
			for (int i = 0; i < m_Size; i++) {
				result.m_Data[i] = m_Data[i];
			}
			result.m_Data[m_Size] = (T)1;
			return result;
		}
		double Length() {
			// Length or magnitude of vector : Sqrt(x^2 + y^2 + ...)
			double variable = 0;
			for (int i = 0; i < m_Size; i++) {
				double value = m_Data[i];
				variable += algo::Power(value, 2);
			}
			return algo::Sqrt(variable);
		}
		bool IsOrthogonal() {
			// Vector is orthogonal if vector length is 1
			if (algo::Absolute(Length() - 1.0) < PRECISION)
				return true;
			return false;
		}
		Vector<T> Normalize() {
			double length = Length();
			Vector<T> normalizedVector(m_Size);
			for (int i = 0; i < m_Size; i++) {
				normalizedVector.m_Data[i] = m_Data[i] / length;
			}
			return normalizedVector;
		}
	};

	enum Type
	{
		Identity = 0,
		RotationX,
		RotationY,
		RotationZ,
		Translation2D,
		Transformation2D,
	};

	template <typename T>
	struct Matrix {
		int m_Rows;
		int m_Cols;
		T* m_Data;

		Matrix(int rows, int cols) : m_Rows(rows), m_Cols(cols) {
			m_Data = new T[m_Rows * m_Cols];
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					int i = x * m_Cols + y;
					m_Data[i] = (T)0;
				}
			}
		}
		Matrix(int rows, int cols, T value) : m_Rows(rows), m_Cols(cols) {
			m_Data = new T[m_Rows * m_Cols];
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					int i = x * m_Cols + y;
					m_Data[i] = value;
				}
			}
		}
		Matrix(int rows, int cols, std::initializer_list<T> values) : m_Rows(rows), m_Cols(cols) {
			assert(m_Rows * m_Cols == values.size());

			m_Data = new T[m_Rows * m_Cols];
			int i = 0;
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					T value = *(values.begin() + i);
					m_Data[x * m_Cols + y] = value;
					i++;
				}
			}
		}
		// Identity matrix
		Matrix(int size, const Type& type) : m_Rows(size), m_Cols(size) {
			m_Data = new T[m_Rows * m_Cols];
			switch (type)
			{
				case Identity:
					assert(m_Rows == m_Cols);
					for (int x = 0; x < m_Rows; x++) {
						for (int y = 0; y < m_Cols; y++) {
							int i = x * m_Cols + y;
							if (x == y)
								m_Data[i] = (T)1;
							else
								m_Data[i] = (T)0;
						}
					}
			}
		}
		// Rotation matrix X, Y, Z
		Matrix(const Type& type, float angleDeg): Matrix(3,3) {
			T angleRad = tg::degToRad(angleDeg);
			switch (type)
			{
			case RotationX:
				m_Data[0] = (T)1;
				m_Data[4] = tg::cos(angleRad);
				m_Data[5] = -tg::sin(angleRad);
				m_Data[7] = tg::sin(angleRad);
				m_Data[8] = tg::cos(angleRad);
				break;

			case RotationY:
				m_Data[0] = tg::cos(angleRad);
				m_Data[2] = tg::sin(angleRad);
				m_Data[4] = (T)1;
				m_Data[6] = -tg::sin(angleRad);
				m_Data[7] = tg::cos(angleRad);
				break;

			case RotationZ:
				m_Data[0] = tg::cos(angleRad);
				m_Data[1] = -tg::sin(angleRad);
				m_Data[3] = tg::sin(angleRad);
				m_Data[4] = tg::cos(angleRad);
				m_Data[8] = (T)1;
				break;
			}
		}
		// Translation matrix
		Matrix(const Type& type, float dx, float dy) : Matrix(3,sm::Identity) {
			switch (type)
			{
			case Translation2D:
				m_Data[2] = dx;
				m_Data[5] = dy;
				break;
			default:
				assert(0);
				break;
			}
		}
		// Transformation matrix (rotation, translation)
		Matrix(const Type& type, float angleDeg, float dx, float dy) : Matrix(sm::RotationZ, angleDeg) {
			switch (type)
			{
			case Transformation2D:
				m_Data[2] = dx;
				m_Data[5] = dy;
				break;
			default:
				assert(0);
				break;
			}
		}
		Matrix(const Matrix& other) : m_Rows(other.m_Rows), m_Cols(other.m_Cols) {
			m_Data = new T[m_Rows * m_Cols];
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
				m_Data = new T[m_Rows * m_Cols];
				for (int x = 0; x < m_Rows; x++) {
					for (int y = 0; y < m_Cols; y++) {
						int i = x * m_Cols + y;
						m_Data[i] = other[i];
					}
				}
			}
			return *this;
		}

		struct RowProxy {
			int m_Size;
			T* row;
			RowProxy(int size, T* rowData) : m_Size(size), row(rowData) {}

			T& operator[](int y) {
				assert(y >= 0 && y < m_Size);
				return row[y];
			}
			const T& operator[](int y) const {
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
		RowProxy& operator=(const RowProxy& other) {
			for (int col = 0; col < m_Rows; col++) {
				this->row[col] = other.row[col];
			}
			return *this;
		}

		Matrix<T> operator+(const Matrix<T>& other) {
			assert(m_Rows == other.m_Rows && m_Cols == other.m_Cols);

			Matrix<T> result(m_Rows, m_Cols);
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					int i = x * m_Cols + y;
					result.m_Data[i] = m_Data[i] + other.m_Data[i];
				}
			}
			return result;
		}
		Matrix<T> operator-(const Matrix<T>& other) {
			assert(m_Rows == other.m_Rows && m_Cols == other.m_Cols);

			Matrix<T> result(m_Rows, m_Cols);
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					int i = x * m_Cols + y;
					result.m_Data[i] = m_Data[i] - other.m_Data[i];
				}
			}
			return result;
		}
		Matrix<T> operator*(const Matrix<T>& other) {
			assert(m_Cols == other.m_Rows);

			int resultRows = m_Rows;
			int resultCols = other.m_Cols;
			Matrix<T> result(resultRows, resultCols);

			for (int rx = 0; rx < resultRows; rx++) {
				for (int ry = 0; ry < resultCols; ry++) {
					int ri = rx * resultCols + ry;
					T midResult = 0;

					for (int j = 0; j < m_Cols; j++) {
						T value = m_Data[rx * m_Cols + j];
						T otherValue = other.m_Data[j * other.m_Cols + ry];
						midResult += value * otherValue;
					}
					result.m_Data[ri] = midResult;
				}
			}
			return result;
		}

		friend std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix) {
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

		Vector<T> operator*(const Vector<T>& vector) {
			assert(m_Cols == vector.m_Size);

			Vector<T> result(m_Rows);
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					int i = x * m_Cols + y;
					result[x] += m_Data[i] * vector[y];
				}
			}
			return result;
		}

		Matrix<T> operator*(const T& value) {
			Matrix<T> result(m_Rows, m_Cols);
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					int i = x * m_Cols + y;
					result.m_Data[i] = m_Data[i] * value;
				}
			}
			return result;
		}
		Matrix<T> operator/(const T& value) {
			Matrix<T> result(m_Rows, m_Cols);
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
				Vector<T> vector(m_Cols);
				for (int y = 0; y < m_Cols; y++) {
					int i = x * m_Cols + y;
					vector.m_Data[y] = m_Data[i];
				}
				if (!vector.IsOrthogonal())
					return false;
			}

			// Check for vertical vectors
			for (int y = 0; y < m_Cols; y++) {
				Vector<T> vector(m_Rows);
				for (int x = 0; x < m_Rows; x++) {
					int i = x * m_Cols + y;
					vector.m_Data[x] = m_Data[i];
				}
				if (!vector.IsOrthogonal())
					return false;
			}
			return true;
		}
		Matrix<T> Homogeneous() {
			assert(m_Rows == m_Cols);

			Matrix<T> result(m_Rows + 1, m_Cols + 1);
			int i = 0;
			for (int x = 0; x < result.m_Rows; x++) {
				for (int y = 0; y < result.m_Cols; y++) {
					int j = x * result.m_Cols + y;
					if (y == m_Cols || x == m_Rows) {
						result.m_Data[j] = (T)0;
					}
					else {
						result.m_Data[j] = m_Data[i];
						i++;
					}
				}
			}
			result.m_Data[result.m_Rows * result.m_Cols - 1] = (T)1;
			return result;
		}
		Matrix<T> SubMatrix(int row, int col) {
			assert(m_Rows == m_Cols);

			Matrix<T> result(m_Rows - 1, m_Cols - 1);
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
		Matrix<T> Transpose() {
			Matrix<T> result(m_Cols, m_Rows);
			for (int x = 0; x < m_Rows; x++) {
				for (int y = 0; y < m_Cols; y++) {
					int i = x * m_Cols + y;
					int newI = y * m_Rows + x;
					result.m_Data[newI] = m_Data[i];
				}
			}
			return result;
		}
		double Determinant() {
			assert(m_Rows == m_Cols);

			double determinant = 0;
			if (m_Rows == 2) {
				determinant = m_Data[0] * m_Data[3] - m_Data[1] * m_Data[2];
			}
			else {
				for (int i = 0; i < m_Cols; i++) {
					Matrix<T> mat = this->SubMatrix(0, i);
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
		Matrix<T> Inverse() {
			if (IsOrthogonal())
				return Transpose();

			assert(m_Rows == m_Cols);
			int determiant = this->Determinant();
			assert(determiant != 0);

			Matrix<T> transposedMatrix = this->Transpose();
			Matrix<T> cofactorMatrix(transposedMatrix.m_Rows, transposedMatrix.m_Cols);

			for (int x = 0; x < cofactorMatrix.m_Rows; x++) {
				for (int y = 0; y < cofactorMatrix.m_Cols; y++) {
					Matrix<T> subMatrix = transposedMatrix.SubMatrix(x, y);
					int subMatrixDeterminant = subMatrix.Determinant();
					int i = x * cofactorMatrix.m_Cols + y;
					if (i % 2 == 0)
						cofactorMatrix.m_Data[i] = subMatrixDeterminant;
					else
						cofactorMatrix.m_Data[i] = -subMatrixDeterminant;
				}
			}
			Matrix<T> inversedMatrix = cofactorMatrix * (1 / determiant);
			return inversedMatrix;
		}
	};

	namespace algo {

		template<typename T>
		Vector<T> Thomas(const Vector<T>& a, const Vector<T>& b, const Vector<T>& c, const Vector<T>& z) {
			int m = b.m_Size;
			sm::Vector<double> U(m);
			sm::Vector<double> L(m-1);
			U[0] = b[0];

			for (int i = 0; i < m - 1; i++) {
				L[i] = a[i] / U[i];
				U[i + 1] = b[i + 1] - L[i] * c[i];
			}

			sm::Vector<double> y = z; //TODO size m
			for (int i = 1; i < m; i++) {
				y[i] = z[i] - L[i - 1] * y[i - 1];
			}

			sm::Vector<double> x = y;
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
		
		template <typename T>
		T Absolute(T value) {
			if (value < 0)
				return -value;
			else
				return value;
		}

		template <typename T>
		T Power(T value, int exponent) {
			T result = 1;
			for (int i = 1; i <= Absolute(exponent); i++)
			{
				result *= value;
			}
			if (exponent < 0)
				return 1 / result;
			else
				return result;
		}

		template <typename T>
		T TaylorSeries(int precision, bool prime, T value) {
			T result = (T)0;
			int j = 0;
			for (int i = 0; i <= precision; i ++) {
				if (!prime && i % 2 == 0)
					continue;
				if (prime && i % 2 == 1)
					continue;

				T upperPart = Power(value, i);
				T lowerPart = Factorial(i);
				
				T midResult = upperPart / lowerPart;

				if (j % 2 == 0)
					result += midResult;
				else
					result -= midResult;
				j++;
			}
			return result;
		}

		template <typename T>
		T NewtonMethod(T value) {
			assert(value >= 0);

			if (value == 0)
				return 0;
			T initialGuess = value / 2.0;
			T newGuess;
			while (true) {
				newGuess = 0.5f * (initialGuess + (value / initialGuess));
				T delta = newGuess - initialGuess;
				if (Absolute(delta) < NEWTON_PRECISION)
					break;
				initialGuess = newGuess;
			}
			return newGuess;
		}

		template <typename T>
		T Sqrt(T value) {
			return NewtonMethod(value);
		}

		template <typename T>
		Vector<T> CrossProduct(const Vector<T>& a, const Vector<T>& b) {
			assert(a.m_Size == 3 && b.m_Size == 3);

			// Make matrix
			Matrix<T> matrix(3, 3);
			for (int x = 0; x < matrix.m_Rows; x++) {
				for (int y = 0; y < matrix.m_Cols; y++) {
					int i = x * matrix.m_Cols + y;
					if (x == 0) {
						if (y % 2 == 0)
							matrix.m_Data[i] = (T)1;
						else
							matrix.m_Data[i] = -(T)1;
					}
					else if (x == 1) {
						matrix.m_Data[i] = a.m_Data[y];
					}
					else if (x == 2) {
						matrix.m_Data[i] = b.m_Data[y];
					}
				}
			}
			
			Vector<T> crossProduct(a.m_Size);
			for (int i = 0; i < crossProduct.m_Size; i++) {
				Matrix<T> subMatrix = matrix.SubMatrix(0, i);
				if (i % 2 == 0) {
					crossProduct.m_Data[i] = subMatrix.Determinant();
				}
				else {
					crossProduct.m_Data[i] = -subMatrix.Determinant();
				}
			}
			return crossProduct;
		}

		template <typename T>
		T DotProduct(const Vector<T>& a, const Vector<T>& b) {
			assert(a.m_Size == b.m_Size);

			T dotProduct = 0;
			for (int i = 0; i < a.m_Size; i++) {
				dotProduct += a.m_Data[i] * b.m_Data[i];
			}
			return dotProduct;
		}
	}
}