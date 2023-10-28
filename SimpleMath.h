#pragma once

#include <iostream>
#include <initializer_list>
#include <cassert>

namespace sm {

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
	};

	enum Type
	{
		Identity = 0,
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
		Matrix(int rows, int cols, const Type& type) : m_Rows(rows), m_Cols(cols) {
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
			return RowProxy(m_Cols, &m_Data[x]);
		}
		const RowProxy operator[](int x) const {
			assert(x >= 0 && x < m_Rows);
			return RowProxy(m_Cols, &m_Data[x]);
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
		int Determinant() {
			assert(m_Rows == m_Cols);

			int determinant = 0;
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
	}
}