#pragma once

#include <iostream>
#include <initializer_list>
#include <cassert>

namespace sm {

	template <typename T, size_t Size>
	struct Vector {
		size_t m_Size = Size;
		T m_Data[Size];

		Vector() {
			for (int i = 0; i < m_Size; i++) {
				m_Data[i] = (T)0;
			}
		}
		Vector(T value) {
			for (int i = 0; i < m_Size; i++) {
				m_Data[i] = value;
			}
		}
		Vector(std::initializer_list<T> values) {
			assert(values.size() == m_Size);

			size_t i = 0;
			for (const T& value : values) {
				m_Data[i++] = value;
			}
		}

		T& operator[](int index){
			assert(index >= 0 && index < m_Size);
			return m_Data[index];
		}
		const T& operator[] (int index) const {
			assert(index >= 0 && index < m_Size);
			return m_Data[index];
		}
		
		Vector<T, Size> operator+(const Vector<T, Size>& other) {
			assert(m_Size == other.m_Size);

			Vector<T, Size> result;
			for (int i = 0; i < Size; i++) {
				result[i] = this->m_Data[i] + other[i];
			}
			return result;
		}
		Vector<T, Size> operator-(const Vector<T, Size>& other) {
			assert(m_Size == other.m_Size);

			Vector<T, Size> result;
			for (int i = 0; i < Size; i++) {
				result[i] = this->m_Data[i] - other[i];
			}
			return result;
		}
		Vector<T, Size> operator*(const Vector<T, Size>& other) {
			assert(m_Size == other.m_Size);

			Vector<T, Size> result;
			for (int i = 0; i < Size; i++) {
				result[i] = this->m_Data[i] * other[i];
			}
			return result;
		}
		Vector<T, Size> operator/(const Vector<T, Size>& other) {
			assert(m_Size == other.m_Size);

			Vector<T, Size> result;
			for (int i = 0; i < Size; i++) {
				result[i] = this->m_Data[i] / other[i];
			}
			return result;
		}

		friend std::ostream& operator<<(std::ostream& os, const Vector<T, Size>& vec) {
			os << "{ ";
			for (int i = 0; i < Size; i++) {
				os << vec[i] << " ";
			}
			os << "}";
			return os;
		}

		Vector<T, Size> operator*(const T& value) {
			Vector<T, Size> result;
			for (int i = 0; i < Size; i++) {
				result[i] = m_Data[i] * value;
			}
			return result;
		}
		Vector<T, Size> operator/(const T& value) {
			Vector<T, Size> result;
			for (int i = 0; i < Size; i++) {
				result[i] = m_Data[i] / value;
			}
			return result;
		}
	};

	template <typename T, size_t Size>
	struct Matrix {
		size_t m_Size = Size;
		T m_Data[Size][Size];

		Matrix() {
			for (int x = 0; x < Size; x++) {
				for (int y = 0; y < Size; y++) {
					m_Data[x][y] = (T)0;
				}
			}
		}
		Matrix(T value) {
			for (int x = 0; x < Size; x++) {
				for (int y = 0; y < Size; y++) {
					m_Data[x][y] = value;
				}
			}
		}
		Matrix(std::initializer_list<T> values) {
			assert(values.size() != Size);
			int i = 0;
			for (int x = 0; x < Size; x++) {
				for (int y = 0; y < Size; y++) {
					m_Data[x][y] = values.begin()[i];
					i++;
				}
			}
		}

		struct RowProxy {
			size_t m_Size;
			T* row;
			RowProxy(T* rowData, size_t size) : row(rowData), m_Size(size) {}

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
			assert(x >= 0 && x < m_Size);
			return RowProxy(m_Data[x], m_Size);
		}
		const RowProxy operator[](int x) const {
			assert(x >= 0 && x < m_Size);
			return RowProxy(m_Data[x], m_Size);
		}

		Matrix<T, Size> operator+(const Matrix<T, Size>& other) {
			assert(Size == other.m_Size);

			Matrix<T, Size> result;
			for (int x = 0; x < Size; x++) {
				for (int y = 0; y < Size; y++) {
					result[x][y] = this->m_Data[x][y] + other.m_Data[x][y];
				}
			}
			return result;
		}
		Matrix<T, Size> operator-(const Matrix<T, Size>& other) {
			assert(Size == other.m_Size);

			Matrix<T, Size> result;
			for (int x = 0; x < Size; x++) {
				for (int y = 0; y < Size; y++) {
					result[x][y] = this->m_Data[x][y] - other.m_Data[x][y];
				}
			}
			return result;
		}
		Matrix<T, Size> operator*(const Matrix<T, Size>& other) {
			assert(Size == other.m_Size);

			Matrix<T, Size> result;
			for (int x = 0; x < Size; x++) {
				for (int y = 0; y < Size; y++) {
					result[x][y] = this->m_Data[x][y] * other.m_Data[x][y];
				}
			}
			return result;
		}
		Matrix<T, Size> operator/(const Matrix<T, Size>& other) {
			assert(Size == other.m_Size);

			Matrix<T, Size> result;
			for (int x = 0; x < Size; x++) {
				for (int y = 0; y < Size; y++) {
					result[x][y] = this->m_Data[x][y] / other.m_Data[x][y];
				}
			}
			return result;
		}

		friend std::ostream& operator<<(std::ostream& os, const Matrix<T, Size>& matrix) {
			os << "{ ";
			for (int x = 0; x < Size; x++) {
				if (x != 0)
					os << "  ";
				for (int y = 0; y < Size; y++) {
					os << matrix.m_Data[x][y] << " ";
				}
				if(x != Size-1)
					os << std::endl;
			}
			os << "}";
			return os;
		}

		Vector<T, Size> operator*(const Vector<T, Size>& vector) {
			assert(Size == vector.m_Size);

			Vector<T, Size> result;
			for (int x = 0; x < Size; x++) {
				for (int y = 0; y < Size; y++) {
					result[x] += m_Data[x][y] * vector[y];
				}
			}
			return result;
		}
		Vector<T, Size> operator/(const Vector<T, Size>& vector) {
			assert(Size == vector.m_Size);

			Vector<T, Size> result;
			for (int x = 0; x < Size; x++) {
				for (int y = 0; y < Size; y++) {
					result[x] += m_Data[x][y] / vector[y];
				}
			}
			return result;
		}

		Matrix<T, Size> operator*(const T& value) {
			Matrix<T, Size> result;
			for (int x = 0; x < Size; x++) {
				for (int y = 0; y < Size; y++) {
					result[x][y] = m_Data[x][y] * value;
				}
			}
			return result;
		}
		Matrix<T, Size> operator/(const T& value) {
			Matrix<T, Size> result;
			for (int x = 0; x < Size; x++) {
				for (int y = 0; y < Size; y++) {
					result[x][y] = m_Data[x][y] / value;
				}
			}
			return result;
		}
	};
}