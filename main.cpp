#include "SimpleMath.h"

double ro = 1200.0; // kg / m3
double cp = 1500.0; // J / kgK
double k = 0.3; // W / mK
double L = 0.06; // m
float Tzac = 80.0 + 273.0; // K
float Tzrak = 20.0 + 273.0; // K
int h = 100; // W / m2K

int dt = 1; // s
const int N = 6; // število vozlišè
int t = 20000; // konèni èas v s
double dx = L / 10;
double diff = k / (ro * cp);
double Fo = diff * dt / (dx* dx);
double Bi = h * dx / k;

int main() {

	sm::Vector<double> a(N-1);
	sm::Vector<double> b(N);
	sm::Vector<double> c(N-1);
	sm::Vector<double> z(N);

	sm::Vector<double> T(N, Tzac);

	// define a
	for (int i = 0; i < a.m_Size; i++) {
		a[i] = -Fo;
	}
	a[a.m_Size - 1] = -2 * Fo;

	// define c
	for (int i = 0; i < c.m_Size; i++) {
		c[i] = -Fo;
	}
	c[0] = -2 * Fo;

	// define b
	for (int i = 0; i < b.m_Size; i++) {
		b[i] = 1 + 2 * Fo;
	}
	b[b.m_Size - 1] = 1 + 2 * Fo + 2 * Fo * Bi;

	int p = t / dt;
	sm::Vector<double> temp_zg(p);
	sm::Vector<double> temp_sp(p);
	for (int i = 0; i < p; i++) {
		for (int j = 0; j < N; j++) {
			z[j] = T[j];
		}
		z[z.m_Size - 1] = T[T.m_Size - 1] + 2 * Fo * Bi * Tzrak;

		auto Tnew = sm::algo::Thomas(a, b, c, z);
		T = Tnew;
		temp_zg[i] = T[T.m_Size - 1] - 273;
		temp_sp[i] = T[0] - 273;
	}
	//std::cout << temp_zg << std::endl;
	
#define RESET "\033[0m"
	int delta = 1;
	float deltaX = (Tzac - Tzrak) / 40;
	float deltaY = (float)p / 170;
	std::cout << "Tzac----------------------------------" << std::endl;
	for (int x = Tzac-273; x > Tzrak-273; x-=deltaX) {
		int i = (Tzac - 273 - x) * deltaX;
		std::cout << temp_zg[i] << "\t";
		for (int y = 0; y < p; y+=deltaY) {
			if (temp_zg[y] > x - delta && temp_zg[y] < x + delta)
				std::cout << "\033[31m"<<"x"<<RESET;
			else if (temp_sp[y] > x - delta && temp_sp[y] < x + delta)
				std::cout << "\033[36m"<<"+"<<RESET;
			else
				std::cout << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "Tzrak----------------------------------" << std::endl;

	sm::Matrix<int> mat1(2, 3, {1,2,3,4,5,6});
	sm::Matrix<int> mat3 = mat1.Transpose();
	std::cout << mat1 << std::endl;
	std::cout << mat3 << std::endl;
	sm::Matrix<float> mat4(3, 3, sm::Identity);
	std::cout << mat4 << std::endl;
	std::cout << mat4.Determinant() << std::endl;
}