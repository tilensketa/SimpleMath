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

	sm::Matrix<int> mat1(3, 3, { 1,2,3,0,1,4,5,6,0 });
	auto inv = mat1.Inverse();

	std::cout << mat1 << std::endl;
	std::cout << inv << std::endl;

	float d = sm::tg::sin(0.523599);
	std::cout << d << std::endl;

	float rad = PI / 6;
	float deg = sm::tg::radToDeg(rad);
	std::cout << deg << std::endl;

	double value = 2;
	std::cout << sm::algo::Power(value,2) << std::endl;
	std::cout << sm::algo::Power(value,-2) << std::endl;

	std::cout << sm::tg::sin(rad) << std::endl;
	std::cout << sm::tg::cos(rad) << std::endl;
	std::cout << sm::tg::tan(rad) << std::endl;
	std::cout << sm::tg::cot(rad) << std::endl;

	sm::Vector<int> vec({ 1,2,3 });
	auto hom = vec.Homogeneous();
	std::cout << vec << std::endl;
	std::cout << hom << std::endl;
	
	sm::Matrix<int> matrix(3,3, { 1,2,3,4,5,6,7,8,9 });
	auto homMatrix = matrix.Homogeneous();
	std::cout << matrix << std::endl;
	std::cout << homMatrix << std::endl;

	sm::Matrix<float> rotMatX(sm::RotationX, 90.0);
	sm::Matrix<float> rotMatY(sm::RotationY, 90.0);
	sm::Matrix<float> rotMatZ(sm::RotationZ, 90.0);
	std::cout << rotMatX << std::endl;
	std::cout << rotMatY << std::endl;
	std::cout << rotMatZ << std::endl;

	std::cout << "-------------------------" << std::endl;
	sm::Matrix<double> rotMat(sm::RotationZ, -30.0f);
	auto transRotMat = rotMat.Transpose();
	sm::Vector<double> point({ 1,3 });
	auto homoPoint = point.Homogeneous();
	auto newPoint = transRotMat * homoPoint;
	std::cout << newPoint << std::endl;

	std::cout << "-------------------------" << std::endl;
	sm::Matrix<double> transMat(sm::Translation2D, 1.0f, 5.0f);
	std::cout << transMat << std::endl;

	std::cout << "-------------------------" << std::endl;
	sm::Matrix<double> transfMat(sm::Transformation2D, 180.0f, 0.0f, -2.0f);
	std::cout << transfMat << std::endl;
	sm::Vector<double> pt({ 1,1 });
	std::cout << pt << std::endl;
	auto newPt = transfMat * pt.Homogeneous();
	std::cout << newPt << std::endl;

	std::cout << "-------------------------------" << std::endl;

	for (int i = 0; i < 101; i++) {
		std::cout << i << ": " << sm::algo::Sqrt((double)i) << std::endl;
	}
	sm::Vector<float> newVec({ 0.86602,0.5,0 });
	std::cout << newVec << std::endl;
	std::cout << newVec.IsOrthogonal() << std::endl;

	sm::Matrix<double> rtMat(sm::RotationZ, 30.0f);
	std::cout << rtMat << std::endl;
	std::cout << rtMat.IsOrthogonal() << std::endl;
	std::cout << rtMat.Inverse() << std::endl;

	std::cout << "-------------------------------------" << std::endl;
	sm::Vector<double> vec1({ 3,2,-1 });
	sm::Vector<double> vec2({ 3,0,5 });
	std::cout << vec1 << std::endl;
	std::cout << vec2 << std::endl;
	std::cout << sm::algo::CrossProduct(vec1, vec2) << std::endl;
	std::cout << vec1.Normalize() << std::endl;
	std::cout << sm::algo::DotProduct(vec1, vec2) << std::endl;
}