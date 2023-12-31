#include "SimpleMath.h"

#include "HeatTransfer.h"

int main() {

	HeatTransfer heat(Example::one);

#pragma region Usage
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
#pragma endregion
#pragma region Trigonometry
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
#pragma endregion
#pragma region Vector
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
	sm::Vector crossProduct = sm::algo::CrossProduct(vector1, vector2);
	double dotProduct = sm::algo::DotProduct(vector1, vector2);
#pragma endregion
#pragma region Matrix
	// Constructors
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
#pragma endregion

	sm::Vector val1({ 0,5,3 });
	sm::Vector val2({ 2,1,2 });
	for (int i = 0; i < 11; i++) {
		Number t = i / 10.0;
		std::cout << t << ": " << sm::algo::Lerp(val1, val2, t) << std::endl;
	}
	

	double nu = 0.4;
	std::cout << "-----------Number: " << nu << std::endl;
	std::cout << "Floor: " << sm::algo::Floor(nu) << std::endl;
	std::cout << "Ceil: " << sm::algo::Ceil(nu) << std::endl;
	std::cout << "Round: " << sm::algo::Round(nu) << std::endl;

	std::cout << sm::algo::Clamp(nu, 0, 5) << std::endl;
	std::cout << sm::algo::Clamp(nu, -1, 1) << std::endl;
	std::cout << sm::algo::Clamp01(nu) << std::endl;

	std::cout << "-----------" << std::endl;
	sm::Matrix coefMat(3, 3, { 3,-1,1,2,6,2,1,2,4 });
	sm::Vector constants({ 2,5,1 });
	std::cout << "Cramers Rule: " << sm::algo::CramersRule(coefMat, constants) << std::endl;
	std::cout << "Gauss Jordan Elimination: " << sm::algo::GaussJordanElimination(coefMat, constants) << std::endl;
	std::cout << "Gauss Seidel: " << sm::algo::GaussSeidel(coefMat, constants) << std::endl;

	sm::Vector v1({ 1,2,3 });
	sm::Vector v2({ 0,1,2 });
	std::cout << (v1 > v2) << std::endl;
	std::cout << sm::algo::Absolute(v1) << std::endl;
}