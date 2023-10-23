#include "SimpleMath.h"

int main() {

	sm::Vector<int, 2> vector({ 1,2 });
	sm::Matrix<int, 2> matrix({ 1,2,3,4 });

	auto result = matrix * vector;
	std::cout << result << std::endl;
}