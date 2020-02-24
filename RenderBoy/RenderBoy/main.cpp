#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;

int main() {

	char a;
	Vector4f v1 = Vector4f(3, -2, 5, 1);
	Vector4f v2 = Vector4f(-2, 3, 1, 0);

	Vector4f finalVec = v1 + v2;

	std::cout << finalVec << " ";

	std::cin >> a;

}