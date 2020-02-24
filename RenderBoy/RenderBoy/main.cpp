#include <Eigen/Dense>
#include <iostream>
#include <fstream>

using namespace Eigen;

int main() {

	int nx = 200;
	int ny = 100;

	//open the file 
	std::ofstream image;
	image.open("render.ppm");
	image << "P3\n" << nx << " " << ny << "\n255\n";
	//draw the ppm image 
	for (int j = ny - 1; j >= 0; j--) {
		for (int i = 0; i < nx; i++) {

			Vector3f col(float(i) / float(nx), float(j) / float(ny), 0.45);

			int ir = int(255.99*col[0]);
			int ig = int(255.99*col[1]);
			int ib = int(255.99*col[2]);
			image << ir << " " << ig << " " << ib << "\n";
		}
	}

}


//this class hold all the ray info 
class ray {
	ray() {}
	ray(const Vector3f &a, const Vector3f &b) {}

	Vector3f origin() const { return A; }
	Vector3f direction() const { return B; }
	Vector3f point_at_T(float t) const { return ( A + t * B); }

	Vector3f A;
	Vector3f B;

};