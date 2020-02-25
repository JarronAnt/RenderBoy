#include <Eigen/Dense>
#include <iostream>
#include <fstream>

using namespace Eigen;

//this class hold all the ray info 
class ray {
public:
	ray(const Vector3f &a, const Vector3f &b) { A = a; B = b; }

	Vector3f A;
	Vector3f B;

	Vector3f origin() const { return A; }
	Vector3f direction() const { return B; }
	Vector3f point_at_T(float t) const { return (A + t * B); }

	

};

float hitSphere(const Vector3f &center, float radius, const ray &r) {

	//get the vector from the ray to the center of the sphere 
	Vector3f OC = r.origin() - center;
	//get the abc components of the qudratic equation 
	float a = r.direction().dot(r.direction());
	float b = 2.0 * OC.dot(r.direction());
	float c = OC.dot(OC) - radius * radius;
	//calculate the discriminate if discrim > 0 we have a real soloution 
	// to the equation of the line hitting the sphere 
	float discrim = b * b - 4 * a*c;

	//return the acutal hit point 
	if (discrim < 0) {
		return (-1.0);
	}
	else {
		return (-b - sqrt(discrim)) / (2.0*a);
	}

}
Vector3f colour( const ray &r) {

	//if the ray hit the sphere at the pixel color it blue
	float t = hitSphere(Vector3f(0, 0, -1), 0.5, r);

	if (t > 0.0) {
		//compute the normal at point t 
		Vector3f N = (r.point_at_T(t) - Vector3f(0, 0, -1)).normalized();
		//generate a colour based on the normal 
		return 0.5*Vector3f(N[0] + 1, N[1] + 1, N[2] + 1);
	}
		

	Vector3f unitDirection = r.direction().normalized();
	//this turns the y vector from -1 to 1 into 0 to 1 
	t = 0.5*(unitDirection[1] + 1.0);
	//this is a linear interprotation (lerp)
	//formula is blendvalue = (1-t)*startColor + t*endColor
	return (1.0 - t)*Vector3f(1.0, 1.0, 1.0) + t * Vector3f(0.5, 0.2, 1.0);
	
}

int main() {

	int nx = 200;
	int ny = 100;

	//open the file 
	std::ofstream image;
	image.open("render.ppm");
	image << "P3\n" << nx << " " << ny << "\n255\n";

	//set the image size info
	Vector3f lowerLeft(-2.0, -1.0, -1.0);
	Vector3f horizontal(4.0, 0.0, 0.0);
	Vector3f vertical(0.0, 2.0, 0.0);
	Vector3f origin(0.0, 0.0, 0.0);

	//draw the ppm image 
	for (int j = ny - 1; j >= 0; j--) {
		for (int i = 0; i < nx; i++) {

			float u = float(i) / float(nx);
			float v = float(j) / float(ny);
			//create a ray to each pixel 
			ray r(origin, lowerLeft + u * horizontal + v * vertical);
			//calculate a color for that pixel 
			Vector3f col = colour(r);
			//set the color 
			int ir = int(255.99*col[0]);
			int ig = int(255.99*col[1]);
			int ib = int(255.99*col[2]);
			image << ir << " " << ig << " " << ib << "\n";
		}
	}

}


