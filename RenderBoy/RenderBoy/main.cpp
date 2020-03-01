/*
READ:
This is my first serious attempt at a raytracer so bare with the code.
Terrible coding practices are in place cause im lazy
using the eigen3 library for all linear algebra operations cause I cant 
be botherd to create my own vector class 
preformance isnt an object of concern in this raytracer so if it compiles slow....that's unfortuante  
this raytracer doesn't leverage the gpu at all its purely cpu based.
*/

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <cstdlib>

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

//hittable object info
struct hitRecord {
	float t;
	Vector3f p;
	Vector3f normal;
};
//base class for all hittable object
class hittable {
public:
	virtual bool hit(const ray &r, float tMin, float tMax, hitRecord &rec) const = 0;
};

//sphere class 
class sphere : public hittable {
public:
	sphere() {}
	sphere(Vector3f cen, float r) :center(cen), radius(r) {};
	virtual bool hit(const ray &r, float tMin, float tMax, hitRecord &rec) const;
	Vector3f center;
	float radius;
};
//hit function for the sphere
bool sphere::hit(const ray &r, float tMin, float tMax, hitRecord &rec) const{
	//get the vector from the ray to the center of the sphere 
	Vector3f OC = r.origin() - center;
	//get the abc components of the qudratic equation 
	float a = r.direction().dot(r.direction());
	float b = OC.dot(r.direction());
	float c = OC.dot(OC) - radius * radius;
	//calculate the discriminate if discrim > 0 we have a real soloution 
	// to the equation of the line hitting the sphere 
	float discrim = b * b -  a*c;

	if (discrim > 0) {
		float temp = (-b - sqrt(discrim)) / a;
		if (temp < tMax && temp > tMin) {
			rec.t = temp;
			rec.p = r.point_at_T(rec.t);
			rec.normal = (rec.p - center) / radius;
			return true;
		}
		temp = (-b + sqrt(discrim)) / a;
		if (temp < tMax && temp > tMin) {
			rec.t = temp;
			rec.p = r.point_at_T(rec.t);
			rec.normal = (rec.p - center) / radius;
			return true;
		}
	}
	return false;
}

//a list of all hittable objects
class hitList : public hittable {
public:
	hitList() {}
	//this is an object that contains a list of hittable objects 
	hitList(hittable **l, int n) { list = l; listSize = n; }
	//hit function for objects in the list
	virtual bool hit(
		const ray& r, float tmin, float tmax, hitRecord& rec) const;
	//double pointer = list 
	hittable **list;
	int listSize;
};

//hit function for the whole list
bool hitList::hit(const ray& r, float t_min, float t_max,hitRecord& rec) const {

	hitRecord temp_rec;
	bool hitAnything = false;
	double closest_so_far = t_max;
	//loop through the list and see if we hit anything 
	for (int i = 0; i < listSize; i++) {
		if (list[i]->hit(r, t_min, closest_so_far, temp_rec)) {
			hitAnything = true;
			closest_so_far = temp_rec.t;
			rec = temp_rec;
		}
	}
	return hitAnything;
}

//color in our scene 
Vector3f colour( const ray &r, hittable *world) {

	hitRecord hr;
	if (world->hit(r, 0.0, FLT_MAX, hr)) {
		return 0.5*Vector3f(hr.normal[0] + 1, hr.normal[1] + 1, hr.normal[2] + 1);
	}
	else {
		Vector3f unitVec = r.direction().normalized();
		float t = 0.5*(unitVec[1] + 1.0);
		return (1.0 - t)*Vector3f(1.0, 1.0, 1.0) + t * Vector3f(0.5, 0.7, 1.0);
	}
	
}

//rand double generator for antialiasing 
double randomDouble() {
	return rand() / (RAND_MAX + 1.0);
}

//abstraction of the camera class 
class camera {
public:
	camera() { 
		lowerLeft = Vector3f(-2.0, -1.0, -1.0);
		horizontal = Vector3f(4.0, 0.0, 0.0);
		vertical = Vector3f(0.0, 2.0, 0.0);
		origin = Vector3f(0.0, 0.0, 0.0);
	}
	ray getRay(float u, float v) {
		return ray(origin, lowerLeft + u * horizontal + v * vertical - origin);
	}

	Vector3f lowerLeft, horizontal, vertical, origin;
};


int main() {
	//screen x,y and sample sizes in terms of pixels
	int nx = 200;
	int ny = 100;
	//the higher the ns value the better the antialiasing but slows down the program 
	int ns = 50;

	//open the file 
	std::ofstream image;
	image.open("render.ppm");
	image << "P3\n" << nx << " " << ny << "\n255\n";

	//create a list of hittable objects 
	hittable *list[2];
	//populate that list 
	list[0] = new sphere(Vector3f(0, 0, -1), 0.5);
	list[1] = new sphere(Vector3f(0, -100.5, -1), 100);
	//create the new hit list 
	hittable *world = new hitList(list, 2);
	camera cam;
	//draw the ppm image 
	for (int j = ny - 1; j >= 0; j--) {
		for (int i = 0; i < nx; i++) {
			//set color to black  for each pixel initally 
			Vector3f col(0, 0, 0);
			//this loop samples the colour in the pixel
			for (int  s = 0; s < ns; s++)
			{
				//get the u,v coords for the ray
				float u = float(i + randomDouble()) / float(nx);
				float v = float(j + randomDouble()) / float(ny);
				//get the ray based on above coords
				ray r = cam.getRay(u, v);
				//add all the colors in the pixel 
				col += colour(r, world);
			}
			//take the average of the pixel color
			col /= ns;
			//set the color 
			int ir = int(255.99*col[0]);
			int ig = int(255.99*col[1]);
			int ib = int(255.99*col[2]);
			image << ir << " " << ig << " " << ib << "\n";
		}
	}

}


