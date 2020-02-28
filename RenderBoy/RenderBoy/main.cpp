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
bool hitList::hit(const ray& r, float t_min, float t_max,
	hitRecord& rec) const {

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

	//create a list of hittable objects 
	hittable *list[2];
	//populate that list 
	list[0] = new sphere(Vector3f(0, 0, -1), 0.5);
	list[1] = new sphere(Vector3f(0, -100.5, -1), 100);
	//create the new hit list 
	hittable *world = new hitList(list, 2);

	//draw the ppm image 
	for (int j = ny - 1; j >= 0; j--) {
		for (int i = 0; i < nx; i++) {

			float u = float(i) / float(nx);
			float v = float(j) / float(ny);

			//create a ray to each pixel 
			ray r(origin, lowerLeft + u * horizontal + v * vertical);

			Vector3f p = r.point_at_T(2.0);

			//calculate a color for that pixel 
			Vector3f col = colour(r,world);
			//set the color 
			int ir = int(255.99*col[0]);
			int ig = int(255.99*col[1]);
			int ib = int(255.99*col[2]);
			image << ir << " " << ig << " " << ib << "\n";
		}
	}

}


