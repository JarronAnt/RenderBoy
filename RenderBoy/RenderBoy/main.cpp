/*
READ:
This is my first serious attempt at a raytracer so bare with the code.
Terrible coding practices are in place cause im lazy
using the eigen3 library for all linear algebra operations cause I cant 
be botherd to create my own vector class (NOTE: I added my own vector class anyways)
preformance isnt an object of concern in this raytracer so if it compiles slow....that's unfortuante  
this raytracer doesn't leverage the gpu at all its purely cpu based.
*/

//#include <Eigen/Dense>
#include "vec3.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

//using namespace Eigen;

class material;

//this class hold all the ray info 
class ray {
public:
	ray() {}
	ray(const vec3 &a, const vec3 &b) { A = a; B = b; }

	vec3 A;
	vec3 B;

	vec3 origin() const { return A; }
	vec3 direction() const { return B; }
	vec3 point_at_T(float t) const { return (A + t * B); }

	

};

//hittable object info
struct hitRecord {
	float t;
	vec3 p;
	vec3 normal;
	material *matPtr;
};

//base class for all hittable object
class hittable {
public:
	virtual bool hit(const ray &r, float tMin, float tMax, hitRecord &rec) const = 0;
};


//base class for all materials
class material {
public:
	virtual bool scatter(const ray &r, const hitRecord &hr, vec3 &attenuation, ray &scattered) const = 0;
};

//sphere class 
class sphere : public hittable {
public:
	sphere() {}
	sphere(vec3 cen, float r, material *m) :center(cen), radius(r), matPtr(m) {};
	virtual bool hit(const ray &r, float tMin, float tMax, hitRecord &rec) const;
	vec3 center;
	float radius;
	material *matPtr;
};
//hit function for the sphere
bool sphere::hit(const ray &r, float tMin, float tMax, hitRecord &rec) const{
	//get the vector from the ray to the center of the sphere 
	vec3 OC = r.origin() - center;
	//get the abc components of the qudratic equation 
	float a = dot(r.direction(), r.direction());
	float b = dot(OC, r.direction());
	float c = dot(OC,OC) - radius * radius;
	//calculate the discriminate if discrim > 0 we have a real soloution 
	// to the equation of the line hitting the sphere 
	float discrim = b * b -  a*c;

	if (discrim > 0) {
		float temp = (-b - sqrt(discrim)) / a;
		if (temp < tMax && temp > tMin) {
			rec.t = temp;
			rec.p = r.point_at_T(rec.t);
			rec.normal = (rec.p - center) / radius;
			rec.matPtr = matPtr;
			return true;
		}
		temp = (-b + sqrt(discrim)) / a;
		if (temp < tMax && temp > tMin) {
			rec.t = temp;
			rec.p = r.point_at_T(rec.t);
			rec.normal = (rec.p - center) / radius;
			rec.matPtr = matPtr;
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

//rand double generator for antialiasing 
double randomDouble() {
	return rand() / (RAND_MAX + 1.0);
}

//this function gives us a random point in the sphere for diffuse lighting
vec3 randInSphere() {
	vec3 p;
	do {

		p = 2.0*vec3(randomDouble(), randomDouble(), randomDouble()) - vec3(1, 1, 1);

	} while ((p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) >= 1.0);
	return p;
}


//this class handles lambertian (diffuse) material
class lambertian : public material {
public:
	lambertian(const vec3 &a) : albedo(a) {}
	virtual bool scatter(const ray &r, const hitRecord &hr, vec3 &attenuation, ray &scattered) const {
		vec3 target = hr.p + hr.normal + randInSphere();
		scattered = ray(hr.p, target - hr.p);
		attenuation = albedo;
		return true;
	}


	vec3 albedo;
};

//this is the function that handles reflections
vec3 reflect(const vec3 &v, const vec3 &n) {
	
	return v-2*dot(v,n)*n;
}


class metal : public material {
public:
	metal(const vec3& a) : albedo(a) {}
	virtual bool scatter(const ray& r_in, const hitRecord& rec,
		vec3& attenuation, ray& scattered) const {
		vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
		scattered = ray(rec.p, reflected);
		attenuation = albedo;
		return (dot(scattered.direction(), rec.normal) > 0);
	}
	vec3 albedo;
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
vec3 colour( const ray &r, hittable *world, int depth) {

	hitRecord hr;
	if (world->hit(r, 0.001, FLT_MAX, hr)) {
		ray scattered;		
		vec3 attenuation;
		if (depth < 50 && hr.matPtr->scatter(r, hr, attenuation, scattered)) {
			return attenuation * colour(scattered, world, depth + 1);
		}
		else {
			return vec3(0, 0, 0);
		}
		
	}
	else {
		vec3 unitVec = unit_vector(r.direction());
		float t = 0.5*(unitVec[1] + 1.0);
		return (1.0 - t)*vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
	}
	
}

//abstraction of the camera class 
class camera {
public:
	camera() { 
		lowerLeft = vec3(-2.0, -1.0, -1.0);
		horizontal = vec3(4.0, 0.0, 0.0);
		vertical = vec3(0.0, 2.0, 0.0);
		origin = vec3(0.0, 0.0, 0.0);
	}
	ray getRay(float u, float v) {
		return ray(origin, lowerLeft + u * horizontal + v * vertical - origin);
	}

	vec3 lowerLeft, horizontal, vertical, origin;
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
	hittable *list[4];
	//populate that list 
	list[0] = new sphere(vec3(0, 0, -1), 0.5, new lambertian(vec3(0.0, 0.5, 0.7)));
	list[1] = new sphere(vec3(0, -100.5, -1), 100, new lambertian(vec3(0.5, 0.5, 0.0)));
	list[2] = new sphere(vec3(1, 0, -1), 0.5, new metal(vec3(0.8, 0.6, 0.2)));
	list[3] = new sphere(vec3(-1, 0, -1), 0.5, new metal(vec3(0.8, 0.8, 0.8)));

	//create the new hit list 
	hittable *world = new hitList(list, 4);
	camera cam;
	//draw the ppm image 
	for (int j = ny - 1; j >= 0; j--) {
		for (int i = 0; i < nx; i++) {
			//set color to black  for each pixel initally 
			vec3 col(0, 0, 0);
			//this loop samples the colour in the pixel
			for (int  s = 0; s < ns; s++)
			{
				//get the u,v coords for the ray
				float u = float(i + randomDouble()) / float(nx);
				float v = float(j + randomDouble()) / float(ny);
				//get the ray based on above coords
				ray r = cam.getRay(u, v);
				//add all the colors in the pixel 
				col += colour(r, world,0);
			}
			//take the average of the pixel color
			col /= ns;
			//gamma correct the image
			col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
			//set the color 
			int ir = int(255.99*col[0]);
			int ig = int(255.99*col[1]);
			int ib = int(255.99*col[2]);
			image << ir << " " << ig << " " << ib << "\n";
		}
	}

}


