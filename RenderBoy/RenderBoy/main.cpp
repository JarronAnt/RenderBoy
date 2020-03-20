/*
READ:
This is my first serious attempt at a raytracer so bare with the code.
Terrible coding practices are in place cause im lazy
using the eigen3 library for all linear algebra operations cause I cant 
be botherd to create my own vector class (NOTE: I added my own vector class anyways)
preformance isnt an object of concern in this raytracer so if it compiles slow....that's unfortuante  
this raytracer doesn't leverage the gpu at all its purely cpu based.
*/

#include "vec3.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "perlin.h"

#define PI 3.14159265358979323846

//using namespace Eigen;

class material;

//this class hold all the ray info 
class ray {
public:
	ray() {}
	ray(const vec3 &a, const vec3 &b, float ti = 0.0) { A = a; B = b; _time = ti; }

	vec3 A;
	vec3 B;
	float _time;

	vec3 origin() const { return A; }
	vec3 direction() const { return B; }
	vec3 point_at_T(float t) const { return (A + t * B); }
	float time() const { return _time; }

	

};

//hittable object info
struct hitRecord {
	float t;
	float u;
	float v;
	vec3 p;
	vec3 normal;
	material *matPtr;
};

//some faster comparison functions
float ffmin(float a, float b) { return a < b ? a : b; }
float ffmax(float a, float b) { return a > b ? a : b; }

//axis aligned bounded box 
class aabb {
public:
	aabb() {}
	aabb(const vec3& a, const vec3& b) { _min = a; _max = b; }

	vec3 min() const { return _min; }
	vec3 max() const { return _max; }

	bool hit(const ray& r, float tmin, float tmax) const {
		for (int a = 0; a < 3; a++) {
			float t0 = ffmin((_min[a] - r.origin()[a]) / r.direction()[a],
				(_max[a] - r.origin()[a]) / r.direction()[a]);
			float t1 = ffmax((_min[a] - r.origin()[a]) / r.direction()[a],
				(_max[a] - r.origin()[a]) / r.direction()[a]);
			tmin = ffmax(t0, tmin);
			tmax = ffmin(t1, tmax);
			if (tmax <= tmin)
				return false;
		}
		return true;
	}

	vec3 _min;
	vec3 _max;
};

//base class for all hittable object
class hittable {
public:
	virtual bool hit(const ray &r, float tMin, float tMax, hitRecord &rec) const = 0;
	virtual bool bounding_box(float t0, float t1, aabb& box) const = 0;
};

 aabb surrounding_box(aabb box0, aabb box1) {
	 vec3 small(ffmin(box0.min().x(), box1.min().x()),
		 ffmin(box0.min().y(), box1.min().y()),
		 ffmin(box0.min().z(), box1.min().z()));
	 vec3 big(ffmax(box0.max().x(), box1.max().x()),
		 ffmax(box0.max().y(), box1.max().y()),
		 ffmax(box0.max().z(), box1.max().z()));
	 return aabb(small, big);
 }

//base class for all materials
class material {
public:
	virtual bool scatter(const ray &r, const hitRecord &hr, vec3 &attenuation, ray &scattered) const = 0;
	virtual vec3 emitted(float u, float v, const vec3 &p) const {
		return vec3(0, 0, 0);
	}
};

//sphere class 
class sphere : public hittable {
public:
	sphere() {}
	sphere(vec3 cen, float r, material *m) :center(cen), radius(r), matPtr(m) {};
	virtual bool hit(const ray &r, float tMin, float tMax, hitRecord &rec) const;
	virtual bool bounding_box(float t0, float t1, aabb& box) const;

	vec3 center;
	float radius;
	material *matPtr;
};

//build the bounding box for the sphere
bool sphere::bounding_box(float t0, float t1, aabb& box) const {
	box = aabb(center - vec3(radius, radius, radius),
		center + vec3(radius, radius, radius));
	return true;
}

//a moving sphere to test motion blur
class movingSphere : public hittable {
public:
	movingSphere() {}
	movingSphere(vec3 cen0, vec3 cen1, float t0, float t1, float r, material *m)
		: center0(cen0), center1(cen1), time0(t0), time1(t1), radius(r), mat_ptr(m)
	{};
	virtual bool hit(const ray& r, float tmin, float tmax, hitRecord& rec) const;
	virtual bool bounding_box(float t0, float t1, aabb& box) const;
	vec3 center(float time) const;
	vec3 center0, center1;
	float time0, time1;
	float radius;
	material *mat_ptr;
};
vec3 movingSphere::center(float time) const {
	return center0 + ((time - time0) / (time1 - time0))*(center1 - center0);
}
bool movingSphere::bounding_box(float t0, float t1, aabb& box) const {
	aabb box0(center(t0) - vec3(radius, radius, radius),
		center(t0) + vec3(radius, radius, radius));
	aabb box1(center(t1) - vec3(radius, radius, radius),
		center(t1) + vec3(radius, radius, radius));
	box = surrounding_box(box0, box1);
	return true;
}

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
bool movingSphere::hit(const ray& r, float t_min, float t_max, hitRecord& rec) const {
	vec3 oc = r.origin() - center(r.time());
	float a = dot(r.direction(), r.direction());
	float b = dot(oc, r.direction());
	float c = dot(oc, oc) - radius * radius;
	float discriminant = b * b - a * c;
	if (discriminant > 0) {
		float temp = (-b - sqrt(discriminant)) / a;
		if (temp < t_max && temp > t_min) {
			rec.t = temp;
			rec.p = r.point_at_T(rec.t);
			rec.normal = (rec.p - center(r.time())) / radius;
			rec.matPtr = mat_ptr;
			return true;
		}
		temp = (-b + sqrt(discriminant)) / a;
		if (temp < t_max && temp > t_min) {
			rec.t = temp;
			rec.p = r.point_at_T(rec.t);
			rec.normal = (rec.p - center(r.time())) / radius;
			rec.matPtr = mat_ptr;
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
	virtual bool bounding_box(float t0, float t1, aabb& box) const;
	//double pointer = list 
	hittable **list;
	int listSize;
};


//this function gives us a random point in the sphere for diffuse lighting
vec3 randInSphere() {
	vec3 p;
	do {

		p = 2.0*vec3(randomDouble(), randomDouble(), randomDouble()) - vec3(1, 1, 1);

	} while ((p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) >= 1.0);
	return p;
}

class texture {
public:
	virtual vec3 value(float u, float v, const vec3& p) const = 0;
};

class constant_texture : public texture {
public:
	constant_texture() {}
	constant_texture(vec3 c) : color(c) {}
	virtual vec3 value(float u, float v, const vec3& p) const {
		return color;
	}
	vec3 color;
};

class checker_texture : public texture {
public:
	checker_texture() {}
	checker_texture(texture *t0, texture *t1) : even(t0), odd(t1) {}
	virtual vec3 value(float u, float v, const vec3& p) const {
		float sines = sin(10 * p.x())*sin(10 * p.y())*sin(10 * p.z());
		if (sines < 0)
			return odd->value(u, v, p);
		else
			return even->value(u, v, p);
	}
	texture *odd;
	texture *even;
};

class diffuseLight : public material {
public:
	diffuseLight(texture *a) : emit(a) {}
	virtual bool scatter(const ray &r, const hitRecord &hr, vec3 &attenuation, ray &scattered) const { return false; }
	virtual vec3 emitted(float u, float v, const vec3 &p) const {
		return emit->value(u, v, p);
	}
	texture *emit;
};


//this class handles lambertian (diffuse) material
class lambertian : public material {
public:
	lambertian(texture *a) : albedo(a) {}
	virtual bool scatter(const ray &r, const hitRecord &hr, vec3 &attenuation, ray &scattered) const {
		vec3 target = hr.p + hr.normal + randInSphere();
		scattered = ray(hr.p, target - hr.p, r.time());
		attenuation = albedo->value(0,0,hr.p);
		return true;
	}


	texture *albedo;
};

//this is the function that handles reflections
vec3 reflect(const vec3 &v, const vec3 &n) {
	
	return v-2*dot(v,n)*n;
}

//this class handles metal material
class metal : public material {
public:
	metal(const vec3& a, float f) : albedo(a) {
		f < 1 ? fuzz = f : fuzz = 1;
	}
	virtual bool scatter(const ray& r_in, const hitRecord& rec,
		vec3& attenuation, ray& scattered) const {
		vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
		scattered = ray(rec.p, reflected + fuzz*randInSphere());
		attenuation = albedo;
		return (dot(scattered.direction(), rec.normal) > 0);
	}
	vec3 albedo;
	float fuzz;
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

bool hitList::bounding_box(float t0, float t1, aabb& box) const {
	if (listSize < 1) return false;
	aabb temp_box;
	bool first_true = list[0]->bounding_box(t0, t1, temp_box);
	if (!first_true)
		return false;
	else
		box = temp_box;
	for (int i = 1; i < listSize; i++) {
		if (list[i]->bounding_box(t0, t1, temp_box)) {
			box = surrounding_box(box, temp_box);
		}
		else
			return false;
	}
	return true;
}

//this function handles refraction of a ray based on snells law
bool refract(const vec3 &v, const vec3 &n, float niOverNt, vec3 &refracted) {
	vec3 uv = unit_vector(v);
	float dt = dot(uv, n);
	float discrim = 1.0 - niOverNt * niOverNt*(1 - dt * dt);
	if (discrim > 0) {
		refracted = niOverNt * (uv - n * dt) - n * sqrt(discrim);
		return true;
	}
	else {
		return false;
	}
}

//this function approximates reflectivity based on angle 
float schlick(float cos, float rfxIdx) {
	float r0 = (1 - rfxIdx) / (1 + rfxIdx);
	r0 = r0 * r0;
	return r0 + (1 - r0)*pow((1 - cos), 5);
}

//this function hold the code for dielectric material
class dielectric : public material {
public:
	dielectric(float ri) :rfxIdx(ri){}
	virtual bool scatter(const ray &r, const hitRecord &hr, vec3 &attenuation, ray &scattered) const {

		vec3 outNormal;
		vec3 reflected = reflect(r.direction(), hr.normal);
		float niOverNt;
		attenuation = vec3(1.0, 1.0 ,1.0);
		vec3 refracted;

		float reflectProb;
		float cos;

		if (dot(r.direction(), hr.normal) > 0) {
			outNormal = -hr.normal;
			niOverNt = rfxIdx;
			cos = rfxIdx * dot(r.direction(), hr.normal) / r.direction().length();
		}
		else {
			outNormal = hr.normal;
			niOverNt = 1.0 / rfxIdx;
			cos =  -dot(r.direction(), hr.normal) / r.direction().length();

		}

		if (refract(r.direction(), outNormal, niOverNt, refracted)) {
			reflectProb = schlick(cos, rfxIdx);
		}
		else {
			reflectProb = 1.0;
		}

		if (randomDouble() < reflectProb) {
			scattered = ray(hr.p, reflected);
		}
		else {
			scattered = ray(hr.p, refracted);
		}

		return true;
	}


	float rfxIdx;
};

//color in our scene 
vec3 colour( const ray &r, hittable *world, int depth) {

	hitRecord hr;
	if (world->hit(r, 0.001, FLT_MAX, hr)) {
		ray scattered;		
		vec3 attenuation;
		vec3 emitted = hr.matPtr->emitted(hr.u, hr.v, hr.p);
		if (depth < 50 && hr.matPtr->scatter(r, hr, attenuation, scattered)) {
			return emitted + attenuation * colour(scattered, world, depth + 1);
		}
		else {
			//return vec3(0, 0, 0);
			return emitted;
		}
		
	}
	else {
		/*vec3 unitVec = unit_vector(r.direction());
		float t = 0.5*(unitVec[1] + 1.0);
		return (1.0 - t)*vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);*/
		return vec3(0, 0, 0);
	}
	
}

//gets a random point in a disk (our camera lens)
vec3 randUnitInDisk() {
	vec3 p;
	do {
		p = 2.0*vec3(randomDouble(), randomDouble(), 0) - vec3(1, 1, 0);
	} while (dot(p, p) >= 1);

	return p;
}

//abstraction of the camera class 
class camera {
public:
	camera(vec3 lookFrom, vec3 lookAt, vec3 vup, float fov, float aspect, float apeture, float focusDist,
		float t0, float t1) {
		time0 = t0;
		time1 = t1;
		lensRadius = apeture / 2;
		float theta = fov * PI / 180;
		float halfHeight = tan(theta / 2);
		float halfWidth = aspect * halfHeight;
		origin = lookFrom;
		w = unit_vector(lookFrom - lookAt);
		u = unit_vector(cross(vup, w));
		v = cross(w, u);

		lowerLeft = origin - halfWidth * focusDist * u - halfHeight * focusDist * v - focusDist * w;
		horizontal = 2 * halfWidth*focusDist*u;
		vertical = 2 * halfHeight*focusDist*v;

	}
	ray getRay(float s, float t) {
		vec3 r = lensRadius * randUnitInDisk();
		vec3 offset = u * r.x() + v * r.y();
		float time = time0 + randomDouble()*(time1 - time0);
		return ray(origin + offset, lowerLeft + s * horizontal + t * vertical - origin - offset,time);
	}

	vec3 lowerLeft, horizontal, vertical, origin;
	float lensRadius;
	vec3 u, v, w;
	float time0, time1;

};

//node of our bvh structure
class bvh_node : public hittable {
public:
	bvh_node() {}
	bvh_node(hittable **l, int n, float time0, float time1);

	virtual bool hit(const ray& r, float tmin, float tmax, hitRecord& rec) const;
	virtual bool bounding_box(float t0, float t1, aabb& box) const;

	hittable *left;
	hittable *right;
	aabb box;
};

bool bvh_node::bounding_box(float t0, float t1, aabb& b) const {
	b = box;
	return true;
}

bool bvh_node::hit(const ray& r, float t_min, float t_max, hitRecord& rec) const {
	if (box.hit(r, t_min, t_max)) {
		hitRecord left_rec, right_rec;
		bool hit_left = left->hit(r, t_min, t_max, left_rec);
		bool hit_right = right->hit(r, t_min, t_max, right_rec);
		if (hit_left && hit_right) {
			if (left_rec.t < right_rec.t)
				rec = left_rec;
			else
				rec = right_rec;
			return true;
		}
		else if (hit_left) {
			rec = left_rec;
			return true;
		}
		else if (hit_right) {
			rec = right_rec;
			return true;
		}
		else
			return false;
	}
	else return false;
}

//these are axis dependent comparioson 
int box_x_compare(const void * a, const void * b) {
	aabb box_left, box_right;
	hittable *ah = *(hittable**)a;
	hittable *bh = *(hittable**)b;
	if (!ah->bounding_box(0, 0, box_left) || !bh->bounding_box(0, 0, box_right))
		std::cerr << "no bounding box in bvh_node constructor\n";
	if (box_left.min().x() - box_right.min().x() < 0.0)
		return -1;
	else
		return 1;
}
int box_y_compare(const void * a, const void * b)
{
	aabb box_left, box_right;
	hittable *ah = *(hittable**)a;
	hittable *bh = *(hittable**)b;
	if (!ah->bounding_box(0, 0, box_left) || !bh->bounding_box(0, 0, box_right))
		std::cerr << "no bounding box in bvh_node constructor\n";
	if (box_left.min().y() - box_right.min().y() < 0.0)
		return -1;
	else
		return 1;
}
int box_z_compare(const void * a, const void * b)
{
	aabb box_left, box_right;
	hittable *ah = *(hittable**)a;
	hittable *bh = *(hittable**)b;
	if (!ah->bounding_box(0, 0, box_left) || !bh->bounding_box(0, 0, box_right))
		std::cerr << "no bounding box in bvh_node constructor\n";
	if (box_left.min().z() - box_right.min().z() < 0.0)
		return -1;
	else
		return 1;
}

//this class constructs our bvh structure
bvh_node::bvh_node(hittable **l, int n, float time0, float time1) {
	int axis = int(3 * randomDouble());

	if (axis == 0)
		qsort(l, n, sizeof(hittable *), box_x_compare);
	else if (axis == 1)
		qsort(l, n, sizeof(hittable *), box_y_compare);
	else
		qsort(l, n, sizeof(hittable *), box_z_compare);

	if (n == 1) {
		left = right = l[0];
	}
	else if (n == 2) {
		left = l[0];
		right = l[1];
	}
	else {
		left = new bvh_node(l, n / 2, time0, time1);
		right = new bvh_node(l + n / 2, n - n / 2, time0, time1);
	}

	aabb box_left, box_right;

	if (!left->bounding_box(time0, time1, box_left) ||
		!right->bounding_box(time0, time1, box_right)) {

		std::cerr << "no bounding box in bvh_node constructor\n";
	}

	box = surrounding_box(box_left, box_right);
}


class noiseTex : public texture {
public:
	noiseTex(){}
	noiseTex(float sc) : scale(sc) {}
	virtual vec3 value(float u, float v, const vec3 &p) const {
		//return vec3(1, 1, 1) * noise.noise(p * scale);
		return vec3(1, 1, 1)*0.5*(1 + sin(scale*p.x() + 5 * noise.turb(scale*p)));
	}
	perlin noise;
	float scale;
};

class xyRect : public hittable {
public:
	xyRect() {}
	xyRect(float _x0, float _x1, float _y0, float _y1, float _k, material *mat)
		: x0(_x0), x1(_x1), y0(_y0), y1(_y1), k(_k), mp(mat) {};
	virtual bool hit(const ray& r, float t0, float t1, hitRecord& hr) const;
	virtual bool bounding_box(float t0, float t1, aabb& box) const {
		box = aabb(vec3(x0, y0, k - 0.0001), vec3(x1, y1, k + 0.0001));
		return true;
	}
	material *mp;
	float x0, x1, y0, y1, k;
};
class xzRect : public hittable {
public:
	xzRect() {}
	xzRect(float _x0, float _x1, float _z0, float _z1, float _k, material *mat)
		: x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k), mp(mat) {};
	virtual bool hit(const ray& r, float t0, float t1, hitRecord& rec) const;
	virtual bool bounding_box(float t0, float t1, aabb& box) const {
		box = aabb(vec3(x0, k - 0.0001, z0), vec3(x1, k + 0.0001, z1));
		return true;
	}
	material *mp;
	float x0, x1, z0, z1, k;
};
class yzRect : public hittable {
public:
	yzRect() {}
	yzRect(float _y0, float _y1, float _z0, float _z1, float _k, material *mat)
		: y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), mp(mat) {};
	virtual bool hit(const ray& r, float t0, float t1, hitRecord& hr) const;
	virtual bool bounding_box(float t0, float t1, aabb& box) const {
		box = aabb(vec3(k - 0.0001, y0, z0), vec3(k + 0.0001, y1, z1));
		return true;
	}
	material  *mp;
	float y0, y1, z0, z1, k;
};



bool xyRect::hit(const ray& r, float t0, float t1, hitRecord& hr) const {
	float t = (k - r.origin().z()) / r.direction().z();
	if (t < t0 || t > t1)
		return false;
	float x = r.origin().x() + t * r.direction().x();
	float y = r.origin().y() + t * r.direction().y();
	if (x < x0 || x > x1 || y < y0 || y > y1)
		return false;
	hr.u = (x - x0) / (x1 - x0);
	hr.v = (y - y0) / (y1 - y0);
	hr.t = t;
	hr.matPtr = mp;
	hr.p = r.point_at_T(t);
	hr.normal = vec3(0, 0, 1);
	return true;
}
bool xzRect::hit(const ray& r, float t0, float t1, hitRecord& hr) const {
	float t = (k - r.origin().y()) / r.direction().y();
	if (t < t0 || t > t1)
		return false;
	float x = r.origin().x() + t * r.direction().x();
	float z = r.origin().z() + t * r.direction().z();
	if (x < x0 || x > x1 || z < z0 || z > z1)
		return false;
	hr.u = (x - x0) / (x1 - x0);
	hr.v = (z - z0) / (z1 - z0);
	hr.t = t;
	hr.matPtr = mp;
	hr.p = r.point_at_T(t);
	hr.normal = vec3(0, 1, 0);
	return true;
}
bool yzRect::hit(const ray& r, float t0, float t1, hitRecord& hr) const {
	float t = (k - r.origin().x()) / r.direction().x();
	if (t < t0 || t > t1)
		return false;
	float y = r.origin().y() + t * r.direction().y();
	float z = r.origin().z() + t * r.direction().z();
	if (y < y0 || y > y1 || z < z0 || z > z1)
		return false;
	hr.u = (y - y0) / (y1 - y0);
	hr.v = (z - z0) / (z1 - z0);
	hr.t = t;
	hr.matPtr = mp;
	hr.p = r.point_at_T(t);
	hr.normal = vec3(1, 0, 0);
	return true;
}


class flipNorms : public hittable {
public:
	flipNorms(hittable *p): ptr(p){}
	virtual bool hit(const ray &r,float tMin, float tMax, hitRecord &hr) const {
		if (ptr->hit(r, tMin, tMax, hr)) {
			hr.normal = -hr.normal;
			return true;
		}
		else {
			return false;
		}
	}

	virtual bool bounding_box(float t0, float t1, aabb &box) const {
		return ptr->bounding_box(t0, t1, box);
	}
	hittable *ptr;
};


//----------------------------------------------Scenes-------------------------------------------------------//

hittable *randScene() {
	int n = 50000;
	hittable **list = new hittable*[n + 1];

	texture *checker = new checker_texture(
		new constant_texture(vec3(0.2, 0.3, 0.1)),
		new constant_texture(vec3(0.9, 0.9, 0.9))
	);

	list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(checker));
	int i = 1;
	for (int a = -10; a < 10; a++) {
		for (int b = -10; b < 10; b++) {
			float choose_mat = randomDouble();
			vec3 center(a + 0.9*randomDouble(), 0.2, b + 0.9*randomDouble());
			if ((center - vec3(4, 0.2, 0)).length() > 0.9) {
				if (choose_mat < 0.8) {  // diffuse
					list[i++] = new movingSphere(
						center,
						center + vec3(0, 0.5*randomDouble(), 0),
						0.0, 1.0, 0.2,
						new lambertian(new constant_texture(
							vec3(randomDouble()*randomDouble(),
								randomDouble()*randomDouble(),
								randomDouble()*randomDouble())
						))
					);
				}
				else if (choose_mat < 0.95) { // metal
					list[i++] = new sphere(
						center, 0.2,
						new metal(
							vec3(0.5*(1 + randomDouble()),
								0.5*(1 + randomDouble()),
								0.5*(1 + randomDouble())),
							0.5*randomDouble()
						)
					);
				}
				else {  // glass
					list[i++] = new sphere(center, 0.2, new dielectric(1.5));
				}
			}
		}
	}

	list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));
	list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new lambertian(new constant_texture(vec3(0.4, 0.2, 0.1))));
	list[i++] = new sphere(vec3(4, 1, 0), 1.0, new lambertian(new constant_texture(vec3(0.7, 0.6, 0.5))));

	return new hitList(list, i);
}
hittable *twoSpheres() {
	texture *checker = new checker_texture(
		new constant_texture(vec3(0.2, 0.3, 0.1)),
		new constant_texture(vec3(0.9, 0.9, 0.9)));
	int n = 50;
	hittable **list = new hittable*[n + 1];
	list[0] = new sphere(vec3(0, -10, 0), 10, new lambertian(checker));
	list[1] = new sphere(vec3(0, 10, 0), 10, new lambertian(checker));
	return new hitList(list, 2);
}
hittable *twoPerlin() {
	texture *pertext = new noiseTex(1);
	hittable **list = new hittable*[2];
	list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(pertext));
	list[1] = new sphere(vec3(0, 2, 0), 2, new lambertian(pertext));
	return new hitList(list, 2);
}
hittable *basicLight() {

	texture *pertext = new noiseTex(4);
	hittable **list = new hittable*[4];
	list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(pertext));
	list[1] = new sphere(vec3(0, 2, 0), 2, new lambertian(pertext));
	list[2] = new sphere(vec3(0, 7, 0), 2,
		new diffuseLight(new constant_texture(vec3(4, 4, 4))));
	list[3] = new xyRect(3, 5, 1, 3, -2,
		new diffuseLight(new constant_texture(vec3(4, 4, 4))));

	return new hitList(list, 4);
}
hittable *cornellBox() {
	hittable **list = new hittable*[6];
	int i = 0;
	material *red = new lambertian(new constant_texture(vec3(0.65, 0.05, 0.05)));
	material *white = new lambertian(new constant_texture(vec3(0.73, 0.73, 0.73)));
	material *green = new lambertian(new constant_texture(vec3(0.12, 0.45, 0.15)));
	material *light = new diffuseLight(new constant_texture(vec3(15, 15, 15)));

	list[i++] = new flipNorms(new yzRect(0, 555, 0, 555, 555, green));
	list[i++] = new yzRect(0, 555, 0, 555, 0, red);
	list[i++] = new xzRect(213, 343, 227, 332, 554, light);
	list[i++] = new flipNorms(new xzRect(0, 555, 0, 555, 555, white));
	list[i++] = new xzRect(0, 555, 0, 555, 0, white);
	list[i++] = new flipNorms(new xzRect(0, 555, 0, 555, 555, white));

	return new hitList(list, i);
}
//------------------------------------------------------------------------------------------------------------//
int main() {
	//screen x,y and sample sizes in terms of pixels
	int nx = 800;
	int ny = 600;
	//the higher the ns value the better the antialiasing but slows down the program 
	int ns = 50;

	//open the file 
	std::ofstream image;
	image.open("render.ppm");
	image << "P3\n" << nx << " " << ny << "\n255\n";

	//create a list of hittable objects 
/*	hittable *list[5];
	//populate that list 
	list[0] = new sphere(vec3(0, 0, -1), 0.5, new lambertian(vec3(0.0, 0.5, 0.7)));
	list[1] = new sphere(vec3(0, -100.5, -1), 100, new lambertian(vec3(0.5, 0.5, 0.0)));
	list[2] = new sphere(vec3(1, 0, -1), 0.5, new metal(vec3(0.8, 0.6, 0.2), 0.7));
	list[3] = new sphere(vec3(-1, 0, -1), 0.5, new dielectric(1.5));
	list[4] = new sphere(vec3(-1, 0, -1), -0.45, new dielectric(1.5));*/


	//camera info
	vec3 lookFrom(278, 278, -800);
	vec3 lookAt(278, 278, 0);
	float distToFocus = 10.0;
	float appeture = 0.0;
	float fov = 40.0;
	camera cam(lookFrom, lookAt, vec3(0, 1, 0), fov	, float(nx) / float(ny),appeture,distToFocus,0.0,1.0);

	//create the new hit list 
	hittable *world = cornellBox();
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


