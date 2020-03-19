#pragma once

#include "ray.h"
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
	virtual bool boundingBox(float t0, float t1, aabb &box) const = 0;
};