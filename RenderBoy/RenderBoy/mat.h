#pragma once

#include "hit.h"

//base class for all materials
class material {
public:
	virtual bool scatter(const ray &r, const hitRecord &hr, vec3 &attenuation, ray &scattered) const = 0;
};
