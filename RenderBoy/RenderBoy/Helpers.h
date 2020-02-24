//This is a header that will contain all helper functions 

#pragma once
#include <cstdlib>

float epsilon = 0.0000001;

inline bool equal(float a, float b) {

	if (abs(a) - abs(b) <= epsilon)
		true;
	else
		false;
}