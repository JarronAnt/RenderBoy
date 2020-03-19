#pragma once
#include "vec3.h"

inline float trilinear_interp(float c[2][2][2], float u, float v, float w) {
	float accum = 0;
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			for (int k = 0; k < 2; k++)
				accum += (i*u + (1 - i)*(1 - u))*
				(j*v + (1 - j)*(1 - v))*
				(k*w + (1 - k)*(1 - w))*c[i][j][k];

	return accum;
}


class perlin {
public:

	float noise(const vec3 &p) const {
		float u = p.x() - floor(p.x());
		float v = p.y() - floor(p.y());
		float w = p.z() - floor(p.z());
		int i = floor(p.x());
		int j = floor(p.y());
		int k = floor(p.z());
		float c[2][2][2];
		for (int di = 0; di < 2; di++)
			for (int dj = 0; dj < 2; dj++)
				for (int dk = 0; dk < 2; dk++)
					c[di][dj][dk] = ranFloat[
						permX[(i + di) & 255] ^
							permY[(j + dj) & 255] ^
							permZ[(k + dk) & 255]
					];
		return trilinear_interp(c, u, v, w);
	}
	static float *ranFloat;
	static int	*permX;
	static int *permY;
	static int *permZ;
};

static float* perlinGen() {
	float *p = new float[256];
	for (int i = 0; i < 256; ++i)
		p[i] = randomDouble();
	return p;
}

void permute(int *p, int n) {
	for (int i = n - 1; i > 0; i--) {
		int target = int(randomDouble()*(i + 1));
		int tmp = p[i];
		p[i] = p[target];
		p[target] = tmp;
	}
	return;
}

static int* perlin_Perm() {
	int * p = new int[256];
	for (int i = 0; i < 256; i++)
		p[i] = i;
	permute(p, 256);
	return p;
}

float *perlin::ranFloat = perlinGen();
int *perlin::permX = perlin_Perm();
int *perlin::permY = perlin_Perm();
int *perlin::permZ = perlin_Perm();