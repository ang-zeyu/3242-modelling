#ifdef _WIN32
#include <Windows.h>
#include "GL\glut.h"
#define M_PI 3.141592654
#elif __APPLE__
#include <OpenGL/gl.h>
#include <GLUT/GLUT.h>
#endif

#include "math.h"
#include <string>
#include <unordered_map>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <queue>
#include <iomanip>
using namespace std;

/*
 Some miscellaneous things
 */

struct pairHash
{
	size_t operator()(pair<int, int> const& p) const noexcept
	{
		size_t h1 = hash<int>{}(p.first);
		size_t h2 = hash<int>{}(p.second);
		return h1 ^ (h2 << 1); // adapted from official c++ docs
	}
};

inline void crossProduct(double result[3], double v1[3], double v2[3])
{
	result[0] = v1[1] * v2[2] - v2[1] * v1[2];
	result[1] = v2[0] * v1[2] - v1[0] * v2[2];
	result[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

inline void subtractV(double result[3], double v1[3], double v2[3])
{
	result[0] = v1[0] - v2[0];
	result[1] = v1[1] - v2[1];
	result[2] = v1[2] - v2[2];
}

inline double mag(double v[3])
{
	return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

inline double dotProduct(double v1[3], double v2[3])
{
	return v1[0] * v2[0]
		+ v1[1] * v2[1]
		+ v1[2] * v2[2];
}
