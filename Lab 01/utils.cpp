#include "mesh.h"

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

struct pairHash
{
	size_t operator()(pair<int, int> const& p) const noexcept
	{
		size_t h1 = hash<int>{}(p.first);
		size_t h2 = hash<int>{}(p.second);
		return h1 ^ (h2 << 1); // adapted from official c++ docs
	}
};