#pragma once

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
#include <unordered_set>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <bitset>
#include <vector>
#include <list>
#include <map>
#include <queue>
#include <iomanip>
using namespace std;

// maximum number of vertices and triangles
#define MAXV 1000000
#define MAXT 1000000

typedef int OrTri;
typedef int tIdx;

inline OrTri makeOrTri(tIdx t, int version)
{
	return (t << 3) | version;
};

inline tIdx idx(OrTri ot)
{
	return ot >> 3;
};

inline int ver(OrTri ot)
{
	return ot & 0b111;
};

inline OrTri enext(OrTri ot)
{
	int v = ver(ot);
	return makeOrTri(idx(ot),
		v < 3 ? ((v + 1) % 3) : (3 + ((v - 1) % 3)));
};

inline OrTri sym(OrTri ot)
{
	int v = ver(ot);
	return v < 3 ? ot + 3 : ot - 3;
};

class myObjType {
	int vcount = 0;
	int tcount = 0;

	double vlist[MAXV][3];   // vertices list
	list<int> vToTList[MAXV];   // vertex : triangle list, for vertex normal calculation
	double vnlist[MAXV][3];  // vertex normal list

	int tlist[MAXT][3];      // triangle list
	int fnlist[MAXT][3];     // fnext list
	double nlist[MAXT][3];   // storing triangle normals
	
	// Lab 2 user selection task
	bitset<MAXT + 1> selectedT;
	bitset<MAXT + 1> visibleT; // for isSelectingFacing mode (ctrl-alt click drag)

	// Laplacian deformation final boss
	int selectedV = 0;

	double lmax[3];          // the maximum coordinates of x,y,z
	double lmin[3];          // the minimum coordinates of x,y,z

	int statMinAngle[18]; // each bucket is  degrees has a 10 degree range from 0 to 180 degree
	int statMaxAngle[18]; 

	int numCc = 0;

public:
	myObjType()
	{
		vcount = 0;
		tcount = 0;
	};
	void readFile(char* filename);  // assumming file contains a manifold
	void writeFile(char* filename);  
	void draw();
    void computeStat();

	void drawOffscreen();            // Lab 2 optional task - user select marquee, for ctrl-alt click drag mode
	void computeSelectedTriangles(); // Lab 2 Optional task - user select marquee

	// Lab 2 Main task
	inline int org(OrTri t)
	{
		/*
		 Versions 0 to 2:
		 // First 3 would just be the version number
		 // Latter 3 takes on the next column (due to sym)
		 0 1 2
		 1 2 0
		 2 0 1
		*/
		int v = ver(t);
		return tlist[idx(t)][v < 3 ? v : (v + 1) % 3];
	};

	// Lab 2 Main task
	inline int dest(OrTri t)
	{
		/*
		 dest(ver 0) = dest(ver 4) = 1 (index of vertex in tlist)
		 dest(ver 1) = dest(ver 5) = 2
		 dest(ver 2) = dest(ver 3) = 0
		*/
		int v = ver(t);
		return tlist[idx(t)][v == 2 ? 0 : ((v + 1) % 4)];
	};

	// Lab 2 Main task
	inline OrTri fnext(OrTri t)
	{
		int v = ver(t);
		return v < 3
			? fnlist[idx(t)][v]
			: sym(fnlist[idx(t)][v - 3]);
	};

	// Final boss(es)
	void subdivide();
	void relax();
	void decimate();

	void displace(double* displacement); // laplacian deformation

private:
	void reset();
	void readObjFile(char* filename);
	void read3dsFile(char* filename); // Optional Task 4
	void postReadFile();

	void computeNormalFor(int t, double result[3]);
	void computeNormals();
	void computeVertexNormals();
	void computeFnlist();
	void computeNumCc();
	bool orientTriangles();

	// Final boss

	// Mesh decimation
	unordered_set<string> linkV(int v);
	unordered_set<string> linkE(OrTri t);
	bool contract(OrTri edge);

	// Laplacian deformation
	void computeSelectedVertex();
};


