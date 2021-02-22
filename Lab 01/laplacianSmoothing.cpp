#include "mesh.h"
#include "utils.h"

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
#include <iterator>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <queue>
#include <iomanip>
using namespace std;

/*
 My functionalities concerning laplacian smoothing
 */

void myObjType::smooth()
{
	unordered_set<int> triangleVertices;
	unordered_set<int> triangleVerticesBlacklist;
	for (int t = 1; t <= tcount; t++)
	{
		if (!selectedT.test(t))
		{
			continue;
		}

		bool isBoundaryEdge[3];
		for (int i = 0; i < 3; i++)
			isBoundaryEdge[i] = fnlist[t][i] == makeOrTri(t, i);

		for (int i = 0; i < 3; i++)
		{
			// boundary vertices should not be processed
			if (isBoundaryEdge[i] || isBoundaryEdge[(i + 2) % 3])
			{
				triangleVerticesBlacklist.insert(tlist[t][i]);
				continue;
			}

			triangleVertices.insert(tlist[t][i]);
		}
	}

	// Post remove blacklisted vertices as it may be blacklisted after it was processed
	unordered_set<int>::iterator blacklistIt = triangleVerticesBlacklist.begin();
	for (; blacklistIt != triangleVerticesBlacklist.end(); blacklistIt++)
	{
		if (triangleVertices.find(*blacklistIt) != triangleVertices.end())
		{
			triangleVertices.erase(*blacklistIt);
		}
	}

	// Laplacian smoothing

	// delay updates so all calculations are based on current state only
	unordered_map<int, double*> vertexUpdates;

	unordered_set<int>::iterator selectedVerticesIt = triangleVertices.begin();
	for (; selectedVerticesIt != triangleVertices.end(); selectedVerticesIt++)
	{
		int v = *selectedVerticesIt;

		unordered_set<int> vertexNeighbours;

		// Find neighbours
		list<int> vertexTriangles = vToTList[v];
		list<int>::iterator vertexTrianglesIt = vertexTriangles.begin();
		for (; vertexTrianglesIt != vertexTriangles.end(); vertexTrianglesIt++)
		{
			for (int i = 0; i < 3; i++)
				vertexNeighbours.insert(tlist[*vertexTrianglesIt][i]);

		}
		vertexNeighbours.erase(v);

		// Sum then average
		double sum[3]; sum[0] = 0; sum[1] = 0; sum[2] = 0;
		unordered_set<int>::iterator neighboursIt = vertexNeighbours.begin();
		for (; neighboursIt != vertexNeighbours.end(); neighboursIt++)
		{
			for (int i = 0; i < 3; i++)
				sum[i] += vlist[*neighboursIt][i];
		}
		for (int i = 0; i < 3; i++)
			sum[i] /= vertexNeighbours.size();

		// Calculate new vertices
		double* update = (double*)malloc(sizeof(double) * 3);
		for (int i = 0; i < 3; i++)
		{
			update[i] = vlist[v][i] + 0.1 * (sum[i] - vlist[v][i]);
			// cout << vlist[v][i] << " " << update[i] << endl;
		}

		vertexUpdates.insert({ v, update });
	}

	// Update vertices to new positions
	unordered_map<int, double*>::iterator updatesIt = vertexUpdates.begin();
	for (; updatesIt != vertexUpdates.end(); updatesIt++)
	{
		for (int i = 0; i < 3; i++)
		{
			vlist[updatesIt->first][i] = updatesIt->second[i];
			// cout << vlist[(*updatesIt).first][i] << " " << (*updatesIt).second[i] << endl;
		}
		free(updatesIt->second);
	}

	computeNormals();
	computeVertexNormals();

	cout << "Smoothed selected portion of mesh" << endl;
}

