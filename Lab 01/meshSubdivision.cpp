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
 My functionalities concerning mesh subdivision and relaxation
 */

inline void average(double result[3], double v1[3], double v2[3])
{
	result[0] = (v1[0] + v2[0]) / 2;
	result[1] = (v1[1] + v2[1]) / 2;
	result[2] = (v1[2] + v2[2]) / 2;
}


void myObjType::subdivide()
{
	/*
	 During subdivision, new vertices may be shared by multiple original triangles
	 Avoid constructing redundant vertices
	 */
	unordered_map<pair<int, int>, int, pairHash> edgeVertexMap;

	int originalTCount = tcount; // fix iteration to original triangles
	for (int t = 1; t <= originalTCount; t++)
	{
		if (!selectedT.test(t))
		{
			continue;
		}

		pair<int, int> edges[3];
		edges[0] = make_pair(min(tlist[t][0], tlist[t][1]), max(tlist[t][0], tlist[t][1]));
		edges[1] = make_pair(min(tlist[t][1], tlist[t][2]), max(tlist[t][1], tlist[t][2]));
		edges[2] = make_pair(min(tlist[t][0], tlist[t][2]), max(tlist[t][0], tlist[t][2]));

		// Add new midpoint vertices as necessary
		int newMidpointVertices[3];
		for (int i = 0; i < 3; i++)
		{
			pair<int, int> edge = edges[i];
			if (edgeVertexMap.find(edge) == edgeVertexMap.end())
			{
				vcount++;
				average(vlist[vcount], vlist[tlist[t][i]], vlist[tlist[t][(i + 1) % 3]]);
				edgeVertexMap.insert({ edge, vcount });
				newMidpointVertices[i] = vcount;

				vToTList[vcount].clear(); // ensure
			}
			else
			{
				newMidpointVertices[i] = edgeVertexMap.at(edge);
			}
		}

		// Construct 3 more subdivided triangles
		for (int i = 0; i < 3; i++)
		{
			tcount++;

			tlist[tcount][0] = tlist[t][i];
			tlist[tcount][1] = newMidpointVertices[i];
			tlist[tcount][2] = newMidpointVertices[(i + 2) % 3];

			// Identical normal
			nlist[tcount][0] = nlist[t][0];
			nlist[tcount][1] = nlist[t][1];
			nlist[tcount][2] = nlist[t][2];

			for (int j = 0; j < 3; j++)
			{
				vToTList[tlist[tcount][j]].push_back(tcount);
			}

			selectedT.set(tcount, true);
		}

		// Readjust current triangle to be central triangle
		for (int i = 0; i < 3; i++)
		{
			vToTList[tlist[t][i]].remove(t);                    // remove old v -> t mapping

			tlist[t][i] = newMidpointVertices[i];
			vToTList[tlist[t][i]].push_back(t); // new v -> t mapping
		}

		// For each version of current triangle
		// If it is a boundary edge (on selected portion), 
		// also need to partially subdivide the adjacent triangle
		for (int v = 0; v < 3; v++)
		{
			OrTri adjacentTriangle = fnext(makeOrTri(t, v)); // triangle adjacent to version i of current triangle
			int adjacentTriangleVer = ver(adjacentTriangle) % 3; // normalized version
			int adjacentTriangleIdx = idx(adjacentTriangle);

			bool isBoundary = !selectedT.test(adjacentTriangleIdx);
			if (isBoundary)
			{
				tcount++;

				// New tlist
				tlist[tcount][0] = tlist[adjacentTriangleIdx][(adjacentTriangleVer + 1) % 3];
				tlist[tcount][1] = tlist[adjacentTriangleIdx][(adjacentTriangleVer + 2) % 3];
				tlist[tcount][2] = newMidpointVertices[v];

				nlist[tcount][0] = nlist[adjacentTriangleIdx][0];
				nlist[tcount][1] = nlist[adjacentTriangleIdx][1];
				nlist[tcount][2] = nlist[adjacentTriangleIdx][2];

				// Add vToTList mapping
				for (int j = 0; j < 3; j++)
					vToTList[tlist[tcount][j]].push_back(tcount);

				// Update old adjacent triangle
				// Update vToTList mapping
				vToTList[newMidpointVertices[v]].push_back(adjacentTriangleIdx);     // new v -> t mapping
				vToTList[tlist[adjacentTriangleIdx][(adjacentTriangleVer + 1) % 3]]
					.remove(adjacentTriangleIdx); // remove old mapping

				// Update tlist - just middle vertex
				tlist[adjacentTriangleIdx][(adjacentTriangleVer + 1) % 3] = newMidpointVertices[v];
			}
		}
	}

	computeFnlist();
	computeVertexNormals();
}

void myObjType::relax()
{
	// in any run, modified triangles shall not be checked for relaxation again
	// unordered_set<int> modifiedTriangles;

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

		/*double* currentTriangleAngles = computeAngles(vlist[tlist[t][0]], vlist[tlist[t][1]], vlist[tlist[t][2]]);
		double currTriMinAngle = min(currentTriangleAngles[0], min(currentTriangleAngles[1], currentTriangleAngles[2]));
		double currTriMaxAngle = max(currentTriangleAngles[0], max(currentTriangleAngles[1], currentTriangleAngles[2]));

		for (int v = 0; v < 3; v++)
		{
			OrTri adjacentTriangle = fnlist[t][v];
			if (adjacentTriangle == makeOrTri(t, v))
			{
				continue; // boundary edge
			}

			int adjacentTriangleIdx = idx(adjacentTriangle);
			int adjacentTriangleVer = ver(adjacentTriangle) % 3; // normalized to 0 - 3, for finding the non-shared vertex
			double* adjacentTriangleAngles = computeAngles(
				vlist[tlist[adjacentTriangleIdx][0]],
				vlist[tlist[adjacentTriangleIdx][1]],
				vlist[tlist[adjacentTriangleIdx][2]]);
			double minAngle = min(currTriMinAngle, min(adjacentTriangleAngles[0], min(adjacentTriangleAngles[1], adjacentTriangleAngles[2])));
			double maxAngle = max(currTriMaxAngle, max(adjacentTriangleAngles[0], max(adjacentTriangleAngles[1], adjacentTriangleAngles[2])));

			double* newTriangleAngles1 = computeAngles(
				vlist[tlist[t][(v + 1) % 3]],
				vlist[tlist[t][(v + 2) % 3]],
				vlist[tlist[adjacentTriangleIdx][(adjacentTriangleVer + 2) % 3]]);
			double* newTriangleAngles2 = computeAngles(
				vlist[tlist[t][v]],
				vlist[tlist[t][(v + 2) % 3]],
				vlist[tlist[adjacentTriangleIdx][(adjacentTriangleVer + 2) % 3]]);

			double newMinAngle = min(
				newTriangleAngles1[0],
				min(newTriangleAngles1[1],
					min(newTriangleAngles1[2],
						min(newTriangleAngles2[0],
							min(newTriangleAngles2[1],
								newTriangleAngles2[2])
						)))
			);

			if (newMinAngle <= minAngle)
			{
				continue;
			}
		}*/
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

	cout << "Relaxed selected portion of mesh" << endl;
}

