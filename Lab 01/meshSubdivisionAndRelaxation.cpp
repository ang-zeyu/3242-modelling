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

extern int relaxationStepSize;

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
			adjacentTriangle = makeOrTri(adjacentTriangleIdx, adjacentTriangleVer);

			bool isBoundary = !selectedT.test(adjacentTriangleIdx);
			if (isBoundary)
			{
				tcount++;

				// New tlist
				tlist[tcount][0] = tlist[adjacentTriangleIdx][(adjacentTriangleVer + 1) % 3]; // given to new triangle
				tlist[tcount][1] = tlist[adjacentTriangleIdx][(adjacentTriangleVer + 2) % 3]; // still common to both
				tlist[tcount][2] = newMidpointVertices[v];                                    // new midpoint common to both

				// Fnlist for the other two triangle versions
				// also needs to be updated manually here, as we may have a "pacman" shaped selection
				// The triangle inside the pacman mouth may be partially subdivided twice or even thrice.

				// defer updating fnlist[tcount][0] to computeFnList, because it will never be subdivided again in the same run
				OrTri otherTriangle1 = fnext(enext(adjacentTriangle));
				int otherTriangle1Idx = idx(otherTriangle1);
				int otherTriangle1Ver = ver(otherTriangle1) % 3; // normalized ver
				OrTri otherTriangle2 = fnext(enext(enext(adjacentTriangle)));
				int otherTriangle2Idx = idx(otherTriangle2);
				int otherTriangle2Ver = ver(otherTriangle2) % 3; // normalized ver

				fnlist[otherTriangle1Idx][otherTriangle1Ver] = makeOrTri(tcount,
					ver(otherTriangle1) < 3 ? 0 : 3);
				fnlist[otherTriangle2Idx][otherTriangle2Ver] = makeOrTri(adjacentTriangleIdx,
					ver(otherTriangle2) < 3 ? (adjacentTriangleVer + 2) % 3 : adjacentTriangleVer + 2);
				
				for (int j = 0; j < 3; j++)
				{
					nlist[tcount][j] = nlist[adjacentTriangleIdx][j]; // identical normal
					vToTList[tlist[tcount][j]].push_back(tcount);     // Add vToTList mapping
				}

				// Update old adjacent triangle
				// Update vToTList mapping
				vToTList[newMidpointVertices[v]].push_back(adjacentTriangleIdx);     // new v -> t mapping
				vToTList[tlist[adjacentTriangleIdx][(adjacentTriangleVer + 1) % 3]].remove(adjacentTriangleIdx); // remove old mapping

				// Update old tlist - just the middle vertex needs to be moved
				tlist[adjacentTriangleIdx][(adjacentTriangleVer + 1) % 3] = newMidpointVertices[v]; // new midpoint common to both
			}
		}
	}

	computeFnlist();
	computeVertexNormals();
}

double relaxationMaxDeviation = 5.0 * M_PI / 180.0;
void myObjType::relax()
{
	auto cmp = [this](const int& e1, const int& e2) // Comparator for edge's vertex degrees
	{
		int org1 = org(e1);
		int dest1 = dest(e1);

		int org2 = org(e2);
		int dest2 = dest(e2);

		return vToTList[org1].size() + vToTList[dest1].size() - vToTList[org2].size() - vToTList[dest2].size();
	};

	// Pick by shortest edge length
	// Exclude boundary edges on selection
	std::priority_queue<int, std::vector<int>, decltype(cmp)> edgeQueue(cmp);

	for (int t = 1; t <= tcount; t++)
	{
		if (!selectedT.test(t))
		{
			continue;
		}

		for (int i = 0; i < 3; i++)
		{
			// Boundary edges on mesh should not be processed
			// Also boundary edges on selected portion of mesh
			OrTri tri = makeOrTri(t, i);
			OrTri adjacentTri = fnext(tri);
			int adjacentTriidx = idx(adjacentTri);
			if (adjacentTriidx == t || !selectedT.test(adjacentTriidx))
			{
				continue;
			}

			edgeQueue.push(tri);
		}
	}

	// cout << edgeQueue.size() << endl;
	// In any run, modified triangles shall not be checked for relaxation again
	unordered_set<int> modifiedTriangles;

	int numRelaxedEdges = 0;
	while (!edgeQueue.empty() && numRelaxedEdges < relaxationStepSize)
	{
		OrTri edge = edgeQueue.top();
		edgeQueue.pop();

		int orgV = org(edge);
		int destV = dest(edge);

		int tri1 = idx(edge);
		int ortri2 = fnext(edge);
		int tri2 = idx(fnext(edge));

		if (modifiedTriangles.find(tri1) != modifiedTriangles.end()
			|| modifiedTriangles.find(tri2) != modifiedTriangles.end())
		{
			continue;
		}

		double normal1[3];
		computeNormalFor(tri1, normal1);
		double normal2[3];
		computeNormalFor(tri2, normal2);

		// dot product formula to get angle
		double planarAngle = acos(dotProduct(normal1, normal2) / (mag(normal1) * mag(normal2)));
		// cout << dotProduct(normal1, normal2) << " " << mag(normal1) * mag(normal2) << " " <<  planarAngle << " " << relaxationMaxDeviation << endl;
		if (planarAngle <= relaxationMaxDeviation)
		{
			// Swap edge
			int ver1 = ver(edge);  // guaranteed to be ver 0 - 2, as pushed earlier into edgeQueue
			int ver2 = ver(ortri2);

			// Triangle quality check - rollback if unsatisfactory
			double* tri1angles = computeTriangleAngles(vlist[tlist[tri1][0]], vlist[tlist[tri1][1]], vlist[tlist[tri1][2]]);
			double* tri2angles = computeTriangleAngles(vlist[tlist[tri2][0]], vlist[tlist[tri2][1]], vlist[tlist[tri2][2]]);

			double minAngleBefore = min(tri1angles[0], min(tri1angles[1], min(tri1angles[2],
				min(tri2angles[0], min(tri2angles[1], tri2angles[2]))))) * 180.0 / M_PI;

			int vertexBackups[2];
			vertexBackups[0] = tlist[tri1][ver1];
			vertexBackups[1] = tlist[tri2][ver2 < 3 ? ((ver2 + 1) % 3) : (ver2 % 3)];

			// Modify triangle vertex mappings
			tlist[tri1][ver1] = tlist[tri2][(ver2 + 2) % 3];
			tlist[tri2][ver2 < 3 ? ((ver2 + 1) % 3) : (ver2 % 3)] = tlist[tri1][(ver1 + 2) % 3];

			// Triangle quality check - rollback if unsatisfactory
			tri1angles = computeTriangleAngles(vlist[tlist[tri1][0]], vlist[tlist[tri1][1]], vlist[tlist[tri1][2]]);
			tri2angles = computeTriangleAngles(vlist[tlist[tri2][0]], vlist[tlist[tri2][1]], vlist[tlist[tri2][2]]);
			double minAngleAfter = min(tri1angles[0], min(tri1angles[1], min(tri1angles[2],
				min(tri2angles[0], min(tri2angles[1], tri2angles[2]))))) * 180.0 / M_PI;

			if (minAngleBefore < minAngleAfter)
			{
				numRelaxedEdges += 1;
				modifiedTriangles.insert(tri1);
				modifiedTriangles.insert(tri2);

				// Update vToTList mappings
				vToTList[vertexBackups[0]].remove(tri1);
				vToTList[vertexBackups[1]].remove(tri2);
				vToTList[tlist[tri1][ver1]].push_back(tri1);
				vToTList[tlist[tri2][ver2 < 3 ? ((ver2 + 1) % 3) : (ver2 % 3)]].push_back(tri2);
			}
			else
			{
				// rollback
				tlist[tri1][ver1] = vertexBackups[0];
				tlist[tri2][ver2 < 3 ? ((ver2 + 1) % 3) : (ver2 % 3)] = vertexBackups[1];
			}
		}
	}

	computeFnlist();
	computeNormals();
	computeVertexNormals();

	if (numRelaxedEdges == relaxationStepSize)
	{
		cout << "Relaxed " << numRelaxedEdges << " edges!" << endl;
	}
	else
	{
		cout << "Relaxed " << numRelaxedEdges << " edges!" << endl;
		cout << "No or not enough edges fufilled relaxation criteria." << endl;
	}
}

