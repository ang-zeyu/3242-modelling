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

extern int decimationStepSize;

/*
 My functionalities concerning mesh decimation
 */

inline unordered_set<string> myObjType::linkV(int v)
{
	list<int> adjacentTriangles = vToTList[v];

	// Star(v)
	// Set of all simplices s.t. each simplex has a face v
	// and
	// Closure(Star(v))
	// Set of all faces of Star(V)...

	unordered_set<string> starV;
	starV.insert(to_string(v)); // vertex itself

	unordered_set<string> closureStarV;

	list<int>::iterator it = adjacentTriangles.begin();
	for (; it != adjacentTriangles.end(); it++)
	{
		int tIdx = *it;
		for (int ver = 0; ver < 3; ver++) // for each triangle version
		{
			// Star(v)

			// The triangle itself
			string triangle = "";
			for (int triVertex = ver; triVertex < (ver + 3); triVertex++)
				triangle += to_string(tlist[tIdx][triVertex % 3]);
			starV.insert(triangle);

			string copyTriangle(triangle);
			std::reverse(copyTriangle.begin(), copyTriangle.end());
			starV.insert(copyTriangle);

			// The adjacent edges
			OrTri adjacentTriangleN = makeOrTri(tIdx, ver);
			int orgV = org(adjacentTriangleN);
			int destV = dest(adjacentTriangleN);
			if (orgV == v || destV == v)
			{
				starV.insert(to_string(orgV) + to_string(destV));
				starV.insert(to_string(destV) + to_string(orgV));
			}

			// Closure(Star(v))

			// Triangles itself
			closureStarV.insert(triangle);
			closureStarV.insert(copyTriangle);

			// ALL edges - regardless of adjacent or not, because it is a face of the triangles
			closureStarV.insert(to_string(orgV) + to_string(destV));
			closureStarV.insert(to_string(destV) + to_string(orgV));
		}

		// All vertices, because it is a face of the triangles
		for (int v = 0; v < 3; v++)
			closureStarV.insert(to_string(tlist[tIdx][v]));
	}

	// Closure(v)
	// Set of all faces of v... is just v =)
	// Then Star(Closure(v)) = Star(v)

	// Link
	// Cl(St(v)) - St(Cl(v))
	unordered_set<string> result;
	unordered_set<string>::iterator clStVIt = closureStarV.begin();
	for (; clStVIt != closureStarV.end(); clStVIt++)
	{
		if (starV.find(*clStVIt) == starV.end())
		{
			result.insert(*clStVIt);
		}
	}

	return result;
}

inline unordered_set<string> myObjType::linkE(OrTri t)
{
	int orgV = org(t);
	int destV = dest(t);

	// Star(e)
	// Set of all simplices s.t. each simplex has a face v
	// and
	// Closure(Star(e))
	// Set of all faces of Star(e)...

	unordered_set<string> starE;
	starE.insert(to_string(orgV) + to_string(destV)); // edge itself
	starE.insert(to_string(destV) + to_string(orgV)); // edge itself

	unordered_set<string> closureStarE;

	list<OrTri> adjacentTriangles({ idx(t), idx(fnext(t)) });

	list<int>::iterator it = adjacentTriangles.begin();
	for (; it != adjacentTriangles.end(); it++)
	{
		int tIdx = *it;
		for (int ver = 0; ver < 3; ver++) // for each triangle version
		{
			// Star(e)

			// The triangle itself
			string triangle = "";
			for (int v = ver; v < (ver + 3); v++)
				triangle += to_string(tlist[tIdx][v % 3]);
			starE.insert(triangle);

			string copyTriangle(triangle);
			std::reverse(copyTriangle.begin(), copyTriangle.end());
			starE.insert(copyTriangle);

			// Closure(Star(e))

			// Triangles itself
			closureStarE.insert(triangle);
			closureStarE.insert(copyTriangle);

			// All edges, because it is a face of the triangles
			OrTri adjacentTriangleN = makeOrTri(tIdx, ver);
			int orgV = org(adjacentTriangleN);
			int destV = dest(adjacentTriangleN);
			closureStarE.insert(to_string(orgV) + to_string(destV));
			closureStarE.insert(to_string(destV) + to_string(orgV));
		}

		// All vertices, because it is a face of the triangles
		for (int v = 0; v < 3; v++)
			closureStarE.insert(to_string(tlist[tIdx][v]));
	}

	// Closure(e)
	// Set of all faces of e... is just orgV, destV and the edge itself
	// So,
	// Star(Closure(e))
	// All simplices that have face of the edge's vertices or the edge itself...
	// = Star(orgV) + Star(destV) + Star(edge)

	unordered_set<string> starClosureE;

	unordered_set<string>::iterator starEIt = starE.begin();
	for (; starEIt != starE.end(); starEIt++)
		starClosureE.insert(*starEIt);

	int closureEVertices[2];
	closureEVertices[0] = orgV;
	closureEVertices[1] = destV;

	for (int i = 0; i < 2; i++)
	{
		int v = closureEVertices[i];
		starClosureE.insert(to_string(v)); // vertex itself

		list<int> vAdjacentTriangles = vToTList[v];

		it = vAdjacentTriangles.begin();
		for (; it != vAdjacentTriangles.end(); it++)
		{
			int tIdx = *it;
			for (int ver = 0; ver < 3; ver++) // for each triangle version
			{
				// Star(v)

				// The triangle itself
				string triangle = "";
				for (int triVertex = ver; triVertex < (ver + 3); triVertex++)
					triangle += to_string(tlist[tIdx][triVertex % 3]);
				starClosureE.insert(triangle);

				string copyTriangle(triangle);
				std::reverse(copyTriangle.begin(), copyTriangle.end());
				starClosureE.insert(copyTriangle);

				// The adjacent edges
				OrTri adjacentTriangleN = makeOrTri(tIdx, ver);
				int adjacentOrgV = org(adjacentTriangleN);
				int adjacentDestV = dest(adjacentTriangleN);
				if (adjacentOrgV == v || adjacentDestV == v)
				{
					starClosureE.insert(to_string(adjacentOrgV) + to_string(adjacentDestV));
					starClosureE.insert(to_string(adjacentDestV) + to_string(adjacentOrgV));
				}
			}
		}
	}

	// Link
	// Cl(St(e)) - St(Cl(e))
	unordered_set<string> result;
	unordered_set<string>::iterator closureStarEIt = closureStarE.begin();
	for (; closureStarEIt != closureStarE.end(); closureStarEIt++)
	{
		// cout << string(*closureStarEIt) << " ";
		if (starClosureE.find(*closureStarEIt) == starClosureE.end())
		{
			result.insert(string(*closureStarEIt));	
		}
	}
	// cout << endl;

	return result;
}

void myObjType::decimate()
{
	auto cmp = [this](const int& e1, const int& e2) // Comparator for edge lengths
	{
		int org1 = org(e1);
		int dest1 = dest(e1);
		double subResult1[3];
		subtractV(subResult1, vlist[org1], vlist[dest1]);
		double edgeLen1 = mag(subResult1);

		int org2 = org(e2);
		int dest2 = dest(e2);
		double subResult2[3];
		subtractV(subResult2, vlist[org2], vlist[dest2]);
		double edgeLen2 = mag(subResult2);

		return edgeLen2 < edgeLen1;
	};

	// Pick by shortest edge length
	// Exclude boundary edges on selection
	std::priority_queue<int, std::vector<int>, decltype(cmp)> edgeQueue(cmp);
	
	// Build boundary triangle exclusion set, which can create holes
	unordered_set<OrTri> boundaryTriangleVertexExclusions;
	for (int t = 1; t <= tcount; t++)
	{
		for (int v = 0; v < 3; v++)
		{
			OrTri currentTriangle = makeOrTri(t, v);
			OrTri adjacentTriangle = fnext(currentTriangle);     // triangle adjacent to version i of current triangle
			if (adjacentTriangle == currentTriangle)
			{
				boundaryTriangleVertexExclusions.insert(org(currentTriangle));
				boundaryTriangleVertexExclusions.insert(dest(currentTriangle));
				break;
			}
		}
	}

	for (int t = 1; t <= tcount; t++)
	{
		if (!selectedT.test(t))
			continue;

		for (int v = 0; v < 3; v++)
		{
			OrTri currentTriangle = makeOrTri(t, v);
			OrTri adjacentTriangle = fnext(currentTriangle);     // triangle adjacent to version i of current triangle
			int adjacentTriangleIdx = idx(adjacentTriangle);			

			// This checks if it is a boundary edge on the **user selection**
			// unlike the next one which checks if it is a boundary edge on the entire mesh
			if (selectedT.test(adjacentTriangleIdx))
			{
				// Use earlier built set to exclude edges touching boundary vertices of the **entire mesh**
				if (boundaryTriangleVertexExclusions.find(org(currentTriangle)) != boundaryTriangleVertexExclusions.end()
					|| boundaryTriangleVertexExclusions.find(dest(currentTriangle)) != boundaryTriangleVertexExclusions.end())
				{
					// cout << "Boundary triangle vertex exclusion fail!" << endl;
					continue;
				}

				edgeQueue.push(currentTriangle);
			}
		}
	}

	int numContractedEdges = 0;
	while (!edgeQueue.empty() && numContractedEdges < decimationStepSize)
	{
		OrTri edge = edgeQueue.top();
		edgeQueue.pop();

		/*cout << idx(edge) << " " << ver(edge) << " " << org(edge) << " " << dest(edge) 
			<< " " << tlist[idx(edge)][0] << " " << tlist[idx(edge)][1] << " " << tlist[idx(edge)][2] << endl;*/
		int orgV = org(edge);
		int destV = dest(edge);

		unordered_set<string> linkOrgV = linkV(orgV);
		unordered_set<string> linkDestV = linkV(destV);
		// cout << linkOrgV.size() << endl;

		unordered_set<string> linkOrgVAndDestV;

		unordered_set<string>::iterator setIt = linkOrgV.begin();
		for (; setIt != linkOrgV.end(); setIt++)
		{
			if (linkDestV.find(*setIt) != linkDestV.end())
			{
				linkOrgVAndDestV.insert(*setIt);
			}
		}

		unordered_set<string> linkEdge = linkE(edge);

		/*cout << linkOrgVAndDestV.size() << " " << linkEdge.size() << endl;;

		setIt = linkOrgVAndDestV.begin();
		for (; setIt != linkOrgVAndDestV.end(); setIt++)
		{
			cout << *setIt << " ";
		}
		cout << endl;

		setIt = linkEdge.begin();
		for (; setIt != linkEdge.end(); setIt++)
		{
			cout << *setIt << " ";
		}
		cout << endl;*/

		if (linkOrgVAndDestV == linkEdge)
		{
			// Contract org to dest
			if (contract(edge))
			{
				numContractedEdges += 1;
			}
			/*else
			{
				cout << "triangle flipped!" << endl;
			}*/
		}
		/*else
		{
			cout << "Link not equal!" << endl;
		}*/
	}

	if (numContractedEdges == decimationStepSize)
	{
		cout << "Contracted " << numContractedEdges << " edges!" << endl;
	}
	else
	{
		cout << "Contracted " << numContractedEdges << " edges!" << endl;
		cout << "No or not enough edges fufilled contraction criteria." << endl;
	}

	computeFnlist();
	computeNormals();
	computeVertexNormals();
}

bool myObjType::contract(OrTri edge)
{
	int orgV = org(edge);
	int destV = dest(edge);

	unordered_set<int> trianglesToDestroy;
	OrTri nextTri = edge;
	do
	{
		trianglesToDestroy.insert(idx(nextTri));
		nextTri = fnext(nextTri);
	} while (idx(nextTri) != idx(edge));

	// Triangles to update are the ones connected to orgV
	unordered_set<int> trianglesToUpdate;
	list<int>::iterator listIt = vToTList[orgV].begin();
	for (; listIt != vToTList[orgV].end(); listIt++)
	{
		if (trianglesToDestroy.find(*listIt) == trianglesToDestroy.end())
		{
			trianglesToUpdate.insert(*listIt);
		}
	}

	unordered_set<int>::iterator setIt = trianglesToUpdate.begin();
	unordered_map<int, int> updatedTriangles; // for undoing if triangle is flipped
	for (; setIt != trianglesToUpdate.end(); setIt++)
	{
		int v = -1;
		for (int i = 0; i < 3; i++)
		{
			if (tlist[*setIt][i] == orgV)
			{
				v = i;
				break;
			}
		}

		if (v == -1)
		{
			cout << "Tried to update triangle not connected to OrgV:" << orgV << " t:" << *setIt << " "
				<< tlist[*setIt][0] << " " << tlist[*setIt][1] << " " << tlist[*setIt][2] << endl;
			continue;
		}

		/*cout << "OrgV:" << orgV << "DestV:" << destV << " t:" << *setIt << " "
			<< tlist[*setIt][0] << " " << tlist[*setIt][1] << " " << tlist[*setIt][2] << endl;*/

		tlist[*setIt][v] = destV;
		updatedTriangles.insert({ *setIt, v });

		/*
		Flipped normal check
		*/

		double newNormal[3];
		computeNormalFor(*setIt, newNormal);
		double dotProd = dotProduct(nlist[*setIt], newNormal);

		if (dotProd < 0)
		{
			// Undo vertex updates
			unordered_map<int, int>::iterator mapIt = updatedTriangles.begin();
			for (; mapIt != updatedTriangles.end(); mapIt++)
			{
				tlist[mapIt->first][mapIt->second] = orgV;
			}

			return false;
		}
	}

	// Successfully contracted, update vToTList mapping
	setIt = trianglesToUpdate.begin();
	for (; setIt != trianglesToUpdate.end(); setIt++)
	{
		vToTList[destV].push_back(*setIt);
	}

	setIt = trianglesToDestroy.begin();
	for (; setIt != trianglesToDestroy.end(); setIt++)
	{
		for (int v = 0; v < 3; v++)
		{
			vToTList[tlist[*setIt][v]].remove(*setIt);
		}
	}

	// Fill tlist holes in, and update some variables
	setIt = trianglesToDestroy.begin();
	int numTrianglesDestroyed = 0;
	for (int t = 1; t <= tcount; t++)
	{
		if (trianglesToDestroy.find(t) != trianglesToDestroy.end())
		{
			numTrianglesDestroyed += 1;
			continue;
		}

		int prev = t - numTrianglesDestroyed;
		if (prev != t)
		{
			for (int v = 0; v < 3; v++)
			{
				tlist[prev][v] = tlist[t][v];
				vToTList[tlist[prev][v]].remove(t);
				vToTList[tlist[prev][v]].push_back(prev);
			}

			selectedT.set(prev, selectedT.test(t));
			selectedT.set(t, false);
		}
	}

	tcount -= trianglesToDestroy.size();

	vToTList[orgV] = list<int>(); // no longer used

	return true;
}
