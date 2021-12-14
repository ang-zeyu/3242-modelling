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

#include <Eigen/Dense>;
using Eigen::MatrixXd;
using Eigen::VectorXd;


/*
 My functionalities concerning laplacian deformation
 */


// Select the vertex common to all selected triangles in the intersection box
// Fail if no such vertex is found (selection box too large)
void myObjType::computeSelectedVertex()
{
	selectedV = 0;

	unordered_map<int, int> vertexDegreeMap;
	int currentMax = 0;
	int numSelectedT = 0;

	for (int t = 1; t <= tcount; t++)
	{
		if (!selectedT.test(t))
		{
			continue;
		}

		numSelectedT++;

		for (int v = 0; v < 3; v++)
		{
			int amount = 0;
			if (vertexDegreeMap.find(tlist[t][v]) == vertexDegreeMap.end())
			{
				amount = 1;
			}
			else
			{
				amount = vertexDegreeMap.at(tlist[t][v]) + 1;
				vertexDegreeMap.erase(tlist[t][v]);
			}

			vertexDegreeMap.insert({ tlist[t][v], amount });
			if (amount > currentMax)
			{
				currentMax = amount;
				selectedV = tlist[t][v];
			}
		}
	}

	if (numSelectedT != currentMax)
	{
		cout << "Selected too many triangles! Select only a single vertex for a more accurate selection." << endl;
	}
}

void myObjType::displace(double* displacement)
{
    if (selectedV == 0)
    {
		cout << "Select some triangles using ctrl or ctrl-alt click drag, then a vertex using alt click drag!" << endl;
        return;
    }

	unordered_set<int> selectedVertices;

	// Bfs - optimization - only compute for vertices of the same component of the selected vertex
	unordered_set<int> processedTriangles;
	queue<int> triangleQueue; // bfs queue
	triangleQueue.push(*(vToTList[selectedV].begin()));
	while (!triangleQueue.empty())
	{
		int currentTriangle = triangleQueue.front();
		triangleQueue.pop();
		processedTriangles.insert(currentTriangle);

		// Add the vertices
		for (int v = 0; v < 3; v++)
		{
			selectedVertices.insert(tlist[currentTriangle][v]);
		}

		OrTri* fnTriangles = fnlist[currentTriangle];
		for (int i = 0; i < 3; i++)
		{
			int fnextTriangleIdx = idx(fnTriangles[i]);
			
			if (
				// Skip self-loops and triangles processed before due to traversal order
				processedTriangles.find(fnextTriangleIdx) == processedTriangles.end()
				// Also, limit traversal to selected triangles from the selected vertex, because in minimising shape,
				// 'Disjoint' selected triangles from the selected vertex but in the same component
				// Would stay the same anyway
				&& selectedT.test(fnextTriangleIdx)
				)
			{
				triangleQueue.push(fnextTriangleIdx);
			}
		}
	}

	// Repeat,
	// add boundary vertices beyond selection as additional position constraints
	unordered_set<int> boundaryVertices;
	int numselected = 0;
	for (int t : processedTriangles)
	{
		if (selectedT.test(t))
		{
			numselected += 1;
			// Boundary vertices beyond selection
			OrTri* fnTriangles = fnlist[t];
			for (int i = 0; i < 3; i++)
			{
				int fnextTriangleIdx = idx(fnTriangles[i]);
				if (!selectedT.test(fnextTriangleIdx))
				{
					// cout << "boundary " << fnextTriangleIdx << endl;
					for (int v = 0; v < 3; v++)
					{
						if (selectedVertices.find(tlist[fnextTriangleIdx][v]) == selectedVertices.end())
						{
							boundaryVertices.insert(tlist[fnextTriangleIdx][v]);
						}
					}
				}
			}
		}
	}
	// cout << numselected << " " << selectedVertices.size() << endl;

	// Structured as - selected vertex -- beyond boundary vertices -- misc vertices
	vector<int> vertexList;
	vertexList.push_back(selectedV);
	for (int v : boundaryVertices)
	{
		vertexList.push_back(v);
	}
	for (int v : selectedVertices)
	{
		if (v != selectedV)
		{
			vertexList.push_back(v);
		}
	}
	// cout << boundaryVertices.size() << endl;

	// cout << "Obtained all " << totalNumVertices << " vertices of selected triangles" << endl;

	/*
	 Form and solve the least squares system.
	 */

	int totalNumVertices = vertexList.size();
	int totalEntries = totalNumVertices * 3;

	/*
	 Lls (linear least squares) system of equal weight laplacian operator.

	 Denote this Ax = b;

	 A is a totalEntries * totalEntries matrix.
	 - For i, j within the same vertex, the 3x3 entries are the identity matrix.
	 - For i, j not within the same vertex,
	   where the vertex represented by j is a neighbour of the vertex represented by i,
	   the 3x3 sub matrix's **diagonal** entries are: -1 / num neighbours
	 - Hence with x and b, this forms a lls system minimizing shape differences. (laplacian operator)

	 x is a totalEntries * 1 vector of *vertex positions* (in x1, y1, z1, x2, y2, z2, ... format)
	 we wish to solve with lls.

	 b is a totalEntries * 1 vector of *laplacian operator results* (similarly in format to x)
	 */
	// 
	MatrixXd A = MatrixXd::Zero(totalEntries, totalEntries);
	VectorXd b = VectorXd::Zero(totalEntries);

	for (int v1 = 0; v1 < totalNumVertices; v1++)
	{
		int currentVertex = vertexList.at(v1);

		// Get neighbours of current vertex
		unordered_set<int> neighbours;
		for (int t : vToTList[currentVertex])
		{
			for (int v = 0; v < 3; v++)
			{
				// cout << "n " << tlist[t][v] << " ";
				if (selectedVertices.find(tlist[t][v]) != selectedVertices.end()
					|| boundaryVertices.find(tlist[t][v]) != boundaryVertices.end())
				{
					neighbours.insert(tlist[t][v]);
				}
			}
		}
		neighbours.erase(currentVertex);
		int numNeighbours = neighbours.size();

		/*if (boundaryVertices.find(currentVertex) != boundaryVertices.end())
		{
			cout << currentVertex << " boundary Num Neighbours: " << numNeighbours << endl;
		}
		else
		{
			cout << currentVertex << " selected Num Neighbours: " << numNeighbours << endl;
		}*/

		// Equal weight factor
		double negativeOneOverN = -1.0 / (double)numNeighbours;

		// Laplacian operator result contribution from current vertex
		int baseRowIdx = v1 * 3;
		for (int row = baseRowIdx; row < baseRowIdx + 3; row++)
		{
			b(row) = vlist[currentVertex][row % 3];
		}

		// Laplacian operator contributions from neighbour vertices
		for (int neighbour : neighbours)
		{
			for (int row = baseRowIdx; row < baseRowIdx + 3; row++)
			{
				b(row) += (negativeOneOverN * vlist[neighbour][row % 3]);
			}
		}

		// Fill A
		for (int v2 = 0; v2 < totalNumVertices; v2++)
		{
			if (v1 == v2)
			{
				// Identity sub matrix
				for (int k = baseRowIdx; k < baseRowIdx + 3; k++)
				{
					A(k, k) = 1;
				}
			}
			else
			{
				// Diagonal - 1 / num neighbour entries
				int neighbourVertex = vertexList.at(v2);
				// cout << neighbourVertex << " ";
				if (neighbours.find(neighbourVertex) != neighbours.end())
				{
					int baseColIdx = v2 * 3;
					for (int k = 0; k < 3; k++)
					{
						A(baseRowIdx + k, baseColIdx + k) = negativeOneOverN;
					}
				}
			}
		}
	}
	// cout << A;

	/*
	 Hard positional constraints are denoted by Cx = d.
	 There are g = 1 (handle vertex) + z (boundaryVertices.size()) positional constraints.

	 C is a g * totalEntries matrix,
	 where the 3x3 sub matrix beginning at (selectedV * 3) is the identity matrix.

	 x is the same as above.

	 d (g * 1 column vector) is simply as described below. (the value of the positional constraints)
	 */
	int numConstraints = 3 + boundaryVertices.size() * 3;
	VectorXd d(numConstraints);
	for (int i = 0; i < 3; i++)
	{
		d(i) = vlist[selectedV][i] + displacement[i];
	}
	for (int i = 3; i < numConstraints; i++)
	{
		d(i) = vlist[vertexList.at(i / 3)][i % 3];
	}

	/*
	 Solve using QR factorization
	 By the special properties of identity matrix, the QR factorization need not be computed.

	 Rather the lls solution to Ax = b satisfying Cx = d is
	 x = Q1 u + Q2 v = [d; v],
	 where v = (A2^T A2)^-1 A2^T (b - A1 d)

	 where A1 (totalEntries * 3) A2 is (totalEntries * (totalEntries - 3)) sub-matrices of A
	 */
	// cout << A.rows() << "," << A.cols() << endl;

	MatrixXd A1 = A.block(0, 0,              totalEntries, numConstraints);
	MatrixXd A2 = A.block(0, numConstraints, totalEntries, totalEntries - numConstraints);

	// Finally,
	MatrixXd A2T = A2.transpose();
	VectorXd v = (A2T * A2).inverse() * A2T * (b - A1 * d);

	VectorXd x(d.size() + v.size()); // solution
	x << d,
		 v;
	// cout << x;

	// Update positions...
	for (int v = 0; v < totalNumVertices; v++)
	{
		int baseRow = v * 3;
		int currentVertex = vertexList.at(v);
		for (int row = baseRow; row < baseRow + 3; row++)
		{
			vlist[currentVertex][row % 3] = x(row);
		}
	}

	// yay!
	free(displacement);
	computeNormals();
	computeVertexNormals();
}

long substrCount(int n, string s) {
	long numSubstrings = 0;

	// a aa aaa ... contribution
	char prevChar = ' ';
	int prevCharCount = 0;
	for (int i = 0; i < s.length(); i++)
	{
		if (s[i] != prevChar && prevChar != ' ')
		{
			numSubstrings += (prevCharCount * (prevCharCount + 1)) / 2;

			prevChar = s[i];
			prevCharCount = 1;
		}
		else
		{
			prevChar = s[i];
			prevCharCount += 1;
		}
	}

	
	numSubstrings += (prevCharCount * (prevCharCount + 1)) / 2;
	

	// ada axa aaxaa ... contributions
	// Normally, such a string can be split into 3 segments,
	// Generating a further n * (n-1) + 1 special substrings.
	// But this is already accounted for by the earlier loop.
	// Here, only take into account the contribution from such "palindromic" strings itself
	char prevprevChar = ' ';
	int prevprevCharCount = 0;
	bool isMatching = false;

	prevChar = ' ';
	prevCharCount = 0;

	for (int i = 0; i < s.length(); i++)
	{
		if (s[i] != prevChar)
		{
			if (s[i] == prevprevChar)
			{
				isMatching = true;
			}
			else
			{
				if (isMatching)
				{
					isMatching = false;
					numSubstrings += min(prevprevCharCount, prevCharCount);
				}

				prevprevChar = prevChar;
				prevprevCharCount = prevCharCount;
			}

			prevChar = s[i];
			prevCharCount = 1;
		}
		else
		{
			prevCharCount += 1;
		}
	}

	if (isMatching)
	{
		numSubstrings += min(prevprevCharCount, prevCharCount);
	}

	return numSubstrings;
}
