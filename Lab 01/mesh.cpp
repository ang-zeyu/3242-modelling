#include "mesh.h"

#include "utils.cpp";

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
#include <vector>
#include <list>
#include <map>
#include <queue>
#include <iomanip>
using namespace std;

extern float mat_diffuse[];
float rmat_diffuse[] = { 1, 0, 0, 1.0f }; // Lab 2 boundary edge visualisation optional task

extern int startx, starty, currX, currY;
extern bool isSelecting;

void myObjType::draw() {

	glEnable(GL_LIGHTING);

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

	glPushMatrix();

	double longestSide = 0.0;
	for (int i = 0; i < 3; i++)
		if ((lmax[i] - lmin[i]) > longestSide)
			longestSide = (lmax[i] - lmin[i]);
	double scaleFactor = 4.0 / longestSide;
	glScalef(scaleFactor, scaleFactor, scaleFactor);
	glTranslated(-(lmin[0] + lmax[0]) / 2.0, -(lmin[1] + lmax[1]) / 2.0, -(lmin[2] + lmax[2]) / 2.0);

	for (int i = 1; i <= tcount; i++)
	{
		glBegin(GL_POLYGON);
// uncomment the following after you computed the normals - done
		// glNormal3dv(nlist[i]);    

		for (int j = 0; j < 3; j++)
		{
			glNormal3dv(vnlist[tlist[i][j]]); // Moved down here, use vnlist for lab 2 optional task (vertex normals)

			// Lab 2 boundary edge visualisation optional task
			/*if (fnlist[i][j] == makeOrTri(i, j))
			{
				glMaterialfv(GL_FRONT, GL_DIFFUSE, rmat_diffuse);
			}
			else
			{
				glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			}*/

			glVertex3dv(vlist[tlist[i][j]]);
		}
			
		glEnd();

		glMaterialfv(GL_FRONT, GL_DIFFUSE, rmat_diffuse);
		glPolygonOffset(1000, 1000);
		glLineWidth(3);
		for (int j = 0; j < 3; j++)
		{
			// Lab 2 boundary edge visualisation optional task
			if (fnlist[i][j] == makeOrTri(i, j))
			{
				glBegin(GL_LINES);

				glVertex3dv(vlist[tlist[i][j]]);
				glVertex3dv(vlist[tlist[i][(j + 1) % 3]]);

				glEnd();
			}
		}
		glLineWidth(1);
		glPolygonOffset(0, 0);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	}
	glDisable(GL_LIGHTING);

	glPopMatrix();
}

void myObjType::writeFile(char* filename)
{
	// Task 3 start

	cout << "Writing to " << filename << endl;
	ofstream outFile;
	outFile.open(filename);
	
	// Vertices
	for (int i = 1; i <= vcount; i++)
	{
		outFile
			<< "v " + to_string(vlist[i][0]) + " " + to_string(vlist[i][1]) + " " + to_string(vlist[i][2])
			<< endl;
	}

	// Faces / triangles
	for (int i = 1; i <= tcount; i++)
	{
		outFile
			<< "f " + to_string(tlist[i][0]) + " " + to_string(tlist[i][1]) + " " + to_string(tlist[i][2])
			<< endl;
	}

	outFile.close();
	cout << "Done writing to " << filename << endl;

	// Task 3 end
}

void myObjType::read3dsFile(char* filename)
{
	cout << "Opening " << filename << endl;
	ifstream inFile;
	inFile.open(filename);
	if (!inFile.is_open()) {
		cout << "We cannot find your file " << filename << endl;
		exit(1);
	}

	string line;
	int i, j;
	bool firstVertex = 1;
	double currCood;

	while (getline(inFile, line))
	{
		if ((line[0] == 'v' || line[0] == 'f') && line[1] == ' ')
		{
			if (line[0] == 'v')
			{
				vcount++;
				i = 1;
				const char* linec = line.data();
				for (int k = 0; k < 3; k++) { // k is 0,1,2 for x,y,z
					while (linec[i] == ' ') i++;
					j = i;
					while (linec[j] != ' ') j++;
					currCood = vlist[vcount][k] = atof(line.substr(i, j - i).c_str());
					if (firstVertex)
						lmin[k] = lmax[k] = currCood;
					else {
						if (lmin[k] > currCood)
							lmin[k] = currCood;
						if (lmax[k] < currCood)
							lmax[k] = currCood;
					}
					i = j;
				}

				firstVertex = 0;
			}

			if (line[0] == 'f')
			{
				tcount++;
				i = 1;
				const char* linec = line.data();
				for (int k = 0; k < 3; k++) {
					while (linec[i] == ' ') i++;
					j = i;
					while (linec[j] != ' ' && linec[j] != '\\') j++;
					tlist[tcount][k] = atof(line.substr(i, j - i).c_str());
					i = j;
					fnlist[tcount][k] = 0;
					while (linec[j] != ' ') j++;
				}
			}
		}
	}

	// Task 1 Normal Computation
	computeNormals();
	// Task 1 End

	cout << "No. of vertices: " << vcount << endl;
	cout << "No. of triangles: " << tcount << endl;
	computeStat();
}

void myObjType::readFile(char* filename)
{
	cout << "Opening " << filename << endl;
	ifstream inFile;
	inFile.open(filename);
	if (!inFile.is_open()) {
		cout << "We cannot find your file " << filename << endl;
		exit(1);
	}

	string line;
	int i, j;
	bool firstVertex = 1;
	double currCood;

	while (getline(inFile, line))
	{
		if ((line[0] == 'v' || line[0] == 'f') && line[1] == ' ')
		{
			if (line[0] == 'v')
			{
				vcount++;
				i = 1;
				const char* linec = line.data();
				for (int k = 0; k < 3; k++) { // k is 0,1,2 for x,y,z
					while (linec[i] == ' ') i++;
					j = i;
					while (linec[j] != ' ') j++;
					currCood = vlist[vcount][k] = atof(line.substr(i, j - i).c_str());
					if (firstVertex) 
						lmin[k] = lmax[k] = currCood;
					else {
						if (lmin[k] > currCood)
							lmin[k] = currCood;
						if (lmax[k] < currCood)
							lmax[k] = currCood;
					}
					i = j;
				}

				firstVertex = 0;
			}

			if (line[0] == 'f')
			{
				tcount++;
				i = 1;
				const char* linec = line.data();
				for (int k = 0; k < 3; k++) {
					while (linec[i] == ' ') i++;
					j = i;
					while (linec[j] != ' ' && linec[j] != '\\') j++;
					tlist[tcount][k] = atof(line.substr(i, j - i).c_str());

					vToTList[tlist[tcount][k]].push_back(tcount); // lab 2 vertex normal calculation optional task
					i = j;
					fnlist[tcount][k] = 0;
					// while (linec[j] != ' ') j++;
				}

				// Support quadrilaterals
				while (linec[i] == ' ') i++;

				if (i < line.size())
				{
					tcount++;
					tlist[tcount][0] = tlist[tcount - 1][0];
					tlist[tcount][1] = tlist[tcount - 1][2];

					while (linec[i] == ' ') i++;
					j = i;
					while (linec[j] != ' ' && linec[j] != '\\') j++;
					tlist[tcount][2] = atof(line.substr(i, j - i).c_str());

					// lab 2 vertex normal calculation optional task
					vToTList[tlist[tcount][0]].push_back(tcount);
					vToTList[tlist[tcount][1]].push_back(tcount);
					vToTList[tlist[tcount][2]].push_back(tcount);
				}
			}
		}
	}

	// We suggest you to compute the normals here

	// Lab 2 fnext task
	computeFnlist();

	// Lab 2 Optional Task - orient triangles
	if (orientTriangles())
	{
		computeFnlist();
	}

	// Lab 1 Task 1 Normal Computation
	computeNormals();
	// Lab 1 Task 1 End

	// Lab 2 Optional Task - Vertex Normal calculation
	computeVertexNormals();

    cout << "No. of vertices: " << vcount << endl;
    cout << "No. of triangles: " << tcount << endl;
    computeStat();
}

// Lab 1 Task 1 Normal Computation
void crossProduct(double result[3], double v1[3], double v2[3])
{
	result[0] = v1[1] * v2[2] - v2[1] * v1[2];
	result[1] = v2[0] * v1[2] - v1[0] * v2[2];
	result[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

void myObjType::computeNormals()
{
	for (int i = 1; i <= tcount; i++)
	{
		int v1Idx = tlist[i][0];
		double* v1 = vlist[v1Idx];
		int v2Idx = tlist[i][1];
		double* v2 = vlist[v2Idx];
		int v3Idx = tlist[i][2];
		double* v3 = vlist[v3Idx];

		double v2MinusV1[3];
		v2MinusV1[0] = v2[0] - v1[0];
		v2MinusV1[1] = v2[1] - v1[1];
		v2MinusV1[2] = v2[2] - v1[2];

		double v3MinusV1[3];
		v3MinusV1[0] = v3[0] - v1[0];
		v3MinusV1[1] = v3[1] - v1[1];
		v3MinusV1[2] = v3[2] - v1[2];

		double crossProd[3];
		crossProduct(crossProd, v2MinusV1, v3MinusV1);
		double magnitude = sqrt(pow(crossProd[0], 2) + pow(crossProd[1], 2) + pow(crossProd[2], 2));

		nlist[i][0] = crossProd[0] / magnitude;
		nlist[i][1] = crossProd[1] / magnitude;
		nlist[i][2] = crossProd[2] / magnitude;
	}
}
// Lab 1 Task 1 Normal Computation End

// Lab 2 Optional Task - Vertex Normals Start
void myObjType::computeVertexNormals()
{
	if (tcount == 0)
	{
		cout << "No triangles to calculate vertex normals for" << endl;
		return;
	}

	// Algo: When reading file, construct vertex -> triangle inverse mapping
	// Read list of triangles for each vertex, compute vertex normal

	unordered_set<int> unprocessedTriangles;
	for (int i = 1; i <= tcount; i++)
		unprocessedTriangles.insert(i);

	for (int i = 1; i <= vcount; i++)
	{
		vnlist[i][0] = 0;
		vnlist[i][1] = 0;
		vnlist[i][2] = 0;

		for (const int idx : vToTList[i])
		{
			double* triangleNormal = nlist[idx];
			vnlist[i][0] += triangleNormal[0];
			vnlist[i][1] += triangleNormal[1];
			vnlist[i][2] += triangleNormal[2];
		}
	}

	/*
	 Attempted alternate algo (too slow):
	 - BFS, treating each triangle like a graph node in bfs
	 - For each triangle's vertices, fan the vertex, getting the normals of its triangles
	 - Average it

	while (!unprocessedTriangles.empty())
	{
		queue<int> triangleQueue; // bfs queue

		int tri = *(unprocessedTriangles.begin()); // any triangle still unprocessed
		triangleQueue.push(tri);
		unprocessedTriangles.erase(tri);

		unordered_set<int> unprocessedVertices;
		for (int i = 1; i <= vcount; i++)
			unprocessedVertices.insert(i);

		while (!triangleQueue.empty())
		{
			int currentTriangle = triangleQueue.front(); triangleQueue.pop();

			// Fanning portion - for each unprocessed vertex of the current triangle,
			// fan it and calculate the vertex normal
			OrTri currStartingTriangle = makeOrTri(currentTriangle, 0);
			for (int i = 0; i < 3; i++)
			{
				int currStartingVertex = tlist[currentTriangle][i];
				if (unprocessedVertices.find(currStartingVertex) == unprocessedVertices.end())
				{
					currStartingTriangle = enext(currStartingTriangle);
					continue;
				}
				unprocessedVertices.insert(currStartingVertex);

				int currStartingTriangleIdx = idx(currStartingTriangle);
				OrTri fanTriangle = currStartingTriangle;
				cout << currStartingTriangleIdx << endl;
				do
				{
					double* fanTriangleNormal = nlist[idx(fanTriangle)];
					vnlist[currStartingVertex][0] += fanTriangleNormal[0];
					vnlist[currStartingVertex][1] += fanTriangleNormal[1];
					vnlist[currStartingVertex][2] += fanTriangleNormal[2];

					// Fan: Fnext -> Sym -> Enext
					fanTriangle = enext(sym(fnext(fanTriangle)));

					cout << idx(fanTriangle) << " ";
				}
				while (idx(fanTriangle) != currStartingTriangleIdx);
				cout << endl;

				currStartingTriangle = enext(currStartingTriangle);
			}

			// Bfs portion - for each of the adjacent triangles...
			OrTri* fnTriangles = fnlist[currentTriangle];
			for (int i = 0; i < 3; i++)
			{
				int fnextTriangleIdx = idx(fnTriangles[i]);

				// Skip self-loops and triangles processed before due to traversal order
				if (unprocessedTriangles.find(fnextTriangleIdx) == unprocessedTriangles.end())
				{
					continue;
				}

				triangleQueue.push(fnextTriangleIdx);
				unprocessedTriangles.erase(fnextTriangleIdx);
			}
		}
	}
	*/

	// Normalize
	for (int i = 1; i <= vcount; i++)
	{
		double magnitude = sqrt(vnlist[i][0] * vnlist[i][0]
			+ vnlist[i][1] * vnlist[i][1]
			+ vnlist[i][2] * vnlist[i][2]);
		vnlist[i][0] /= magnitude;
		vnlist[i][1] /= magnitude;
		vnlist[i][2] /= magnitude;
	}

	cout << "Vertex normals calculated" << endl;
}
// Lab 2 Optional Task - Vertex Normals End

// Lab 2 Main task - fnext fnlist building
void myObjType::computeFnlist()
{
	// Build edge triangle map
	unordered_map<pair<int, int>, OrTri*, pairHash> edgeTriangleMap;
	for (int i = 1; i <= tcount; i++)
	{
		for (int v = 0; v < 6; v++)
		{
			pair<int, int> currentEdge = v < 3
				? make_pair(tlist[i][v], tlist[i][(v + 1) % 3])
				: make_pair(tlist[i][(v + 1) % 3], tlist[i][v % 3]);

			if (edgeTriangleMap.find(currentEdge) == edgeTriangleMap.end())
			{
				OrTri* initialArrayPair = (OrTri*)malloc(sizeof(OrTri) * 2);
				initialArrayPair[0] = 0; initialArrayPair[1] = 0;
				edgeTriangleMap.insert({ currentEdge, initialArrayPair });
			}

			OrTri* orTriPair = edgeTriangleMap.at(currentEdge);
			orTriPair[orTriPair[0] == 0 ? 0 : 1] = makeOrTri(i, v);
		}
	}

	// Build fnlist...
	for (int i = 1; i <= tcount; i++)
	{
		// 012, 120, 201
		for (int v = 0; v < 3; v++)
		{
			pair<int, int> currentEdge = make_pair(tlist[i][v], tlist[i][(v + 1) % 3]);

			if (edgeTriangleMap.find(currentEdge) == edgeTriangleMap.end())
			{
				throw exception("Tried to build fnlist with unprocessed edge");
			}

			OrTri* orTriPair = edgeTriangleMap.at(currentEdge);
			OrTri current = makeOrTri(i, v);
			if (orTriPair[0] != 0 && orTriPair[0] != current)
			{
				fnlist[i][v] = orTriPair[0];
			}
			else if (orTriPair[1] != 0 && orTriPair[1] != current)
			{
				fnlist[i][v] = orTriPair[1];
			}
			else // back to itself
			{
				fnlist[i][v] = makeOrTri(i, v);
			}
		}
	}
}

// Lab 2 Main task
inline OrTri myObjType::fnext(OrTri t)
{
	int v = ver(t);
	return v < 3
		? fnlist[idx(t)][v]
		: sym(fnlist[idx(t)][v - 3]);
}

// Lab 2 optional task
// Called in computeStat()
void myObjType::computeNumCc()
{
	numCc = 0;

	if (tcount == 0)
	{
		cout << "Num Connected Components = " << numCc << endl;
		return;
	}

	/*
	 My algo:
	 BFS, treating each triangle like a graph node in bfs
	 */

	unordered_set<int> unprocessedTriangles;
	for (int i = 1; i <= tcount; i++)
		unprocessedTriangles.insert(i);

	while (!unprocessedTriangles.empty())
	{
		queue<int> triangleQueue; // bfs queue

		int tri = *(unprocessedTriangles.begin()); // any triangle still unprocessed
		triangleQueue.push(tri);
		unprocessedTriangles.erase(tri);

		while (!triangleQueue.empty())
		{
			int currentTriangle = triangleQueue.front(); triangleQueue.pop();

			OrTri* fnTriangles = fnlist[currentTriangle];
			for (int i = 0; i < 3; i++)
			{
				int fnextTriangleIdx = idx(fnTriangles[i]);

				// Skip self-loops and triangles processed before due to traversal order
				if (unprocessedTriangles.find(fnextTriangleIdx) == unprocessedTriangles.end())
				{
					continue;
				}

				triangleQueue.push(fnextTriangleIdx);
				unprocessedTriangles.erase(fnextTriangleIdx);
			}
		}

		numCc += 1;
	}

	cout << "Num Connected Components = " << numCc << endl;
}

// Lab 2 optional task
// Called in readFile()
bool myObjType::orientTriangles()
{
	if (tcount == 0)
	{
		cout << "0 triangles, skipping orientTriangles" << endl;
		return false;
	}

	/*
	 My algo:
	 BFS, treating each triangle like a graph node in bfs
	 */
	bool successful = true;

	// Clone to reset if unsuccessful
	int** tlistClone = (int**)malloc(MAXT * 3 * sizeof(int));
	memcpy(tlistClone, tlist, MAXT * 3 * sizeof(int));

	unordered_set<int> unprocessedTriangles;
	for (int i = 1; i <= tcount; i++)
		unprocessedTriangles.insert(i);

	// Orient triangles for the whole mesh using bfs
	while (!unprocessedTriangles.empty())
	{
		queue<int> triangleQueue; // bfs queue

		int tri = *(unprocessedTriangles.begin()); // any triangle still unprocessed
		triangleQueue.push(tri);
		unprocessedTriangles.erase(tri);

		// Orient triangles for a surface
		while (!triangleQueue.empty())
		{
			int currentTriangle = triangleQueue.front(); triangleQueue.pop();

			int* currTriangleTlist = tlist[currentTriangle];
			unordered_set<pair<int, int>, pairHash> currTriangleEdges{
				{currTriangleTlist[0], currTriangleTlist[1]},
				{currTriangleTlist[1], currTriangleTlist[2]},
				{currTriangleTlist[2], currTriangleTlist[0]}
			};

			OrTri* fnTriangles = fnlist[currentTriangle];

			for (int i = 0; i < 3; i++)
			{
				int fnextTriangleIdx = idx(fnTriangles[i]);

				// Always skip self-loops
				if (fnextTriangleIdx == currentTriangle)
					continue;

				/*
				 Direction check:
				 - Two **different** faces are facing the same direction if
				   they DON'T share the same edge by right hand rule.
				 - Check if the version 0 to 2 of the current triangle
				   shares the edge (currTriangleEdges set) with the version 0 of the neighbouring triangle.
				 */
				int* fnextTriangleTlist = tlist[fnextTriangleIdx];
				bool sameDirection = currTriangleEdges.find({ fnextTriangleTlist[0], fnextTriangleTlist[1] }) == currTriangleEdges.end()
					&& currTriangleEdges.find({ fnextTriangleTlist[1], fnextTriangleTlist[2] }) == currTriangleEdges.end()
					&& currTriangleEdges.find({ fnextTriangleTlist[2], fnextTriangleTlist[0] }) == currTriangleEdges.end();

				// Skip triangles processed before due to traversal order
				if (unprocessedTriangles.find(fnextTriangleIdx) == unprocessedTriangles.end())
				{
					/*
					 But before skipping, conduct a non well-orientable check.
					 That is, if we hit a processed neighbour again, then upon doing the
					 direction check it is not the same direction, the mesh must be non well-orientable.
					 If it is non well-orientable, immediately terminate and goto end.
					 */
					if (!sameDirection)
					{
						successful = false;
						goto end;
					}

					continue;
				}

				// For newly encountered triangles,

				// Update the bfs queue and visited set
				triangleQueue.push(fnextTriangleIdx);
				unprocessedTriangles.erase(fnextTriangleIdx);
				
				// Invert the vertices if it is not the same direction
				if (!sameDirection)
				{
					// Swap tlist here only, fnlist is recomputed entirely if the whole swap is successful later (in readFile)
					int temp = fnextTriangleTlist[2];
					fnextTriangleTlist[2] = fnextTriangleTlist[0];
					fnextTriangleTlist[0] = temp;
				}
			}
		}
	}

end:
	if (!successful)
	{
		memcpy(tlist, tlistClone, MAXT * 3 * sizeof(int));
		cout << "Orient Triangles is unsuccessful" << endl;
	}
	else
	{
		cout << "Orient Triangles is successful" << endl;
	}

	free(tlistClone);

	return successful;
}

void myObjType::computeStat()
{
	int i;
    double minAngle = 0;
    double maxAngle = 0;

	// Lab 1 Task 2 start
	for (int i = 1; i <= tcount; i++)
	{
		int* t = tlist[i];
		double* v1 = vlist[t[0]];
		double* v2 = vlist[t[1]];
		double* v3 = vlist[t[2]];

		double v1ToV2Len = sqrt(abs(pow(v2[0] - v1[0], 2) + pow(v2[1] - v1[1], 2) + pow(v2[2] - v1[2], 2)));
		double v1ToV3Len = sqrt(abs(pow(v3[0] - v1[0], 2) + pow(v3[1] - v1[1], 2) + pow(v3[2] - v1[2], 2)));
		double v2ToV3Len = sqrt(abs(pow(v3[0] - v2[0], 2) + pow(v3[1] - v2[1], 2) + pow(v3[2] - v2[2], 2)));

		// cosine rule
		double v1Angle = acos(
			(pow(v1ToV2Len, 2) + pow(v1ToV3Len, 2) - pow(v2ToV3Len, 2)) / (2 * v1ToV2Len * v1ToV3Len));
		double v2Angle = acos(
			(pow(v1ToV2Len, 2) + pow(v2ToV3Len, 2) - pow(v1ToV3Len, 2)) / (2 * v1ToV2Len * v2ToV3Len));
		double v3Angle = acos(
			(pow(v2ToV3Len, 2) + pow(v1ToV3Len, 2) - pow(v1ToV2Len, 2)) / (2 * v2ToV3Len * v1ToV3Len));

		double currTriMinAngle = min(v1Angle, min(v2Angle, v3Angle)) * 180.0 / M_PI;
		double currTriMaxAngle = max(v1Angle, max(v2Angle, v3Angle)) * 180.0 / M_PI;
		// printf("%f %f %f %f %f\n", v1ToV2Len, v1ToV3Len, v2ToV3Len, min(v1Angle, min(v2Angle, v3Angle)), max(v1Angle, max(v2Angle, v3Angle)));

		int minBucketIdx = currTriMinAngle / 10;
		int maxBucketIdx = currTriMaxAngle / 10;
		statMinAngle[minBucketIdx] += 1;
		statMaxAngle[maxBucketIdx] += 1;

		minAngle = min(minAngle, currTriMinAngle);
		maxAngle = max(maxAngle, currTriMaxAngle);
	}
	// Lab 1 Task 2 end
    
    cout << "Min. angle = " << minAngle << endl;
    cout << "Max. angle = " << maxAngle << endl;

	cout << "Statistics for Maximum Angles" << endl;
	for (i = 0; i < 18; i++)
		cout << statMaxAngle[i] << " ";
	cout << endl;
	cout << "Statistics for Minimum Angles" << endl;
	for (i = 0; i < 18; i++)
		cout << statMinAngle[i] << " ";
	cout << endl;

	computeNumCc();
}

// Lab 2 Main task
inline int myObjType::org(OrTri t)
{
	/*
	 Versions 0 to 2:
	 // First 3 would just be the version number
	 // Latter 3 takes on the middle column (due to sym)
	 0 1 2
	 1 2 0
	 2 0 1
	*/
	int v = ver(t);
	return tlist
		[idx(t)]
		[v < 3 ? v : (v + 1) % 3];
}

// Lab 2 Main task
inline int myObjType::dest(OrTri t)
{
	/*
	 dest(ver 0) = dest(ver 3)
	 dest(ver 1) = dest(ver 4)
	 dest(ver 2) = dest(ver 5)
	*/
	int v = ver(t);
	return tlist
		[idx(t)]
		[(v + 2) % 3];
}
