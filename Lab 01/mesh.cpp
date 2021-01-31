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
#include <unordered_map>;
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "mesh.h"
#include <map>
#include <queue>
#include <iomanip>
using namespace std;



void myObjType::draw() {

	glEnable(GL_LIGHTING);

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

	glPushMatrix();
	double longestSide = 0.0;
	for (int i = 0; i < 3; i++)
		if ((lmax[i] - lmin[i]) > longestSide)
			longestSide = (lmax[i] - lmin[i]);
	glScalef(4.0 / longestSide, 4.0 / longestSide, 4.0 / longestSide);
	glTranslated(-(lmin[0] + lmax[0]) / 2.0, -(lmin[1] + lmax[1]) / 2.0, -(lmin[2] + lmax[2]) / 2.0);
	for (int i = 1; i <= tcount; i++)
	{
		glBegin(GL_POLYGON);
// uncomment the following after you computed the normals - done
		glNormal3dv(nlist[i]);    
		for (int j = 0; j < 3; j++)
			glVertex3dv(vlist[tlist[i][j]]);
		glEnd();
	
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
					i = j;
					fnlist[tcount][k] = 0;
					while (linec[j] != ' ') j++;
				}
			}
		}
	}

	// We suggest you to compute the normals here

	// Task 1 Normal Computation
	computeNormals();
	// Task 1 End

	// Lab 2 fnext task
	computeFnlist();

    cout << "No. of vertices: " << vcount << endl;
    cout << "No. of triangles: " << tcount << endl;
    computeStat();
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

		double v2MinusV1x = v2[0] - v1[0];
		double v2MinusV1y = v2[1] - v1[1];
		double v2MinusV1z = v2[2] - v1[2];

		double v3MinusV1x = v3[0] - v1[0];
		double v3MinusV1y = v3[1] - v1[1];
		double v3MinusV1z = v3[2] - v1[2];

		double normalX = (v2MinusV1y * v3MinusV1z) - (v3MinusV1y * v2MinusV1z);
		double normalY = (v3MinusV1x * v2MinusV1z) - (v2MinusV1x * v3MinusV1z);
		double normalZ = (v2MinusV1x * v3MinusV1y) - (v3MinusV1x * v2MinusV1y);
		double magnitude = sqrt(pow(normalX, 2) + pow(normalY, 2) + pow(normalZ, 2));

		nlist[i][0] = normalX;
		nlist[i][1] = normalY;
		nlist[i][2] = normalZ;
	}
}

void myObjType::computeFnlist()
{
	// Build edge triangle map
	unordered_map<string, OrTri*> edgeTriangleMap;
	for (int i = 1; i <= tcount; i++)
	{
		for (int v = 0; v < 6; v++)
		{
			string currentEdge = v < 3
				? tlist[i][v] + tlist[i][(v + 1) % 3] + ""
				: tlist[i][(v + 1) % 3] + tlist[i][v % 3] + ""; // TODO probably not needed

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
			string currentEdge = tlist[i][v] + tlist[i][(v + 1) % 3] + "";

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

void myObjType::computeStat()
{
	int i;
    double minAngle = 0;
    double maxAngle = 0;

	// Task 2 start

	for (int i = 1; i <= tcount; i++)
	{
		int* t = tlist[i];
		double* v1 = vlist[t[0]];
		double* v2 = vlist[t[1]];
		double* v3 = vlist[t[2]];

		double v1ToV2Len = abs(v2[0] - v1[0] + v2[1] - v1[1] + v2[2] - v1[2]);
		double v1ToV3Len = abs(v3[0] - v1[0] + v3[1] - v1[1] + v3[2] - v1[2]);
		double v2ToV3Len = abs(v3[0] - v2[0] + v3[1] - v2[1] + v3[2] - v2[2]);

		// cosine rule
		double v1Angle = acos(
			(pow(v1ToV2Len, 2) + pow(v1ToV3Len, 2) - pow(v2ToV3Len, 2)) / (2 * v1ToV2Len * v1ToV3Len));
		double v2Angle = acos(
			(pow(v1ToV2Len, 2) + pow(v2ToV3Len, 2) - pow(v1ToV3Len, 2)) / (2 * v1ToV2Len * v2ToV3Len));
		double v3Angle = acos(
			(pow(v2ToV3Len, 2) + pow(v1ToV3Len, 2) - pow(v1ToV2Len, 2)) / (2 * v2ToV3Len * v1ToV3Len));

		double currTriMinAngle = min(v1Angle, min(v2Angle, v3Angle)) * 180.0 / M_PI;
		double currTriMaxAngle = max(v1Angle, max(v2Angle, v3Angle)) * 180.0 / M_PI;

		int minBucketIdx = currTriMinAngle / 10;
		int maxBucketIdx = currTriMaxAngle / 10;
		statMinAngle[minBucketIdx] += 1;
		statMaxAngle[maxBucketIdx] += 1;

		minAngle = min(minAngle, currTriMinAngle);
		maxAngle = max(minAngle, currTriMaxAngle);
	}

	// Task 2 end
    
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
}

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

inline OrTri myObjType::fnext(OrTri t)
{
	int v = ver(t);
	return v < 3
		? fnlist[idx(t)][v]
		: sym(fnlist[idx(t)][v % 3]);
}
