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

// Lab 2 user marquee selection task
extern float projectionMatrix[];
float modelViewMatrix[16]; // store opengl's mv matrix at time of selection
float cyanmat_diffuse[] = { 0, 1, 1, 1.0f };   // selected triangles highlight
float greenmat_diffuse[] = { 0, 1, 0.5, 1.0f }; // selected triangles border highlight
extern double selectBoxCoords[];
extern bool isSelecting, isDeselecting;
// for isSelectingFacing mode (ctrl-alt click drag)
extern bool isSelectingFacing, m_Smooth;


void myObjType::draw() {

	glEnable(GL_LIGHTING);

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

	double longestSide = 0.0;
	for (int i = 0; i < 3; i++)
		if ((lmax[i] - lmin[i]) > longestSide)
			longestSide = (lmax[i] - lmin[i]);
	double scaleFactor = 4.0 / longestSide;
	glScalef(scaleFactor, scaleFactor, scaleFactor);
	glTranslated(-(lmin[0] + lmax[0]) / 2.0, -(lmin[1] + lmax[1]) / 2.0, -(lmin[2] + lmax[2]) / 2.0);

	// Lab 2 user marquee selection task
	// Store mv matrix at this point (after glScalef, glTranslated, ...) for getting transformed vertex pos later
	if (isSelecting)
	{
		glGetFloatv(GL_MODELVIEW_MATRIX, modelViewMatrix);
	}

	// Normal render run
	for (int i = 1; i <= tcount; i++)
	{
		glBegin(GL_POLYGON);

		// uncomment the following after you computed the normals - done
		// recommented for lab 2 vertex normals optional task (moved to below glNormal3dv)
		// glNormal3dv(nlist[i]);    

		// Lab 2 user marquee selection task - highlight triangle
		if (selectedT.test(i))
			glMaterialfv(GL_FRONT, GL_DIFFUSE, cyanmat_diffuse);

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

		// Lab 2 user marquee selection task - highlight triangle - reset colour
		if (selectedT.test(i))
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			
		glEnd();

		glMaterialfv(GL_FRONT, GL_DIFFUSE, rmat_diffuse);
		glPolygonOffset(10000, 10000);
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
			// Lab 2 user marquee selection task - highlight triangle border colour
			else if (selectedT.test(i))
			{
				glMaterialfv(GL_FRONT, GL_DIFFUSE, greenmat_diffuse);
				glBegin(GL_LINES);
					glVertex3dv(vlist[tlist[i][j]]);
					glVertex3dv(vlist[tlist[i][(j + 1) % 3]]);
				glEnd();
				glMaterialfv(GL_FRONT, GL_DIFFUSE, rmat_diffuse);
			}
		}
		glLineWidth(1);
		glPolygonOffset(0, 0);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	}

	glDisable(GL_LIGHTING);
}

void myObjType::drawOffscreen()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glShadeModel(GL_FLAT);

	unsigned int r = 1; // start from 1, as background is black
	unsigned int g = 0;
	unsigned int b = 0;

	unordered_map<unsigned int, int> colourTriangleMap;

	for (int i = 1; i <= tcount; i++)
	{
		glColor3ub((byte)r, (byte)g, (byte)b);
		colourTriangleMap.insert({
			((r << 24) | (g << 16) | (b << 8)) & 0xFFFFFF00, // ensure last 8 bits are 0
			i
			});

		glBegin(GL_POLYGON);

		for (int j = 0; j < 3; j++)
		{
			glVertex3dv(vlist[tlist[i][j]]);
		}

		glEnd();

		b += 1;
		if (b == 256)
		{
			if (g == 255)
			{
				r += 1;
				g = 0;
				b = 0;
			}
			else
			{
				g += 1;
				b = 0;
			}
		}
	}

	glReadBuffer(GL_BACK);

	int width = glutGet(GLUT_WINDOW_WIDTH);
	int height = glutGet(GLUT_WINDOW_HEIGHT);
	unsigned int size = width * height * 3;

	byte* readColours = (byte*)malloc(sizeof(byte) * size);
	glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, readColours);

	visibleT.reset();
	for (int i = 0; i < size; i += 3)
	{
		unsigned int fullColour = ((readColours[i] << 24) | (readColours[i + 1] << 16) | (readColours[i + 2] << 8)) & 0xFFFFFF00;
		if (colourTriangleMap.find(fullColour) != colourTriangleMap.end())
		{
			visibleT.set(colourTriangleMap.at(fullColour));
		}
	}

	free(readColours);

	if (m_Smooth)
		glShadeModel(GL_SMOOTH);
}

void multMatrix(float transformMatrix[16], double v[3], double vTransformed[3])
{
	double homogenousFactor = transformMatrix[3] * v[0]
		+ transformMatrix[7] * v[1]
		+ transformMatrix[11] * v[2]
		+ transformMatrix[15] * 1;

	vTransformed[0] = (transformMatrix[0] * v[0]
		+ transformMatrix[4] * v[1]
		+ transformMatrix[8] * v[2]
		+ transformMatrix[12] * 1) / homogenousFactor; // 1 because vlist is in non-homogenous coordinates

	vTransformed[1] = (transformMatrix[1] * v[0]
		+ transformMatrix[5] * v[1]
		+ transformMatrix[9] * v[2]
		+ transformMatrix[13] * 1) / homogenousFactor;

	// Can disregard z coordinate as object is rotated, not the camera
	vTransformed[2] = (transformMatrix[2] * v[0]
		+ transformMatrix[6] * v[1]
		+ transformMatrix[10] * v[2]
		+ transformMatrix[14] * 1) / homogenousFactor;
}

// Lab 2 Optional task - user select marquee
void myObjType::computeSelectedTriangles()
{
	/*cout << "Select box coords x1 y1 x2 y2: " 
		<< selectBoxCoords[0] << " "
		<< selectBoxCoords[1] << " "
		<< selectBoxCoords[3] << " "
		<< selectBoxCoords[4] << " " << endl;*/

	double selectBoxTopleft[3];
	double selectBoxBottomRight[3];
	multMatrix(projectionMatrix, selectBoxCoords, selectBoxTopleft);
	multMatrix(projectionMatrix, selectBoxCoords + 3, selectBoxBottomRight);

	/*cout << "Select box coords x1 y1 x2 y2 new: "
		<< selectBoxTopleft[0] << " "
		<< selectBoxTopleft[1] << " "
		<< selectBoxBottomRight[0] << " "
		<< selectBoxBottomRight[1] << " " << endl;*/

	// Update bitset
	for (int i = 1; i <= tcount; i++)
	{
		if (isSelectingFacing && !visibleT.test(i)) // "alt" mode
		{
			continue;
		}

		double vTransformed[3][3];

		// Transform each triangle's vertices using opengl's model view matrix.
		for (int j = 0; j < 3; j++)
		{
			double vTransformedTemp[3];
			multMatrix(modelViewMatrix, vlist[tlist[i][j]], vTransformedTemp);
			multMatrix(projectionMatrix, vTransformedTemp, vTransformed[j]);

			// cout << j << ": " << vTransformed[j][0] << " " << vTransformed[j][1] << " " << vTransformed[j][2] << endl;
		}

		// Then perform cohen-sutherland intersection (simplified as no need to clip)
		// with the selection box **for each edge** in the triangle.
		// If any edge tests true, the triangle should be selected
		for (int j = 0; j < 3; j++)
		{
			double v1[2];
			v1[0] = vTransformed[j][0];
			v1[1] = vTransformed[j][1];
			double v2[2];
			v2[0] = vTransformed[(j + 1) % 3][0];
			v2[1] = vTransformed[(j + 1) % 3][1];

		restart: // cohen sutherland repeat point for case 4 - reset outcodes but keep shortened line
			unsigned int v1Outcode = 0;
			unsigned int v2Outcode = 0;

			// v1
			if (v1[1] > selectBoxTopleft[1])
				v1Outcode |= 0b00001000; // b0
			else if (v1[1] < selectBoxBottomRight[1])
				v1Outcode |= 0b00000100; // b1

			if (v1[0] > selectBoxBottomRight[0])
				v1Outcode |= 0b00000010; // b2
			else if (v1[0] < selectBoxTopleft[0])
				v1Outcode |= 0b00000001; // b3

			// v2
			if (v2[1] > selectBoxTopleft[1])
				v2Outcode |= 0b00001000; // b0
			else if (v2[1] < selectBoxBottomRight[1])
				v2Outcode |= 0b00000100; // b1

			if (v2[0] > selectBoxBottomRight[0])
				v2Outcode |= 0b00000010; // b2
			else if (v2[0] < selectBoxTopleft[0])
				v2Outcode |= 0b00000001; // b3

			if (v1Outcode == 0 && v2Outcode == 0)
			{
				// cout << "Both zero " << v1[0] << " " << v1[1] << " " << (int)v1Outcode << " " << v2[0] << " " << v2[1] << endl;

				// line within - accept
				selectedT.set(i, !isDeselecting);
				break;
			}

			if ((v1Outcode & v2Outcode) != 0)
			{
				// both lines outside on the same side - reject, continue checking other edges
				// cout << "Both on same side" << endl;
				continue;
			}

			if ((v1Outcode == 0 && v2Outcode != 0) || (v1Outcode != 0 && v2Outcode == 0))
			{
				// cout << "One zero " << v1[0] << " " << v1[1] << " " << v2[0] << " " << v2[1] << endl;

				// one inside, one outside - normally intersection needs to be calculated here
				// but we are not doing clipping - accept
				selectedT.set(i, !isDeselecting);
				break;
			}

			// Last case - still tricky - neither 0, but bitwise AND yields 0
			// Try shortening v2 into the window by intersecting with one side

			double v1MinusV2[2]; // delta
			v1MinusV2[0] = v1[0] - v2[0];
			v1MinusV2[1] = v1[1] - v2[1];

			/*cout << "Cohen sutherland last case " << j << " v2Outcode: " << v2Outcode
				<< " v2: " << v2[0] << " " << v2[1] << " "
				<< " v1: " << v1[0] << " " << v1[1] << " "
				<< " v1MinusV2: " << v1MinusV2[0] << " " << v1MinusV2[1] << endl;*/

			if ((v2Outcode & 0b00001000) != 0) // top
			{
				v2[0] += (v2[1] - selectBoxTopleft[1]) * v1MinusV2[0] / abs(v1MinusV2[1]);
				v2[1] = selectBoxTopleft[1];
			}
			else if ((v2Outcode & 0b00000100) != 0) // down
			{
				v2[0] += (selectBoxBottomRight[1] - v2[1]) * v1MinusV2[0] / abs(v1MinusV2[1]);
				v2[1] = selectBoxBottomRight[1];
			}
			else if ((v2Outcode & 0b00000010) != 0) // right
			{
				v2[1] += (v2[0] - selectBoxBottomRight[0]) * v1MinusV2[1] / abs(v1MinusV2[0]);
				v2[0] = selectBoxBottomRight[0];
			}
			else if ((v2Outcode & 0b00000001) != 0) // left
			{
				v2[1] += (selectBoxTopleft[0] - v2[0]) * v1MinusV2[1] / abs(v1MinusV2[0]);
				v2[0] = selectBoxTopleft[0];
			}

			// cout << "                           new " << v2[0] << " " << v2[1] << endl;

			goto restart; // rerun outcode check
		}

		// Last edge case check cohen sutherland cannot find
		// Selection box is entirely within the polygon
		// Logic referenced from
		// https://math.stackexchange.com/questions/51326/determining-if-an-arbitrary-point-lies-inside-a-triangle-defined-by-three-points
		// Test the first point of the selection box
		if (isDeselecting ? selectedT.test(i) : !selectedT.test(i))
		{
			double AB[2];
			AB[0] = vTransformed[1][0] - vTransformed[0][0];
			AB[1] = vTransformed[1][1] - vTransformed[0][1];
			double BC[2];
			BC[0] = vTransformed[2][0] - vTransformed[1][0];
			BC[1] = vTransformed[2][1] - vTransformed[1][1];
			double CA[2];
			CA[0] = vTransformed[0][0] - vTransformed[2][0];
			CA[1] = vTransformed[0][1] - vTransformed[2][1];

			double AP[2];
			AP[0] = selectBoxTopleft[0] - vTransformed[0][0];
			AP[1] = selectBoxTopleft[1] - vTransformed[0][1];
			double BP[2];
			BP[0] = selectBoxTopleft[0] - vTransformed[1][0];
			BP[1] = selectBoxTopleft[1] - vTransformed[1][1];
			double CP[2];
			CP[0] = selectBoxTopleft[0] - vTransformed[2][0];
			CP[1] = selectBoxTopleft[1] - vTransformed[2][1];

			// Just get third term as z position does not matter after projection for checking this
			double ABAP = AB[0] * AP[1] - AP[0] * AB[1];
			double BABP = BC[0] * BP[1] - BP[0] * BC[1];
			double CACP = CA[0] * CP[1] - CP[0] * CA[1];

			if ((ABAP >= 0 && BABP >= 0 && CACP >= 0)
				|| (ABAP < 0 && BABP < 0 && CACP < 0))
			{
				selectedT.set(i, !isDeselecting);
			}
		}
	}
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

void myObjType::readFile(char* filename)
{
	reset();

	int len = strlen(filename);
	if (filename[len - 3] == '3' && filename[len - 2] == 'd' && filename[len - 1] == 's')
	{
		read3dsFile(filename);
	}
	else
	{
		readObjFile(filename);
	}
}

void myObjType::reset()
{
	for (int i = 0; i <= vcount; i++)
		vToTList[i] = list<int>();
	for (int i = 0; i < 18; i++)
	{
		statMinAngle[i] = 0;
		statMaxAngle[i] = 0;
	}

	tcount = vcount = 0;
	numCc = 0;
	selectedT.reset();
	visibleT.reset();
}

void myObjType::read3dsFile(char* filename)
{
	cout << "Opening 3ds file: " << filename << endl;
	FILE* fp;
	fopen_s(&fp, filename, "rb");
	if (fp == NULL) {
		cout << "We cannot find your file " << filename << endl;
		exit(1);
	}

	unsigned short id = 0;
	unsigned int len = 0;

	while (ftell(fp) != -1)
	{
		id = 0; len = 0;
		fread((char*)&id, 1, 2, fp);
		fread((char*)&len, 1, 4, fp);

		if (id == 0x4D4D // main chunk
			|| id == 0x3D3D // 3d editor chunk
			|| id == 0x4100) // triangle mesh
		{
			// cout << "Entering chunk " << id << endl;
			continue;
		}
		else if (id == 0x4000)
		{
			// object chunk has an ASCII name to skip
			// http://paulbourke.net/dataformats/3ds/
			char name[64];
			fread(name, 1, 64, fp);

			int endIdx = 0;
			for (; endIdx < 64; endIdx++)
			{
				if (name[endIdx] == '\0')
				{
					endIdx++;
					break;
				}
			}

			// cout << "Object chunk name was n characters long: " << endIdx << " ftell " << ftell(fp) - 64 << endl;

			fseek(fp, endIdx - 64, SEEK_CUR);
		}
		else if (id == 0x4110) // vertices list
		{
			bool firstVertex = 1;

			unsigned short numVertices = 0;
			fread((char*)&numVertices, 1, 2, fp);

			float x = 0;
			float y = 0;
			float z = 0;

			for (short i = 0; i < numVertices; i++)
			{
				vcount++;

				fread((char*)&x, 1, 4, fp);
				vlist[vcount][0] = x;

				fread((char*)&y, 1, 4, fp);
				vlist[vcount][1] = y;

				fread((char*)&z, 1, 4, fp);
				vlist[vcount][2] = z;

				/*cout << "Reading triangle " << i << " x y z "
					<< vlist[vcount][0] << " " << vlist[vcount][1] << " " << vlist[vcount][2] << endl;*/

				for (int k = 0; k < 3; k++) {
					if (firstVertex)
					{
						lmin[k] = lmax[k] = vlist[vcount][k];
					}
					else
					{
						if (lmin[k] > vlist[vcount][k])
							lmin[k] = vlist[vcount][k];
						if (lmax[k] < vlist[vcount][k])
							lmax[k] = vlist[vcount][k];
					}
				}

				firstVertex = 0;
			}
		}
		else if (id == 0x4120) // face list
		{
			unsigned short numPoly = 0;
			fread((char*)&numPoly, 1, 2, fp);

			unsigned short v1Num, v2Num, v3Num;

			for (int i = 0; i < numPoly; i++)
			{
				tcount++;

				// cout << "Face number " << tcount << endl;

				fread((char*)&v1Num, 1, 2, fp);
				tlist[tcount][0] = ((int)v1Num) + 1; // + 1 as 3ds vertices are 0 indexed, cast to int first in case of short overflow

				fread((char*)&v2Num, 1, 2, fp);
				tlist[tcount][1] = ((int)v2Num) + 1;

				fread((char*)&v3Num, 1, 2, fp);
				tlist[tcount][2] = ((int)v3Num) + 1;

				fseek(fp, 2, SEEK_CUR);

				// lab 2 vertex normal calculation optional task
				vToTList[tlist[tcount][0]].push_back(tcount);
				vToTList[tlist[tcount][1]].push_back(tcount);
				vToTList[tlist[tcount][2]].push_back(tcount);
			}

			break;
		}
		else
		{
			// cout << "skipping id " << id << " len " << len << endl;

			long seekAmt = len - 6;
			// cout << seekAmt << " ftell: " << ftell(fp) << endl;
			fseek(fp, seekAmt, SEEK_CUR);
		}
	}

	fclose(fp);

	postReadFile();
}

void myObjType::readObjFile(char* filename)
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
					vlist[vcount][k] = atof(line.substr(i, j - i).c_str());
					if (firstVertex) 
						lmin[k] = lmax[k] = vlist[vcount][k];
					else {
						if (lmin[k] > vlist[vcount][k])
							lmin[k] = vlist[vcount][k];
						if (lmax[k] < vlist[vcount][k])
							lmax[k] = vlist[vcount][k];
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

	postReadFile();
}

void myObjType::postReadFile()
{
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
	
	unordered_set<int> unprocessedTriangles;
	for (int i = 1; i <= tcount; i++)
		unprocessedTriangles.insert(i);

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
	for (int t = 1; t <= tcount; t++)
	{
		for (int v = 0; v < 6; v++)
		{
			pair<int, int> currentEdge = v < 3
				? make_pair(tlist[t][v], tlist[t][(v + 1) % 3])
				: make_pair(tlist[t][(v + 1) % 3], tlist[t][v % 3]);

			if (edgeTriangleMap.find(currentEdge) == edgeTriangleMap.end())
			{
				OrTri* initialArrayPair = (OrTri*)malloc(sizeof(OrTri) * 2);
				initialArrayPair[0] = 0; initialArrayPair[1] = 0;
				edgeTriangleMap.insert({ currentEdge, initialArrayPair });
			}

			OrTri* orTriPair = edgeTriangleMap.at(currentEdge);
			orTriPair[orTriPair[0] == 0 ? 0 : 1] = makeOrTri(t, v);
		}
	}

	// Build fnlist...
	for (int t = 1; t <= tcount; t++)
	{
		// 012, 120, 201
		for (int v = 0; v < 3; v++)
		{
			pair<int, int> currentEdge = make_pair(tlist[t][v], tlist[t][(v + 1) % 3]);

			if (edgeTriangleMap.find(currentEdge) == edgeTriangleMap.end())
			{
				throw exception("Tried to build fnlist with unprocessed edge");
			}

			OrTri* orTriPair = edgeTriangleMap.at(currentEdge);
			OrTri current = makeOrTri(t, v);
			if (orTriPair[0] != 0 && orTriPair[0] != current)
			{
				fnlist[t][v] = orTriPair[0];
			}
			else if (orTriPair[1] != 0 && orTriPair[1] != current)
			{
				fnlist[t][v] = orTriPair[1];
			}
			else // back to itself
			{
				fnlist[t][v] = makeOrTri(t, v);
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

			nlist[tcount][0] = nlist[t][0];
			nlist[tcount][1] = nlist[t][1];
			nlist[tcount][2] = nlist[t][2];

			for (int j = 0; j < 3; j++)
				vToTList[tlist[tcount][j]].push_back(tcount);

			selectedT.set(tcount, true);
		}

		// Readjust current triangle to be central triangle
		for (int i = 0; i < 3; i++)
		{
			vToTList[tlist[t][i]].remove(t);                    // remove old v -> t mapping
			vToTList[newMidpointVertices[i]].push_back(tcount); // new v -> t mapping
		}

		tlist[t][0] = newMidpointVertices[0];
		tlist[t][1] = newMidpointVertices[1];
		tlist[t][2] = newMidpointVertices[2];

		// For each version of current triangle
		// If it is a boundary edge (on selected portion), 
		// also need to partially subdivide the adjacent triangle
		for (int v = 0; v < 3; v++)
		{
			OrTri adjacentTriangle = fnext(makeOrTri(t, v)); // triangle adjacent to version i of current triangle
			int adjacentTriangleVer = ver(adjacentTriangle) % 3; // normalized version
			int adjacentTriangleIdx =  idx(adjacentTriangle);

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
		for (;  vertexTrianglesIt != vertexTriangles.end(); vertexTrianglesIt++)
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

inline double* computeAngles(double* v1, double* v2, double* v3)
{
	double v1ToV2Len = sqrt(abs(pow(v2[0] - v1[0], 2) + pow(v2[1] - v1[1], 2) + pow(v2[2] - v1[2], 2)));
	double v1ToV3Len = sqrt(abs(pow(v3[0] - v1[0], 2) + pow(v3[1] - v1[1], 2) + pow(v3[2] - v1[2], 2)));
	double v2ToV3Len = sqrt(abs(pow(v3[0] - v2[0], 2) + pow(v3[1] - v2[1], 2) + pow(v3[2] - v2[2], 2)));

	// cosine rule
	double angles[3];
	angles[0] = acos(
		(pow(v1ToV2Len, 2) + pow(v1ToV3Len, 2) - pow(v2ToV3Len, 2)) / (2 * v1ToV2Len * v1ToV3Len));
	angles[1] = acos(
		(pow(v1ToV2Len, 2) + pow(v2ToV3Len, 2) - pow(v1ToV3Len, 2)) / (2 * v1ToV2Len * v2ToV3Len));
	angles[2] = acos(
		(pow(v2ToV3Len, 2) + pow(v1ToV3Len, 2) - pow(v1ToV2Len, 2)) / (2 * v2ToV3Len * v1ToV3Len));

	return angles;
}

void myObjType::computeStat()
{
	int i;
    double minAngle = 181;
    double maxAngle = 0;

	// Lab 1 Task 2 start
	for (int i = 1; i <= tcount; i++)
	{
		int* t = tlist[i];
		double* v1 = vlist[t[0]];
		double* v2 = vlist[t[1]];
		double* v3 = vlist[t[2]];

		double* angles = computeAngles(v1, v2, v3);

		double currTriMinAngle = min(angles[0], min(angles[1], angles[2])) * 180.0 / M_PI;
		double currTriMaxAngle = max(angles[0], max(angles[1], angles[2])) * 180.0 / M_PI;
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
