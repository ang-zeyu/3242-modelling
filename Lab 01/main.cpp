// CS3241Lab1.cpp : Defines the entry point for the console application.


#include <iostream>
#include "mesh.h"

#ifdef _WIN32
#include <Windows.h>
#include "GL\glut.h"
#define M_PI 3.141592654
#elif __APPLE__
#include <OpenGL/gl.h>
#include <GLUT/GLUT.h>
#endif

using namespace std;




//#define M_PI 3.141592654
#define DEFAULT_WINDOW_LENGTH_WIDTH 600.0

myObjType myObj;


// global variable

bool m_Smooth = TRUE;
bool m_Highlight = FALSE;
GLfloat angle = 0;   /* in degrees */
GLfloat angle2 = 0;   /* in degrees */
GLfloat zoom = 1.0;
int mouseButton = 0;
int moving, startx, starty;

// For lab 2 optional task - user select marquee
int currX, currY;
double selectBoxCoords[6];
float projectionMatrix[16];
bool isSelecting = false;
bool isSelectingFacing = false;
bool triggerOffscreenDraw = false;
bool isDeselecting = false;

// Final boss - laplacian deformation
double laplacianDeformStepSize = 0.1;
bool isSelectingVertex = false;

float mat_diffuse[] = { 0.1f, 0.5f, 0.8f, 1.0f };

#define NO_OBJECT 4;
int current_object = 0;

using namespace std;

void setupLighting()
{
	glShadeModel(GL_SMOOTH);
	glEnable(GL_NORMALIZE);

	// Lights, material properties
    GLfloat	ambientProperties[]  = {0.7f, 0.7f, 0.7f, 1.0f};
	GLfloat	diffuseProperties[]  = {0.8f, 0.8f, 0.8f, 1.0f};
    GLfloat	specularProperties[] = {1.0f, 1.0f, 1.0f, 1.0f};
	GLfloat lightPosition[] = {-100.0f,100.0f,100.0f,1.0f};
	
    glClearDepth( 1.0 );

	glLightfv( GL_LIGHT0, GL_POSITION, lightPosition);
	
    glLightfv( GL_LIGHT0, GL_AMBIENT, ambientProperties);
    glLightfv( GL_LIGHT0, GL_DIFFUSE, diffuseProperties);
    glLightfv( GL_LIGHT0, GL_SPECULAR, specularProperties);
    glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 0.0);

	// Default : lighting
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);

}



void display(void)
{

	float mat_specular[] = { 0.3f, 0.3f, 0.3f, 1.0f };
	float mat_ambient[] = { 0.3f, 0.3f, 0.3f, 1.0f };
	float mat_ambient_color[] = { 0.8f, 0.8f, 0.2f, 1.0f };
	// float mat_diffuse[] = { 0.1f, 0.5f, 0.8f, 1.0f }; - moved to mesh.h
	float shininess = 20;
	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialf(GL_FRONT, GL_SHININESS, shininess);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glPushMatrix();
		gluLookAt(0, 0, 10, 0, 0, 0, 0, 1, 0);

		// Lab 2 Optional task - draw user select marquee, store projection matrix
		if (isSelecting)
		{
			glGetFloatv(GL_PROJECTION_MATRIX, projectionMatrix);
			glPushMatrix();
				if (isSelectingVertex)
					glColor4d(0, 200, 0, 0.3);
				else if (isDeselecting)
					glColor4d(200, 30, 0, 0.3);
				else
					glColor4d(0, 30, 200, 0.3);
				glTranslated(0, 0, 5.0);
				int windowWidth = glutGet(GLUT_WINDOW_WIDTH);
				int windowHeight = glutGet(GLUT_WINDOW_HEIGHT);
				double nearPlaneLength = tan(20.0 / 180.0 * M_PI) * 5 * 2;
				selectBoxCoords[0] = (startx - windowWidth / 2.0) / windowWidth * nearPlaneLength;
				selectBoxCoords[1] = ((windowHeight - starty) - windowHeight / 2.0) / windowHeight * nearPlaneLength;
				selectBoxCoords[2] = -5;
				selectBoxCoords[3] = (currX - windowWidth / 2.0) / windowWidth * nearPlaneLength;
				selectBoxCoords[4] = ((windowHeight - currY) - windowHeight / 2.0) / windowHeight * nearPlaneLength;
				selectBoxCoords[5] = -5;
				glRectdv(selectBoxCoords, selectBoxCoords + 3);
			glPopMatrix();
		}

		glRotatef(angle2, 1.0, 0.0, 0.0);
		glRotatef(angle, 0.0, 1.0, 0.0);
		glScalef(zoom, zoom, zoom);

		myObj.draw();
		glutSwapBuffers(); // moved up here for lab 2 selection marquee optional task

		// Lab 2 user marquee selection task
		// Render to a offscreen buffer with unique colours for all polygons
		// For finding 'visible' polygons
		if (triggerOffscreenDraw)
		{
			myObj.drawOffscreen();
			myObj.computeSelectedTriangles();
			triggerOffscreenDraw = false;
			isDeselecting = false;
			isSelectingFacing = false;
			glutPostRedisplay();
		}
	glPopMatrix();
}




void keyboard (unsigned char key, int x, int y)
{
	char filename[256];
	switch (key) {
	case 'p':
	case 'P':
		glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		break;			
	case 'w':
	case 'W':
		glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
		break;			
	case 'v':
	case 'V':
		glPolygonMode(GL_FRONT_AND_BACK,GL_POINT);
		break;			
	case 's':
	case 'S':
		cout << "Switching shading mode" << endl;
		m_Smooth = !m_Smooth;
		glShadeModel(m_Smooth ? GL_SMOOTH : GL_FLAT); // Lab 2 optional task
		break;
	case 'h':
	case 'H':
		m_Highlight = !m_Highlight;
		break;
	// Optional task 4
	case 'r':
	case 'R':
		cout << "Enter the 3ds or obj filename you want to read:";
		cin >> filename;
		myObj.readFile(filename); // Lab 1 optional task
		break;
	case 'o':
	case 'O':
		cout << "Enter the filename you want to write:";
		cin >> filename;
		myObj.writeFile(filename);
		break;
	// Final boss
	case 'b':
	case 'B':
		myObj.subdivide();
		break;
	case 'e':
	case 'E':
		myObj.relax();
		break;
	case 't':
	case 'T':
		break;
	case 'd':
	case 'D':
		myObj.decimate();
		break;
	case '1':
	case '2':
	case '3':
	case '4':
		current_object = key - '1';
		break;

	case 'Q':
	case 'q':
		exit(0);
	break;

	default:
	break;
	}

	glutPostRedisplay();
}

void arrows(int key, int x, int y)
{
	double* displacement = (double*)malloc(sizeof(double) * 3);
	displacement[0] = 0;
	displacement[1] = 0;
	displacement[2] = 0;

	switch (key)
	{
	case GLUT_KEY_UP:
		displacement[1] += laplacianDeformStepSize;
		break;
	case GLUT_KEY_DOWN:
		displacement[1] -= laplacianDeformStepSize;
		break;
	case GLUT_KEY_LEFT:
		displacement[0] -= laplacianDeformStepSize;
		break;
	case GLUT_KEY_RIGHT:
		displacement[0] += laplacianDeformStepSize;
		break;
	case GLUT_KEY_PAGE_UP:
		displacement[2] += laplacianDeformStepSize;
		break;
	case GLUT_KEY_PAGE_DOWN:
		displacement[2] -= laplacianDeformStepSize;
		break;
	case GLUT_KEY_HOME:
		laplacianDeformStepSize -= 0.1;
		cout << "Laplacian deform step size: " << laplacianDeformStepSize << endl;
		return;
	case GLUT_KEY_END:
		laplacianDeformStepSize += 0.1;
		cout << "Laplacian deform step size: " << laplacianDeformStepSize << endl;
		return;
	default:
		return;
	}

	myObj.displace(displacement);
	glutPostRedisplay();
}

void mouse(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN)
	{
		mouseButton = button;

		// Lab 2 Optional task - user select marquee
		if ((glutGetModifiers() & GLUT_ACTIVE_CTRL) != 0)
		{
			isSelecting = true;
			if ((glutGetModifiers() & GLUT_ACTIVE_SHIFT) != 0)
				isDeselecting = true;
			if ((glutGetModifiers() & GLUT_ACTIVE_ALT) != 0)
				isSelectingFacing = true;
		}
		else if ((glutGetModifiers() & GLUT_ACTIVE_ALT) != 0)
		{
			isSelecting = true;
			isSelectingVertex = true;
			isSelectingFacing = true;
		}
		else
		{
			moving = 1;
		}

		startx = x;
		starty = y;
    }
	
	// Lab 2 Optional task - user select marquee
	if (isSelecting)
	{
		currX = x;
		currY = y;
	}
	
	if (state == GLUT_UP)
	{
		mouseButton = button;
		moving = 0;

		// Lab 2 Optional task - user select marquee
		if (isSelecting)
		{
			// If using ctrl-alt click drag mode, delay computation to display callback,
			// where a offscreen draw is triggered to populate triangles with unique colours
			// (determine which primitives are on screen)
			triggerOffscreenDraw = isSelectingFacing;
			if (!triggerOffscreenDraw)
			{
				myObj.computeSelectedTriangles();
				isDeselecting = false;
				isSelectingFacing = false;
			}
			isSelecting = false;
		}
    }

	glutPostRedisplay();
	// cout << "Mouse action detected - x y: " << x << " " << y << endl;
}

void motion(int x, int y)
{
    if (moving)
	{
		if (mouseButton == GLUT_LEFT_BUTTON)
		{
			angle = angle + (x - startx);
			angle2 = angle2 + (y - starty);
		}
		else
		{
			zoom += ((y - starty) * 0.001);
		}

		startx = x;
		starty = y;
		glutPostRedisplay();
    }
	else if (isSelecting)
	{
		// Lab 2 Optional task - user select marquee
		currX = x;
		currY = y;
		glutPostRedisplay();
	}
}

int main(int argc, char **argv)
{
	char filename[255];
	cout<< "CS3242" << endl << endl;

	cout << "Enter the 3ds or obj filename you want to open:";
	cin >> filename;
	myObj.readFile(filename);

	//cout << "1-4: Draw different objects" << endl;
	cout << "R: Read another file (3ds or obj)" << endl; // Lab 1 optional task - 3ds
	cout << "O: Save object to file" << endl;
	cout << "S: Toggle Smooth Shading"<<endl;
	cout << "H: Toggle Highlight"<<endl;
	cout << "W: Draw Wireframe"<<endl;
	cout << "P: Draw Polygon"<<endl;
	cout << "V: Draw Vertices"<<endl;
	cout << "Q: Quit" <<endl<< endl;

	cout << "Left mouse click and drag: rotate the object"<<endl;
	cout << "Right mouse click and drag: zooming" << endl;

	// Lab 2 optional task - user selection
	cout << endl << "Triangle selection:" << endl;
	cout << "Left mouse ctrl-click and drag: selection box" << endl; 
	cout << "Left mouse ctrl-shift-click and drag: deselection box" << endl;
	cout << "Combine the above two with \"alt\" to select / deselect only visible triangles" << endl;

	// Final boss
	cout << "D: Decimate / simplify selected portion of mesh" << endl;
	cout << "B: Barycentric subdivide selected portion of mesh" << endl;
	cout << "E: Edge swapping to relax mesh" << endl;

	cout << "Alt (without ctrl) click drag to select a vertex for laplacian deformation," << endl
		<< "  with some selected triangles first (see above instructions)" << endl;
	cout << "Use left (-x) right(+x) up(+y) down(-y) pgup(+z) pgdown(-z) keys to laplacian-deform" << endl;
	cout << "Use home (decrease) and end (increase) keys to increase / decrease laplacian deform step size" << endl;

	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize (DEFAULT_WINDOW_LENGTH_WIDTH, DEFAULT_WINDOW_LENGTH_WIDTH);
	glutInitWindowPosition (50, 50);
	glutCreateWindow ("CS3242 Assignment");
	glClearColor (1.0,1.0,1.0, 1.0);
	glutDisplayFunc(display);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(arrows); // laplacian deformation final boss
	setupLighting();
	glDisable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST); 
	glDepthMask(GL_TRUE);

	glEnable(GL_POLYGON_OFFSET_FILL); // lab 2 visualize boundary edges task

    glMatrixMode(GL_PROJECTION);
    gluPerspective( /* field of view in degree */ 40.0,
  /* aspect ratio */ 1.0,
    /* Z near */ 1.0, /* Z far */ 80.0);
	glMatrixMode(GL_MODELVIEW);
	glutMainLoop();

	return 0;
}
