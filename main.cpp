#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <windows.h>
#include <glut.h>

#include "1605078_classes.h"
#include "bitmap_image.hpp"

#define pi (2*acos(0.0))
#define UP 1
#define RIGHT 2
#define LOOK 3
#define GUN 4

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;

double rotate_angle;

struct point pos;
struct point u, r, l;

int recursionLevel, imageDimension;
double floorWidth = 1000;
double tileWidth = 20;

vector<Object*> objects;
vector<Light*> lights;

void loadData(const string& sceneFileName) {
    string objectType;
    ifstream inputFile;

    Color color;
    Vector3D referencePoint;

    int shine;
	int totalObjects, totalLights;
    double length;
    double ambient, diffuse, specular, recursiveReflection;
    
    inputFile.open(sceneFileName);

    inputFile >> recursionLevel;
    inputFile >> imageDimension;
    inputFile >> totalObjects;

    objects.reserve(totalObjects+1);
    for (int i = 0; i < totalObjects; ++i) {
        inputFile >> objectType;
        //cout << "objectType: " << objectType << endl;
        if(objectType == "sphere") {
            inputFile >> referencePoint.x;
            inputFile >> referencePoint.y;
            inputFile >> referencePoint.z;
            inputFile >> length;

            inputFile >> color.r;
            inputFile >> color.g;
            inputFile >> color.b;

            inputFile >> ambient;
            inputFile >> diffuse;
            inputFile >> specular;
            inputFile >> recursiveReflection;

            inputFile >> shine;
            auto* obj = new Sphere(referencePoint, length);
            obj->setType("sphere");
            obj->setColor(color);
            obj->setCoEfficients(ambient, diffuse, specular, recursiveReflection);
            obj->setShine(shine);

            //cout << obj->toString() << endl;
            objects.push_back(obj);
        }
        else if(objectType == "triangle") {
            Vector3D a, b, c;
            inputFile >> a.x;
            inputFile >> a.y;
            inputFile >> a.z;
            inputFile >> b.x;
            inputFile >> b.y;
            inputFile >> b.z;
            inputFile >> c.x;
            inputFile >> c.y;
            inputFile >> c.z;

            inputFile >> color.r;
            inputFile >> color.g;
            inputFile >> color.b;

            inputFile >> ambient;
            inputFile >> diffuse;
            inputFile >> specular;
            inputFile >> recursiveReflection;

            inputFile >> shine;

            auto* triangle = new Triangle(a, b, c);
            triangle->setColor(color);
            triangle->setCoEfficients(ambient, diffuse, specular, recursiveReflection);
            triangle->setShine(shine);
            triangle->setType("triangle");

            //cout << triangle->toString() << endl;
            objects.push_back(triangle);
        }
        else if(objectType == "general") {
            auto* obj = new GeneralQuadraticSurface();
            
            inputFile >> obj->A;
            inputFile >> obj->B;
            inputFile >> obj->C;
            inputFile >> obj->D;
            inputFile >> obj->E;
            inputFile >> obj->F;
            inputFile >> obj->G;
            inputFile >> obj->H;
            inputFile >> obj->I;
            inputFile >> obj->J;

            inputFile >> obj->referencePoint.x;
            inputFile >> obj->referencePoint.y;
            inputFile >> obj->referencePoint.z;
            inputFile >> obj->length;
            inputFile >> obj->width;
            inputFile >> obj->height;

            inputFile >> color.r;
            inputFile >> color.g;
            inputFile >> color.b;

            inputFile >> ambient;
            inputFile >> diffuse;
            inputFile >> specular;
            inputFile >> recursiveReflection;

            inputFile >> shine;

            obj->setType("general");
            obj->setColor(color);
            obj->setCoEfficients(ambient, diffuse, specular, recursiveReflection);
            obj->setShine(shine);

            //cout << obj->toString() << endl;
            objects.push_back(obj);
        }
    }

    inputFile >> totalLights;
    lights.reserve(totalLights);
    for (int i = 0; i < totalLights; ++i) {
        inputFile >> referencePoint.x;
        inputFile >> referencePoint.y;
        inputFile >> referencePoint.z;

        inputFile >> color.r;
        inputFile >> color.g;
        inputFile >> color.b;

        lights.push_back(new Light(referencePoint, color));
    }

    inputFile.close();

    Floor* floor = new Floor(floorWidth, tileWidth);
    floor->setType("floor");
    floor->setShine(10);
    floor->setCoEfficients(0.4, 0.2, 0.2, 0.2);
	objects.push_back(floor);
}

void drawStuff() {
	for (auto & object : objects) {
		object->draw();
	}

	for (auto & light : lights) {
		light->draw();
	}
}

void capture() {
	bitmap_image image(imageDimension, imageDimension);
    for(int i=0;i<imageDimension;i++){
        for(int j=0;j<imageDimension;j++){
            image.set_pixel(i,j, 0, 0, 0);
        }
    }

	double windowHeight = 500;
	double windowWidth = 500;
	double viewAngle = 80*(pi/180);
	double planeDistance = (windowHeight/2.0) / tan(viewAngle/2.0);
	
	Vector3D eye(pos.x, pos.y, pos.z);
	Vector3D look(l.x, l.y, l.z);
	Vector3D right(r.x, r.y, r.z);
	Vector3D up(u.x, u.y, u.z);
	
	Vector3D topLeft = eye + look*planeDistance - right*(windowWidth/2.0) + up*(windowHeight/2.0);
	
	double du = windowWidth/imageDimension;
	double dv = windowHeight/imageDimension;

	topLeft = topLeft + right*(du/2) - up*(dv/2);
	//cout << "topLeft: " << topLeft.toString() << endl;

	double t, t_min = 100000;
	for(int i=0;i<imageDimension;i++){
        for(int j=0;j<imageDimension;j++){
            Vector3D curPixel = topLeft + right*(du*i) - up*(dv*j);
			Ray ray(eye, curPixel-eye);
			
			t_min = 10000000;
			Color color;
			Object* nearestObj = nullptr;
			for (auto & object : objects) {
				t = object->intersect(ray, color, 0);
				if(t > 0 && t < t_min) {
					t_min = t;
					nearestObj = object;
					//if(object->type != "floor") cout << i << "," << j << " " << t << " " << object->type << "\n";
				}
			}
			
			if(nearestObj) {
				//cout << "nearestObj: " << nearestObj->toString() << endl;
				t_min = nearestObj->intersect(ray, color, 1);
				image.set_pixel(i,j, color.r*255, color.g*255, color.b*255);
			}
        }
    }
	//cout << "out of loop" << endl;
	image.save_image("out.bmp");
						   
}

void reset_pos() {
    pos = {100, 100, 25};
	u = {0, 0, 1};
	r = {-1.0/sqrt(2), 1.0/sqrt(2), 0};
	l = {-1.0/sqrt(2), -1.0/sqrt(2), 0};
}

void rotate_vector(struct point a, struct point b, double theta, int dir){
    struct point p, c;

    p.x = (b.y*a.z - b.z*a.y);
    p.y = (b.z*a.x - b.x*a.z);
    p.z = (b.x*a.y - b.y*a.x);
    //printf("perpendicular to a and b: p.x: %f, p.y: %f, p.z: %f\n", p.x, p.y, p.z);

	double dot = a.x*b.x + a.y*b.y + a.z*b.z;

	c.x = a.x * cos(theta) + p.x * sin(theta) + b.x * dot * (1 - cos(theta));
	c.y = a.y * cos(theta) + p.y * sin(theta) + b.y * dot * (1 - cos(theta));
	c.z = a.z * cos(theta) + p.z * sin(theta) + b.z * dot * (1 - cos(theta));

	double unit = sqrt(c.x*c.x + c.y*c.y + c.z*c.z);
	//unit = 1;

	if(dir == LOOK){
        l.x = c.x/unit;
        l.y = c.y/unit;
        l.z = c.z/unit;
        //printf("l: l.x: %f, l.y: %f, l.z: %f\n", l.x, l.y, l.z);
	}
	else if (dir == RIGHT) {
        r.x = c.x/unit;
        r.y = c.y/unit;
        r.z = c.z/unit;
        //printf("r: r.x: %f, r.y: %f, r.z: %f\n", r.x, r.y, r.z);
	}
	else if (dir == UP) {
        u.x = c.x/unit;
        u.y = c.y/unit;
        u.z = c.z/unit;
        //printf("u: u.x: %f, u.y: %f, u.z: %f\n", u.x, u.y, u.z);
	}
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '1':
			//look left
			//rotate counterclockwise l, r w.r.t u
			rotate_vector(l, u, rotate_angle, LOOK);
			rotate_vector(r, u, rotate_angle, RIGHT);
			break;
		case '2':
			//look right
			//rotate clockwise l, r w.r.t u
			rotate_vector(l, u, -rotate_angle, LOOK);
			rotate_vector(r, u, -rotate_angle, RIGHT);
			break;
		case '3':
			//look up
			//rotate counterclockwise l, u w.r.t r
			rotate_vector(l, r, rotate_angle, LOOK);
			rotate_vector(u, r, rotate_angle, UP);
			break;
		case '4':
			//look down
			//rotate clockwise l, u w.r.t r
			rotate_vector(l, r, -rotate_angle, LOOK);
			rotate_vector(u, r, -rotate_angle, UP);
			break;
		case '5':
			//tilt clockwise
			//rotate clockwise r, u w.r.t l
			rotate_vector(r, l, -rotate_angle, RIGHT);
			rotate_vector(u, l, -rotate_angle, UP);
			break;
		case '6':
			//tilt counterclockwise
			//rotate counterclockwise r, u w.r.t l
			rotate_vector(r, l, rotate_angle, RIGHT);
			rotate_vector(u, l, rotate_angle, UP);
			break;

		case '0':
            //reset_pos();
            capture();
			break;

		default:
			break;
	}
}

void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			pos.x += -l.x;
			pos.y += -l.y;
			pos.z += -l.z;
			break;
		case GLUT_KEY_UP:		//up arrow key
			pos.x += l.x;
			pos.y += l.y;
			pos.z += l.z;
			break;

		case GLUT_KEY_RIGHT:	//right arrow key
			pos.x += r.x;
			pos.y += r.y;
			pos.z += r.z;
			break;
		case GLUT_KEY_LEFT:		//left arrow key
			pos.x += -r.x;
			pos.y += -r.y;
			pos.z += -r.z;
			break;

		case GLUT_KEY_PAGE_UP:
			pos.x += u.x;
			pos.y += u.y;
			pos.z += u.z;
			break;
		case GLUT_KEY_PAGE_DOWN:
			pos.x += -u.x;
			pos.y += -u.y;
			pos.z += -u.z;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}

void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				//shoot();
			}
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}

void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	//gluLookAt(0,0,200,	0,0,0,	0,1,0);
	gluLookAt(pos.x, pos.y, pos.z, pos.x+l.x, pos.y+l.y, pos.z+l.z, u.x, u.y, u.z);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	//drawAxes();
	//drawGrid();

	drawStuff();
	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}

void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=0;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;
	rotate_angle=pi/10.0;

	pos = {90, -90, 50};
	u = {0, 0, 1};
	r = {1.0/sqrt(2), 1.0/sqrt(2), 0};
	l = {-1.0/sqrt(2), 1.0/sqrt(2), 0};

	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("Ray Tracing Field");

	string dir = "sample/";
    string sceneFileName = dir+"scene.txt";
	loadData(sceneFileName);
	//toString();
	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
