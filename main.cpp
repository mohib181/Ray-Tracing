#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <windows.h>
#include <glut.h>

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
double q_angle;
double e_angle;
double a_angle;
double d_angle;
double angle_inc;

double plane_dist;
double plane_len;

double big_rad;
double small_rad;
double cylinder_len;
double header_len;

struct point
{
	double x,y,z;
};

struct point pos;
struct point u, r, l;
struct point gun_l;

int bullet_count, max_bullets;
struct point bullets[100];


void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 700,0,0);
			glVertex3f(-700,0,0);

			glVertex3f(0,-700,0);
			glVertex3f(0, 700,0);

			glVertex3f(0,0, 700);
			glVertex3f(0,0,-700);
		}glEnd();
	}
}

void drawGrid()
{
	int i;
	if(drawgrid==1)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -90, 0);
				glVertex3f(i*10,  90, 0);

				//lines parallel to X-axis
				glVertex3f(-90, i*10, 0);
				glVertex3f( 90, i*10, 0);
			}
		}glEnd();
	}
}

void drawSquare(double a)
{
    //glColor3f(1.0,0.0,0.0);
	glBegin(GL_QUADS);{
		glVertex3f( a, 0, a);
		glVertex3f( a, 0, -a);
		glVertex3f(-a, 0, -a);
		glVertex3f(-a, 0, a);
	}glEnd();
}

void drawCircle(double radius,int segments)
{
    int i;
    struct point points[100];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawCone(double radius,double height,int segments)
{
    int i;
    double shade;
    struct point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw triangles using generated points
    for(i=0;i<segments;i++)
    {
        //create shading effect
        if(i<segments/2)shade=2*(double)i/(double)segments;
        else shade=2*(1.0-(double)i/(double)segments);
        glColor3f(shade,shade,shade);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(0,0,height);
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			//points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			//points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			//points[i][j].z=h;
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].z=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].y=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);

		for(j=0;j<slices;j++)
		{
		    if(j%2 == 0) glColor3f(0, 0, 0);
            else glColor3f(1, 1, 1);

			glBegin(GL_QUADS);{
			    //upper hemisphere
			    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x,-points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,-points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,-points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,-points[i+1][j].y,points[i+1][j].z);
			}glEnd();
		}
	}
}

void drawHalfSphere(double radius, int slices, int stacks, bool dir)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].z=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].y=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        for(j=0;j<slices;j++)
		{
		    if(j%2 == 0) glColor3f(0, 0, 0);
            else glColor3f(1, 1, 1);

			glBegin(GL_QUADS);{
			    if(dir) {
					//upper hemisphere
			    	glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
					glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
					glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
					glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
				}else {
					//lower hemisphere
					glVertex3f(points[i][j].x,-points[i][j].y,points[i][j].z);
					glVertex3f(points[i][j+1].x,-points[i][j+1].y,points[i][j+1].z);
					glVertex3f(points[i+1][j+1].x,-points[i+1][j+1].y,points[i+1][j+1].z);
					glVertex3f(points[i+1][j].x,-points[i+1][j].y,points[i+1][j].z);
				}
			}glEnd();
		}
	}
}

void drawCylinder(double radius, int slices, int stacks)
{
	struct point points[150][150];
	int i,j;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=radius*cos(((double)j/(double)slices)*2*pi);
			points[i][j].z=radius*sin(((double)j/(double)slices)*2*pi);
			points[i][j].y=i;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
		for(j=0;j<slices;j++)
		{
		    if(j%2 == 0) glColor3f(0, 0, 0);
            else glColor3f(1, 1, 1);
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
			}glEnd();
		}
	}
}

void drawCanonHead(double radius, int slices, int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		r=2*radius - radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].z=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].y=i;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        for(j=0;j<slices;j++)
		{
		    if(j%2 == 0) glColor3f(0, 0, 0);
            else glColor3f(1, 1, 1);

			glBegin(GL_QUADS);{
			    //upper hemisphere
			    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
            }glEnd();
		}
	}
}


void drawShootingLine() {
	glPushMatrix();
	glColor3f(1, 0, 0);
	glBegin(GL_LINES);{
		glVertex3f(gun_l.x, gun_l.y, gun_l.z);
		glVertex3f(gun_l.x*700, gun_l.y*700, gun_l.z*700);
	}glEnd();
	glPopMatrix();
}

void shoot() {
	double t = plane_dist/gun_l.y;
	struct point c = {gun_l.x*t, gun_l.y*t, gun_l.z*t};

	if((c.x < plane_len && c.x > -plane_len) && (c.z < plane_len && c.z > -plane_len)) {
		printf("hit\n");
		//printf("c: %f, %f, %f\n", c.x, c.y, c.z);
		bullets[(bullet_count++)%max_bullets] = {c.x, c.y-5, c.z};
	}
	else {
		printf("miss\n");
	}
}

void drawBulletShot() {
	for (int i = 0; i < bullet_count; i++)
	{
		glPushMatrix();

		glTranslatef(bullets[i].x, bullets[i].y, bullets[i].z);
    	glColor3f(1, 0, 0);
    	drawSquare(5);

		glPopMatrix();
	}

}

void drawStuff() {
    glPushMatrix();

    glRotatef(q_angle, 0, 0, 1);
    drawHalfSphere(big_rad, 40, 40, false);

	glRotatef(e_angle, 1, 0, 0);
	drawHalfSphere(big_rad, 40, 40, true);
    //drawSphere(big_rad, 40, 40);

	glTranslatef(0, big_rad, 0);
    glRotatef(a_angle, 1, 0, 0);
    glRotatef(d_angle, 0, 1, 0);

	glTranslatef(0, small_rad, 0);
	drawHalfSphere(small_rad, 40, 40, false);
    drawCylinder(small_rad, 40, cylinder_len);

    glTranslatef(0, cylinder_len, 0);
    drawCanonHead(small_rad, 40, header_len);

    glPopMatrix();

	//drawShootingLine();

	glPushMatrix();

    glTranslatef(0, plane_dist, 0);
    glColor3f(0.5, 0.5, 0.5);
    drawSquare(plane_len);

    glPopMatrix();

	drawBulletShot();
}

void reset_pos() {
    pos = {100, 100, 25};
	u = {0, 0, 1};
	r = {-1.0/sqrt(2), 1.0/sqrt(2), 0};
	l = {-1.0/sqrt(2), -1.0/sqrt(2), 0};

	gun_l = {0, -1, 0};
	q_angle = 0;
	e_angle = 0;
	a_angle = 0;
	d_angle = 0;
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
	else if (dir == GUN) {
		gun_l.x = c.x/unit;
        gun_l.y = c.y/unit;
        gun_l.z = c.z/unit;
		//printf("gun_l: gun_l.x: %f, gun_l.y: %f, gun_l.z: %f\n", gun_l.x, gun_l.y, gun_l.z);
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

		case 'q':
			if(q_angle+angle_inc < 60) {
				q_angle += angle_inc;
				rotate_vector(gun_l, {0, 0, 1}, angle_inc*(pi/180), GUN);
			}
            break;
        case 'w':
			if (q_angle-angle_inc > -60) {
				q_angle -= angle_inc;
				rotate_vector(gun_l, {0, 0, 1}, -angle_inc*(pi/180), GUN);
			}
            break;
        case 'e':
            if (e_angle+angle_inc < 50) {
				e_angle += angle_inc;
				rotate_vector(gun_l, {1, 0, 0}, angle_inc*(pi/180), GUN);
			}
			break;
        case 'r':
            if (e_angle-angle_inc > -50) {
				e_angle -= angle_inc;
				rotate_vector(gun_l, {1, 0, 0}, -angle_inc*(pi/180), GUN);
			}
			break;
        case 'a':
            if (a_angle+angle_inc < 50) {
				a_angle += angle_inc;
				rotate_vector(gun_l, {1, 0, 0}, angle_inc*(pi/180), GUN);
			}
			break;
        case 's':
            if (a_angle-angle_inc > -50) {
				a_angle -= angle_inc;
				rotate_vector(gun_l, {1, 0, 0}, -angle_inc*(pi/180), GUN);
			}
			break;
        case 'd':
            d_angle = d_angle+angle_inc > 60 ? d_angle : d_angle+angle_inc;
            break;
        case 'f':
            d_angle = d_angle-angle_inc < -60 ? d_angle : d_angle-angle_inc;
            break;

		case '0':
            reset_pos();
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
				shoot();
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

	drawAxes();
	drawGrid();

	drawStuff();

    //glColor3f(1,0,0);
    //drawSquare(10);

    //drawSS();

    //drawCircle(30,24);

    //drawCone(20,50,24);

	//drawSphere(30,24,20);




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
	angle_inc = 2;

	pos = {200, 150, 50};
	u = {0, 0, 1};
	r = {-1.0/sqrt(2), 1.0/sqrt(2), 0};
	l = {-1.0/sqrt(2), -1.0/sqrt(2), 0};

	gun_l = {0, 1, 0};

	q_angle = 0;
	e_angle = 0;
	a_angle = 0;
	d_angle = 0;

	bullet_count = 0;
	max_bullets = 100;

	plane_dist = 500;
	plane_len = 250;

	big_rad = 40;
	small_rad = 10;
	cylinder_len = 100;
	header_len = 20;

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
