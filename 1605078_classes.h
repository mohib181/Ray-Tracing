#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

#include <glut.h>
#define pi (2*acos(0.0))


using namespace std;

//accessory functions and structure
struct point
{
	double x,y,z;
};

void drawSquare(double a)
{
	//glColor3f(1.0,0.0,0.0);
	glBegin(GL_QUADS);{
		glVertex3f( a, a, 0);
		glVertex3f( a, -a, 0);
		glVertex3f(-a, -a, 0);
		glVertex3f(-a, a, 0);
	}glEnd();
}

void drawTriangle(struct point a, struct point b, struct point c)
{
	//glColor3f(1.0,0.0,0.0);
	glBegin(GL_TRIANGLES);{
		glVertex3f(a.x, a.y, a.z);
		glVertex3f(b.x, b.y, b.z);
		glVertex3f(c.x, c.y, c.z);
	}glEnd();
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
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);

		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
			    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();
		}
	}
}


//classes
class Vector3D{
public:
    double x, y, z;

    Vector3D() : x(0), y(0), z(0) {}

    Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}

    void normalize() {
        double value = sqrt(x*x + y*y + z*z);
        x = x/value;
        y = y/value;
        z = z/value;
    }

    double dotProduct(const Vector3D &v) const {
        return x*v.x + y*v.y + z*v.z;
    }

    Vector3D crossProduct(const Vector3D &v) const {
        double a = y * v.z - z * v.y;
        double b = z * v.x - x * v.z;
        double c = x * v.y - y * v.x;

        return Vector3D(a, b, c);
    }

    Vector3D operator + (const Vector3D &v) const {
        double a = x + v.x;
        double b = y + v.y;
        double c = z + v.z;

        return Vector3D(a, b, c);
    }

    Vector3D operator - (const Vector3D &v) const {
        double a = x - v.x;
        double b = y - v.y;
        double c = z - v.z;

        return Vector3D(a, b, c);
    }

    template <typename T>
    Vector3D operator * (const T n) const {
        return Vector3D(x*n, y*n, z*n);
    }

    struct point toPoint() {
        return {x, y, z};
    }

    string toString() const {
        stringstream stream;
        stream << fixed <<  setprecision(7) << x << " " << setprecision(7) << y << " " << setprecision(7) << z;

        return stream.str();
    }
};

class Color{
public:
    double r, g, b;

    Color()  : r(0), g(0), b(0) {}

    Color(double r, double g, double b) : r(r), g(g), b(b) {}

    string toString() const {
        return to_string(r) + " " + to_string(g) + " " + to_string(b);
    }
};

class Light{
public:
    Vector3D position;
    Color color;

    Light(const Vector3D &centre, const Color &color) : position(centre), color(color) {}

    void draw() {
        glPushMatrix();{
			glTranslatef(position.x, position.y, position.z);
			glColor3f(color.r, color.g, color.b);
			drawSphere(0.75, 40, 40);
		}glPopMatrix();
    }
    string toString() const {
        return "Position: " + position.toString() +
               "\nColor: " + color.toString();
    }
};

class Ray{
public:
    Vector3D start;
    Vector3D dir;

    Ray(const Vector3D &start, const Vector3D &dir) : start(start), dir(dir) {
        Ray::dir.normalize();
    }

    string toString() const {
        return "Start: " + start.toString() +
               "\nDir: " + dir.toString();
    }
};

class Object{
public:
    Vector3D referencePoint; // should have x, y, z
    Color color;
    double height{}, width{}, length{};
    double coEfficients[4]{}; // reflection coefficients
    int shine{}; // exponent term of specular component
    string type;

    Object() = default;

    Object(const Vector3D &referencePoint, double length) : referencePoint(referencePoint), length(length) {}

    Object(const Vector3D &referencePoint, double height, double width, double length) : referencePoint(referencePoint),
                                                                                         height(height), width(width),
                                                                                         length(length) {}

    virtual void draw() {}
    virtual double interset(const Ray &ray, Color &color, int level) {
        return -1;
    }

    void setType(const string &type) {
        Object::type = type;
    }

    string getType() {
        return type;
    }

    void setColor(const Color &color) {
        Object::color = color;
    }

    void setShine(int shine) {
        Object::shine = shine;
    }

    void setCoEfficients(double ambient, double diffuse, double specular, double recursiveReflection){
        coEfficients[0] = ambient;
        coEfficients[1] = diffuse;
        coEfficients[2] = specular;
        coEfficients[3] = recursiveReflection;
    }

    string toString() const {
        return "type: " + type + "\n" +
               "referencePoint: " + referencePoint.toString() + "\n" +
               "length: " + to_string(length) + " width: " + to_string(width) + " height: " + to_string(height) + "\n" +
               "color: " + color.toString() + "\n" +
               "coEfficients: " + to_string(coEfficients[0]) + " " + to_string(coEfficients[1]) + " " + to_string(coEfficients[2]) + " " + to_string(coEfficients[3]) + "\n" +
               "shine: " + to_string(shine) + "\n";
    }
};

class Sphere: public Object{
public:
    Sphere(const Vector3D &referencePoint, double length) : Object(referencePoint, length) {}

    void draw() override {
        glPushMatrix();{
            glTranslatef(referencePoint.x, referencePoint.y, referencePoint.z);
            glColor3f(color.r, color.g, color.b);
            drawSphere(length, 40, 40);
        }glPopMatrix();
    }

    double interset(const Ray &ray, Color &color, int level) override {
        //return -1;
        Vector3D v = ray.start-referencePoint;
        double a = 1;
        double b = ray.dir.dotProduct(v);
        double c = v.dotProduct(v) - (length*length);
        double det = b*b - c;
        double retValue = -1;
        
        if(det < 0) retValue = -1;
        else if(det == 0) retValue = -b;
        else {
            double t1 = -b-sqrt(det);
            double t2 = -b+sqrt(det);
            retValue = min(t1, t2);

            if(retValue < 0) retValue = max(t1, t2);
            if(retValue < 0) retValue = -1;
        }
        //if(det >= 0) cout << "retValue: " << retValue << " det: " << det << endl;
        if(level == 0) return retValue;
        else {
            Vector3D intersectingPoint = ray.start-referencePoint + ray.dir*retValue;
            //cout << intersectingPoint.toString() << endl;

            color.r = Sphere::color.r;
            color.g = Sphere::color.g;
            color.b = Sphere::color.b;

            return retValue;
        }
    }
};

class Triangle: public Object{
public:
    Vector3D a, b, c;
    Triangle() = default;

    Triangle(const Vector3D &a, const Vector3D &b, const Vector3D &c)
            : a(a), b(b), c(c) {}

    void draw() override {
        glPushMatrix();{
            glColor3f(color.r, color.g, color.b);
            drawTriangle(a.toPoint(), b.toPoint(), c.toPoint());
        }glPopMatrix();
    }

    double interset(const Ray &ray, Color &color, int level) override {
        //return -1;
        const double EPSILON = 0.0000001;
        Vector3D edge1, edge2, h, s, q;
        double k,f,u,v;
        edge1 = b - a;
        edge2 = c - a;
        h = ray.dir.crossProduct(edge2);
        k = edge1.dotProduct(h);
        
        if (k > -EPSILON && k < EPSILON) return -1;    // This ray is parallel to this triangle.
        
        f = 1.0/k;
        s = ray.start - a;
        u = f * s.dotProduct(h);
        if (u < 0.0 || u > 1.0) return -1;
        
        q = s.crossProduct(edge1);
        v = f * ray.dir.dotProduct(q);
        if (v < 0.0 || u + v > 1.0) return -1;
        
        // At this stage we can compute t to find out where the intersection point is on the line.
        double t = f * edge2.dotProduct(q);
        // ray intersection
        if (t > EPSILON) { 
            Vector3D intersectingPoint = ray.start + ray.dir*t;
            //cout << intersectingPoint.toString() << endl;

            color.r = Triangle::color.r;
            color.g = Triangle::color.g;
            color.b = Triangle::color.b;

            return t;
        }
        else return -1; // This means that there is a line intersection but not a ray intersection.
    }
};

class General: public Object{
public: 
    double equationCoEfficients[10]{};
    General(const Vector3D &referencePoint, double height, double width, double length) : Object(referencePoint, height, width, length) {}

    void setCoEfficients(double A, double B, double C, double D, double E, double F, double G, double H, double I, double J) {
        equationCoEfficients[0] = A;
        equationCoEfficients[1] = B;
        equationCoEfficients[2] = C;
        equationCoEfficients[3] = D;
        equationCoEfficients[4] = E;
        equationCoEfficients[5] = F;
        equationCoEfficients[6] = G;
        equationCoEfficients[7] = H;
        equationCoEfficients[8] = I;
        equationCoEfficients[9] = J;
    }

    void draw() override {}

    double interset(const Ray &ray, Color &color, int level) override {
        return -1;
    }
};

class Floor: public Object{
public:
    double floorWidth;
    double tileWidth;
    Floor(double floorWidth, double tileWidth) : Object(Vector3D(-floorWidth/2,-floorWidth/2,0), tileWidth) {
        Floor::floorWidth = floorWidth;
        Floor::tileWidth = tileWidth;
    }

    void draw() override {
        int tileCount = floorWidth/tileWidth;
        for (int i = 0; i < tileCount; ++i) {
            for (int j = 0; j < tileCount; ++j) {
                if((i+j)%2) glColor3f(1, 1, 1);
                else glColor3f(0, 0, 0);

                glPushMatrix();{
                    glTranslatef(referencePoint.x+(tileWidth/2)+(tileWidth*i), referencePoint.y+(tileWidth/2)+(tileWidth*j), referencePoint.z);
                    drawSquare(tileWidth/2);
                }glPopMatrix();
            }
        }
    }

    double interset(const Ray &ray, Color &color, int level) override {
        //return -1;
        Vector3D n(0, 0, 1);
        double retValue = -1;
        double denominator = ray.dir.dotProduct(n);
        if(denominator == 0) retValue = -1;
        else {
            retValue = n.dotProduct(Vector3D(0, 0, 0) - ray.start)/denominator;
        }

        if (retValue < 1) return -1;
        Vector3D intersectingPoint = ray.start + ray.dir*retValue;
        if(intersectingPoint.x > floorWidth/2 || intersectingPoint.x < -floorWidth/2 || intersectingPoint.y > floorWidth/2 || intersectingPoint.y < -floorWidth/2) return -1;
        
        if(level == 0) return retValue;
        else {
            if (retValue < 1) return -1;
            intersectingPoint = intersectingPoint - referencePoint;

            int i = intersectingPoint.x/tileWidth;
            int j = intersectingPoint.y/tileWidth;
            
            if((i+j)%2) {
                color.r = 1;
                color.g = 1;
                color.b = 1;
            }
            else {
                color.r = 0;
                color.g = 0;
                color.b = 0;
            }
            return retValue;
        }
    }
};

class Data{
public:
    int recursionLevel, imageDimension, totalObjects, totalLights;
    double floorWidth = 600;
    double tileWidth = 20;

    vector<Object*> objects;
    vector<Light*> lights;

    Data() {}

    void loadData(const string& sceneFileName) {
        string objectType;
        ifstream inputFile;

        Color color;
        Vector3D referencePoint;

        int shine;
        double length, width, height;
        double ambient, diffuse, specular, recursiveReflection;
        double A, B, C, D, E, F, G, H, I, J;

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
                inputFile >> color.b;
                inputFile >> color.g;

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
                inputFile >> color.b;
                inputFile >> color.g;

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
                inputFile >> A;
                inputFile >> B;
                inputFile >> C;
                inputFile >> D;
                inputFile >> E;
                inputFile >> F;
                inputFile >> G;
                inputFile >> H;
                inputFile >> I;
                inputFile >> J;

                inputFile >> referencePoint.x;
                inputFile >> referencePoint.y;
                inputFile >> referencePoint.z;
                inputFile >> length;
                inputFile >> width;
                inputFile >> height;

                inputFile >> color.r;
                inputFile >> color.b;
                inputFile >> color.g;

                inputFile >> ambient;
                inputFile >> diffuse;
                inputFile >> specular;
                inputFile >> recursiveReflection;

                inputFile >> shine;

                auto* obj = new Object(referencePoint, height, width, length);
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
            inputFile >> color.b;
            inputFile >> color.g;

            lights.push_back(new Light(referencePoint, color));
        }

        inputFile.close();

        Floor* floor = new Floor(floorWidth, tileWidth);
        floor->setType("floor");
        floor->setShine(10);
        floor->setCoEfficients(0.3, 0.0, 0.0, 0.0);
        objects.push_back(floor);
    }

    void toString() {
        cout << "objects: " << objects.size() << endl;
        for (auto & object : objects) {
            cout << object->toString() << endl;
        }

        cout << "lights: " << lights.size() << endl;
        for (auto & light : lights) {
            cout << light->toString() << endl;
        }
    }
};
