#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

#include <glut.h>
#define pi (2*acos(0.0))
#define AMB 0
#define DIFF 1
#define SPEC 2
#define REF 3

using namespace std;

//accessory functions and variables
extern int recursionLevel;
extern double floorWidth;
extern double tileWidth;

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

    double calculateDistance(const Vector3D &v) const {
        return sqrt((x-v.x)*(x-v.x) + (y-v.y)*(y-v.y) + (z-v.z)*(z-v.z));
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

    Color operator + (const Color &c) const {
        double red = min(this->r + c.r, 1.0);
        double green = min(this->g + c.g, 1.0);
        double blue = min(this->b + c.b, 1.0);

        return Color(red, green, blue);
    }

    Color operator * (const Color &c) const {
        double red = min(this->r*c.r, 1.0);
        double green = min(this->g*c.g, 1.0);
        double blue = min(this->b*c.b, 1.0);

        return Color(red, green, blue);
    }

    template <typename T>
    Color operator * (const T n) const {
        return Color(r*n, g*n, b*n);
    }

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
    virtual double intersect(const Ray &ray, Color &color, int level) {
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

extern vector<Object*> objects;
extern vector<Light*> lights;

void getRayColor(const Light* light, const Ray & ray, const Vector3D &intersectingPoint, const Vector3D &normal, double diffConstant, double specConstant, double shine, const Color &objectColor, Color &color) {
    Vector3D dir = light->position-intersectingPoint;
    dir.normalize();
    
    Vector3D start = intersectingPoint + dir*1;
    Ray lightRay(start, dir);
    
    double distance = start.calculateDistance(light->position);
    double t;
    bool inShadow = false;
    Color dummyColor;
    for (auto & object : objects) {
        t = object->intersect(lightRay, dummyColor, 0);
        if(t > 0 && t < distance) {
            inShadow = true;
            break;
        }
    }

    //cout << intersectingPoint.toString() << "->" << inShadow << endl;
    if(!inShadow) {
        Vector3D reflectedRay = ray.dir - normal*(normal.dotProduct(ray.dir)*2);
        reflectedRay.normalize();

        //cout << normal.toString() << "\n" << reflectedRay.toString() << endl;

        double lambertValue = max(normal.dotProduct(lightRay.dir), 0.0);
        double phongValue = max(reflectedRay.dotProduct(ray.dir), 0.0);
        //cout << lambertValue << " " << phongValue << " diff " << diffConstant << " spec: " << specConstant << endl;
        
        Color lambertColor = light->color*objectColor*lambertValue*diffConstant;
        Color phongColor = light->color*objectColor*pow(phongValue, shine)*specConstant;

        //cout << lambertColor.toString() << " " << phongColor.toString() << endl;

        color = color + lambertColor + phongColor; 
    }
}

void getRecursionColor(const Ray &ray, const Vector3D &intersectionPoint, const Vector3D &normal, double refConstant, int level, Color &color) {
    if(level < recursionLevel) {
        Vector3D reflectionDir = ray.dir - normal*(normal.dotProduct(ray.dir)*2);
        reflectionDir.normalize();
        
        Vector3D start = intersectionPoint + reflectionDir*1;
        Ray reflectionRay(start, reflectionDir);

        double t, t_min = 10000000;
        Color reflectedColor;
        Object* nearestObj = nullptr;
        for (auto & object : objects) {
            t = object->intersect(reflectionRay, reflectedColor, 0);
            if(t > 0 && t < t_min) {
                t_min = t;
                nearestObj = object;
            }
        }

        if(nearestObj) {
            //cout << "nearestObj: " << nearestObj->toString() << endl;
            t_min = nearestObj->intersect(reflectionRay, reflectedColor, level+1);
            color = color + reflectedColor*refConstant;
        }
    }
}

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

    double intersect(const Ray &ray, Color &color, int level) override {
        //return -1;
        Vector3D v = ray.start-referencePoint;
        double a = 1;
        double b = ray.dir.dotProduct(v);
        double c = v.dotProduct(v) - (length*length);
        double det = b*b - c;
        double retValue = -1;
        
        if(det < 0) return -1;
        else if(det == 0) retValue = -b;
        else {
            double t1 = -b-sqrt(det);
            double t2 = -b+sqrt(det);
            retValue = min(t1, t2);

            if(retValue < 0) retValue = max(t1, t2);
            if(retValue < 0) return -1;
        }
        //if(det >= 0) cout << "retValue: " << retValue << " det: " << det << endl;
        if(level == 0) return retValue;
        else {
            Vector3D intersectingPoint = ray.start + ray.dir*retValue;
            //cout << intersectingPoint.toString() << endl;
            Vector3D normal = intersectingPoint - referencePoint;
            normal.normalize();

            color = this->color*coEfficients[AMB];
            //cout << "before ray: " << color.toString() << endl;
            for(auto& light: lights) {
                getRayColor(light, ray, intersectingPoint, normal, coEfficients[DIFF], coEfficients[SPEC], shine, this->color, color);
            }
            //cout << "after ray: " << color.toString() << endl;
            getRecursionColor(ray, intersectingPoint, normal, coEfficients[REF], level, color);
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

    double intersect(const Ray &ray, Color &color, int level) override {
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
            if(level == 0) return t;

            Vector3D intersectingPoint = ray.start + ray.dir*t;
            Vector3D normal = edge1.crossProduct(edge2);
            normal.normalize();
            //cout << intersectingPoint.toString() << endl;

            color = this->color*coEfficients[AMB];

            for(auto& light: lights) {
                getRayColor(light, ray, intersectingPoint, normal, coEfficients[DIFF], coEfficients[SPEC], shine,this->color, color);
            }

            getRecursionColor(ray, intersectingPoint, normal, coEfficients[REF], level, color);

            return t;
        }
        else return -1; // This means that there is a line intersection but not a ray intersection.
    }
};

class GeneralQuadraticSurface: public Object{
public: 
    double A, B, C, D, E, F, G, H, I, J;
    GeneralQuadraticSurface() {}

    bool boundaryCheck(const Ray &ray, double t) {
        Vector3D intersectingPoint = ray.start + ray.dir*t;

        if(length > 0 && (intersectingPoint.x < referencePoint.x || intersectingPoint.x > referencePoint.x+length)) return false;
        if(width > 0 && (intersectingPoint.y < referencePoint.y || intersectingPoint.y > referencePoint.y+width)) return false;
        if(height > 0 && (intersectingPoint.z < referencePoint.z || intersectingPoint.z > referencePoint.z+height)) return false;

        return true;
    }

    void draw() override {}

    Vector3D getNormal(const Vector3D& intersectingPoint) {
        double x = 2*A*intersectingPoint.x+
                    D*intersectingPoint.y+
                    E*intersectingPoint.z+
                    G;
        double y = 2*B*intersectingPoint.y+
                    D*intersectingPoint.x+
                    F*intersectingPoint.z+
                    H;
        double z = 2*C*intersectingPoint.z+
                    E*intersectingPoint.x+
                    F*intersectingPoint.y+
                    I;

        return Vector3D(x, y, z);                    
    }

    double intersect(const Ray &ray, Color &color, int level) override {
        //return -1;
        double retValue = -1;
        double a = A*ray.dir.x*ray.dir.x + 
                    B*ray.dir.y*ray.dir.y + 
                    C*ray.dir.z*ray.dir.z + 
                    D*ray.dir.x*ray.dir.y + 
                    E*ray.dir.x*ray.dir.z + 
                    F*ray.dir.y*ray.dir.z;
        double b = 2*A*ray.start.x*ray.dir.x + 2*B*ray.start.y*ray.dir.y + 2*C*ray.start.z*ray.dir.z +
                    D*ray.start.x*ray.dir.y + 
                    D*ray.start.y*ray.dir.x + 
                    E*ray.start.x*ray.dir.z + 
                    E*ray.start.z*ray.dir.x + 
                    F*ray.start.y*ray.dir.z + 
                    F*ray.start.z*ray.dir.y + 
                    G*ray.dir.x + 
                    H*ray.dir.y + 
                    I*ray.dir.z;
        double c = A*ray.start.x*ray.start.x + 
                    B*ray.start.y*ray.start.y + 
                    C*ray.start.z*ray.start.z + 
                    D*ray.start.x*ray.start.y + 
                    E*ray.start.x*ray.start.z + 
                    F*ray.start.y*ray.start.z + 
                    G*ray.start.x + 
                    H*ray.start.y + 
                    I*ray.start.z +
                    J;

        double det = b*b-4*a*c;

        if(det < 0) return -1;
        else if(det == 0) retValue = -b/(2*a);
        else {
            double t1 = (-b-sqrt(det))/(2*a);
            double t2 = (-b+sqrt(det))/(2*a);

            if(boundaryCheck(ray, min(t1, t2))) retValue = min(t1, t2);
            if(retValue < 0 && boundaryCheck(ray, max(t1, t2))) retValue = max(t1, t2);
            if(retValue < 0) return -1;
        }

        if(level == 0) return retValue;
        else {
            Vector3D intersectingPoint = ray.start + ray.dir*retValue;
            Vector3D normal = getNormal(intersectingPoint);
            normal.normalize();
            
            color = this->color*coEfficients[AMB];

            for(auto& light: lights) {
                getRayColor(light, ray, intersectingPoint, normal, coEfficients[DIFF], coEfficients[SPEC], shine, this->color, color);
            }

            getRecursionColor(ray, intersectingPoint, normal, coEfficients[REF], level, color);

            return retValue;
        }
    }
};

class Floor: public Object{
public:
    Floor(double floorWidth, double tileWidth) : Object(Vector3D(-floorWidth/2,-floorWidth/2,0), tileWidth) {}

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

    double intersect(const Ray &ray, Color &color, int level) override {
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
            Vector3D floorPoint = intersectingPoint - referencePoint;

            int i = floorPoint.x/tileWidth;
            int j = floorPoint.y/tileWidth;
            Color floorColor;
            
            if((i+j)%2) {
                floorColor.r = 1;
                floorColor.g = 1;
                floorColor.b = 1;
            }
            else {
                floorColor.r = 0;
                floorColor.g = 0;
                floorColor.b = 0;
            }
            color = floorColor*coEfficients[AMB];

            Vector3D normal(0, 0, 1);
            for(auto& light: lights) {
                getRayColor(light, ray, intersectingPoint, normal, coEfficients[DIFF], coEfficients[SPEC], shine, floorColor, color);
            }
            getRecursionColor(ray, intersectingPoint, normal, coEfficients[REF], level, color);
            
            return retValue;
        }
    }
};


