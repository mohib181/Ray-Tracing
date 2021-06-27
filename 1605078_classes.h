#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

using namespace std;

struct point
{
	double x,y,z;
};

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

    string toString() const {
        return "Position: " + position.toString() +
               "\nColor: " + color.toString();
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

    virtual void draw(){}

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

    }
};

class Triangle: public Object{
public:
    Vector3D a, b, c;
    Triangle() = default;

    Triangle(const Vector3D &a, const Vector3D &b, const Vector3D &c)
            : a(a), b(b), c(c) {}

    void draw() override {

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

    }
};

class Data{
public:
    int recursionLevel, imageDimension, totalObjects, totalLights;
    double floorWidth = 1000;
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

        auto* floor = new Floor(floorWidth, tileWidth);
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
