#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <vector>

typedef unsigned char RGB[3];

using namespace parser;
typedef struct{
    Vec3f origin;
    Vec3f direction;
} Ray;

using namespace parser;
Vec3f normalize(const Vec3f& vec) {
    float magnitude = std::sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
    Vec3f zero = {0,0,0};
    Vec3f normalize = {vec.x / magnitude, vec.y / magnitude, vec.z / magnitude};
    if (magnitude == 0) {
        // Handle zero magnitude vector, could return zero vector or throw an error
        return zero;
    }
    return normalize;
}

float dot(const Vec3f& u, const Vec3f& v) {
    return (u.x*v.x + u.y*v.y + u.z*v.z);
}

Vec3f cross(const Vec3f& u, const Vec3f& v) {
    return {u.y*v.z - u.z*v.y, u.z*v.x - u.x*v.z, u.x*v.y - u.y*v.x};
}

Vec3f minus(const Vec3f& u, const Vec3f& v) {
    return {u.x - v.x, u.y - v.y, u.z - v.z};
}

using namespace parser;
Vec3f raySphereIntersection(Ray ray, Vec3f center, float radius){
    Vec3f oc = minus(ray.origin, center);
    float a = dot(ray.direction, ray.direction);
    float b = 2.0 * dot(oc, ray.direction);
    float c = dot(oc, oc) - radius * radius;
    float discriminant = b*b - 4*a*c;
    if (discriminant < 0) {
        return {0,0,0};
    }
    float t = (-b - std::sqrt(discriminant)) / (2.0 * a);
    if (t < 0) {
        return {0,0,0};
    }
    Vec3f intersection = {ray.origin.x + t * ray.direction.x, ray.origin.y + t * ray.direction.y, ray.origin.z + t * ray.direction.z};
    return intersection;
}

using namespace parser;
Vec3f rayTriangleIntersection(Ray ray, Vec3f triA, Vec3f triB, Vec3f triC){
    Vec3f edge1 = minus(triB, triA);
    Vec3f edge2 = minus(triC, triA);
    Vec3f h = cross(ray.direction, edge2);
    float a = dot(edge1, h);
    if (a > -0.00001 && a < 0.00001) {
        return {0,0,0};
    }
    float f = 1.0 / a;
    Vec3f s = minus(ray.origin, triA);
    float u = f * dot(s, h);
    if (u < 0.0 || u > 1.0) {
        return {0,0,0};
    }
    Vec3f q = cross(s, edge1);
    float v = f * dot(ray.direction, q);
    if (v < 0.0 || u + v > 1.0) {
        return {0,0,0};
    }
    float t = f * dot(edge2, q);
    if (t > 0.00001) {
        Vec3f intersection = {ray.origin.x + t * ray.direction.x, ray.origin.y + t * ray.direction.y, ray.origin.z + t * ray.direction.z};
        return intersection;
    }
    return {0,0,0};
}

using namespace parser;
Vec3f rayIntersection(Ray ray) {}

using namespace parser;
Ray generateRay(const Camera& cam, int i, int j) {
    int z_sign;
    if(cam.gaze.z < 0) {z_sign = -1;}
    else {z_sign = 1;}

    float left = cam.near_plane.x;
    float right = cam.near_plane.y;
    float bottom = cam.near_plane.z;
    float top = cam.near_plane.w;

    float pixel_width = abs((right - left) / cam.image_width);
    float pixel_height = abs((top - bottom) / cam.image_height);

    float target_x = left + ((i+.5)*pixel_width);
    float target_y = top - ((j+.5)*pixel_height);
    float target_z = cam.position.z + (z_sign * cam.near_distance);

    Vec3f unit_ray;
    unit_ray.x = target_x - cam.position.x;
    unit_ray.y = target_y - cam.position.y;
    unit_ray.z = target_z - cam.position.z;

    unit_ray = normalize(unit_ray);

    return {cam.position, unit_ray};
}

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.
    //
    // Normally, you would be running your ray tracing
    // code here to produce the desired image.

    const RGB BAR_COLOR[8] =
    {
        { 255, 255, 255 },  // 100% White
        { 255, 255,   0 },  // Yellow
        {   0, 255, 255 },  // Cyan
        {   0, 255,   0 },  // Green
        { 255,   0, 255 },  // Magenta
        { 255,   0,   0 },  // Red
        {   0,   0, 255 },  // Blue
        {   0,   0,   0 },  // Black
    };

    int width = 640, height = 480;
    int columnWidth = width / 8;

    unsigned char* image = new unsigned char [width * height * 3];

    int i = 0;
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            //int colIdx = x / columnWidth;
            
        }
    }

    write_ppm("test.ppm", image, width, height);

}
