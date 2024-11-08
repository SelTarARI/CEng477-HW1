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

Vec3f plus(const Vec3f& u, const Vec3f& v) {
    return {u.x + v.x, u.y + v.y, u.z + v.z};
}

Vec3f timesScalar(const Vec3f& u, float t) {
    return {u.x * t, u.y * t, u.z * t};
}

Vec3f negate(const Vec3f& u) {
    return {-u.x, -u.y, -u.z};
}

using namespace parser;
Vec3f raySphereIntersection(Ray ray, Vec3f center, float radius, bool& hit) {
    Vec3f oc = minus(ray.origin, center);
    float a = dot(ray.direction, ray.direction);
    float b = 2.0 * dot(oc, ray.direction);
    float c = dot(oc, oc) - radius * radius;
    float discriminant = b*b - 4*a*c;
    if (discriminant < 0) {
        hit = false;
        return {0,0,0};
    }
    float t = (-b - std::sqrt(discriminant)) / (2.0 * a);
    if (t < 0) {
        hit = false;
        return {0,0,0};
    }
    Vec3f intersection = plus(ray.origin, timesScalar(ray.direction, t));
    hit = true;
    return intersection;
}

using namespace parser;
Vec3f rayTriangleIntersection(Ray ray, Vec3f triA, Vec3f triB, Vec3f triC, bool& hit){
    Vec3f edge1 = minus(triB, triA);
    Vec3f edge2 = minus(triC, triA);
    Vec3f n = cross(edge1, edge2);

    float nDotRayDirection = dot(n, ray.direction);
    if (nDotRayDirection == 0) {
        hit = false;
        return {0,0,0};
    }

    float d = dot(n, triA);
    float t = -(dot(n, ray.origin) + d) / nDotRayDirection;
    if (t < 0) {
        hit = false;
        return {0,0,0};
    }

    Vec3f p = plus(ray.origin, timesScalar(ray.direction, t));
    Vec3f c;
    Vec3f edge3 = minus(triC, triB);
    Vec3f p_triB = minus(p, triB);
    c = cross(edge3, p_triB);
    if (dot(n, c) < 0) {
        hit = false;
        return {0,0,0};
    }

    Vec3f edge2Minus = negate(edge2);
    Vec3f p_triC = minus(p, triC);
    c = cross(edge2Minus, p_triC);
    if (dot(n, c) < 0) {
        hit = false;
        return {0,0,0};
    }

    hit = true;
    return p;
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

    for (int j =0; j < scene.cameras.size(); j++) {
        std::cout << scene.cameras[0].image_name << std::endl;
    }

    int i = 0;
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            //int colIdx = x / columnWidth;
            bool hit;
            Ray ray = generateRay(scene.cameras[0], x, y);
            Vec3f coordinate = raySphereIntersection(ray, scene.vertex_data[scene.spheres[0].center_vertex_id], scene.spheres[0].radius, hit);
            if (coordinate.x != 0 && coordinate.y != 0 && coordinate.z != 0) {
                image[i++] = 255;
                image[i++] = 255;
                image[i++] = 255;
            } else {
                image[i++] = 0;
                image[i++] = 0;
                image[i++] = 255;
            }
        }
    }

    write_ppm("test.ppm", image, width, height);

}
