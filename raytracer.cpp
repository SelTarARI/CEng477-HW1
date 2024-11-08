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

//using namespace parser;
//Vec3f rayIntersection(Ray ray) {}

using namespace parser;
Ray generateRay(Camera camera, int x, int y) {
    Vec3f w = normalize(minus(camera.position, camera.gaze));
    Vec3f u = normalize(cross(camera.up, w));
    Vec3f v = cross(w, u);
    float aspectRatio = (float)camera.image_width / camera.image_height;
    float top = camera.near_plane.w;
    float right = top * aspectRatio;
    float bottom = -top;
    float left = -right;
    float uCoord = left + (right - left) * (x + 0.5) / camera.image_width;
    float vCoord = bottom + (top - bottom) * (y + 0.5) / camera.image_height;
    Vec3f direction = plus(timesScalar(w, -camera.near_distance), plus(timesScalar(u, uCoord), timesScalar(v, vCoord)));
    return {camera.position, normalize(direction)};
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



    for (int j =0; j < scene.cameras.size(); j++) {
        int width = scene.cameras[j].image_width;
        int height = scene.cameras[j].image_height;
        unsigned char* image = new unsigned char[width * height * 3];
        std::cout << (scene.cameras[j].image_name) << std::endl;

        int i = 0;
        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                bool hit=false;

                Ray ray = generateRay(scene.cameras[0], x, y);
                //printf("Ray Origin: %f %f %f\n", ray.origin.x, ray.origin.y, ray.origin.z);
                //printf("Ray Direction: %f %f %f\n", ray.direction.x, ray.direction.y, ray.direction.z);
                Vec3f coordinate = raySphereIntersection(ray, scene.vertex_data[scene.spheres[0].center_vertex_id], scene.spheres[0].radius, hit);
                //Vec3f triA = scene.vertex_data[scene.meshes[0].faces[0].v0_id];
                //Vec3f triB = scene.vertex_data[scene.meshes[0].faces[0].v1_id];
                //Vec3f triC = scene.vertex_data[scene.meshes[0].faces[0].v2_id];
                //Vec3f coordinateTri = rayTriangleIntersection(ray, triA, triB, triC, hit);

                if (hit) { //red
                    printf("Hit\n");
                    image[i++] = 255;
                    image[i++] = 0;
                    image[i++] = 0;
                } else { //blue
                    image[i++] = 0;
                    image[i++] = 0;
                    image[i++] = 255;
                }
            }
        }
        write_ppm("test.ppm", image, width, height);
    }


}
