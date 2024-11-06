#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <vector>

typedef unsigned char RGB[3];

typedef struct{
    Vec3f origin;
    Vec3f direction;
} Ray;

Vec3f normalize(const Vec3f& vec) {
    float magnitude = std::sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
    
    if (magnitude == 0) {
        // Handle zero magnitude vector, could return zero vector or throw an error
        return Vec3f(0, 0, 0);
    }
    
    return Vec3f(vec.x / magnitude, vec.y / magnitude, vec.z / magnitude);
}

Vec3f raySphereIntersection(Ray ray, Vec3f center, float radius){
    
}

Vec3f rayIntersection(Ray ray) {

}

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

    return unit_ray;
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
            int colIdx = x / columnWidth;
            image[i++] = 255;
            image[i++] = 0;
            image[i++] = 255;
        }
    }

    write_ppm("test.ppm", image, width, height);

}
