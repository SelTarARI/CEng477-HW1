#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <vector>

using namespace parser;

typedef struct
{
    Vec3f origin;
    Vec3f direction;
    int depth = 0;
} Ray;

typedef struct
{
    Vec3f hitpoint;
    Vec3f normal;
    Material material;
    bool hitHappened = false;
} Hit;

Vec3f colorFinder(Scene scene, Ray ray);

Vec3f normalize(const Vec3f &vec)
{
    float magnitude = std::sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
    Vec3f zero = {0, 0, 0};
    if (magnitude == 0)
    {
        // Handle zero magnitude vector, could return zero vector or throw an error
        return zero;
    }
    Vec3f normalize = {vec.x / magnitude, vec.y / magnitude, vec.z / magnitude};

    return normalize;
}

float dot(const Vec3f &u, const Vec3f &v)
{
    return (u.x * v.x + u.y * v.y + u.z * v.z);
}

Vec3f cross(const Vec3f &u, const Vec3f &v)
{
    return {u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
}

Vec3f minus(const Vec3f &u, const Vec3f &v)
{
    return {u.x - v.x, u.y - v.y, u.z - v.z};
}

Vec3f plus(const Vec3f &u, const Vec3f &v)
{
    return {u.x + v.x, u.y + v.y, u.z + v.z};
}

Vec3f timesScalar(const Vec3f &u, float t)
{
    return {u.x * t, u.y * t, u.z * t};
}

Vec3f timesForColor(const Vec3f &u, const Vec3f &v)
{
    return {u.x * v.x, u.y * v.y, u.z * v.z};
}

Vec3f negate(const Vec3f &u)
{
    return {-u.x, -u.y, -u.z};
}

float distance(Vec3f a, Vec3f b)
{
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2));
}

Vec3f raySphereIntersection(Ray ray, Vec3f center, float radius, bool &hit)
{
    Vec3f oc = minus(ray.origin, center);
    float a = dot(ray.direction, ray.direction);
    float b = 2.0 * dot(oc, ray.direction);
    float c = dot(oc, oc) - radius * radius;
    float discriminant = b * b - 4 * a * c;
    if (discriminant < 0)
    {
        hit = false;
        return {0, 0, 0};
    }
    float t = (-b - std::sqrt(discriminant)) / (2.0 * a);
    if (t < 0)
    {
        hit = false;
        return {0, 0, 0};
    }
    Vec3f intersection = plus(ray.origin, timesScalar(ray.direction, t));
    hit = true;
    return intersection;
}

Vec3f rayTriangleIntersection(Ray ray, Vec3f triA, Vec3f triB, Vec3f triC, bool &hit)
{
    Vec3f edge1 = minus(triB, triA);
    Vec3f edge2 = minus(triC, triA);
    Vec3f h = cross(ray.direction, edge2);
    float a = dot(edge1, h);
    if (a > -0.00001 && a < 0.00001)
    {
        hit = false;
        return {0, 0, 0};
    }
    float f = 1.0 / a;
    Vec3f s = minus(ray.origin, triA);
    float u = f * dot(s, h);
    if (u < 0 || u > 1)
    {
        hit = false;
        return {0, 0, 0};
    }
    Vec3f q = cross(s, edge1);
    float v = f * dot(ray.direction, q);
    if (v < 0 || u + v > 1)
    {
        hit = false;
        return {0, 0, 0};
    }
    float t = f * dot(edge2, q);
    if (t > 0.00001)
    {
        hit = true;
        return plus(ray.origin, timesScalar(ray.direction, t));
    }
    hit = false;
    return {0, 0, 0};
}

Vec3f rayMeshIntersection(Ray ray, Scene scene, int meshIndex, bool &hit)
{
    std::vector<Vec3f> hitList;
    Mesh mesh = scene.meshes[meshIndex];
    Vec3f intersection = {0, 0, 0};
    for (int i = 0; i < mesh.faces.size(); i++)
    {
        Face face = mesh.faces[i];
        Vec3f v0 = scene.vertex_data[face.v0_id - 1];
        Vec3f v1 = scene.vertex_data[face.v1_id - 1];
        Vec3f v2 = scene.vertex_data[face.v2_id - 1];
        Vec3f coordinate = rayTriangleIntersection(ray, v0, v1, v2, hit);
        if (hit)
        {
            hitList.push_back(coordinate);
        }
    }
    if (hitList.size() == 0)
    {
        hit = false;
        return {0, 0, 0};
    }
    hit = true;
    for (int i = 0; i < hitList.size(); i++)
    {
        if (i == 0)
        {
            intersection = hitList[i];
        }
        else
        {
            if (distance(ray.origin, hitList[i]) < distance(ray.origin, intersection))
            {
                intersection = hitList[i];
            }
        }
    }
    return intersection;
}

Ray generateRay(Camera camera, int x, int y)
{
    Vec3f w = normalize(minus(camera.position, camera.gaze));
    Vec3f u = normalize(cross(camera.up, w));
    Vec3f v = cross(w, u);
    float aspectRatio = (float)camera.image_width / camera.image_height;
    float top = camera.near_plane.w;
    float right = top * aspectRatio;
    float bottom = -top;
    float left = -right;
    float uCoord = left + (right - left) * (x + 0.5) / camera.image_width;
    float vCoord = top - (top - bottom) * (y + 0.5) / camera.image_height;
    Vec3f direction = plus(timesScalar(w, -camera.near_distance), plus(timesScalar(u, uCoord), timesScalar(v, vCoord)));
    Ray out;
    out.origin = camera.position;
    out.direction = normalize(direction);
    return out;
}

Vec3f reflectionColor(Scene scene, Ray ray, Hit hit)
{
    Vec3f output = {0, 0, 0};
    if (hit.material.is_mirror)
    {
        Ray reflectionRay;
        reflectionRay.origin = hit.hitpoint;
        reflectionRay.direction = normalize(plus(ray.direction, timesScalar(hit.normal, 2 * dot(negate(ray.direction), hit.normal))));
        reflectionRay.depth = ray.depth + 1;
        output = timesForColor(colorFinder(scene, reflectionRay), hit.material.mirror);
    }
    return output;
}

Vec3f shadeComputation(Scene scene, Ray ray, Hit hit)
{
    Vec3f output = timesForColor(hit.material.ambient, scene.ambient_light);

    if (hit.material.is_mirror)
    {
        output = plus(output, reflectionColor(scene, ray, hit));
    }

    for (int k = 0; k < scene.point_lights.size(); k++)
    {
        Vec3f lightPosition = scene.point_lights[k].position;
        Vec3f lightDirection = normalize(minus(lightPosition, hit.hitpoint));

        // shadow ray
        Ray shadowRay;
        shadowRay.origin = hit.hitpoint;
        shadowRay.direction = lightDirection;
        bool shadowHit = false;
        bool controller = false;
        for (int l = 0; l < scene.spheres.size(); l++)
        {
            Vec3f center = scene.vertex_data[scene.spheres[l].center_vertex_id - 1];
            float radius = scene.spheres[l].radius;
            Vec3f coordinate = raySphereIntersection(shadowRay, center, radius, shadowHit);
            if (shadowHit && distance(hit.hitpoint, coordinate) < distance(hit.hitpoint, lightPosition) && distance(hit.hitpoint, coordinate) > 0.001)
            {
                controller = true;
                break;
            }
        }
        if (controller)
        {
            continue;
        }
        for (int l = 0; l < scene.triangles.size(); l++)
        {
            Vec3f v0 = scene.vertex_data[scene.triangles[l].indices.v0_id - 1];
            Vec3f v1 = scene.vertex_data[scene.triangles[l].indices.v1_id - 1];
            Vec3f v2 = scene.vertex_data[scene.triangles[l].indices.v2_id - 1];
            Vec3f coordinate = rayTriangleIntersection(shadowRay, v0, v1, v2, shadowHit);
            if (shadowHit) //&& distance(hit.hitpoint, coordinate) < distance(hit.hitpoint, lightPosition) && distance(hit.hitpoint, coordinate) > 0.001)
            {
                controller = true;
                break;
            }
        }
        if (controller)
        {
            continue;
        }
        for (int l = 0; l < scene.meshes.size(); l++)
        {
            Vec3f v0 = scene.vertex_data[scene.meshes[l].faces[0].v0_id - 1];
            Vec3f v1 = scene.vertex_data[scene.meshes[l].faces[0].v1_id - 1];
            Vec3f v2 = scene.vertex_data[scene.meshes[l].faces[0].v2_id - 1];
            Vec3f coordinate = rayMeshIntersection(shadowRay, scene, l, shadowHit);
            if (shadowHit) //&& distance(hit.hitpoint, coordinate) < distance(hit.hitpoint, lightPosition) && distance(hit.hitpoint, coordinate) > 0.001)
            {
                controller = true;
                break;
            }
        }
        if (controller)
        {
            continue;
        }

        // point light diffuse
        Vec3f lightIntensity = scene.point_lights[k].intensity;
        Vec3f diffuseCoefficient = hit.material.diffuse;
        Vec3f normal = hit.normal;

        float normalDotLight = dot(normal, lightDirection);
        float cosTheta = std::max(0.0f, normalDotLight);
        float oneOverDistanceSquare = 1 / pow(distance(lightPosition, hit.hitpoint), 2);
        Vec3f radiance = timesForColor(timesScalar(lightIntensity, oneOverDistanceSquare * cosTheta), diffuseCoefficient);
        output = plus(output, radiance);

        // specular reflection
        Vec3f viewDirection = normalize(minus(ray.origin, hit.hitpoint));
        Vec3f halfVector = normalize(plus(lightDirection, viewDirection));
        float normalDotHalf = dot(normal, halfVector);
        float cosAlpha = std::max(0.0f, normalDotHalf);
        Vec3f specularCoefficient = hit.material.specular;
        float phongExponent = hit.material.phong_exponent;
        Vec3f specularRadiance = timesForColor(timesScalar(lightIntensity, oneOverDistanceSquare * pow(cosAlpha, phongExponent)), specularCoefficient);
        output = plus(output, specularRadiance);
    }
    return output;
}

Vec3f colorFinder(Scene scene, Ray ray)
{
    if (ray.depth > scene.max_recursion_depth)
    {
        return {0, 0, 0};
    }

    Hit hitSphere;
    Hit hitSphereFinal;
    bool sphereHitFirst = true;

    Hit hitTriangle;
    Hit hitTriangleFinal;
    bool triangleHitFirst = true;

    Hit hitMesh;
    Hit hitMeshFinal;
    bool meshHitFirst = true;

    for (int k = 0; k < scene.spheres.size(); k++)
    {
        Vec3f center = scene.vertex_data[scene.spheres[k].center_vertex_id - 1];
        float radius = scene.spheres[k].radius;
        Vec3f coordinate = raySphereIntersection(ray, center, radius, hitSphere.hitHappened);
        if (hitSphere.hitHappened)
        {
            hitSphere.hitpoint = coordinate;
            hitSphere.normal = normalize(minus(coordinate, center));
            hitSphere.material = scene.materials[scene.spheres[k].material_id - 1];
            if (sphereHitFirst)
            {
                hitSphereFinal.hitpoint = hitSphere.hitpoint;
                hitSphereFinal.normal = hitSphere.normal;
                hitSphereFinal.material = hitSphere.material;
                hitSphereFinal.hitHappened = true;
                sphereHitFirst = false;
            }
            else
            {
                if (distance(ray.origin, hitSphere.hitpoint) < distance(ray.origin, hitSphereFinal.hitpoint))
                {
                    hitSphereFinal.hitpoint = hitSphere.hitpoint;
                    hitSphereFinal.normal = hitSphere.normal;
                    hitSphereFinal.material = hitSphere.material;
                    hitSphereFinal.hitHappened = true;
                }
            }
        }
    }

    for (int k = 0; k < scene.triangles.size(); k++)
    {
        Vec3f v0 = scene.vertex_data[scene.triangles[k].indices.v0_id - 1];
        Vec3f v1 = scene.vertex_data[scene.triangles[k].indices.v1_id - 1];
        Vec3f v2 = scene.vertex_data[scene.triangles[k].indices.v2_id - 1];
        Vec3f coordinate = rayTriangleIntersection(ray, v0, v1, v2, hitTriangle.hitHappened);
        if (hitTriangle.hitHappened)
        {
            hitTriangle.hitpoint = coordinate;
            hitTriangle.normal = normalize(cross(minus(v1, v0), minus(v2, v0)));
            hitTriangle.material = scene.materials[scene.triangles[k].material_id - 1];
            if (triangleHitFirst)
            {
                hitTriangleFinal.hitpoint = hitTriangle.hitpoint;
                hitTriangleFinal.normal = hitTriangle.normal;
                hitTriangleFinal.material = hitTriangle.material;
                hitTriangleFinal.hitHappened = true;
                triangleHitFirst = false;
            }
            else
            {
                if (distance(ray.origin, hitTriangle.hitpoint) < distance(ray.origin, hitTriangleFinal.hitpoint))
                {
                    hitTriangleFinal.hitpoint = hitTriangle.hitpoint;
                    hitTriangleFinal.normal = hitTriangle.normal;
                    hitTriangleFinal.material = hitTriangle.material;
                    hitTriangleFinal.hitHappened = true;
                }
            }
        }
    }

    for (int k = 0; k < scene.meshes.size(); k++)
    {
        Vec3f v0 = scene.vertex_data[scene.meshes[k].faces[0].v0_id - 1];
        Vec3f v1 = scene.vertex_data[scene.meshes[k].faces[0].v1_id - 1];
        Vec3f v2 = scene.vertex_data[scene.meshes[k].faces[0].v2_id - 1];
        Vec3f coordinate = rayMeshIntersection(ray, scene, k, hitMesh.hitHappened);
        if (hitMesh.hitHappened)
        {
            hitMesh.hitpoint = coordinate;
            hitMesh.normal = normalize(cross(minus(v1, v0), minus(v2, v0)));
            hitMesh.material = scene.materials[scene.meshes[k].material_id - 1];
            if (meshHitFirst)
            {
                hitMeshFinal.hitpoint = hitMesh.hitpoint;
                hitMeshFinal.normal = hitMesh.normal;
                hitMeshFinal.material = hitMesh.material;
                hitMeshFinal.hitHappened = true;
                meshHitFirst = false;
            }
            else
            {
                if (distance(ray.origin, hitMesh.hitpoint) < distance(ray.origin, hitMeshFinal.hitpoint))
                {
                    hitMeshFinal.hitpoint = hitMesh.hitpoint;
                    hitMeshFinal.normal = hitMesh.normal;
                    hitMeshFinal.material = hitMesh.material;
                    hitMeshFinal.hitHappened = true;
                }
            }
        }
    }

    Hit finalHit = hitSphereFinal;

    if (hitTriangleFinal.hitHappened & (!finalHit.hitHappened || distance(ray.origin, hitTriangleFinal.hitpoint) < distance(ray.origin, finalHit.hitpoint)))
    {
        finalHit = hitTriangleFinal;
    }
    if (hitMeshFinal.hitHappened & (!finalHit.hitHappened || distance(ray.origin, hitMeshFinal.hitpoint) < distance(ray.origin, finalHit.hitpoint)))
    {
        finalHit = hitMeshFinal;
    }

    if (finalHit.hitHappened)
    {
        return shadeComputation(scene, ray, finalHit);
    }

    else if (ray.depth == 0)
    {
        Vec3i color = scene.background_color;
        return {float(color.x), float(color.y), float(color.z)};
    }
    else
    {
        return {0, 0, 0};
    }
}

int main(int argc, char *argv[])
{
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    for (int j = 0; j < scene.cameras.size(); j++)
    {
        int width = scene.cameras[j].image_width;
        int height = scene.cameras[j].image_height;
        unsigned char *image = new unsigned char[width * height * 3];

        int i = 0;
        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                Ray ray = generateRay(scene.cameras[j], x, y);
                Vec3f color = colorFinder(scene, ray);
                image[i++] = color.x > 255 ? 255 : color.x;
                image[i++] = color.y > 255 ? 255 : color.y;
                image[i++] = color.z > 255 ? 255 : color.z;
            }
        }
        write_ppm(scene.cameras[j].image_name.c_str(), image, width, height);
    }
}
