#include <iostream>
#include <vector>
#include "CompFab.h"
#include "Mesh.h"

__attribute__((always_inline))
void intersect_tri_1(CompFab::Ray &ray, CompFab::Triangle &triangle,
                     int &out, double &t, double &u, double &v) {
    // From https://www.graphics.cornell.edu/pubs/1997/MT97.pdf
    CompFab::Vec3 E1, E2, T, P, Q;
    double det, inv_det;

    out = 0;

    E1 = triangle.m_v2 - triangle.m_v1;
    E2 = triangle.m_v3 - triangle.m_v1;

    P = ray.m_direction % E2;

    det = E1 * P;

    if (det < EPSILON)
        return;

    T = ray.m_origin - triangle.m_v1;

    u = T * P;
    if (u < 0 || u > det)
        return;

    Q = T % E1;

    v = ray.m_direction * Q;
    if (v < 0 || u + v > det)
        return;

    t = E2 * Q;
    inv_det = 1.0 / det;

    t *= inv_det;
    u *= inv_det;
    v *= inv_det;

    out = 1;
}

__attribute__((always_inline))
void intersect_tri_2(CompFab::Ray &ray, CompFab::Triangle &triangle,
                     int &out, double &t, double &u, double &v) {
    // From https://www.graphics.cornell.edu/pubs/1997/MT97.pdf
    CompFab::Vec3 E1, E2, T, P, Q;
    double det, inv_det;

    out = 0;

    E1 = triangle.m_v2 - triangle.m_v1;
    E2 = triangle.m_v3 - triangle.m_v1;

    P = ray.m_direction % E2;

    det = E1 * P;

    if (det > -EPSILON && det < EPSILON)
        return;

    inv_det = 1.0 / det;

    T = ray.m_origin - triangle.m_v1;

    u = (T * P) * inv_det;
    if (u < 0.0 || u > 1.0)
        return;

    Q = T % E1;

    v = (ray.m_direction * Q) * inv_det;
    if (v < 0.0 || u + v > 1.0)
        return;

    t = (E2 * Q) * inv_det;

    if (t < 0)
        return;

    out = 1;
}

//Ray-Triangle intersection
//Returns 1 if triangle and ray intersect, 0 otherwise
int rayTriangleIntersection(CompFab::Ray &ray, CompFab::Triangle &triangle) {
    /********* ASSIGNMENT *********/
    /* Ray-Triangle intersection test: Return 1 if ray intersects triangle,
     * 0 otherwise */

    int out;
    double t, u, v;

    intersect_tri_2(ray, triangle, out, t, u, v);

    return out;
}

int main(int argc, char **argv) {
    using CompFab::Vec3;
    using CompFab::Triangle;
    using CompFab::Ray;

    Vec3 v1(-1, 1, 0);
    Vec3 v2(-1, -1, 0);
    Vec3 v3(1, 0, 0);

    Triangle tri(v1, v2, v3);

    for (double i = -2; i < 2;) {
        for (double j = -2; j < 2;) {
            Vec3 o(i, j, 2);
            Vec3 d(0, 0, -0.00001);
            Ray ray(o, d);

            std::cout << rayTriangleIntersection(ray, tri) << " ";

            j += 0.2;
        }

        std::cout << "\n";

        i += 0.2;
    }
    return 0;
}
