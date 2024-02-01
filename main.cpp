//Computational Fabrication Assignment #1
// By David Levin 2014
#include <iostream>
#include <vector>
#include <ctime>
#include "CompFab.h"
#include "Mesh.h"

//Triangle list (global)
typedef std::vector<CompFab::Triangle> TriangleList;

TriangleList g_triangleList;
CompFab::VoxelGrid *g_voxelGrid;

int rayTriangleIntersection(CompFab::Ray &ray, CompFab::Triangle &triangle);

unsigned int numSurfaceIntersections(CompFab::Vec3 &voxelPos, CompFab::Vec3 &dir);

bool loadMesh(char *filename, unsigned int dim);

void saveVoxelsToObj(const char *outfile);

int main(int argc, char **argv) {
    // Load OBJ
    if (argc < 5) {
        std::cout << "Usage: voxelizer <input_mesh_file> <output_mesh_file> <resolution> <multidirection?>\n";
        std::cout << "<input_mesh_file>\tInput mesh file name\n";
        std::cout << "<output_mesh_file>\tOutput mesh file name\n";
        std::cout << "<resolution>\t\tVoxelizer grid resolution, e.g., 16, 32, 64\n";
        std::cout << "<multidirection?>\t0 for uni-directional ray casting\n";
        std::cout << "\t\t\t1 for tri-directional ray casting\n";
        exit(1);
    }

    unsigned int dim = strtoul(argv[3], nullptr, 10);
    unsigned int use_multicast = strtoul(argv[4], nullptr, 10);

    std::cout << "Load Mesh : " << argv[1] << "\n";
    loadMesh(argv[1], dim);


    // Cast ray, check if voxel is inside or outside
    // even number of surface intersections = outside (OUT then IN then OUT)
    // odd number = inside (IN then OUT)
    CompFab::Vec3 voxelPos;
    CompFab::Vec3 direction_mono(1.0, 0.0, 0.0);
    CompFab::Vec3 directions[3] = {
            {1., 0., 0.},
            {0., 1., 0.},
            {0., 0., 1.},
    };

    clock_t start, end;
    double cpu_time_used;

    start = clock();
    /********* ASSIGNMENT *********/
    /* Iterate over all voxels in g_voxelGrid and test whether they are inside our outside of the
     * surface defined by the triangles in g_triangleList */

//    for (CompFab::Triangle &triangle: g_triangleList) {
//        std::cout << "(" << triangle.m_v1.m_x << "," << triangle.m_v1.m_y << "," << triangle.m_v1.m_z << ")\n";
//        std::cout << "(" << triangle.m_v2.m_x << "," << triangle.m_v2.m_y << "," << triangle.m_v2.m_z << ")\n";
//        std::cout << "(" << triangle.m_v3.m_x << "," << triangle.m_v3.m_y << "," << triangle.m_v3.m_z << ")\n";
//        std::cout << "\n";
//    }

    unsigned int nx = g_voxelGrid->m_dimX;
    unsigned int ny = g_voxelGrid->m_dimY;
    unsigned int nz = g_voxelGrid->m_dimZ;

    unsigned int num_hits;

    if (use_multicast) {
        for (CompFab::Vec3 &direction: directions) {
            for (unsigned int i = 0; i < nx; ++i) {
                for (unsigned int j = 0; j < ny; ++j) {
                    for (unsigned int k = 0; k < nz; ++k) {
                        voxelPos = CompFab::Vec3(
                                static_cast<double>(i) / static_cast<double>(nx),
                                static_cast<double>(j) / static_cast<double>(ny),
                                static_cast<double>(k) / static_cast<double>(nz)
                        );

                        num_hits = numSurfaceIntersections(voxelPos, direction);
                        (g_voxelGrid->m_insideArray)[k * (nx * ny) + j * ny + i] |= (num_hits % 2);
                    }
                }
            }
        }
    } else {
        for (unsigned int i = 0; i < nx; ++i) {
            for (unsigned int j = 0; j < ny; ++j) {
                for (unsigned int k = 0; k < nz; ++k) {
                    voxelPos = CompFab::Vec3(
                            static_cast<double>(i) / static_cast<double>(nx),
                            static_cast<double>(j) / static_cast<double>(ny),
                            static_cast<double>(k) / static_cast<double>(nz)
                    );

                    num_hits = numSurfaceIntersections(voxelPos, direction_mono);
                    (g_voxelGrid->m_insideArray)[k * (nx * ny) + j * ny + i] = (num_hits % 2);
                }
            }
        }
    }

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("Voxelizer (Plain) uses %f seconds for %u grid.\n", cpu_time_used, dim);

    // Write out voxel data as obj
    saveVoxelsToObj(argv[2]);
    delete g_voxelGrid;

    return 0;
}

void intersect_tri_1(CompFab::Ray &ray, CompFab::Triangle &triangle,
                     int &out, double &t, double &u, double &v) {
    // From https://www.graphics.cornell.edu/pubs/1997/MT97.pdf
    // Culling to one-sided face
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

void intersect_tri_2(CompFab::Ray &ray, CompFab::Triangle &triangle,
                     int &out, double &t, double &u, double &v) {
    // From https://www.graphics.cornell.edu/pubs/1997/MT97.pdf
    // Two-sided face (with two direction option)
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

    // Comment out for bidirectional ray
    if (t < 0)
        return;
    // End comment

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

//    static size_t cnt = 0;
//    if (out)
//        printf("%zu (%f,%f,%f): %f\n", cnt, ray.m_origin.m_x, ray.m_origin.m_y, ray.m_origin.m_z, t);
//
//    ++cnt;

    return out;
}

//Number of intersections with surface made by a ray originating at voxel and cast in direction.
unsigned int numSurfaceIntersections(CompFab::Vec3 &voxelPos, CompFab::Vec3 &dir) {
    unsigned int numHits = 0;

    /********* ASSIGNMENT *********/
    /* Check and return the number of times a ray cast in direction dir,
     * from voxel center voxelPos intersects the surface */

    CompFab::Ray ray(voxelPos, dir);

    for (CompFab::Triangle &triangle: g_triangleList) {
        numHits += rayTriangleIntersection(ray, triangle);
    }

    return numHits;
}

bool loadMesh(char *filename, unsigned int dim) {
    g_triangleList.clear();

    Mesh *tempMesh = new Mesh(filename, true);

    CompFab::Vec3 v1, v2, v3;

    // Addition to speed up
    g_triangleList.reserve(tempMesh->t.size());

    //copy triangles to global list
    for (unsigned int tri = 0; tri < tempMesh->t.size(); ++tri) {
        v1 = tempMesh->v[tempMesh->t[tri][0]];
        v2 = tempMesh->v[tempMesh->t[tri][1]];
        v3 = tempMesh->v[tempMesh->t[tri][2]];
        g_triangleList.emplace_back(v1, v2, v3);
    }

    //Create Voxel Grid
    CompFab::Vec3 bbMax, bbMin;
    BBox(*tempMesh, bbMin, bbMax);

    //Build Voxel Grid
    double bbX = bbMax[0] - bbMin[0];
    double bbY = bbMax[1] - bbMin[1];
    double bbZ = bbMax[2] - bbMin[2];
    double spacing;

    if (bbX > bbY && bbX > bbZ) {
        spacing = bbX / (double) (dim - 2);
    } else if (bbY > bbX && bbY > bbZ) {
        spacing = bbY / (double) (dim - 2);
    } else {
        spacing = bbZ / (double) (dim - 2);
    }

    CompFab::Vec3 hspacing(0.5 * spacing, 0.5 * spacing, 0.5 * spacing);

    g_voxelGrid = new CompFab::VoxelGrid(bbMin - hspacing, dim, dim, dim, spacing);

    delete tempMesh;

    return true;
}

void saveVoxelsToObj(const char *outfile) {
    Mesh box;
    Mesh mout;
    unsigned int nx = g_voxelGrid->m_dimX;
    unsigned int ny = g_voxelGrid->m_dimY;
    unsigned int nz = g_voxelGrid->m_dimZ;
    double spacing = g_voxelGrid->m_spacing;

    CompFab::Vec3 hspacing(0.5 * spacing, 0.5 * spacing, 0.5 * spacing);

    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                if (!g_voxelGrid->isInside(ii, jj, kk)) {
                    continue;
                }
                CompFab::Vec3 coord(0.5f + ((double) ii) * spacing, 0.5f + ((double) jj) * spacing,
                                    0.5f + ((double) kk) * spacing);
                CompFab::Vec3 box0 = coord - hspacing;
                CompFab::Vec3 box1 = coord + hspacing;
                makeCube(box, box0, box1);
                mout.append(box);
            }
        }
    }

    mout.save_obj(outfile);
}
