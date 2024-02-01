#include <iostream>
#include <vector>
#include <ctime>
#include "omp.h"
#include "CompFab.h"
#include "Mesh.h"
#include "bvh/bvh.h"

using TriangleList = std::vector<CompFab::Triangle>;

TriangleList g_triangleList;
CompFab::VoxelGrid *g_voxelGrid;

bool loadMesh(char *filename, unsigned int dim);

void saveVoxelsToObj(const char *outfile);

int main(int argc, char **argv) {
    if (argc < 4) {
        std::cout << "Usage: voxelizer_cpu <input_mesh_file> <output_mesh_file> <resolution> <multidirection?>\n";
        std::cout << "<input_mesh_file>\tInput mesh file name\n";
        std::cout << "<output_mesh_file>\tOutput mesh file name\n";
        std::cout << "<resolution>\t\tVoxelizer grid resolution, e.g., 16, 32, 64\n";
        exit(1);
    }

    unsigned int dim = strtoul(argv[3], nullptr, 10);

    std::cout << "Load Mesh : " << argv[1] << "\n";
    loadMesh(argv[1], dim);
    size_t nx = dim, ny = dim, nz = dim;

    double start, end;
    double cpu_time_used;

    start = omp_get_wtime();
    // ************************** BEGIN LOGIC *******************************
    auto bvh = vt::bvh_t<64>::from_triangles(
            g_triangleList.data(),
            g_triangleList.size()
    );

    vt::vec3_t direction = {1., 0., 0.};
    size_t num;
    size_t i, j, k;

#pragma omp parallel default(none) shared(i, j, k, nx, ny, nz, bvh, g_voxelGrid, direction) private(num)
    {
#pragma omp for nowait collapse(3)
        for (i = 0; i < nx; ++i) {
            for (j = 0; j < ny; ++j) {
                for (k = 0; k < nz; ++k) {
                    vt::ray_t ray;
                    vt::utils::vec::assn(ray.direction, direction);

                    ray.origin[0] = static_cast<vt::real_t>(i) / static_cast<vt::real_t>(nx);
                    ray.origin[1] = static_cast<vt::real_t>(j) / static_cast<vt::real_t>(ny);
                    ray.origin[2] = static_cast<vt::real_t>(k) / static_cast<vt::real_t>(nz);
                    ray.t = vt::numeric::infinity;
                    num = bvh.intersect_ray(ray);

                    (g_voxelGrid->m_insideArray)[k * (nx * ny) + j * ny + i] = (num % 2 == 1);
//                printf("(%f,%f,%f): %f, %zu\n", ray.origin[0], ray.origin[1], ray.origin[2], ray.t, num);
                }
            }
        }
    }
    // **************************  END LOGIC  *******************************
    end = omp_get_wtime();
    cpu_time_used = ((double) (end - start));

    printf("Voxelizer (BVH+AVX512+OPENMP) uses %f seconds for %u grid.\n", cpu_time_used, dim);

    saveVoxelsToObj(argv[2]);
    delete g_voxelGrid;

    return 0;
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
