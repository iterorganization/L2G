#ifndef ACCELL_TLAS_H
#define ACCELL_TLAS_H

#if defined(_WIN32)
    #ifdef flt_STATIC
        #define TLAS_API
    #elif flt_EXPORTS
        #define TLAS_API __declspec(dllexport)
    #else
        #define TLAS_API __declspec(dllimport)
    #endif
#else
    #ifdef flt_EXPORTS
        #define TLAS_API __attribute__ ((visibility ("default")))
    #else
        #define TLAS_API
    #endif
#endif

#include <tiny_bvh.h>

struct BLAS {
#ifdef BVH_USEAVX2
    std::unique_ptr<tinybvh::BVH8_CPU> bvh;
#else
    #ifdef BVH_USESSE
    std::unique_ptr<tinybvh::BVH4_CPU> bvh;
    #else
    std::unique_ptr<tinybvh::BVH> bvh;
    #endif
#endif
    std::vector<tinybvh::bvhvec4> triangles;
};


/// TLAS
class TLAS_API TLAS
{
private:
    /// Top Level Accelerated structure
    tinybvh::BVH m_tlas;
    std::vector<BLAS> m_blases;

    TLAS(const TLAS&) = delete;
    TLAS& operator=(const TLAS&) = delete;

    TLAS(TLAS&&) = default;
    TLAS& operator=(TLAS&&) = default;

    /// Not used for storage but for building TLAS.
    std::vector<tinybvh::BLASInstance> m_instances;
    std::vector<tinybvh::BVHBase*> m_bvhs_ptrs; // For C-stype **ptr

    // Ray struct, contains the information for hit, primitive hit and
    // geometry hit
    tinybvh::Ray m_ray;
    bool m_last_hit=false;

    void rebuildTLAS();

public:
    TLAS() = default;
    ~TLAS() = default;

    /// Function for loading a triangular mesh into tinybvh. Float precision.
    /// @returns geom_id ID of the loaded geometry. Used for identification
    ///          or removal of the geometry.
    /// @param[in] vertices 1D array of points (p1.x, p1.y, p1.z, p2.x, p2.y,...).
    /// @param[in] n_vertices Number of vertices (len(vertices)/3).
    /// @param[in] triangles 1D array containing all triangles (t1.p1, t1.p2,
    ///                      t1.p3, t2.p1, ...).
    /// @param[in] n_triangles Number of triangles (len(triangles)/3).
    unsigned int commitMesh(float* vertices, long int n_vertices,
                            unsigned* triangles, long int n_triangles);

    /// Deletes the loaded geometry with the assigned id.
    /// @returns successful Returns True if the geometry was successfully
    ///                     removed.
    /// @param[in] geom_id Id of the geometry that we wish to remove.
    bool deleteMesh(unsigned int geom_id);

    /// OpenMP friendly castRay function as the ray tracing structs are
    /// provided by the calling function
    ///
    /// @param[in, out] ray Contains the information of the ray and the hit
    void castRay(tinybvh::Ray *ray);

    /// Function for performing ray tracing. After calling this function you
    /// call checkIfHit to get confirmation if hit happened.
    /// @param[in] ox X component of origin point
    /// @param[in] oy Y component of origin point
    /// @param[in] oz Z component of origin point
    /// @param[in] dx X component of direction vector
    /// @param[in] dy Y component of direction vector
    /// @param[in] dz Z component of direction vector
    /// @param[in] tfar  Ending length along the direction vector to check
    ///                  for intersection
    void castRay(float ox, float oy, float oz, float dx, float dy, float dz,
                 float tfar);

    /// Function that tells you if there is an intersection hit. Call this
    /// after caling castRay.
    /// @returns True if a geometry was hit, false otherwise.
    bool checkIfHit() { return m_last_hit; };

    /// This function returns the geometry ID if an intersection hit occured.
    /// @returns geomId of the geometry that was hit.
    int returnGeomId() { return m_ray.hit.inst; };

    /// This function returns the geometry ID if an intersection hit occured.
    /// @returns primId of the triangle that belongs to the geometry that was
    ///                 hit.
    int returnPrimId() { return m_ray.hit.prim; };

    bool isEmpty(){ return m_blases.empty(); };
};

#endif /*ACCELL_TLAS_H*/