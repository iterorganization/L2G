#include <climits>
#include <memory>
#define TINYBVH_IMPLEMENTATION
#include <tlas.hpp>


unsigned int TLAS::commitMesh(float* vertices, long int n_vertices,
                              unsigned* triangles, long int n_triangles){

    // Create the bvec4 list
    BLAS mesh;
    mesh.triangles.resize(n_triangles*3);

    int tri_p_idx;
    for (int i = 0; i < n_triangles; ++i)
    {
        // i = tri_idx
        tri_p_idx = triangles[3*i];

        // Single triangle
        mesh.triangles[3*i].x = vertices[3*tri_p_idx];
        mesh.triangles[3*i].y = vertices[3*tri_p_idx+1];
        mesh.triangles[3*i].z = vertices[3*tri_p_idx+2];

        tri_p_idx = triangles[3*i+1];

        mesh.triangles[3*i+1].x = vertices[3*tri_p_idx];
        mesh.triangles[3*i+1].y = vertices[3*tri_p_idx+1];
        mesh.triangles[3*i+1].z = vertices[3*tri_p_idx+2];

        tri_p_idx = triangles[3*i+2];

        mesh.triangles[3*i+2].x = vertices[3*tri_p_idx];
        mesh.triangles[3*i+2].y = vertices[3*tri_p_idx+1];
        mesh.triangles[3*i+2].z = vertices[3*tri_p_idx+2];
    }

#ifdef BVH_USEAVX2
    mesh.bvh = std::make_unique<tinybvh::BVH8_CPU>();
    mesh.bvh->Build(mesh.triangles.data(), n_triangles);
#else
    #ifdef BVH_USESSE
    mesh.bvh = std::make_unique<tinybvh::BVH4_CPU>();
    mesh.bvh->BuildAVX(mesh.triangles.data(), n_triangles);
    #else
    mesh.bvh = std::make_unique<tinybvh::BVH>();
    mesh.bvh->Build(mesh.triangles.data(), n_triangles);
    #endif
#endif

    // Add a instance to the instances
    unsigned int idx = m_blases.size();
    m_blases.push_back(std::move(mesh));
    rebuildTLAS();
    return idx;
}

bool TLAS::deleteMesh(unsigned int geom_id){
    if (geom_id>=m_blases.size()){
        return false;
    }
    m_blases.erase(m_blases.begin()+geom_id);
    // Rebuild the TLAS
    rebuildTLAS();
    return true;
}

void TLAS::rebuildTLAS(){
    m_instances.clear();
    m_bvhs_ptrs.clear();

    if (m_blases.empty()) return;

    for (unsigned i = 0; i < m_blases.size(); ++i)
    {
        m_instances.push_back(tinybvh::BLASInstance(i));
        m_bvhs_ptrs.push_back(m_blases[i].bvh.get());
    }

    m_tlas.Build( m_instances.data(), m_instances.size(), m_bvhs_ptrs.data(), m_blases.size());
}

void TLAS::castRay(tinybvh::Ray *ray){
    m_tlas.Intersect(*ray);
}

void TLAS::castRay(float ox, float oy, float oz, float dx, float dy, float dz,
                   float t){
    m_last_hit=false;
    m_ray.O.x = ox;
    m_ray.O.y = oy;
    m_ray.O.z = oz;
    m_ray.D.x = dx;
    m_ray.D.y = dy;
    m_ray.D.z = dz;
    m_ray.hit.t = t;
    m_ray.hit.inst = UINT_MAX;
    m_ray.hit.prim = UINT_MAX;

    castRay(&m_ray);
    m_last_hit = m_ray.hit.inst != UINT_MAX;
}

