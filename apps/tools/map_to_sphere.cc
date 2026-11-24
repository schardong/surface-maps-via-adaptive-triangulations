/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Maps a single genus-0 triangle mesh to the unit sphere.
 * The output is the mesh with the same connectivity, whose vertices lie on the sphere.
 * The map is piecewise linear on the faces (each face is mapped onto a spherical triangle).
 *
 * Usage: map_to_sphere <mesh> [output.obj]
 */

#include <SurfaceMaps/MultiRes/MultiResSphereEmbedding.hh>
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeSphereEmbedding.hh>
#include <SurfaceMaps/Utils/Genus.hh>
#include <SurfaceMaps/Utils/IO.hh>
#include <SurfaceMaps/Utils/MeshNormalization.hh>

#include <iostream>

using namespace SurfaceMaps;

int main(int argc, char** argv)
{
    if (argc < 2 || argc > 3)
    {
        std::cerr << "usage: " << argv[0] << " <mesh> [output.obj]" << std::endl;
        return 1;
    }

    const fs::path mesh_path = argv[1];
    const fs::path default_out_path = fs::path(mesh_path).replace_extension("").string() + "_sphere.obj";
    const fs::path out_path = (argc == 3) ? fs::path(argv[2]) : default_out_path;

    if (fs::exists(out_path))
    {
        std::cerr << "Output file " << out_path << " already exists. Remove it or choose another path." << std::endl;
        return 1;
    }

    TriMesh mesh = read_mesh(mesh_path);
    if (mesh.n_vertices() == 0)
    {
        std::cerr << "Mesh " << mesh_path << " has no vertices." << std::endl;
        return 1;
    }

    // The sphere embedding only exists for closed genus-0 surfaces
    if (!closed_surface(mesh))
    {
        std::cerr << "Mesh " << mesh_path << " is not closed." << std::endl;
        return 1;
    }
    if (genus(mesh) != 0)
    {
        std::cerr << "Mesh " << mesh_path << " has genus " << genus(mesh) << " (expected 0)." << std::endl;
        return 1;
    }

    // Same normalization as init_map (the mapping is invariant to it)
    center_mesh(mesh);
    normalize_surface_area(mesh);

    // Bijective sphere embedding (no landmarks needed for a single mesh)
    ExternalProperty<VH, Vec3d> embedding = multi_res_sphere_embedding(mesh);
    if (!sphere_embedding_bijective(mesh, embedding))
    {
        std::cerr << "Could not compute a bijective sphere embedding for " << mesh_path << "." << std::endl;
        return 1;
    }

    write_embedding(mesh, embedding, out_path);
    std::cout << "Wrote sphere map of " << mesh_path << " to " << out_path << std::endl;
    return 0;
}
