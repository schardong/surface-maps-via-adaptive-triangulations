/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Creates landmarks from the extreme vertices of a mesh along its principal axes.
 *   e1: direction of largest spread (largest covariance eigenvalue)
 *   e2: direction of second largest spread, orthogonal to e1
 *   e3: e1 x e2 (right-handed frame)
 * For each axis the vertex with minimal and maximal projection is selected,
 * giving 6 landmarks in the fixed order:
 *   e1min, e1max, e2min, e2max, e3min, e3max
 * The order matters: landmark i of one mesh is matched to landmark i of the other.
 *
 * Eigenvectors are only defined up to sign, which would swap min and max between
 * meshes. Each axis is therefore oriented so that the projected point cloud has
 * positive skewness (its longer tail points in the positive direction).
 *
 * Usage: create_landmarks <mesh> [output.pinned]
 */

#include <SurfaceMaps/Utils/IO.hh>

#include <Eigen/Dense>

#include <iostream>
#include <set>

using namespace SurfaceMaps;

namespace
{

Eigen::Vector3d to_eigen(const TriMesh::Point& _p)
{
    return Eigen::Vector3d(_p[0], _p[1], _p[2]);
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2 || argc > 3)
    {
        std::cerr << "usage: " << argv[0] << " <mesh> [output.pinned]" << std::endl;
        return 1;
    }

    const fs::path mesh_path = argv[1];
    const fs::path out_path = (argc == 3) ? fs::path(argv[2]) : fs::path(mesh_path).replace_extension(".pinned");

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

    // Centroid and covariance of the vertex positions
    std::vector<Eigen::Vector3d> points;
    points.reserve(mesh.n_vertices());
    Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
    for (const VH vh : mesh.vertices())
    {
        points.push_back(to_eigen(mesh.point(vh)));
        centroid += points.back();
    }
    centroid /= double(points.size());

    Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
    for (const Eigen::Vector3d& p : points)
    {
        const Eigen::Vector3d d = p - centroid;
        covariance += d * d.transpose();
    }
    covariance /= double(points.size());

    // Eigenvalues are sorted ascending: column 2 is the largest spread, column 1 the second largest
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covariance);
    if (solver.info() != Eigen::Success)
    {
        std::cerr << "Eigen decomposition of the covariance matrix failed." << std::endl;
        return 1;
    }
    Eigen::Vector3d e1 = solver.eigenvectors().col(2);
    Eigen::Vector3d e2 = solver.eigenvectors().col(1);

    // Orient e1 and e2 by positive skewness of the projections
    const auto skewness = [&](const Eigen::Vector3d& _axis) {
        double s = 0.0;
        for (const Eigen::Vector3d& p : points)
        {
            const double t = (p - centroid).dot(_axis);
            s += t * t * t;
        }
        return s;
    };
    if (skewness(e1) < 0.0)
        e1 = -e1;
    if (skewness(e2) < 0.0)
        e2 = -e2;

    // Third axis completes a right-handed frame
    const Eigen::Vector3d e3 = e1.cross(e2).normalized();

    const Eigen::Vector3d axes[3] = { e1, e2, e3 };
    const char* axis_names[3] = { "e1", "e2", "e3" };

    std::cout << "Principal axes (centered at " << centroid.transpose() << "):" << std::endl;
    for (int a = 0; a < 3; ++a)
        std::cout << "  " << axis_names[a] << " = (" << axes[a].transpose() << ")" << std::endl;

    // Extreme vertices along each principal axis
    std::vector<VH> landmarks;
    for (int a = 0; a < 3; ++a)
    {
        VH vmin = *mesh.vertices_begin();
        VH vmax = vmin;
        double tmin = std::numeric_limits<double>::max();
        double tmax = std::numeric_limits<double>::lowest();
        for (const VH vh : mesh.vertices())
        {
            const double t = (to_eigen(mesh.point(vh)) - centroid).dot(axes[a]);
            if (t < tmin)
            {
                tmin = t;
                vmin = vh;
            }
            if (t > tmax)
            {
                tmax = t;
                vmax = vh;
            }
        }

        std::cout << axis_names[a] << "min: vertex " << vmin.idx() << " (projection " << tmin << ")" << std::endl;
        std::cout << axis_names[a] << "max: vertex " << vmax.idx() << " (projection " << tmax << ")" << std::endl;
        landmarks.push_back(vmin);
        landmarks.push_back(vmax);
    }

    // Warn about vertices that are extreme along several axes: the landmark list then contains duplicates
    std::set<int> unique_ids;
    for (const VH vh : landmarks)
        unique_ids.insert(vh.idx());
    if (unique_ids.size() != landmarks.size())
        std::cerr << "Warning: " << landmarks.size() - unique_ids.size() << " duplicated landmark(s). "
                  << "Duplicates may break the map." << std::endl;

    // Keep the landmarks in all optimization stages ("lk")
    const std::vector<LandmarkType> types(landmarks.size(), LandmarkType::Keep);
    write_landmarks(landmarks, out_path, types, false);

    std::cout << "Wrote " << landmarks.size() << " landmarks to " << out_path << std::endl;
    return 0;
}
