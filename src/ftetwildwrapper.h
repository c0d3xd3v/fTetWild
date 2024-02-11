#ifndef FTETWILDWRAPPER_H
#define FTETWILDWRAPPER_H

#include <tbb/task_scheduler_init.h>
#include <thread>

#include <floattetwild/AABBWrapper.h>
#include <floattetwild/FloatTetDelaunay.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/Simplification.h>
#include <floattetwild/Statistics.h>
#include <floattetwild/TriangleInsertion.h>
#include <floattetwild/CSGTreeParser.hpp>
#include <floattetwild/Mesh.hpp>
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/Logger.hpp>

#include <Eigen/Dense>

#include <igl/Timer.h>
#include <igl/write_triangle_mesh.h>
#include <igl/remove_unreferenced.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/mesh/mesh.h>
#include <bitset>

#include <geogram/basic/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/basic/numeric.h>
#include <floattetwild/Predicates.hpp>

#include <floattetwild/MshLoader.h>
#include <geogram/mesh/mesh_AABB.h>

#include <floattetwild/FloatTetDelaunay.h>
#include <floattetwild/MeshImprovement.h>

class FTetWildWrapper
{
private:
    GEO::Mesh* input_mesh;
    floatTetWild::Mesh* mesh;
    floatTetWild::AABBWrapper* tree;
    std::vector<floatTetWild::Vector3>  points;
    std::vector<floatTetWild::Vector3i> faces;
    std::vector<int> input_tags;

    bool skip_simplify;
    bool nobinary;
    bool nocolor;
    bool export_raw;
    double stop_energy;
    double ideal_edge_length_rel;
    double eps_rel;

    static bool is_initialized;

    tbb::task_scheduler_init* scheduler;

    void init();

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    FTetWildWrapper(double stop_energy = 10.0, double ideal_edge_length_rel = 0.05, double eps_rel = 0.001);
    ~FTetWildWrapper();
    void loadMeshGeometry(Eigen::MatrixXf &nodes, Eigen::MatrixXi &tris);
    void tetrahedralize();
    void save();
    void getSurfaceIndices(Eigen::MatrixXi &tris, Eigen::MatrixXi &tets, Eigen::MatrixXf &nodes);
};

#endif // FTETWILDWRAPPER_H
