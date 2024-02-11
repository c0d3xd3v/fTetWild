#include "ftetwildwrapper.h"

#include <memory>

bool FTetWildWrapper::is_initialized = false;


FTetWildWrapper::FTetWildWrapper(double stop_energy, double ideal_edge_length_rel, double eps_rel) :
    skip_simplify(false),
    nobinary(false),
    nocolor(false),
    export_raw(false),
    mesh(nullptr),
    scheduler(nullptr),
    input_mesh(nullptr),
    tree(nullptr),
    stop_energy(stop_energy),
    ideal_edge_length_rel(ideal_edge_length_rel),
    eps_rel(eps_rel)
{
    this->init();
}

FTetWildWrapper::~FTetWildWrapper()
{
    if(mesh != nullptr) delete mesh;
    if(scheduler != nullptr) delete scheduler;
    if(input_mesh != nullptr) delete input_mesh;
    if(tree != nullptr) delete tree;
}

void FTetWildWrapper::init()
{
    if(!is_initialized)
    {
        GEO::initialize();

        std::vector<int> indices(20);
        std::iota(std::begin(indices), std::end(indices), 0);
        floatTetWild::Random::shuffle(indices);
        for (int a : indices)
            std::cout << a << " ";
        std::cout << std::endl;

        // Import standard command line arguments, and custom ones
        GEO::CmdLine::import_arg_group("standard");
        GEO::CmdLine::import_arg_group("pre");
        GEO::CmdLine::import_arg_group("algo");

        is_initialized = true;
    }

//TODO: remove hack
#define FLOAT_TETWILD_USE_TBB

#ifndef WIN32
        setenv("GEO_NO_SIGNAL_HANDLER", "1", 1);
#endif

#ifdef FLOAT_TETWILD_USE_TBB
        const size_t MB          = 1024 * 1024;
        const size_t stack_size  = 64 * MB;
        unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
        unsigned int max_threads = std::numeric_limits<unsigned int>::max();
        num_threads              = std::min(max_threads, num_threads);
        //params.num_threads       = num_threads;
        std::cout << "TBB threads " << num_threads << std::endl;
        scheduler = new tbb::task_scheduler_init(num_threads, stack_size);
#endif
}

void FTetWildWrapper::loadMeshGeometry(Eigen::MatrixXf &nodes, Eigen::MatrixXi &tris)
{
    if(mesh != nullptr) delete mesh;
    mesh = new floatTetWild::Mesh();
    mesh->params.stop_energy = stop_energy;
    mesh->params.ideal_edge_length_rel = ideal_edge_length_rel;
    mesh->params.eps_rel = eps_rel;

    for(int i = 0; i < nodes.rows(); i++)
    {
        floatTetWild::Vector3 p = nodes.row(i).cast<floatTetWild::Scalar>();
        points.push_back(p);
    }

    for(int i = 0; i < tris.rows(); i++)
    {
        floatTetWild::Vector3i f = tris.row(i);
        faces.push_back(f);
    }

    if(input_mesh != nullptr) delete input_mesh;
    input_mesh = new GEO::Mesh();
    std::vector<int> flags;
    floatTetWild::MeshIO::load_mesh(points, faces, *input_mesh, flags);

    if(tree != nullptr) delete tree;
    tree = new floatTetWild::AABBWrapper(*input_mesh);
    floatTetWild::Parameters& params = mesh->params;

    if (!params.init(tree->get_sf_diag()))
    {
        std::cout << "initialization failed ..." << std::endl;
    }

    input_tags.resize(faces.size());
    std::fill(input_tags.begin(), input_tags.end(), 0);

    floatTetWild::simplify(points, faces, input_tags, *tree, params, skip_simplify);
    tree->init_b_mesh_and_tree(points, faces, *mesh);
}

void FTetWildWrapper::tetrahedralize()
{
    std::vector<bool> is_face_inserted(faces.size(), false);
    floatTetWild::FloatTetDelaunay::tetrahedralize(points, faces, *tree, *mesh, is_face_inserted);
    floatTetWild::insert_triangles(points, faces, input_tags, *mesh, is_face_inserted, *tree, false);
    floatTetWild::optimization(points, faces, input_tags, is_face_inserted, *mesh, *tree, {{1, 1, 1, 1}});
    floatTetWild::correct_tracked_surface_orientation(*mesh, *tree);

    floatTetWild::Parameters& params = mesh->params;
    if (params.smooth_open_boundary) {
        floatTetWild::smooth_open_boundary(*mesh, *tree);
        for (floatTetWild::MeshTet& t : mesh->tets) {
            if (t.is_outside)
                t.is_removed = true;
        }
    }
    else {
        if (!params.disable_filtering) {
            if (params.use_floodfill) {
                floatTetWild::filter_outside_floodfill(*mesh);
            }
            else if (params.use_input_for_wn) {
                floatTetWild::filter_outside(*mesh, points, faces);
            }
            else
                floatTetWild::filter_outside(*mesh);
        }
    }
}

void FTetWildWrapper::save()
{
    Eigen::MatrixXd V_sf;
    Eigen::MatrixXi F_sf;
    floatTetWild::Parameters& params = mesh->params;
    if (params.manifold_surface) {
        floatTetWild::manifold_surface(*mesh, V_sf, F_sf);
    }
    else {
        floatTetWild::get_surface(*mesh, V_sf, F_sf);
    }
    std::cout << "write : " << (params.output_path + "_" + params.postfix + ".obj") << std::endl;
    igl::write_triangle_mesh(params.output_path + "_" + params.postfix + ".obj", V_sf, F_sf);
}

void FTetWildWrapper::getSurfaceIndices(Eigen::MatrixXi &tris, Eigen::MatrixXi &tets, Eigen::MatrixXf &nodes)
{
    Eigen::MatrixXd flags;
    Eigen::MatrixXd V;
    Eigen::MatrixXi T;

    Eigen::MatrixXd Vout;
    Eigen::MatrixXi Fout;

    const auto skip_tet = [this](const int i)
    { return mesh->tets[i].is_removed; };
    const auto skip_vertex = [this](const int i)
    { return mesh->tet_vertices[i].is_removed; };
    std::vector<int> t_ids(mesh->tets.size());
    std::iota(std::begin(t_ids), std::end(t_ids), 0);

    int cnt_v = 0;
    std::map<int, int> old_2_new;
    for (int i = 0; i < mesh->tet_vertices.size(); i++)
    {
        if (!skip_vertex(i))
        {
            old_2_new[i] = cnt_v;
            cnt_v++;
        }
    }
    int cnt_t = 0;
    for (const int i : t_ids)
    {
        if (!skip_tet(i))
            cnt_t++;
    }

    V.resize(cnt_v, 3);
    int index = 0;
    for (size_t i = 0; i < mesh->tet_vertices.size(); i++)
    {
        if (skip_vertex(i))
            continue;
        V.row(index++) << mesh->tet_vertices[i][0], mesh->tet_vertices[i][1], mesh->tet_vertices[i][2];
    }

    T.resize(cnt_t, 4);
    flags.resize(cnt_t, 1);
    index = 0;

    const std::array<int, 4> new_indices = {{0, 1, 3, 2}};

    for (const int i : t_ids)
    {
        if (skip_tet(i))
            continue;
        for (int j = 0; j < 4; j++)
        {
            T(index, j) = old_2_new[mesh->tets[i][new_indices[j]]];
        }
        flags(index) = mesh->tets[i].scalar;
        index++;
    }

    Eigen::MatrixXi I;
    Eigen::MatrixXd Vs;
    Eigen::MatrixXi Ts;
    igl::remove_unreferenced(V, T, Vs, Ts, I);

    std::cout << V.rows() << ", " << V.cols() << std::endl;
    std::cout << Vs.rows() << ", " << Vs.cols() << std::endl;


    std::cout << T.rows() << ", " << T.cols() << std::endl;
    std::cout << Ts.rows() << ", " << Ts.cols() << std::endl;

    floatTetWild::Mesh _mesh;
    for(int i = 0; i < V.rows(); i++)
        _mesh.tet_vertices.push_back(floatTetWild::MeshVertex(V.row(i)));
    for(int i = 0; i < T.rows(); i++)
        _mesh.tets.push_back(floatTetWild::MeshTet(T.row(i)));

    nodes.resize(_mesh.tet_vertices.size(), 3);
    for(unsigned int i = 0; i < nodes.rows(); i++)
        nodes.row(i) = _mesh.tet_vertices[i].pos.cast<float>();

    tets.resize(_mesh.tets.size(), 4);
    for(unsigned int i = 0; i < tets.rows(); i++)
        tets.row(i) = _mesh.tets[i].indices;

    floatTetWild::get_boundary_surface_indices(_mesh, tris);
}
