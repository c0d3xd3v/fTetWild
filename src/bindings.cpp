#include "Utils.hpp"

//#include "triangulate_data.hpp"
//#include "triangulate.hpp"
//#include "tetrahedralize.hpp"
#include "ftetwildwrapper.h"

PYBIND11_MODULE(pytetwild, m)
{
    //wildmeshing_binding::tetrahedralize(m);
    py::class_<FTetWildWrapper>(m, "FTetWildWrapper")
            .def(py::init<
                     double, double, double>(),
                     py::arg("stop_energy") = 10,
                     py::arg("ideal_edge_length_rel") = 0.05,
                     py::arg("eps_rel") = 0.001
                )
            .def(
                "loadMeshGeometry", [](FTetWildWrapper &t, Eigen::MatrixXf &V, Eigen::MatrixXi &F)
                {
                    t.loadMeshGeometry(V, F);
                },
                "set mesh",  py::arg("V"), py::arg("F"))
            .def(
                "tetrahedralize", [](FTetWildWrapper &t)
                {
                    t.tetrahedralize();
                },
                "create tet mesh")
            .def(
                "getSurfaceIndices", [](FTetWildWrapper &t)
                {
                    Eigen::MatrixXf nodes;
                    Eigen::MatrixXi tris;
                    Eigen::MatrixXi tets;
                    t.getSurfaceIndices(tris, tets, nodes);
                    return py::make_tuple(tris, tets, nodes);
                },
                "get mesh")
                .def(
                    "save", [](FTetWildWrapper &t, const std::string &path)
                    {
                        t.save(path);
                    },
                    "save mesh", py::arg("path"))
            ;
}
