#include <QApplication>
#include <QPushButton>

#include <igl/read_triangle_mesh.h>

#include <floattetwild/ftetwildwrapper.h>
#include <floattetwild/tetrahedralize.hpp>

class Test {
private:
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    void test()
    {
        const std::string oldLocale = std::setlocale(LC_NUMERIC, nullptr);
        std::setlocale(LC_NUMERIC, "C");

        std::string path = "../../../data/ridex-Body004.obj";

        /*
        std::vector<double> epsr_tags;
        wildmeshing_binding::Tetrahedralizer* tetra = new wildmeshing_binding::Tetrahedralizer();
        tetra->load_mesh(path, "", epsr_tags);
        tetra->tetrahedralize();
        tetra->save("ridex-Body004.msh", false, true, false, true, false, false, false);
        delete tetra;
        */

        Eigen::MatrixXf V;
        Eigen::MatrixXi F;
        igl::read_triangle_mesh(path, V, F);

        FTetWildWrapper* ftetwildWrapper = new FTetWildWrapper();
        ftetwildWrapper->loadMeshGeometry(V, F);
        ftetwildWrapper->tetrahedralize();
        ftetwildWrapper->save();

        Eigen::MatrixXi tris;
        Eigen::MatrixXi tets;
        Eigen::MatrixXf nodes;
        ftetwildWrapper->getSurfaceIndices(tris, tets, nodes);

        std::setlocale(LC_NUMERIC, oldLocale.c_str());
    }
};

int main(int argc, char** argv)
{
    std::cout << "EIGEN_WORLD_VERSION : " << EIGEN_WORLD_VERSION << std::endl;
    std::cout << "EIGEN_MAJOR_VERSION : " << EIGEN_MAJOR_VERSION << std::endl;
    std::cout << "EIGEN_WORLD_VERSION : " << EIGEN_MINOR_VERSION << std::endl;

    QApplication app(argc, argv);

    QPushButton btn;
    btn.show();
    btn.connect(&btn, &QPushButton::clicked, [=](){Test test; test.test();});
    return app.exec();
}
