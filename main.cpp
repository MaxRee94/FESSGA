#include <igl/opengl/glfw/Viewer.h>
#include "io.h"
#include "gui.h"

using namespace Eigen;
using namespace std;
using namespace mvis;


void transform(
    igl::opengl::glfw::Viewer& viewer, MatrixXd& Vhom_orig, MatrixXd& Vhom, Matrix4d& T,
    MatrixXd& V, MatrixXi F, Vector3d pos_offset, float rot_y_offset
) {
    // Update translation component of T matrix
    T(3, 0) = pos_offset(0);
    T(3, 1) = pos_offset(1);
    T(3, 2) = pos_offset(2);

    // Update rotation component of T matrix
    if (rot_y_offset != 0) {
        float cos_ang = cos(rot_y_offset), sin_ang = sin(rot_y_offset);
        T(0, 0) = cos_ang;
        T(2, 2) = cos_ang;
        T(0, 2) = sin_ang;
        T(2, 0) = -sin_ang;
    }

    // Transform mesh
    Vhom = Vhom_orig * T;
    V = Vhom.rowwise().hnormalized();

    // Update viewer
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
    viewer.data().set_face_based(true);
}

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
    //std::cout << "Key: " << key << " " << (unsigned int)key << std::endl;
    return false;
}

int main(int argc, char *argv[])
{
    // Load objects
    int numObjs = 10;
    vector<MatrixXd> V_list;
    vector<MatrixXi> F_list;
    vector<MatrixXd> Vhom_orig_list;
    vector<MatrixXd> Vhom_list;
    MatrixXd V;
    MatrixXi F;
    IO::ReadMesh("E:/Development/MeshVisTempl/assets/test_objects/cube.obj", V, F);
    MatrixXd Vhom_orig = V.rowwise().homogeneous();
    MatrixXd Vhom = Vhom_orig;
    double interobj_offset = 0.5;
    for (int i = 0; i < numObjs; i++) {
        MatrixXd T = Eigen::MatrixXd::Zero(V.rows(), 3);
        MatrixXd Thom = Eigen::MatrixXd::Zero(V.rows(), 4);
        VectorXd offs(V.rows());
        for (int j = 0; j < V.rows(); j++) offs(j) = i * interobj_offset;
        T.col(0) + offs;
        Thom.col(0) + offs;

        V_list.push_back(V + T);
        F_list.push_back(F);
        Vhom_orig_list.push_back(Vhom_orig + Thom);
        Vhom_list.push_back(Vhom + Thom);
    }

    // Create transform matrix
    Vector3d pos_offset = Vector3d(0, 0, 0);
    Vector3d transl = Vector3d(0, 0, 0);
    float rot_y_offset = 0;
    float rot_y_speed = 0.008f;
    Matrix4d T = Matrix4d::Identity(4, 4);
    T = T.transpose();

    // Initialize viewer
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V_list[0], F_list[0]);
    viewer.data().set_face_based(true);
    viewer.callback_key_down = &key_down;
    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&)->bool
    {
        if (viewer.core().is_animating) {
            pos_offset += transl;
            rot_y_offset += rot_y_speed;
            for (int i = 0; i < numObjs; i++) transform(
                viewer, Vhom_orig_list[i], Vhom_list[i], T, V_list[i], F_list[i], pos_offset, rot_y_offset
            );
        }
        return false;
    };
    viewer.launch();

    //GUI(V_list, F_list);
}
