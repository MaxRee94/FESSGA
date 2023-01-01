#include <igl/opengl/glfw/Viewer.h>
#include "io.h"
#include "gui.h"

using namespace Eigen;
using namespace std;
using namespace mvis;


void transform(
    igl::opengl::glfw::Viewer& viewer, MatrixXd& Vhom_orig, MatrixXd& Vhom, Matrix4d& T,
    MatrixXd& V, MatrixXi F, Vector3d& pos_offset, Vector3d& transl, float& rot_y_offset, float rot_y_speed
) {
    // Update translation component of T matrix
    pos_offset += transl;
    T(3, 0) = pos_offset(0);
    T(3, 1) = pos_offset(1);
    T(3, 2) = pos_offset(2);

    // Update rotation component of T matrix
    if (rot_y_speed != 0) {
        rot_y_offset += rot_y_speed;
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
    vector<MatrixXd> V_list;
    vector<MatrixXi> F_list;
    MatrixXd V;
    MatrixXi F;
    IO::ReadMesh("E:/Development/MeshVisTempl/assets/test_objects/cube.obj", V, F);
    for (int i = 0; i < 1; i++) {
        V_list.push_back(V);
        F_list.push_back(F);
    }
    MatrixXd Vhom_orig = V_list[0].rowwise().homogeneous();
    MatrixXd Vhom = Vhom_orig;

    // Create transform matrix
    Vector3d pos_offset = Vector3d(0, 0, 0);
    Vector3d transl = Vector3d(0, 0, 0);
    float rot_y_offset = 0;
    float rot_y_speed = 0.008f;
    Matrix4d T = (Matrix4d(4, 4) <<
        1,          0,        0,        0,
        0,          1,        0,        0,
        0,          0,        1,        0,
        0,          0,        0,        1).finished();
    T = T.transpose();

    // Initialize viewer
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V_list[0], F_list[0]);
    viewer.data().set_face_based(true);

    viewer.callback_key_down = &key_down;
    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&)->bool
    {
        if (viewer.core().is_animating) transform(
            viewer, Vhom_orig, Vhom, T, V_list[0], F_list[0], pos_offset, transl, rot_y_offset, rot_y_speed
        );
        return false;
    };
    viewer.launch();
}
