//#include <igl/opengl/glfw/Viewer.h>
#include "gui.h"
#include "io.h"

using namespace std;

fessga::GUI::GUI(vector<MatrixXd> _V_list, vector<MatrixXi> _F_list) {
    V_list = _V_list;
    F_list = _F_list;
    //viewer.data().set_face_based(true);
}

// This function is called every time a keyboard button is pressed
//bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
//{
//    //std::cout << "Key: " << key << " " << (unsigned int)key << std::endl;
//    return false;
//}

void fessga::GUI::load_example(MatrixXd* V, MatrixXi* F) {
    V_list.clear();
    F_list.clear();
    V_list.push_back(*V);
    F_list.push_back(*F);
}

//void fessga::GUI::transform(
//    igl::opengl::glfw::Viewer& viewer, MatrixXd& Vhom_orig, MatrixXd& Vhom, Matrix4d& T,
//    MatrixXd& V, MatrixXi F, Vector3d pos_offset, float rot_y_offset
//) {
//    // Update translation component of T matrix
//    T(3, 0) = pos_offset(0);
//    T(3, 1) = pos_offset(1);
//    T(3, 2) = pos_offset(2);
//
//    // Update rotation component of T matrix
//    if (rot_y_offset != 0) {
//        float cos_ang = cos(rot_y_offset), sin_ang = sin(rot_y_offset);
//        T(0, 0) = cos_ang;
//        T(2, 2) = cos_ang;
//        T(0, 2) = sin_ang;
//        T(2, 0) = -sin_ang;
//    }
//
//    // Transform mesh
//    Vhom = Vhom_orig * T;
//    V = Vhom.rowwise().hnormalized();
//
//    // Update viewer
//    viewer.data().clear();
//    viewer.data().set_mesh(V, F);
//    viewer.data().set_face_based(true);
//}

//void fessga::GUI::show() {
//    // Update viewer with current mesh lists
//    for (int i = 0; i < V_list.size(); i++)
//        viewer.data().set_mesh(V_list.at(0), F_list.at(0));
//    viewer.launch();
//}
