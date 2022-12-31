#include <igl/opengl/glfw/Viewer.h>

using namespace Eigen;
using namespace std;

int main(int argc, char *argv[])
{
  // Inline mesh of a cube
  MatrixXd V = (MatrixXd(8,3)<<
    0.0,0.0,0.0,
    0.0,0.0,1.0,
    0.0,1.0,0.0,
    0.0,1.0,1.0,
    1.0,0.0,0.0,
    1.0,0.0,1.0,
    1.0,1.0,0.0,
    1.0,1.0,1.0).finished();
  MatrixXd V_homog_orig = V.rowwise().homogeneous();
  MatrixXd V_homog = V_homog_orig;
  const MatrixXi F = (MatrixXi(12,3)<<
    0,6,4,
    0,2,6,
    0,3,2,
    0,1,3,
    2,7,6,
    2,3,7,
    4,6,7,
    4,7,5,
    0,4,5,
    0,5,1,
    1,5,7,
    1,7,3).finished();

  // Create transform matrix
  Vector3d pos_offset = Vector3d(0, 0, 0);
  Vector3d transl = Vector3d(0.001f, 0, 0);
  Matrix4d T = (Matrix4d(4, 4) <<
      1.0,          0.0,        0.0,        0.0,
      0.0,          1.0,        0.0,        0.0,
      0.0,          0.0,        1.0,        0.0,
      transl(0),    transl(1),  transl(2),  1.0).finished();
  T = T.transpose();

  // Initialize viewer
  igl::opengl::glfw::Viewer viewer;
  viewer.core().is_animating = true;
  viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&)->bool
  {
      // Update transform matrix
      pos_offset += transl;
      T(3, 0) = pos_offset(0);
      T(3, 1) = pos_offset(1);
      T(3, 2) = pos_offset(2);

      // Transform mesh
      V_homog = V_homog_orig * T;
      V = V_homog.rowwise().hnormalized();

      // Update viewer
      viewer.data().clear();
      viewer.data().set_mesh(V, F);
      viewer.data().set_face_based(true);

      return false;
  };
  viewer.launch();
}
