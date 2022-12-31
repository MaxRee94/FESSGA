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

  igl::opengl::glfw::Viewer viewer;
  
  Vector3d transform = Vector3d(0.001f, 0, 0);
  Matrix4d T = (Matrix4d(4, 4) <<
      1.0, 0.0, 0.0, 0.0,
      0.0, 1.0, 0.0, 0.0,
      0.0, 0.0, 1.0, 0.0,
      transform(0), transform(1), transform(2), 1.0).finished();
  T = T.transpose();

  // Set initial cube mesh
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);

  viewer.core().is_animating = true;

  MatrixXd V_homog_prev;
  Vector3d d_transform = transform;
  viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&)->bool
  {
      // Update transform matrix
      transform += d_transform;
      T(3, 0) = transform(0);
      T(3, 1) = transform(1);
      T(3, 2) = transform(2);

      viewer.data().clear();
      viewer.data().set_mesh(V, F);
      viewer.data().set_face_based(true);
      V_homog_prev = V_homog;
      V_homog = V_homog_orig * T;
      V = V_homog.rowwise().hnormalized();
      cout << "V homologous: \n" << V_homog << "\n" << endl;
      cout << "V: \n" << V << "\n" << endl;
      cout << "Change: \n" << (V_homog - V_homog_prev) << "\n" << endl;

      return false;
  };
  viewer.launch();
}
