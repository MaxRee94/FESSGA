#pragma once
#include <iostream>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>

using namespace Eigen;
using namespace std;

namespace fessga {
	class GUI {
	public:
		GUI() = default;
		GUI(vector<MatrixXd> _V_list, vector<MatrixXi> _F_list);
		void transform(igl::opengl::glfw::Viewer& viewer, MatrixXd& Vhom_orig, MatrixXd& Vhom, Matrix4d& T,
			MatrixXd& V, MatrixXi F, Vector3d pos_offset, float rot_y_offset);
		void load_example();
		void show();
		igl::opengl::glfw::Viewer viewer;
		vector<MatrixXd> V_list;
		vector<MatrixXi> F_list;
	};

};
