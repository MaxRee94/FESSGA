#pragma once
#include <iostream>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/readOBJ.h>
#include <igl/writeOFF.h>
#include <igl/writeOBJ.h>
#include <igl/writePLY.h>

class IO {
public:
	static void ReadMesh(std::string fpath, Eigen::MatrixXd& V, Eigen::MatrixXi& F, bool suppress_output = false);
};

