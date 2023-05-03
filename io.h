#pragma once
#include <iostream>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/readOBJ.h>
#include <igl/writeOFF.h>
#include <igl/writeOBJ.h>
#include <igl/writePLY.h>

using namespace std;

namespace mvis {
	class IO {
	public:
		static void ReadMesh(string fpath, Eigen::MatrixXd& V, Eigen::MatrixXi& F, bool suppress_output = false);

		static bool FileExists(string fpath);

		static void RenameFile(string source_path, string target_path, bool verbose = false);

		static void WriteToCSV(string csv_path, vector<string> data, string headers);

		static string GetUniquePath(string fpath);

		static string GetLatestPath(string templ);
	};
};
