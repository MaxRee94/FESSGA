#pragma once
#include <iostream>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/readOBJ.h>
#include <igl/writeOFF.h>
#include <igl/writeOBJ.h>
#include <igl/writePLY.h>


namespace mvis {
	class IO {
	public:
		static void ReadMesh(std::string fpath, Eigen::MatrixXd& V, Eigen::MatrixXi& F, bool suppress_output = false);

		static bool FileExists(std::string fpath);

		static void RenameFile(std::string source_path, std::string target_path, bool verbose = false);

		static void WriteToCSV(std::string csv_path, std::vector<std::string> data, std::string headers);

		static std::string GetUniquePath(std::string fpath);

		static std::string GetLatestPath(std::string templ);

		static void write_text_to_file(std::string text, std::string path);
	};
};
