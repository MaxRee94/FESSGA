#pragma once
#include <iostream>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/readOBJ.h>
#include <igl/writeOFF.h>
#include <igl/writeOBJ.h>
#include <igl/writePLY.h>


namespace fessga {
	class IO {
	public:
		static void read_mesh(std::string fpath, Eigen::MatrixXd& V, Eigen::MatrixXi& F, bool suppress_output = false);

		static void read_file_content(std::string fpath, std::vector<std::string>& content);

		static bool file_exists(std::string fpath);

		static void copy_file(std::string _source, std::string _target, bool verbose = false);

		static void rename_file(std::string source_path, std::string target_path, bool verbose = false);

		static void write_to_csv(std::string csv_path, std::vector<std::string> data, std::string headers);

		static std::string get_unique_file_path(std::string fpath);

		static std::string get_unique_path(std::string templ);
		
		static std::string get_latest_path(std::string templ);

		static void write_text_to_file(std::string text, std::string path);

		static std::string create_folder_if_not_exists(std::string folder_path);

		static bool is_empty(std::string folder);

		static std::string get_fullpath(std::string relative_path);

		static void append_to_file(std::string fpath, std::string text);

		static void get_files_in_directory(std::vector<std::string>& filepaths, std::string source_dir);

		static void remove_files(std::vector<std::string>* files);

		static void remove_directory_incl_contents(std::string dir);
	};
};
