#pragma once
#include <Eigen/Core>
#include <vector>
#include <map>

using namespace std;

namespace mvis {
	class help
	{
	public:
		static void NormalizeElementWise(Eigen::Matrix3d* M);

		static vector<size_t> FindAll(string basestring, string target);

		static string GetUniquePath(string fpath);

		static string GetLatestPath(string templ);

		static bool FileExists(string fpath);

		static void RenameFile(string source_path, string target_path, bool verbose = false);

		static std::string ReplaceOccurrences(string basestring, string toReplace, string replaceWith);

		static void WriteToCSV(string csv_path, vector<string> data, string headers);

		static double fisqrt(float n);

		static void increment_key(std::map<std::string, int>* map, std::string key);

		static float get_value(std::map<std::string, float>* map, std::string key);

		static int get_value(std::map<std::string, int>* map, std::string key);
	};
};
