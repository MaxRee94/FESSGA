#pragma once
#include <Eigen/Core>
#include <vector>
#include <map>

using namespace std;

namespace mvis {
	class help
	{
	public:
		static int get_key(std::map<int, int>* _map, int value);

		/*
		Return whether the given vector <vec> contains the integer <item>
		 */
		static bool is_in(std::vector<int>* vec, int item);

		static void print_map(std::map<int, int>* map);

		static void print_vector(std::vector<int>* vec);

		static void NormalizeElementWise(Eigen::Matrix3d* M);

		static vector<size_t> FindAll(string basestring, string target);

		static std::string ReplaceOccurrences(string basestring, string toReplace, string replaceWith);

		static double fisqrt(float n);

		static void increment_key(std::map<std::string, int>* map, std::string key);

		static float get_value(std::map<std::string, float>* map, std::string key);

		static int get_value(std::map<std::string, int>* map, std::string key);

		static int get_value(std::map<int, int>* map, int key);

	};
};
