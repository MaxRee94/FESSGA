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

		static std::string ReplaceOccurrences(string basestring, string toReplace, string replaceWith);

		static double fisqrt(float n);

		static void increment_key(std::map<std::string, int>* map, std::string key);

		static float get_value(std::map<std::string, float>* map, std::string key);

		static int get_value(std::map<std::string, int>* map, std::string key);
	};
};
