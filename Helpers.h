#pragma once
#include <Eigen/Core>
#include <vector>
#include <map>
#include <cstdlib>
#include <set>

#define VERBOSE false

using namespace std;

typedef unsigned int uint;

// Comparison function for sorting the
// set by increasing order of its pair's
// second value
struct _comparator {
	template <typename T>

	// Comparator function
	bool operator()(const T& l, const T& r) const
	{
		if (l.second != r.second) {
			return l.second < r.second;
		}
		return l.first < r.first;
	}
};

typedef std::set < std::pair<int, double>, _comparator> PairSet;

namespace fessga {

	class help
	{
	public:

		static double max(double val1, double val2);

		static void sort(std::map<int, double>& _map, PairSet& _set);

		static int get_key(std::map<int, int>* _map, int value);

		static void populate_with_zeroes(double* _array, int dim_x, int dim_y);
		
		static void populate_with_zeroes(uint* _array, int dim_x, int dim_y);

		static void split(string basestring, string separator, vector<string>& substrings);

		/*
		Return whether the given vector <vec> contains the integer <item>
		 */
		static bool is_in(std::vector<int>* vec, int item);

		static void print_map(std::map<int, int>* map);

		static void print_vector(std::vector<int>* vec);

		static void print_pairs(std::vector<pair<int, int>>* vec);

		static void print(std::string);

		static void NormalizeElementWise(Eigen::Matrix3d* M);

		static vector<size_t> FindAll(string basestring, string target);

		static bool is_in(string basestring, string target);

		static std::string replace_occurrences(string basestring, string toReplace, string replaceWith);

		static double fisqrt(float n);

		static void increment_key(std::map<std::string, int>* map, std::string key);

		static float get_value(std::map<std::string, float>* map, std::string key);

		static int get_value(std::map<std::string, int>* map, std::string key);

		static int get_value(std::map<int, int>* map, int key);

		static int get_value(std::map<uint32_t, uint32_t>* map, uint32_t key);

		static float get_rand_float(float min, float max);

		static uint get_rand_uint(float min, float max);

		static std::string add_padding(std::string basestring, int version);

		static string join_as_string(vector<int> numbers, string separator);

		static string join_as_string(vector<pair<int, int>> numbers, string separator);

	};
};
