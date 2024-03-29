#include "helpers.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <filesystem.>
#include <time.h>


using namespace Eigen;
using namespace std;

float INV_RAND_MAX = 1.0 / (float)RAND_MAX;

void fessga::help::init_RNG() {
    srand(time(NULL));
    int x = rand();
    int z = x * 3;
}

float fessga::help::get_rand_float(float min, float max) {
    return min + (float)rand() * INV_RAND_MAX * (max - min);
}

uint fessga::help::get_rand_uint(float min, float max) {
    float float_rand_range = (float)rand() * INV_RAND_MAX * (max - min);
    return round(min + float_rand_range);
}

bool fessga::help::is_in(std::vector<int>* vec, int item) {
    return find(vec->begin(), vec->end(), item) != vec->end();
}

// Add padding as suffix to given basestring
std::string fessga::help::add_padding(std::string basestring, int version) {
    int padding = 4;
    if (version > 9) {
        if (version > 99) {
            if (version > 999) {
                if (version > 9999) {
                    padding = 0;
                }
                else padding = 1;
            }
            else padding = 2;
        }
        else padding = 3;
    }
    string pad = "";
    for (int i = 0; i < padding; i++) {
        pad += "0";
    }
    return basestring + pad;
}

void fessga::help::print_map(std::map<int, int>* map) {
    int i = 0;
    for (auto const& [key, val] : (*map))
    {
        if (i > 0) {
            std::cout << ", ";
        }
        std::cout << key        // string (key)
            << ':'
            << val;        // string's value
        i++;
    }
    if (i == 0) {
        cout << "<empty map>" << endl;
    }
    else {
        cout << endl;
    }
}

string fessga::help::join_as_string(vector<int> numbers, string separator) {
    string result = "";
    for (auto number : numbers) {
        result += to_string(number) + separator;
    }
    return result;
}

string fessga::help::join_as_string(vector<float> numbers, string separator) {
    string result = "";
    for (auto number : numbers) {
        result += to_string(number) + separator;
    }
    return result;
}

string fessga::help::join_as_string(vector<pair<int, int>> numbers, string separator) {
    string result = "";
    for (auto pair : numbers) {
        result += "(" + to_string(pair.first) + ", " + to_string(pair.second) + ")" + separator;
    }
    return result;
}

string fessga::help::join(vector<string>* strings, string separator) {
    string result = "";
    for (int i = 0; i < strings->size(); i++) {
        result += strings->at(i);
        if (i < strings->size() - 1) result += separator;
    }
    return result;
}

string fessga::help::join(list<string>* _strings, string separator) {
    vector<string> strings;
    for (auto const& string : *_strings) strings.push_back(string);
    return join(&strings, separator);
}

void fessga::help::print_vector(std::vector<int>* vec) {
    for (int i = 0; i < vec->size(); i++) {
        if (i > 0) cout << ", ";
        cout << vec->at(i);
    }
    cout << endl;
}

void fessga::help::print(std::string str) {
#if VERBOSE
    cout << str;
#endif
}

void fessga::help::print_pairs(std::vector<pair<int, int>>* pairs) {
    for (int i = 0; i < pairs->size(); i++) {
        if (i > 0) cout << " ";
        cout << "(";
        cout << to_string(pairs->at(i).first / 6) + ", " + to_string(pairs->at(i).first % 6);
        //cout << to_string(pairs->at(i).first) + ", " + to_string(pairs->at(i).second);
        cout << ")";
    }
    cout << endl;
}

void fessga::help::NormalizeElementWise(Matrix3d* M) {
    for (int i = 0; i < M->cols(); i++) {
        for (int j = 0; j < M->rows(); j++) {
            if ((*M)(i, j) != 0) {
                (*M)(i, j) /= abs((*M)(i, j));
            }
        }
    };
}

string fessga::help::replace_occurrences(string basestring, string toReplace, string replaceWith) {
    string newstring = basestring;
    int pos = newstring.find(toReplace);
    while (pos != string::npos) {
        newstring.replace(pos, toReplace.size(), replaceWith);
        pos = newstring.find(toReplace);
    }

    return newstring;
}


vector<size_t> fessga::help::FindAll(string basestring, string target) {
    vector<size_t> occurrences;
    size_t found = 0;
    while (true) {
        found = basestring.find(target, found);
        if (found != string::npos) {
            occurrences.push_back(found);
        }
        else {
            break;
        }
        found++;
    }

    return occurrences;
}

bool fessga::help::is_in(string basestring, string target) {
    size_t found = 0;
    found = basestring.find(target, found);
    if (found == string::npos) {
        return false;
    }
    else return true;
}


double fessga::help::fisqrt(float n)
{
    float y = n;
    long i = *(long*)&y;
    i = 0x5f3759df - (i >> 1);
    y = *(float*)&i;

    return y * (1.5f - ((n * 0.5f) * y * y));
}


void fessga::help::increment_key(std::map<std::string, int>* _map, std::string key) {
    if (_map->find(key) == _map->end()) {
        (*_map)[key] = 1;
    }
    else {
        (*_map)[key]++;
    }
}

float fessga::help::get_value(std::map<std::string, float>* _map, std::string key) {
    if (_map == 0) return 0;
    if (_map->find(key) == _map->end()) {
        return 0;
    }
    else {
        return _map->at(key);
    }
}

int fessga::help::get_value(std::map<std::string, int>* _map, std::string key) {
    if (_map == 0) return 0;
    if (_map->find(key) == _map->end()) {
        return 0;
    }
    else {
        return _map->at(key);
    }
}

int fessga::help::get_value(std::map<int, int>* _map, int key) {
    if (_map == 0) return -1;
    map<int, int>::iterator it = _map->find(key);
    if (it == _map->end()) {
        return -1;
    }
    else {
        return it->second;
    }
}

int fessga::help::get_value(std::map<int, double>* _map, int key) {
    if (_map == 0) return -1;
    map<int, double>::iterator it = _map->find(key);
    if (it == _map->end()) {
        return -1;
    }
    else {
        return it->second;
    }
}

int fessga::help::get_value(std::map<uint32_t, uint32_t>* _map, uint32_t key) {
    if (_map == 0) return -1;
    std::map<uint32_t, uint32_t>::iterator it = _map->find(key);
    if (it == _map->end()) {
        return -1;
    }
    else {
        return it->second;
    }
}

int fessga::help::get_key(std::map<int, int>* _map, int value) {
    if (_map == 0) return -1;
    for (auto const& [key, val] : (*_map))
    {
        if (val == value) {
            return key;
        }
    }
    return -1;
}

void fessga::help::populate_with_zeroes(double* _array, int dim_x, int dim_y) {
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            _array[x * dim_y + y] = 0.0;
        }
    }
}

void fessga::help::populate_with_zeroes(uint* _array, int dim_x, int dim_y) {
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            _array[x * dim_y + y] = 0.0;
        }
    }
}

// Function to sort the map according
// to value in a (key-value) pairs
void fessga::help::sort(std::map<int, double>& _map, PairSet& _set)
{
    _set = PairSet(_map.begin(), _map.end());
}

double fessga::help::max(double val1, double val2) {
    if (val1 > val2) return val1;
    else return val2;
}

void fessga::help::split(string basestring, string separator, vector<string>& substrings) {
    vector<size_t> occurrences = FindAll(basestring, separator);
    if (occurrences.size() == 0) {
        substrings = { basestring };
        return;
    }
    substrings.push_back(basestring.substr(0, occurrences[0]));
    for (int i = 0; i < occurrences.size() - 1; i++) {
        int name_length = occurrences[i + 1] - (occurrences[i] + separator.size());
        string name = basestring.substr(occurrences[i] + separator.size(), name_length);
        substrings.push_back(name);
    }
    string end = basestring.substr(occurrences[occurrences.size() - 1] + separator.size(), string::npos);
    if (end != "") {
        substrings.push_back(end);
    }
}

// Remove the largest item from the given vector
void fessga::help::remove_largest_vector(vector<vector<int>>* vectors, int& max_size) {
    max_size = 0;
    int largest_item_idx = -1;
    for (int i = 0; i < vectors->size(); i++) {
        if (vectors->at(i).size() > max_size) {
            max_size = vectors->at(i).size();
            largest_item_idx = i;
        }
    }
    vectors->erase(vectors->begin() + largest_item_idx);
}

void fessga::help::remove(vector<int>* vec, int item) {
    auto position = find(vec->begin(), vec->end(), item);
    if (position != vec->end()) {
        int idx = position - vec->begin();
        vec->erase(vec->begin() + idx);
    }
    else throw("Error: Cannot remove item from vector because it is not present\n");
}

bool fessga::help::ends_with(string full_string, string ending) {
    if (full_string.length() >= ending.length()) {
        return (0 == full_string.compare(full_string.length() - ending.length(), ending.length(), ending));
    }
    else {
        return false;
    }
}

bool fessga::help::have_overlap(vector<int>* larger_vector, vector<int>* smaller_vector) {
    for (auto& item : *smaller_vector) {
        if (help::is_in(larger_vector, item)) return true;
    }
    return false;
}

void fessga::help::append_vector(vector<int>& result, vector<int>* vec2) {
    for (auto& item : *vec2) result.push_back(item);
}

void fessga::help::append_vector(vector<int>& result, vector<int> vec2) {
    append_vector(result, &vec2);
}

void fessga::help::append_vector(vector<pair<int, int>>& result, vector<pair<int, int>>* vec2) {
    for (auto& item : *vec2) result.push_back(item);
}

void fessga::help::append_vector(vector<pair<int, int>>& result, vector<pair<int, int>> vec2) {
    append_vector(result, &vec2);
}

void fessga::help::append_vector(vector<string>& result, vector<string>* vec2) {
    for (auto& item : *vec2) result.push_back(item);
}

vector<float> fessga::help::get_free_memory() {
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    unsigned long long RAM_bytes = status.ullAvailPhys;
    unsigned long long VM_bytes = status.ullAvailVirtual;
    unsigned long long Pagefile_bytes = status.ullAvailPageFile;
    auto percent_memory = status.dwMemoryLoad;

    float RAM_gigabytes = (float)RAM_bytes / (float)(1 << 30);
    float VM_gigabytes = (float)VM_bytes / (float)(1 << 30);
    float Pagefile_gigabytes = (float)Pagefile_bytes / (float)(1 << 30);
    return { RAM_gigabytes, VM_gigabytes, Pagefile_gigabytes, (float)percent_memory };
}

double fessga::help::get_mean(vector<double>* distribution) {
    double sum = 0;
    for (double& sample : *distribution) {
        sum += sample;
    }
    return sum / distribution->size();
}

double fessga::help::get_stdev(vector<double>* distribution, double mean) {
    if (mean == -999999) {
        mean = get_mean(distribution);
    }
    double sos = 0;
    for (auto& sample : *distribution) sos += (sample - mean) * (sample - mean);
    double variance = sos / (distribution->size() - 1);
    return sqrt(variance);
}


double fessga::help::get_max(vector<double>* distribution) {
    double max = -INFINITY;
    for (double value : *distribution) {
        if (value > max) max = value;
    }
    return max;
}


double fessga::help::get_min(vector<double>* distribution) {
    double min = INFINITY;
    for (double value : *distribution) {
        if (value < min) min = value;
    }
    return min;
}

