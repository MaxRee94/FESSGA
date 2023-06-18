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

float fessga::help::get_rand_float(float min, float max) {
    return min + (float)rand() * INV_RAND_MAX * (max - min);
}

uint fessga::help::get_rand_uint(float min, float max) {
    float float_rand_range = (float)rand() * INV_RAND_MAX * (max - min);
    return (int)(min + float_rand_range);
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

string fessga::help::join_as_string(vector<pair<int, int>> numbers, string separator) {
    string result = "";
    for (auto pair : numbers) {
        result += "(" + to_string(pair.first) + ", " + to_string(pair.second) + ")" + separator;
    }
    return result;
}

void fessga::help::print_vector(std::vector<int>* vec) {
    for (int i = 0; i < vec->size(); i++) {
        if (i > 0) cout << ", ";
        cout << vec->at(i);
    }
    cout << endl;
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
    int pos = 0;
    string newstring = basestring;
    while (basestring.find(toReplace, pos) != string::npos) {
        pos = basestring.find(toReplace);
        newstring = basestring.replace(pos, toReplace.size(), replaceWith);
        pos++;
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
    // Declare set of pairs and insert
    // pairs according to the comparator
    // function comp()
    auto beginning = _map.end();
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
