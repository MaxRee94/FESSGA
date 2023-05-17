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

void fessga::help::print_vector(std::vector<int>* vec) {
    for (int i = 0; i < vec->size(); i++) {
        if (i > 0) cout << ", ";
        cout << vec->at(i);
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

string fessga::help::ReplaceOccurrences(string basestring, string toReplace, string replaceWith) {
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


double fessga::help::fisqrt(float n)
{
    float y = n;
    long i = *(long*)&y;
    i = 0x5f3759df - (i >> 1);
    y = *(float*)&i;

    return y * (1.5f - ((n * 0.5f) * y * y));
}


void fessga::help::increment_key(std::map<std::string, int>* map, std::string key) {
    if (map->find(key) == map->end()) {
        (*map)[key] = 1;
    }
    else {
        (*map)[key]++;
    }
}

float fessga::help::get_value(std::map<std::string, float>* map, std::string key) {
    if (map->find(key) == map->end()) {
        return 0;
    }
    else {
        return map->at(key);
    }
}

int fessga::help::get_value(std::map<std::string, int>* map, std::string key) {
    if (map->find(key) == map->end()) {
        return 0;
    }
    else {
        return map->at(key);
    }
}

int fessga::help::get_value(std::map<int, int>* _map, int key) {
    map<int, int>::iterator it = _map->find(key);
    if (it == _map->end()) {
        return -1;
    }
    else {
        return it->second;
    }
}

int fessga::help::get_value(std::map<uint32_t, uint32_t>* _map, uint32_t key) {
    std::map<uint32_t, uint32_t>::iterator it = _map->find(key);
    if (it == _map->end()) {
        return -1;
    }
    else {
        return it->second;
    }
}

int fessga::help::get_key(std::map<int, int>* _map, int value) {
    for (auto const& [key, val] : (*_map))
    {
        if (val == value) {
            return key;
        }
    }
    return -1;
}
