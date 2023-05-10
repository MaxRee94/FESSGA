#include "Helpers.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <filesystem.>

using namespace Eigen;
using namespace std;



void mvis::help::print_map(std::map<int, int>* map) {
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

void mvis::help::NormalizeElementWise(Matrix3d* M) {
    for (int i = 0; i < M->cols(); i++) {
        for (int j = 0; j < M->rows(); j++) {
            if ((*M)(i, j) != 0) {
                (*M)(i, j) /= abs((*M)(i, j));
            }
        }
    };
}

string mvis::help::ReplaceOccurrences(string basestring, string toReplace, string replaceWith) {
    int pos = 0;
    string newstring = basestring;
    while (basestring.find(toReplace, pos) != string::npos) {
        pos = basestring.find(toReplace);
        newstring = basestring.replace(pos, toReplace.size(), replaceWith);
        pos++;
    }

    return newstring;
}


vector<size_t> mvis::help::FindAll(string basestring, string target) {
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


double mvis::help::fisqrt(float n)
{
    float y = n;
    long i = *(long*)&y;
    i = 0x5f3759df - (i >> 1);
    y = *(float*)&i;

    return y * (1.5f - ((n * 0.5f) * y * y));
}


void mvis::help::increment_key(std::map<std::string, int>* map, std::string key) {
    if (map->find(key) == map->end()) {
        (*map)[key] = 1;
    }
    else {
        (*map)[key]++;
    }
}

float mvis::help::get_value(std::map<std::string, float>* map, std::string key) {
    if (map->find(key) == map->end()) {
        return 0;
    }
    else {
        return map->at(key);
    }
}

int mvis::help::get_value(std::map<std::string, int>* map, std::string key) {
    if (map->find(key) == map->end()) {
        return 0;
    }
    else {
        return map->at(key);
    }
}
