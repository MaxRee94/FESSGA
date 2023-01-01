#include "Helpers.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <filesystem.>

using namespace Eigen;
using namespace std;


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

void mvis::help::RenameFile(string _source, string _target, bool verbose) {
    if (verbose) { cout << "Renaming file " << _source << " to " << _target << "..." << endl; }

    std::ifstream  src(_source, std::ios::binary);
    std::ofstream  dst(_target, std::ios::binary);

    dst << src.rdbuf();

    dst.close();
    src.close();

    const char* source = _source.c_str();
    if (remove(source) != 0)
        perror("Error deleting file");
    else
        if (verbose) { puts("File successfully deleted"); }
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

bool mvis::help::FileExists(std::string fpath) {
    struct stat buffer;
    return (stat(fpath.c_str(), &buffer) == 0);
}

string mvis::help::GetUniquePath(string fpath) {
    int vcount = 2;
    int padding = 2;
    while (help::FileExists(fpath)) {
        if (vcount > 9) {
            if (vcount > 99) {
                padding = 0;
            }
            else {
                padding = 1;
            }
        }
        string pad = "";
        for (int i = 0; i < padding; i++) {
            pad += "0";
        }
        if (fpath.find("_v") == string::npos) {
            fpath = fpath.replace(fpath.find_last_of("."), 1, "_v" + pad + to_string(vcount) + ".");
        }
        else {
            fpath = fpath.replace(fpath.find_last_of(".") - 5, 5, "_v" + pad + to_string(vcount));
        }
        vcount++;
    }
    return fpath;
}

string mvis::help::GetLatestPath(string templ) {
    std::string path = help::ReplaceOccurrences(templ, "#", "_v002");
    std::string _path = path;
    int vcount = 2;
    int padding = 2;
    while (help::FileExists(_path)) {
        path = _path;
        if (vcount > 9) {
            if (vcount > 99) {
                padding = 0;
            }
            else {
                padding = 1;
            }
        }
        std::string suffix = "_v";
        string pad = "";
        for (int i = 0; i < padding; i++) {
            pad += "0";
        }
        suffix += pad;
        suffix += to_string(vcount);
        _path = help::ReplaceOccurrences(templ, "#", suffix);
        vcount++;
    }

    return path;
}

double mvis::help::fisqrt(float n)
{
    float y = n;
    long i = *(long*)&y;
    i = 0x5f3759df - (i >> 1);
    y = *(float*)&i;

    return y * (1.5f - ((n * 0.5f) * y * y));
}

void mvis::help::WriteToCSV(string csv_path, vector<string> data, string headers) {
    std::ofstream fileStream;
    fileStream.open(csv_path, std::fstream::out);
    fileStream << headers << "\n";

    for (int i = 0; i < data.size(); i++) {
        string line = data[i];
        fileStream << line + "\n";
    }

    fileStream.close();
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
