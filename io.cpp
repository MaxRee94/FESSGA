#include <iostream>
#include <filesystem>
#include "io.h"
#include "Helpers.h"

using namespace std;

void fessga::IO::read_mesh(std::string fpath, Eigen::MatrixXd& V, Eigen::MatrixXi& F, bool suppress_output)
{
    if (!suppress_output) cout << "Reading mesh " << fpath << endl;

    // Check if file path exists
    struct stat buffer;
    if (stat(fpath.c_str(), &buffer) != 0) {
        cout << "File '" + fpath + "' does not exist." << endl;
        throw std::runtime_error("File '" + fpath + "' does not exist.");
    }

    // Read file format matching extension
    std::string extension = fpath.substr(fpath.find_last_of(".") + 1);
    if (extension == "off") {
        igl::readOFF(fpath, V, F);
    }
    else if (extension == "ply") {
        igl::readPLY(fpath, V, F);
    }
    else if (extension == "obj") {
        igl::readOBJ(fpath, V, F);
    }
    else {
        throw std::runtime_error("File type '." + extension + "' not supported.");
    }

    cout << "Finished reading mesh." << endl;
}

void fessga::IO::create_folder_if_not_exists(std::string folder_path) {
    std::filesystem::path dir_path = std::filesystem::path(folder_path);
    if (!std::filesystem::is_directory(dir_path)) {
        cout << "Creating directory " << folder_path << endl;
        std::filesystem::create_directories(folder_path);
    }
}

bool fessga::IO::file_exists(std::string fpath) {
    struct stat buffer;
    return (stat(fpath.c_str(), &buffer) == 0);
}

void fessga::IO::copy_file(std::string _source, std::string _target, bool verbose) {
    if (verbose) { cout << "Renaming file " << _source << " to " << _target << "..." << endl; }

    std::ifstream  src(_source, std::ios::binary);
    std::ofstream  dst(_target, std::ios::binary);

    dst << src.rdbuf();

    dst.close();
    src.close();
}

void fessga::IO::rename_file(string _source, string _target, bool verbose) {
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

void fessga::IO::write_to_csv(string csv_path, vector<string> data, string headers) {
    std::ofstream fileStream;
    fileStream.open(csv_path, std::fstream::out);
    fileStream << headers << "\n";

    for (int i = 0; i < data.size(); i++) {
        string line = data[i];
        fileStream << line + "\n";
    }

    fileStream.close();
}

void fessga::IO::write_text_to_file(string text, string path) {
    std::ofstream fileStream;
    fileStream.open(path, std::fstream::out);
    fileStream << text << "\n";
    fileStream.close();
}

string fessga::IO::get_unique_path(string fpath) {
    int vcount = 2;
    int padding = 2;
    while (IO::file_exists(fpath)) {
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

string fessga::IO::get_latest_path(string templ) {
    std::string path = help::replace_occurrences(templ, "#", "_v002");
    std::string _path = path;
    int vcount = 2;
    int padding = 2;
    while (IO::file_exists(_path)) {
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
        _path = help::replace_occurrences(templ, "#", suffix);
        vcount++;
    }

    return path;
}

bool fessga::IO::is_empty(string folder) {
    for (const auto& _ : std::filesystem::directory_iterator(folder)) {
        return true;
    }
    return false;
}

std::string fessga::IO::get_fullpath(string relative_path) {
    char* _relative_path = (char*)relative_path.c_str();
    char* buffer = 0;
    char* _absolute_path = _fullpath(buffer, _relative_path, 1024);
    string absolute_path(_absolute_path);
    absolute_path = help::replace_occurrences(absolute_path, "\\", "/");

    return absolute_path;
}

void fessga::IO::read_file_content(std::string fpath, std::vector<std::string>& content) {
    // Check if file path exists
    struct stat buffer;
    if (stat(fpath.c_str(), &buffer) != 0) {
        cout << "ERROR: File '" + fpath + "' does not exist." << endl;
        throw std::runtime_error("File '" + fpath + "' does not exist.");
    }

    ifstream fstream;
    string line;
    fstream.open(fpath);
    if (fstream.is_open()) {
        while (getline(fstream, line)) {
            help::replace_occurrences(line, "\n", "");
            content.push_back(line);
        }
    }
}
