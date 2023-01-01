#include <iostream>
#include "io.h"

using namespace std;

void IO::ReadMesh(std::string fpath, Eigen::MatrixXd& V, Eigen::MatrixXi& F, bool suppress_output)
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
}
