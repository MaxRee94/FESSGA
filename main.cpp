#include <igl/opengl/glfw/Viewer.h>
#include "io.h"
#include "controller.h"
#include "fess.h"
#include "images.h"

using namespace Eigen;
using namespace std;
using namespace fessga;



// Parse cli args
void parse_args(
    int argc, char* argv[], Input& input, string& output_folder, string& action,
    int& dim_x, int& dim_y, int& dim_z, float& size
    ) {
    if (argc > 1) action = argv[1];
    if (argc > 2) output_folder = "E:/Development/FESSGA/data/" + string(argv[2]);
    if (argc > 3) {
        string relative_path = string(argv[3]);
        input.path = "../data/objects/" + string(argv[3]);
        input.name = string(argv[3]);
        if (help::ends_with(string(argv[3]), ".jpg")) {
            input.type = "image";
            input.path = "../data/images/" + string(argv[3]);
            size = atof(argv[7]);
        }
        else if (string(argv[3]) == "distribution") {
            input.type = "distribution";
            input.path = output_folder + "/" + string(argv[3]) + ".dens";
            size = atof(argv[7]);
        }
        cout << "Input type: " << input.type << endl;
    }

    if (argc > 4) dim_x = stoi(string(argv[4]));
    else dim_x = 40;
    if (argc > 5) dim_y = stoi(string(argv[5]));
    else dim_y = 40;
    if (argc > 6) dim_z = stoi(string(argv[6]));
    else dim_z = 40;
    if (argc > 8) input.name = string(argv[8]);
}


int main(int argc, char* argv[])
{
    // Parse arguments
    string output_folder, action;
    Input input;
    bool load_distribution;
    int dim_x, dim_y, dim_z;
    float size;
    parse_args(argc, argv, input, output_folder, action, dim_x, dim_y, dim_z, size);

    Controller(input, output_folder, action, dim_x, dim_y, dim_z, size);
}
