#include <igl/opengl/glfw/Viewer.h>
#include "controller.h"
#include "tests.h"

using namespace Eigen;
using namespace std;
using namespace fessga;



// Parse cli args
void parse_args(
    int argc, char* argv[], Input& input, string& output_folder, string& action,
    int& dim_x, float& size
    ) {
    action = argv[1];
    output_folder = "E:/Development/FESSGA/data/" + string(argv[2]);
    string relative_path = string(argv[3]);
    if (help::ends_with(string(argv[3]), ".obj")) {
        input.path = "../data/objects/" + string(argv[3]);
        input.name = string(argv[3]);
    }
    else if (help::ends_with(string(argv[3]), ".jpg")) {
        input.type = "image";
        input.path = "../data/images/" + string(argv[3]);
        size = atof(argv[5]);
    }
    else if (string(argv[3]) == "distribution2d") {
        input.type = "distribution2d";
        input.path = output_folder + "/distribution2d" + ".dens";
        size = atof(argv[5]);
    }
    else if (string(argv[3]) == "distribution3d") {
        input.type = "distribution3d";
        input.path = output_folder + "/distribution3d" + ".dens";
        size = atof(argv[5]);
    }
    cout << "Input type: " << input.type << endl;
    dim_x = stoi(string(argv[4]));
    input.name = string(argv[6]);
}


void run_tests(Controller* controller) {
    Tester tester(controller);
    tester.run_tests();
}


int main(int argc, char* argv[])
{
    // Parse arguments
    string output_folder, action;
    Input input;
    bool load_distribution;
    int dim_x;
    float size;
    parse_args(argc, argv, input, output_folder, action, dim_x, size);

    Controller controller = Controller(input, output_folder, action, dim_x, size);
    if (action == "test") {
        run_tests(&controller);
    }
}
