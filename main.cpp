// debug_malloc.cpp
// compile by using: cl /EHsc /W4 /D_DEBUG /MDd debug_malloc.cpp
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include "crtdbg.h"

//#include <igl/opengl/glfw/Viewer.h>
#define _HAS_STD_BYTE 0
#include "controller.h"
#include "tests.h"


// Parse cli args
void parse_args(
    int argc, char* argv[], Input& input, string& base_folder, string& action,
    int& dim_x
    ) {
    action = argv[1];
    base_folder = "E:/Development/FESSGA/data/" + string(argv[2]);
    string relative_path = string(argv[3]);
    if (help::ends_with(string(argv[3]), ".obj")) {
        input.path = base_folder + "/" + string(argv[3]);
        input.name = string(argv[3]);
    }
    else if (help::ends_with(string(argv[3]), ".jpg")) {
        input.type = "image";
        input.path = base_folder + "/" + string(argv[3]);
        input.size = atof(argv[5]);
    }
    else if (string(argv[3]) == "distribution2d") {
        input.type = "distribution2d";
        input.path = base_folder + "/distribution2d" + ".dens";
        input.size = atof(argv[5]);
    }
    else if (string(argv[3]) == "distribution3d") {
        input.type = "distribution3d";
        input.path = base_folder + "/distribution3d" + ".dens";
        input.size = atof(argv[5]);
    }
    cout << "Input type: " << input.type << endl;
    dim_x = stoi(string(argv[4]));
    input.name = string(argv[6]);
    input.max_stress = atof(argv[7]);
    input.max_iterations = atoi(argv[8]);
}


void run_tests(Controller* controller) {
    Tester tester(controller);
    tester.run_tests();
}


int main(int argc, char* argv[])
{
    _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);

    // Parse arguments
    string base_folder, action;
    Input input;
    bool load_distribution;
    int dim_x;
    parse_args(argc, argv, input, base_folder, action, dim_x);

    cout << "size (input): " << input.size << endl;
    Controller controller = Controller(input, base_folder, action, dim_x);
    if (action == "test") {
        run_tests(&controller);
    }
}
