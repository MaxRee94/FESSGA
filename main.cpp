#include <igl/opengl/glfw/Viewer.h>
#include "io.h"
#include "gui.h"
#include "meshing.h"
#include "physics.h"
#include "controller.h"
#include "fess.h"
#include "images.h"

using namespace Eigen;
using namespace std;
using namespace fessga;


struct Input {
    string type = "object";
    string path;
};


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
        if (help::ends_with(string(argv[3]), ".jpg")) {
            input.type = "image";
            input.path = "../data/images/" + string(argv[3]);
            size = atof(argv[7]);
        }
        else if (argv[3] == "distribution") {
            input.type = "distribution";
            input.path = output_folder + string(argv[3]);
        }
    }

    if (argc > 4) dim_x = stoi(string(argv[4]));
    else dim_x = 40;
    if (argc > 5) dim_y = stoi(string(argv[5]));
    else dim_y = 40;
    if (argc > 6) dim_z = stoi(string(argv[6]));
    else dim_z = 40;
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

    // Initialize mesh lists
    vector<MatrixXd> V_list;
    vector<MatrixXi> F_list;

    // Initialize gui
    GUI gui = GUI(V_list, F_list);

    // Load mesh
    MatrixXd V;
    MatrixXi F;
    mesher::SurfaceMesh surface_mesh;
    if (input.type == "image" || input.type == "distribution") {
        Vector3d object_size(size, size, size);
        surface_mesh = mesher::SurfaceMesh(object_size);
    }
    else {
        fessga::IO::read_mesh(input.path, V, F);
        gui.load_example(&V, &F);
        surface_mesh = mesher::SurfaceMesh(&V, &F);
    }

    // Create 3d grid
    mesher::Grid3D grid = mesher::create_grid3d(100, 100, dim_z, surface_mesh.diagonal);

    // Create output folder
    IO::create_folder_if_not_exists(output_folder);

    // Load distribution from file or generate it on the fly
    uint* slice_2d = new uint[grid.x * grid.y];
    if (input.type == "distribution") {
        cout << "loading densities from .dens file...\n";
        mesher::import_densities(input.path, slice_2d);
    }
    else if (input.type == "image") {
        cout << "Loading densities from image...\n";
        char img_path[200];
        strcpy(img_path, input.path.c_str());
        img::load_distribution_from_image(img_path, slice_2d, grid);
    }
    else {
        // Generate grid-based binary density distribution based on the given (unstructured) mesh file
        uint* densities = new uint[grid.x * grid.y * grid.z];
        mesher::generate_3d_density_distribution(grid, surface_mesh, &gui.V_list[0], &gui.F_list[0], densities);

        // Create slice from 3d binary density distribution for 2d surface generation
        int z = grid.x / 2;
        mesher::create_2d_slice(densities, slice_2d, grid, z);
        mesher::filter_2d_density_distrib(slice_2d, grid.x, grid.y);

        // Export density distribution
        string densities_file = mesher::export_density_distrib(output_folder, slice_2d, grid.x, grid.y);
        cout << "FESS: Exported current density distribution to file: " << densities_file << endl;
    }

    // Execute action
    if (action == "export_mesh") {
        mesher::print_density_distrib(slice_2d, grid.x, grid.y);

        // Obtain a grid-based FE representation based on the chosen mesh
        mesher::FEMesh2D fe_mesh;
        mesher::generate_FE_mesh(grid, surface_mesh, slice_2d, fe_mesh);

        // Export the FE mesh in .msh format
        mesher::export_as_msh_file(&fe_mesh, output_folder);

        // Export the FE mesh in Elmer format (in .header, .nodes, .elments, .boundaries files)
        mesher::export_as_elmer_files(&fe_mesh, output_folder);

        cout << "Finished." << endl;
    }
    else if (action == "evolve") {
        // Do crossover test
        Controller controller = Controller();
        controller.test_2d_crossover();
    }
    else if (action == "fess") {
        // Do FESS
        Controller controller = Controller();
        controller.run_fess();
    }
}
