#include <igl/opengl/glfw/Viewer.h>
#include "io.h"
#include "gui.h"
#include "meshing.h"
#include "physics.h"
#include "controller.h"
#include "fess.h"

using namespace Eigen;
using namespace std;
using namespace fessga;


int main(int argc, char* argv[])
{
    // Initialize mesh lists
    vector<MatrixXd> V_list;
    vector<MatrixXi> F_list;

    // Initialize gui
    GUI gui = GUI(V_list, F_list);

    // Load example object
    MatrixXd V;
    MatrixXi F;
    fessga::IO::read_mesh("../data/test_objects/teapot.obj", V, F);
    gui.load_example(&V, &F);
    mesher::SurfaceMesh surface_mesh = mesher::create_surface_mesh(&V, &F);

    // Create 3d Grid
    mesher::Grid3D grid = mesher::create_grid3d(40, 40, 40, surface_mesh.diagonal);

    // Set output folder
    string output_folder = "E:/Development/FESSGA/data/msh_output/FESS_lowres_test";
    IO::create_folder_if_not_exists(output_folder);

    uint* slice_2d = new uint[grid.x * grid.y];
#if 0:
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
#elif 0
    string densities_path = output_folder + "/distribution.dens";
    mesher::import_densities(densities_path, slice_2d);
#endif

#if 0:
    mesher::print_density_distrib(slice_2d, grid.x, grid.y);

    // Obtain a grid-based FE representation based on the chosen mesh
    mesher::FEMesh2D fe_mesh;
    mesher::generate_FE_mesh(grid, surface_mesh, slice_2d, fe_mesh);

    // Export the FE mesh in .msh format
    mesher::export_as_msh_file(&fe_mesh, output_folder);

    // Export the FE mesh in Elmer format (in .header, .nodes, .elments, .boundaries files)
    mesher::export_as_elmer_files(&fe_mesh, output_folder);

    cout << "Finished." << endl;

    //gui.show();
#elif 0
    float domain_size = 2.0;
    Vector2d inv_cell_size = Vector2d(1.0 / grid.cell_size(0), 1.0 / grid.cell_size(1));
    Vector2d offset = Vector2d(surface_mesh.offset(0), surface_mesh.offset(1));
    double* vonmises = new double[(grid.x + 1) * (grid.y + 1)]; // Nodes grid has +1 width along each dimension
    string filename = "../data/msh_output/test/7_element_project/fessga_generated_elmer_files/case0001.vtk";
    cout << "started reading physics data..." << endl;
    physics::load_2d_physics_data(filename, vonmises, grid.x + 1, grid.y + 1, offset, inv_cell_size);
    cout << "finished reading physics data." << endl;
#elif 0
    // Do crossover test
    Controller controller = Controller();
    controller.test_2d_crossover();
#elif 1
    // Do FESS
    Controller controller = Controller();
    controller.run_fess();
#endif
}
