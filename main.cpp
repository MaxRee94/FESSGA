#include <igl/opengl/glfw/Viewer.h>
#include "io.h"
#include "gui.h"
#include "meshing.h"
#include "physics.h"
#include "tests.h"
#include "fess.h"

using namespace Eigen;
using namespace std;
using namespace fessga;


int main(int argc, char* argv[])
{
    // Grid parameters
    mesher::Grid3D grid3d;
    grid3d.x = 5;
    grid3d.y = 5;
    grid3d.z = 5;
    float domain_size = 2.0; // TODO: replace

    // Initialize mesh lists
    vector<MatrixXd> V_list;
    vector<MatrixXi> F_list;

    // Initialize gui
    GUI gui = GUI(V_list, F_list);

    // Load example object
    MatrixXd V;
    MatrixXi F;
    fessga::IO::ReadMesh("../data/test_objects/teapot.obj", V, F);
    gui.load_example(&V, &F);
    mesher::SurfaceMesh surface_mesh = mesher::create_surface_mesh(&V, &F);

    // -- Normalize mesh --
    // Align barycenter to world origin
    // Get bounding box (min and max for x,y,z)
    // Get cell size along each dimension
    double cell_size = domain_size / (double)grid3d.x;
    // Get offset along each dimension
    Vector3d offset = -cell_size * 0.5 * Vector3d((double)grid3d.x, (double)grid3d.y, (double)grid3d.z);

    string output_folder = "E:/Development/FESSGA/data/msh_output/test/7_element_project";

#if 0:
#elif 1:

    // Generate grid-based binary density distribution based on the given (unstructured) mesh file
    uint32_t* densities = new uint32_t[grid3d.x * grid3d.y * grid3d.z];
    mesher::generate_3d_density_distribution(grid3d, offset, cell_size, &gui.V_list[0], &gui.F_list[0], densities);

    // Create slice from 3d binary density distribution for 2d test
    int z = grid3d.x / 2;
    uint* slice_2d = new uint[grid3d.x * grid3d.y];
    mesher::create_2d_slice(densities, slice_2d, grid3d, z);
    mesher::filter_2d_density_distrib(slice_2d, grid3d.x, grid3d.y);
    mesher::print_2d_density_distrib(slice_2d, grid3d.x, grid3d.y);

    // Obtain a grid-based FE representation based on the chosen mesh
    mesher::FEMesh2D fe_mesh;
    mesher::generate_FE_mesh(grid3d.x, grid3d.y, offset, cell_size, slice_2d, fe_mesh);

    // Encode the FE mesh in a .msh format
    string msh_output_path = output_folder + "/mesh.msh";
    string msh_description;
    mesher::generate_msh_description(&fe_mesh, msh_description);

    // Export the .msh description to a .msh file
    IO::write_text_to_file(msh_description, msh_output_path);

    // Export the FE mesh in Elmer format (in .header, .nodes, .elments, .boundaries files)
    mesher::export_as_elmer_files(&fe_mesh, output_folder);

    cout << "Finished." << endl;

    //gui.show();
#elif 0
    float domain_size = 2.0;
    float cell_size = domain_size / (float)dim_x;
    float inv_cell_size = 1.0 / cell_size;
    Vector2d offset = -cell_size * 0.5 * Vector2d((double)dim_x, (double)dim_y);
    double* vonmises = new double[(dim_x + 1) * (dim_y + 1)]; // Nodes grid has +1 width along each dimension
    string filename = "../data/msh_output/case0001.vtk";
    cout << "started reading physics data..." << endl;
    load_2d_physics_data(filename, vonmises, dim_x + 1, dim_y + 1, offset, inv_cell_size);
    cout << "finished reading physics data." << endl;
#elif 0
    // Do crossover test
    Tester tester = Tester();
    tester.test_2d_crossover();
#elif 0
    // Do FESS test
    Tester tester = Tester();
    tester.test_fess();
#endif
}
