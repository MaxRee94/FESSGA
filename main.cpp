#include <igl/opengl/glfw/Viewer.h>
#include "io.h"
#include "gui.h"
#include "meshing.h"
#include "physics.h"
#include "tests.h"

using namespace Eigen;
using namespace std;
using namespace fessga;


int main(int argc, char* argv[])
{
    int dim_x = 100;
    int dim_y = 100;
    int dim_z = 100;
#if 1:
    // Initialize mesh lists
    vector<MatrixXd> V_list;
    vector<MatrixXi> F_list;

    // Initialize gui
    GUI gui = GUI(V_list, F_list);

    // Load example object
    gui.load_example();

    // Obtain a grid-based FE representation based on the chosen mesh, encoded in a .msh format
    map<uint32_t, uint32_t> line_bounds;
    string mesh_description = "";
    float domain_size = 2.0;
    
    float cell_size = domain_size / (float)dim_x;
    Vector3d offset = -cell_size * 0.5 * Vector3d((double)dim_x, (double)dim_y, (double)dim_z);
    uint32_t* densities = new uint32_t[dim_x * dim_y * dim_z];
    mesher::generate_msh(
        dim_x, dim_y, dim_z, cell_size, offset, &gui.V_list[0], &gui.F_list[0],
        &line_bounds, densities, mesh_description
    );
    //cout << mesh_description << endl;

    // Write mesh description to .msh file
    string output_path = "../data/msh_output/test.msh";
    IO::write_text_to_file(mesh_description, output_path);

    gui.show();
#elif 0
    double* vonmises = new double[dim_x * dim_y * dim_z];
    string filename = "D:/OneDrive/Documenten/CSYST_Project/geometry/gmesh/test_2elements/case0001.vtk";
    cout << "started reading physics data..." << endl;
    load_physics_data(filename, vonmises, dim_x, dim_y, dim_z);
    cout << "finished reading physics data." << endl;
#elif 0
    // Do crossover test
    Tester tester = Tester();
    tester.test_2d_crossover();

#endif
}
