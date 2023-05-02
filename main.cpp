#include <igl/opengl/glfw/Viewer.h>
#include "io.h"
#include "gui.h"

using namespace Eigen;
using namespace std;
using namespace mvis;


/* Generate a binary density distribution on the grid based on the given mesh file
Input:
Returns (by reference):
    DensityDistrib (MatrixXi*): Matrix which contains a binary density value for each cell in the grid
*/
void generate_density_distribution(
    int dim_x, int dim_y, int dim_z, float cell_size,
    MatrixXd* V, MatrixXi* F, MatrixXi* DensityDistrib,
    int* no_nodes, vector<vector<int>>* lines, vector<vector<int>>* quads)
{
    no_nodes = 0;
}


/* Generate a grid-based description of a FE mesh that can be output as a .msh file
Input: 
    dim_x, dim_y, dim_z (int):  Number of cells along each dimension of the grid
    csize (float): Size of a single cell (size=width=height=depth)
    V (MatrixXd*): Pointer to the matrix of vertex positions for the given mesh
    F (MatrixXi*): Pointer to the matrix of faces for the given mesh
    msh (string*): Pointer to the msh file description string to be generated
*/
void generate_msh(
    int dim_x, int dim_y, int dim_z, float cell_size, MatrixXd* V, MatrixXi* F, string* msh
) {
    // Generate a binary density distribution on the grid based on the given mesh file
    MatrixXi DensityDistrib;
    int no_nodes;
    vector<vector<int>> lines;  // List of lines; one line is a vector {<tag>, <node1>, <node2>}
    vector<vector<int>> quads;  // List of quads; one quad is a vector {<tag>, <node1>, <node2>, <node3>, <node4>}
    vector<vector<int>> voxels; // List of voxels; one voxel is a vector {<tag>, <node1>, <node2>, ..., <node8>}
    generate_density_distribution(
        dim_x, dim_y, dim_z, cell_size, V, F, &DensityDistrib, &no_nodes, &lines, &quads
    );

    // Encode mesh data into .msh-description
    // -- Format section
    string format_str = {
        "$MeshFormat \n"
        "2.0 0 8 \n"
        "$EndMeshFormat"
    };

    // -- Nodes section
    string entities_str = "$Entities\n";
    entities_str += to_string(no_nodes) + "\n";
    entities_str += "$EndEntities\n";

    // -- Elements section
    string elements_str = "$Elements\n";
    int num_elements = lines.size() + quads.size() + voxels.size();
    elements_str += to_string(num_elements) + "\n";
    
    // Add lines
    

    // Add quads

    // End elements section
    elements_str += "$EndElements";
}


int main(int argc, char *argv[])
{
    vector<MatrixXd> V_list;
    vector<MatrixXi> F_list;

    GUI gui = GUI(V_list, F_list);
    gui.load_example();

    string mesh_description;
    generate_msh(5, 5, 0, 1, &gui.V_list->at(0), &gui.F_list->at(0), &mesh_description);

    gui.show();
}
