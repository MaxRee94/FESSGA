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
    int dim_x, int dim_y, int dim_z, float cell_size, MatrixXd* V, MatrixXi* F, MatrixXi* DensityDistrib,
    MatrixXd* Nodes, vector<vector<int>>* elements, vector<long long int>* bounds
    )
{
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
    int dim_x, int dim_y, int dim_z, float cell_size, MatrixXd* V, MatrixXi* F,
    vector<long long int>* bounds, string* msh
) {
    // Generate a binary density distribution on the grid based on the given mesh file
    MatrixXi DensityDistrib;
    MatrixXd Nodes;
    vector<vector<int>> elements;  // List of elements. One element is {<number>, <type>, <tag>, <node_1>, ..., <node_n>}
    generate_density_distribution(
        dim_x, dim_y, dim_z, cell_size, V, F, &DensityDistrib, &Nodes, &elements, bounds
    );

    // Encode mesh data into .msh-description
    // -- Format section
    *msh += {
        "$MeshFormat\n"
        "2.0 0 8\n"
        "$EndMeshFormat\n"
    };

    // -- Nodes section
    *msh += "$Nodes\n";
    int no_nodes = Nodes.rows();
    *msh += to_string(no_nodes) + "\n";         // Number of nodes
    for (int i = 0; i < no_nodes; i++) {        // List of nodes, with each node encoded as <node_idx> <x> <y> <z>
        *msh += to_string(i + 1) + " ";                // 1-based node index
        *msh += to_string(Nodes.row(i)[0]) + " ";      // x
        *msh += to_string(Nodes.row(i)[1]) + " ";      // y
        *msh += to_string(Nodes.row(i)[2]) + "\n";     // z
    }
    *msh += "$EndNodes\n";

    // -- Elements section
    *msh += "$Elements\n";
    int no_elements = elements.size();
    *msh += to_string(no_elements) + "\n";
    
    // Iterate over all elements and add to description string
    // Each element is encoded as <elm-number> <elm-type> <number-of-tags> <tag> <node_1> ... <node_n>
    // Element types:
    //      1 : 2-node line
    //      2 : 3-node triangle
    //      3 : 4-node quad
    //      4 : 4-node tetrahedron
    //      5 : 8-node hexahedron (a cube is a regular hexahedron)
    //      15: 1-node point
    // For example: 3 1 2 0 1 1 4 encodes a line from node 1 to node 4 with boundary tag 1
    for (int i = 0; i < elements.size(); i++) {
        vector<int> element = elements[i];
        *msh += to_string(element[0]) + " ";        // Element number
        *msh += to_string(element[1]) + " ";        // Element type
        *msh += "2 ";                               // Number of tags  TODO: Find out if this has any effect
        *msh += to_string(element[2]) + " ";        // Element tag (indicating which boundary it belongs to, if any)
        for (int j = 2; j < element.size(); j++) {  // Node list
            *msh += to_string(element[j]) + " ";
        }
    }

    // End elements section
    *msh += "$EndElements\n";
}


int main(int argc, char *argv[])
{
    vector<MatrixXd> V_list;
    vector<MatrixXi> F_list;

    GUI gui = GUI(V_list, F_list);
    gui.load_example();

    string mesh_description;
    vector<long long int> bounds;
    generate_msh(5, 5, 0, 1, &gui.V_list->at(0), &gui.F_list->at(0), &bounds, &mesh_description);

    gui.show();
}
