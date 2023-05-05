#include <igl/opengl/glfw/Viewer.h>
#include "io.h"
#include "gui.h"
#include "voxelize.h"

using namespace Eigen;
using namespace std;
using namespace mvis;


/* Generate a grid-based description of a FE mesh that can be output as a .msh file
Input: 
    dim_x, dim_y, dim_z (int):  Number of cells along each dimension of the grid
    csize (float): Size of a single cell (size=width=height=depth)
    V (MatrixXd*): Pointer to the matrix of vertex positions for the given mesh
    F (MatrixXi*): Pointer to the matrix of faces for the given mesh
    msh (string): The generated description string in .msh-format
*/
void generate_msh(
    const int dim_x, const int dim_y, const int dim_z, const float cell_size, MatrixXd* V, MatrixXi* F,
    vector<uint64_t>* bounds, uint32_t* densities, string& msh
) {
    // Generate grid-based FE mesh based on the given (unstructured) mesh file
    vector<vector<float>> nodes;
    vector<vector<int>> elements;  // List of elements. One element is {<number>, <type>, <tag>, <node_1>, ..., <node_n>}
    generate_density_distribution(dim_x, dim_y, dim_z, cell_size, V, F, densities, nodes, elements, bounds);
    print_density_distrib(densities);
    //generate_FE_mesh();

#if 1 // Test with hardcoded nodes- and elements lists
    nodes.push_back({ 0.0, 0.0, 0.0 });
    nodes.push_back({ 1.0, 0.0, 0.0 });
    nodes.push_back({ 1.0, 1.0, 0.0 });
    nodes.push_back({ 0.0, 1.0, 0.0 });
    nodes.push_back({ 2.0, 0.0, 0.0 });
    nodes.push_back({ 2.0, 1.0, 0.0 });
    elements.push_back({ 3, 1, 0, 1, 1, 4 });
    elements.push_back({ 4, 1, 0, 2, 4, 3 });
    elements.push_back({ 5, 1, 0, 2, 3, 6 });
    elements.push_back({ 6, 1, 0, 3, 6, 5 });
    elements.push_back({ 7, 1, 0, 4, 5, 2 });
    elements.push_back({ 8, 1, 0, 4, 2, 1 });
    elements.push_back({ 1, 3, 0, 2, 1, 2, 3, 4 });
    elements.push_back({ 2, 3, 0, 2, 2, 5, 6, 3 });
#endif

    // Encode mesh data into .msh-description
    // -- Format section
    msh = {
        "$MeshFormat\n"
        "2.0 0 8\n"
        "$EndMeshFormat\n"
    };

    // -- Nodes section
    msh += "$Nodes\n";
    msh += to_string(nodes.size()) + "\n";         // Number of nodes
    for (int i = 0; i < nodes.size(); i++) {        // List of nodes, with each node encoded as <node_idx> <x> <y> <z>
        msh += to_string(i + 1) + " ";             // 1-based node index
        msh += to_string(nodes[i][0]) + " ";       // x
        msh += to_string(nodes[i][1]) + " ";       // y
        msh += to_string(nodes[i][2]) + "\n";      // z
    }
    msh += "$EndNodes\n";

    // -- Elements section
    msh += "$Elements\n";
    int no_elements = elements.size();
    msh += to_string(no_elements) + "\n";

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
        msh += to_string(element[0]) + " ";        // Element number
        msh += to_string(element[1]) + " ";        // Element type
        msh += "2 ";                               // Number of tags  TODO: Find out if this has any effect
        msh += to_string(element[2]) + " ";        // Element tag (indicating which boundary it belongs to, if any)
        for (int j = 3; j < element.size(); j++) {  // Node list
            msh += to_string(element[j]) + " ";
        }
        msh += "\n";
    }

    // End elements section
    msh += "$EndElements\n";
}


int main(int argc, char* argv[])
{
    // Initialize mesh lists
    vector<MatrixXd> V_list;
    vector<MatrixXi> F_list;

    // Initialize gui
    GUI gui = GUI(V_list, F_list);

    // Load example cube
    gui.load_example();

    // Obtain a grid-based FE representation based on the chosen mesh, encoded in a .msh format
    vector<uint64_t> bounds = {};
    string mesh_description = "";
    int x_dim = 6;
    int y_dim = 6;
    int z_dim = 6;
    uint32_t* densities = new uint32_t[x_dim * y_dim * z_dim];
    generate_msh(x_dim, y_dim, z_dim, 1, &gui.V_list[0], &gui.F_list[0], &bounds, densities, mesh_description);
    cout << mesh_description << endl;

    // Write mesh description to .msh file
    string output_path = "../data/msh_output/test.msh";
    IO::write_text_to_file(mesh_description, output_path);

    gui.show();
}
