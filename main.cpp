#include <igl/opengl/glfw/Viewer.h>
#include "io.h"
#include "gui.h"
#include "meshing.h"
#include "physics.h"

using namespace Eigen;
using namespace std;
using namespace mvis;


#define TEST2D true;
#define TEST3D false;

/* Generate a grid-based description of a FE mesh that can be output as a .msh file
Input: 
    dim_x, dim_y, dim_z (int):  Number of cells along each dimension of the grid
    csize (float): Size of a single cell (size=width=height=depth)
    V (MatrixXd*): Pointer to the matrix of vertex positions for the given mesh
    F (MatrixXi*): Pointer to the matrix of faces for the given mesh
    msh (string): The generated description string in .msh-format
*/
void generate_msh(
    const int dim_x, const int dim_y, const int dim_z, const float cell_size, Vector3d offset, MatrixXd* V, MatrixXi* F,
    map<uint32_t, uint32_t>* line_bounds, uint32_t* densities, string& msh
) {
    // Generate grid-based binary density distribution based on the given (unstructured) mesh file
    vector<string> nodes;
    vector<vector<int>> elements;  // List of elements. One element is {<number>, <type>, <tag>, <node_1>, ..., <node_n>}
    generate_3d_density_distribution(dim_x, dim_y, dim_z, offset, cell_size, V, F, densities);
    print_density_distrib(densities, dim_x);

    // Create slice from 3d binary density distribution for 2d test
    int z = dim_x / 2;
    uint32_t* slice_2d = new uint32_t[dim_x * dim_y];
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            slice_2d[x * dim_x + y] = densities[z * dim_x * dim_y + x * dim_x + y];
        }
    }

    // Generate Finite Element mesh from binary density distribution
    generate_2d_FE_mesh(dim_x, dim_y, offset, cell_size, slice_2d, nodes, elements, line_bounds);

    // Encode mesh data into .msh-description
    // -- Format section
    msh = {
        "$MeshFormat\n"
        "2.0 0 8\n"
        "$EndMeshFormat\n"
    };

    // -- Nodes section
    msh += "$Nodes\n";
    msh += to_string(nodes.size()) + "\n";          // Number of nodes
    for (int i = 0; i < nodes.size(); i++) {        // List of nodes, with each node encoded as <node_idx> <x> <y> <z>
        msh += nodes[i];
    }
    msh += "$EndNodes\n";

    // -- Elements section
    msh += "$Elements\n";
    int no_elements = elements.size();
    msh += to_string(no_elements) + "\n";

    // Iterate over all elements and add to description string
    // Each element is encoded as <elm-number> <elm-type> <number-of-tags> <physical entity> <tag> <node_1> ... <node_n>
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
        for (int j = 0; j < element.size(); j++) {
            msh += to_string(element[j]) + " ";
        }
        msh += "\n";
    }

    // End elements section
    msh += "$EndElements\n";
}


int main(int argc, char* argv[])
{
#if 1:
    // Initialize mesh lists
    vector<MatrixXd> V_list;
    vector<MatrixXi> F_list;

    // Initialize gui
    GUI gui = GUI(V_list, F_list);

    // Load example cube
    gui.load_example();

    // Obtain a grid-based FE representation based on the chosen mesh, encoded in a .msh format
    map<uint32_t, uint32_t> line_bounds;
    string mesh_description = "";
    float domain_size = 2.0;
    int x_dim = 30;
    int y_dim = 30;
    int z_dim = 30;
    float cell_size = domain_size / (float)x_dim;
    Vector3d offset = -cell_size * 0.5 * Vector3d((double)x_dim, (double)y_dim, (double)z_dim);
    uint32_t* densities = new uint32_t[x_dim * y_dim * z_dim];
    generate_msh(
        x_dim, y_dim, z_dim, cell_size, offset, &gui.V_list[0], &gui.F_list[0],
        &line_bounds, densities, mesh_description
    );
    cout << mesh_description << endl;

    // Write mesh description to .msh file
    string output_path = "../data/msh_output/test.msh";
    IO::write_text_to_file(mesh_description, output_path);

    gui.show();
#elif 0
    MatrixXd stress;
    string filename = "D:/OneDrive/Documenten/CSYST_Project/geometry/gmesh/test_2elements/case_t0001.vtu";
    load_physics_data(filename, stress);
#elif 0
    std::string filename = "D:/OneDrive/Documenten/CSYST_Project/geometry/gmesh/test_2elements/case_t0001.vtu";

    // read all the data from the file
    vtkNew<vtkXMLUnstructuredGridReader> reader;
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkNew<vtkNamedColors> colors;

    // Create a mapper and actor
    vtkNew<vtkDataSetMapper> mapper;
    mapper->SetInputConnection(reader->GetOutputPort());
    mapper->ScalarVisibilityOff();

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->EdgeVisibilityOn();
    actor->GetProperty()->SetLineWidth(2.0);
    actor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());

    vtkNew<vtkProperty> backFace;
    backFace->SetColor(colors->GetColor3d("Tomato").GetData());
    actor->SetBackfaceProperty(backFace);

    // Create a renderer, render window, and interactor
    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->AddRenderer(renderer);
    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // Add the actor to the scene
    renderer->AddActor(actor);
    renderer->SetBackground(colors->GetColor3d("Wheat").GetData());

    // Render and interact
    renderWindow->SetSize(640, 480);

    renderWindow->Render();
    renderWindowInteractor->Start();

#endif
}
