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
    int x_dim = 30;
    int y_dim = 30;
    int z_dim = 30;
#if 0:
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
    
    float cell_size = domain_size / (float)x_dim;
    Vector3d offset = -cell_size * 0.5 * Vector3d((double)x_dim, (double)y_dim, (double)z_dim);
    uint32_t* densities = new uint32_t[x_dim * y_dim * z_dim];
    mesher::generate_msh(
        x_dim, y_dim, z_dim, cell_size, offset, &gui.V_list[0], &gui.F_list[0],
        &line_bounds, densities, mesh_description
    );
    cout << mesh_description << endl;

    // Write mesh description to .msh file
    string output_path = "../data/msh_output/test.msh";
    IO::write_text_to_file(mesh_description, output_path);

    gui.show();
#elif 1
    SparseMatrix<double> VonmisesStress;
    string filename = "D:/OneDrive/Documenten/CSYST_Project/geometry/gmesh/test_2elements/case0001.vtk";
    //string filename = "D:/OneDrive/Documenten/CSYST_Project/geometry/gmesh/test_2elements/case_t0001.vtu";
    load_physics_data(filename, &VonmisesStress);
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

#elif 0
    // Do tests
    Tester tester = Tester();
    tester.test_2d_crossover();

#endif
}
