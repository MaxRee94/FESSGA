#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkAutoInit.h>
#include <vtkArrayCalculator.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2)

void load_physics_data(string filename, MatrixXd stress)
{
    // Read data from file
    vtkNew<vtkXMLUnstructuredGridReader> reader;
    reader->SetFileName(filename.c_str());
    reader->Update();
    vtkUnstructuredGrid* output = reader->GetOutput();

    // Obtain vector- and scalar fields
    //aa->SetInput(output);
    //aa->Assign("displacement", vtkDataSetAttributes::VECTORS, vtkAssignAttribute::POINT_DATA);
    vtkArrayCalculator* calc = vtkArrayCalculator::New();
    //calc->SetInputData(output);

}