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
#include <vtkTypedDataArray.h>
#include <vtkTypedDataArrayIterator.h>
#include <vtkArrayCalculator.h>
#include <vtkTypeInt64Array.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include "vtkUnstructuredGridAlgorithm.h"

VTK_MODULE_INIT(vtkRenderingOpenGL2)


void load_physics_data(string filename, SparseMatrix<double>* Vonmises)
{
    // Read data from file
    vtkNew<vtkXMLUnstructuredGridReader> reader;
    reader->SetFileName(filename.c_str());
    reader->Update();
    vtkUnstructuredGrid* output = reader->GetOutput();
    vtkPointData* pointData = output->GetPointData();

    // Obtain Von Mises stress array
    vtkDataArray* vonmises_array = pointData->GetArray("vonmises");
    vtkDoubleArray* vonmises_doubles = dynamic_cast<vtkDoubleArray*>(vonmises_array);

    // Reformat data as sparse matrix which contains an entry for each cell in the grid
    // (including cells not present in the FE mesh, these will have the value 0. Hence the sparse matrix representation).
    // TODO: do this by first obtaining the node indices (with output->getCellData?) and then using these to infer which values
    // in the stress array correspond to which cells in the grid.
    for (int i = 0; i < vonmises_doubles->GetNumberOfValues(); i++) {
        cout << vonmises_doubles->GetValue(i) << endl;
    }
}