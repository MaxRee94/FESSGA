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
#include <vtkUnstructuredGridReader.h>
#include <vtkAutoInit.h>
#include <vtkTypedDataArray.h>
#include <vtkTypedDataArrayIterator.h>
#include <vtkArrayCalculator.h>
#include <vtkTypeInt64Array.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include "vtkUnstructuredGridAlgorithm.h"

VTK_MODULE_INIT(vtkRenderingOpenGL2)


void load_physics_data(string filename, SparseMatrix<double>* VonmisesStress)
{
    // Read data from file
    //vtkNew<vtkXMLUnstructuredGridReader> reader;
    vtkNew<vtkUnstructuredGridReader> reader;
    reader->SetFileName(filename.c_str());
    reader->ReadAllScalarsOn();
    reader->Update();
    vtkUnstructuredGrid* output = reader->GetOutput();

    // Get node indices
    vtkPoints* points = output->GetPoints();
    for (int i = 0; i < points->GetNumberOfPoints(); i++) {
        double* point = points->GetData()->GetTuple(i);
        //cout << "point " << i + 1 << ": " << point[0] << ", " << point[1] << endl;
    }
    cout << endl;

    // Get point ids
    /*vtkIdList* point_ids;
    output->GetCellPoints(8, point_ids);
    for (int i = 0; i < point_ids->GetNumberOfIds(); i++) {
        cout << point_ids->GetId(0) << endl;
    }
    cout << endl;*/
    
    // Get point data (this object contains the physics data)
    vtkPointData* point_data = output->GetPointData();

    cout << "number of scalars: " << reader->GetNumberOfScalarsInFile() << endl;

    // Obtain Von Mises stress array
    vtkDataArray* vonmises_data = point_data->GetScalars("Vonmises");
    cout << "vonmises data array size: " << vonmises_data->GetDataSize() << endl;
    vtkDoubleArray* vonmises_array = dynamic_cast<vtkDoubleArray*>(vonmises_data);

    cout << "num values: " << vonmises_array->GetNumberOfValues() << endl;

    // Reformat data as sparse matrix which contains an entry for each cell in the grid
    // (including cells not present in the FE mesh, these will have the value 0. Hence the sparse matrix representation).
    // TODO: do this by first obtaining the node indices (with output->getCellData?) and then using these to infer which values
    // in the stress array correspond to which cells in the grid.
    for (int i = 0; i < vonmises_array->GetNumberOfValues(); i++) {
        cout << vonmises_array->GetValue(i) << endl;
    }
}