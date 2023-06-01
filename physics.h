#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkProperty.h>
#include <vtkSmartPointer.h>
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
#include <vtkPoints.h>
#include "vtkUnstructuredGridAlgorithm.h"

VTK_MODULE_INIT(vtkRenderingOpenGL2)


void load_physics_data(string filename, double* vonmises, int dim_x, int dim_y, int dim_z)
{
    // Read data from file
    vtkNew<vtkUnstructuredGridReader> reader;
    reader->SetFileName(filename.c_str());
    reader->ReadAllScalarsOn();
    reader->Update();
    vtkUnstructuredGrid* output = reader->GetOutput();

    // Initially, populate vonmises array with zeroes (nodes on the grid which are part of the FE mesh will have their
    // corresponding values in the vonmises array overwritten later)
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            for (int z = 0; z < dim_z; z++) {
                vonmises[x * dim_z * dim_y + y * dim_z + z] = 0;
            }
        }
    }

    // Get point data (this object contains the physics data)
    vtkPointData* point_data = output->GetPointData();

    // Obtain Von Mises stress array
    vtkDoubleArray* vonmises_array = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars("Vonmises"));

    // Overwrite grid values with values from vonmises array (only for nodes with coordinates that lie within the FE mesh)
    vtkPoints* points = output->GetPoints();
    for (int i = 0; i < points->GetNumberOfPoints(); i++) {
        double* point = points->GetData()->GetTuple(i);
        int coord = (int)(point[0] * dim_x * dim_y + point[1] * dim_x + point[2]);
        vonmises[coord] = vonmises_array->GetValue(i);
    }
}