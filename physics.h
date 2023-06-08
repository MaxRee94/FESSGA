#include <iostream>
#include <stdio.h>
#include <stdlib.h>
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
#define VTK_DEBUG_LEAKS true


namespace fessga {
    class physics {
    public:
        static void call_elmer(string bat_file) {
            std::string command = bat_file;
            std::array<char, 80> buffer;
            FILE* pipe = _popen(command.c_str(), "r");
            while (fgets(buffer.data(), 80, pipe) != NULL) {
                std::cout << buffer.data();
            }
            _pclose(pipe);
        }

        static void load_2d_physics_data(
            string filename, double* vonmises, int dim_x, int dim_y, Vector2d offset, Vector2d inv_cell_size)
        {
            // Read data from file
            vtkUnstructuredGridReader* reader = vtkUnstructuredGridReader::New();
            reader->SetFileName(filename.c_str());
            reader->ReadAllScalarsOn();
            reader->Update();
            vtkUnstructuredGrid* output = reader->GetOutput();

            // Initially, populate vonmises array with zeroes (nodes on the grid which are part of the FE mesh will have their
            // corresponding values in the vonmises array overwritten later)
            for (int x = 0; x < dim_x; x++) {
                for (int y = 0; y < dim_y; y++) {
                    vonmises[x * dim_y + y] = 0.0;
                }
            }

            // Get point data (this object contains the physics data)
            vtkPointData* point_data = output->GetPointData();

            // Obtain Von Mises stress array
            vtkDoubleArray* vonmises_array = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars("Vonmises"));

            // Overwrite grid values with values from vonmises array (only for nodes with coordinates that lie within the FE mesh)
            vtkPoints* points = output->GetPoints();
            double* point = new double[3];
            vector<int> coords = {};
            for (int i = 0; i < points->GetNumberOfPoints(); i++) {
                point = points->GetData()->GetTuple(i);
                Vector2d origin_aligned_coord = Vector2d(point[0], point[1]) - offset;
                Vector2d gridscale_coord = inv_cell_size.cwiseProduct(origin_aligned_coord);
                int coord = (int)(gridscale_coord[0] * dim_y + gridscale_coord[1]);
                int x = coord / dim_y;
                int y = coord % dim_y;
                coords.push_back(coord);
                vonmises[coord] = (double)vonmises_array->GetValue(i);
            }
        }
    };
}



