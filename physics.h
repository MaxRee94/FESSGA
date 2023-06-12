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
#include "meshing.h"

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

        static void populate_with_zeroes(double* _array, int dim_x, int dim_y) {
            for (int x = 0; x < dim_x; x++) {
                for (int y = 0; y < dim_y; y++) {
                    _array[x * dim_y + y] = 0.0;
                }
            }
        }

        static void load_2d_physics_data(
            string filename, double* results_cellwise, mesher::Grid3D grid, Vector2d offset, char* data_type)
        {
            // Read data from file
            vtkUnstructuredGridReader* reader = vtkUnstructuredGridReader::New();
            reader->SetFileName(filename.c_str());
            reader->ReadAllScalarsOn();
            reader->Update();
            vtkUnstructuredGrid* output = reader->GetOutput();

            // Initially, populate results array with zeroes (nodes on the grid which are part of the FE mesh will have their
            // corresponding values in the results array overwritten later)
            double* results_nodewise = new double[(grid.x + 1) * (grid.y + 1)]; // Nodes grid has +1 width along each dimension
            populate_with_zeroes(results_nodewise, grid.x + 1, grid.y + 1);
            populate_with_zeroes(results_cellwise, grid.x, grid.y);

            // Get point data (this object contains the physics data)
            vtkPointData* point_data = output->GetPointData();

            // Obtain Von Mises stress array
            vtkDoubleArray* results_array = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars(data_type));

            // Overwrite grid values with values from results array (only for nodes with coordinates that lie within the FE mesh)
            vtkPoints* points = output->GetPoints();
            double* point = new double[3];
            vector<int> coords = {};
            Vector2d inv_cell_size = Vector2d(1.0 / grid.cell_size(0), 1.0 / grid.cell_size(1));
            for (int i = 0; i < points->GetNumberOfPoints(); i++) {
                point = points->GetData()->GetTuple(i);
                Vector2d origin_aligned_coord = Vector2d(point[0], point[1]) - offset;
                Vector2d gridscale_coord = inv_cell_size.cwiseProduct(origin_aligned_coord);
                int coord = (int)(gridscale_coord[0] * grid.y + gridscale_coord[1]);
                int x = coord / grid.y;
                int y = coord % grid.y;
                coords.push_back(coord);
                results_nodewise[coord] = (double)results_array->GetValue(i);
            }

            // Create cellwise results distribution by averaging all groups of 4 corners of a cell
            for (int i = 0; i < grid.size2d; i++) {
                int x = i / grid.y;
                int y = i % grid.y;
                double neighbors_sum = (
                    results_nodewise[i] + results_nodewise[(x + 1) * grid.x + y] + results_nodewise[(x + 1) * grid.x + (y + 1)]
                    + results_nodewise[x * grid.x + y + 1]
                );
                results_cellwise[i] = neighbors_sum / 4.0;
            }
        }
    };
}



