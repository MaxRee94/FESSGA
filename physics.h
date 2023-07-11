#pragma once
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
#include "helpers.h"
#include "io.h"

//VTK_MODULE_INIT(vtkRenderingOpenGL2)
//#define VTK_DEBUG_LEAKS true


using namespace Eigen;


namespace fessga {

    class phys {
    public:

        class FEACase {
        public:
            FEACase() = default;
            FEACase(string _path) : path(_path) {}
            string path;
            vector<string> names, sections;
            vector<int> cells_to_keep;
            vector<int> cutout_cells;
            string content;
            double max_stress_threshold = INFINITY;
            bool maintain_boundary_connection = true;
        };

        class FEAResults2D {
        public:
            FEAResults2D() = default;
            FEAResults2D(int dim_x, int dim_y) { x = dim_x, y = dim_y; }
            PairSet data;
            map<int, double> data_map;
            int x, y;
            string type;
            double min, max;
        };

        static void call_elmer(string bat_file, vector<FILE*>* pipes = 0, bool wait = true, bool verbose = false) {
            std::string command = bat_file;
            if (wait) {
                std::array<char, 80> buffer;
                FILE* pipe = _popen(command.c_str(), "r");
                while (fgets(buffer.data(), 80, pipe) != NULL) {
                    if (verbose) std::cout << buffer.data();
                }
                _pclose(pipe);
            }
            else {
                // If a pipe is provided as a parameter, don't wait for it to finish
                FILE* pipe = _popen(command.c_str(), "r");
                pipes->push_back(pipe);
            }
        }

        static bool load_2d_physics_data(
            string filename, FEAResults2D* results, int dim_x, int dim_y, Vector2d cell_size, Vector3d _offset, char* data_type)
        {
            // Read data from file
            vtkUnstructuredGridReader* reader = vtkUnstructuredGridReader::New();
            reader->SetFileName(filename.c_str());
            reader->ReadAllScalarsOn();
            reader->Update();
            vtkUnstructuredGrid* output = reader->GetOutput();

            // Initially, populate results array with zeroes (nodes on the grid which are part of the FE mesh will have their
            // corresponding values in the results array overwritten later)
            double* results_nodewise = new double[(dim_x + 1) * (dim_y + 1)]; // Nodes grid has +1 width along each dimension
            help::populate_with_zeroes(results_nodewise, dim_x + 1, dim_y + 1);

            // Get point data (this object contains the physics data)
            vtkPointData* point_data = output->GetPointData();

            // Obtain Von Mises stress array
            vtkDoubleArray* results_array = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars(data_type));
            if (results_array->GetSize() == 0) return false; // If the array is empty, there is no physics data to load.

            // Overwrite grid values with values from results array (only for nodes with coordinates that lie within the FE mesh)
            vtkPoints* points = output->GetPoints();
            double* point = new double[3];
            vector<int> coords = {};
            Vector2d inv_cell_size = Vector2d(1.0 / cell_size(0), 1.0 / cell_size(1));
            Vector2d offset = Vector2d(_offset(0), _offset(1));
            //cout << "inv cell size: " << inv_cell_size.transpose() << endl;
            //cout << "offset: " << offset.transpose() << endl;
            for (int i = 0; i < points->GetNumberOfPoints(); i++) {
                point = points->GetData()->GetTuple(i);
                //cout << "\nvtk coord: " << point[0] << ", " << point[1] << endl;
                Vector2d origin_aligned_coord = Vector2d(point[0], point[1]) - offset;
                //cout << "origin aligned coord: " << origin_aligned_coord.transpose() << endl;
                Vector2d gridscale_coord = inv_cell_size.cwiseProduct(origin_aligned_coord);
                //cout << "gridscale coord: " << gridscale_coord.transpose() << endl;
                int coord = (round(gridscale_coord[0]) * (dim_y + 1) + round(gridscale_coord[1]));
                int x = coord / (dim_y + 1);
                int y = coord % (dim_y + 1);
                coords.push_back(coord);
                results_nodewise[coord] = (double)results_array->GetValue(i);
            }

            // Create density distribution indicating which cells have nonzero physics values (FOR DEBUGGING PURPOSES ONLY)
            /*uint* nonzero_physics = new uint[(dim_x) * (dim_y)];
            help::populate_with_zeroes(nonzero_physics, dim_x, dim_y);*/

            // Create cellwise results distribution by averaging all groups of 4 corners of a cell
            double min_stress = 1e30;
            double max_stress = 0;
            for (int i = 0; i < coords.size(); i++) {
                int coord = coords[i];
                int x = coord / (dim_y + 1);
                int y = coord % (dim_y + 1);
                if (x == dim_x || y == dim_y) continue; // Skip coordinates outside cell domain
                Vector4d neighbors;
                neighbors[0] = results_nodewise[coord];
                neighbors[1] = results_nodewise[(x + 1) * (dim_y + 1) + y];
                neighbors[2] = results_nodewise[(x + 1) * (dim_y + 1) + (y + 1)];
                neighbors[3] = results_nodewise[x * (dim_y + 1) + (y + 1)];
                if (neighbors.minCoeff() == 0) continue; // Skip cells with corners that have stress value 0
                //double cell_stress = neighbors.sum() / 4.0;
                double cell_stress = neighbors.maxCoeff();
                if (cell_stress > max_stress) max_stress = cell_stress;
                if (cell_stress < min_stress) min_stress = cell_stress;
                int cell_coord = x * dim_y + y;
                results->data_map.insert(pair(cell_coord, cell_stress));
            }
            help::sort(results->data_map, results->data);
            results->min = min_stress;
            results->max = max_stress;

            //cout << "\nNonzero physics values: " << endl;

            return true;
        }
    };
}



