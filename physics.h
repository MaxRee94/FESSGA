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
#include "helpers.h"

VTK_MODULE_INIT(vtkRenderingOpenGL2)
#define VTK_DEBUG_LEAKS true



namespace fessga {
    class physics {
    public:

        struct FEResults2D {
            FEResults2D(mesher::Grid3D grid) { x = grid.x; y = grid.y; }
            PairSet data;
            map<int, double> data_map;
            int x, y;
            string type;
            double min, max;
        };

        static void call_elmer(string bat_file) {
            std::string command = bat_file;
            std::array<char, 80> buffer;
            FILE* pipe = _popen(command.c_str(), "r");
            while (fgets(buffer.data(), 80, pipe) != NULL) {
                //std::cout << buffer.data();
            }
            _pclose(pipe);
        }

        static int remove_low_stress_cells(
            PairSet* fe_results, uint* densities, mesher::Case* fe_case, mesher::Grid3D grid, int no_cells_to_remove,
            vector<int>& removed_cells, double max_stress_threshold, mesher::Piece* smaller_piece = 0, int total_no_cells = -1
        ) {
            cout << "total no cells (in physics): " << total_no_cells << endl;

            int count = 0;
            int cell_from_smaller_piece;
            for (auto& item : (*fe_results)) {
                int cell_coord = item.first;
                double cell_stress = item.second;

                // If the cell's stress exceeds the maximum, break the loop (all subsequent cells will also exceed the maximum since the list is ordered).
                if (cell_stress > max_stress_threshold) continue;

                // If the cell has a line on which a boundary condition was applied, skip deletion
                if (help::is_in(&fe_case->boundary_cells, cell_coord)) {
                    continue;
                }

                // If the cell was previously whitelisted, skip deletion
                if (help::is_in(&fe_case->whitelisted_cells, cell_coord)) {
                    continue;
                }

                // If the cell was already empty, skip deletion (more importantly: don't count this as a deletion)
                if (!densities[cell_coord]) continue; 

                // If a 'smaller piece' vector was provided, perform cell removal in 'careful mode'. This means: check whether cell deletion results in 
                // multiple pieces. If so, undo the deletion and whitelist the cell.
                if (smaller_piece != 0) {
                    vector<int> _removed_cell = { cell_coord };
                    cout << "\n total no cells: " << total_no_cells - count << endl;
                    cout << "is single piece before element removal: " << mesher::is_single_piece(
                        densities, grid, fe_case, total_no_cells - count, &_removed_cell, smaller_piece->cells[0]) << endl;
                    densities[cell_coord] = 0;
                    cell_from_smaller_piece = smaller_piece->cells[0];
                    bool is_single_piece = mesher::is_single_piece(densities, grid, fe_case, total_no_cells - count, &_removed_cell, cell_from_smaller_piece);
                    cout << "is single piece? " << is_single_piece << endl;
                    if (!is_single_piece) {
                        // Cell removal resulted in multiple pieces. Restore cell and whitelist it.
                        densities[cell_coord] = 1;
                        cout << "restored cell." << endl;
                        fe_case->whitelisted_cells.push_back(cell_coord);
                        continue;
                    }
                    count++;
                }

                // If cell deletion leads to infeasibility, skip deletion
                int no_deleted_neighbors = 0;
                if (!mesher::cell_is_safe_to_delete(densities, grid, cell_coord, no_deleted_neighbors, &fe_case->boundary_cells)) {
                    densities[cell_coord] = 1;
                    cout << "cell is not safe to delete" << endl;
                    continue;
                }
                count += no_deleted_neighbors;

                if (smaller_piece != 0) {
                    vector<int> _removed_cell = { cell_coord };
                    cout << "total no cells: " << total_no_cells - count << endl;
                    cout << "is single piece after neighbor checking / deletion: " << mesher::is_single_piece(
                        densities, grid, fe_case, total_no_cells - count, &_removed_cell, smaller_piece->cells[0]) << endl;
                }

                // Set cell to zero, making it empty
                densities[cell_coord] = 0;
                if (!smaller_piece) count++;
                removed_cells.push_back(cell_coord);
                if (count >= no_cells_to_remove) break;
            }
            return count;
        }

        // Remove a floating piece (if possible)
        static bool remove_floating_piece(
            uint* densities, mesher::Grid3D grid, mesher::Piece* piece, double max_stress_threshold, mesher::Case* fe_case, FEResults2D* fe_results
        ) {
            vector<int> removed_cells;
            for (auto& cell : piece->cells) {
                if (fe_results->data_map[cell] > max_stress_threshold) {
                    // If cell's stress exceeds maximum, restore the already deleted cells and exit this function.
                    mesher::restore_removed_cells(densities, grid, &removed_cells);
                    piece->is_removable = false;
                    return false;
                }
                removed_cells.push_back(cell);
                densities[cell] = 0;
            }
            return true;
        }

        // Remove all pieces smaller than the largest piece from the mesh, if possible
        static vector<int> remove_smaller_pieces(
            uint* densities, mesher::Grid3D grid, mesher::Case* fe_case, int total_no_cells, vector<mesher::Piece>* pieces, vector<int>* removed_cells,
            double max_stress_threshold, FEResults2D* fe_results
        ) {
            vector<int> unremoved_piece_indices;
            int size_largest_piece;
            mesher::remove_largest_piece(pieces, size_largest_piece);
            for (int i = 0; i < pieces->size(); i++) {
                mesher::Piece piece = pieces->at(i);
                if (!piece.is_removable) { unremoved_piece_indices.push_back(i); continue; }
                bool piece_removed = remove_floating_piece(densities, grid, &piece, max_stress_threshold, fe_case, fe_results);
                if (!piece_removed) unremoved_piece_indices.push_back(i);
            }

            return unremoved_piece_indices;
        }

        static bool load_2d_physics_data(
            string filename, FEResults2D& results, mesher::Grid3D grid, Vector3d _offset, char* data_type)
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
            help::populate_with_zeroes(results_nodewise, grid.x + 1, grid.y + 1);

            // Get point data (this object contains the physics data)
            vtkPointData* point_data = output->GetPointData();

            // Obtain Von Mises stress array
            vtkDoubleArray* results_array = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars(data_type));
            if (results_array->GetSize() == 0) return false; // If the array is empty, there is no physics data to load.

            // Overwrite grid values with values from results array (only for nodes with coordinates that lie within the FE mesh)
            vtkPoints* points = output->GetPoints();
            double* point = new double[3];
            vector<int> coords = {};
            Vector2d inv_cell_size = Vector2d(1.0 / grid.cell_size(0), 1.0 / grid.cell_size(1));
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
                int coord = (round(gridscale_coord[0]) * (grid.y + 1) + round(gridscale_coord[1]));
                int x = coord / (grid.y + 1);
                int y = coord % (grid.y + 1);
                coords.push_back(coord);
                results_nodewise[coord] = (double)results_array->GetValue(i);
            }

            // Create density distribution indicating which cells have nonzero physics values (FOR DEBUGGING PURPOSES ONLY)
            /*uint* nonzero_physics = new uint[(grid.x) * (grid.y)];
            help::populate_with_zeroes(nonzero_physics, grid.x, grid.y);*/

            // Create cellwise results distribution by averaging all groups of 4 corners of a cell
            double min_stress = 1e30;
            double max_stress = 0;
            for (int i = 0; i < coords.size(); i++) {
                int coord = coords[i];
                int x = coord / (grid.y + 1);
                int y = coord % (grid.y + 1);
                if (x == grid.x || y == grid.y) continue; // Skip coordinates outside cell domain
                Vector4d neighbors;
                neighbors[0] = results_nodewise[coord];
                neighbors[1] = results_nodewise[(x + 1) * (grid.y + 1) + y];
                neighbors[2] = results_nodewise[(x + 1) * (grid.y + 1) + (y + 1)];
                neighbors[3] = results_nodewise[x * (grid.y + 1) + (y + 1)];
                if (neighbors.minCoeff() == 0) continue; // Skip cells with corners that have stress value 0
                //double cell_stress = neighbors.sum() / 4.0;
                double cell_stress = neighbors.maxCoeff();
                if (cell_stress > max_stress) max_stress = cell_stress;
                if (cell_stress < min_stress) min_stress = cell_stress;
                int cell_coord = x * grid.y + y;
                results.data_map[cell_coord] = cell_stress;
            }
            help::sort(results.data_map, results.data);
            results.min = min_stress;
            results.max = max_stress;

            //cout << "\nNonzero physics values: " << endl;
            //mesher::print_density_distrib(nonzero_physics, grid.x, grid.y);

            return true;
        }
    };
}



