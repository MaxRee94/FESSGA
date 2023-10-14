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
#include <vtkUnstructuredGridWriter.h>
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
            FEACase(string _path, int _dim_x, int _dim_y, double _max_stress_threshold) :
                path(_path), dim_x(_dim_x), dim_y(_dim_y), max_stress_threshold(_max_stress_threshold) {
            }
            Vector2d compute_cell_barycenter(vector<int>* cells);
            void compute_node_barycenters();
            void compute_cell_barycenters();
            void compute_barycenters();
            Vector2d get_node_coords(int idx);

            string path;
            string content;
            string name;
            map<string, vector<int>> bound_cond_nodes;
            map<string, vector<pair<int, int>>> bound_cond_lines;
            map<string, map<string, vector<int>>> bound_cond_cells;
            map<string, Vector2d> node_barycenters;
            map<string, map<string, Vector2d>> cell_barycenters;
            vector<string> names, sections;
            int dim_x, dim_y;
            double max_stress_threshold;
        };
        
        class FEACaseManager {
        public:
            FEACaseManager() = default;
            FEACaseManager(
                vector<phys::FEACase> _sources, vector<phys::FEACase> _targets, double _max_stress_threshold,
                int _dim_x, int _dim_y
            ) {
                sources = _sources;
                targets = _targets;
                active_cases = { _sources };
                dim_x = _dim_x;
                dim_y = _dim_y;
                dynamic = true;
            }
            FEACaseManager(
                vector<phys::FEACase> _active_cases, double _max_stress_threshold, int _dim_x, int _dim_y
            ) {
                active_cases = _active_cases;
                max_stress_threshold = _max_stress_threshold;
                dim_x = _dim_x;
                dim_y = _dim_y;
                dynamic = false;
            }
            Vector2d get_node_coords(int idx);
            Vector2d get_vector2nn(Vector2d node, vector<int>& target_nodes);
            void order_nodes_by_distance(vector<int>* source_nodes, Vector2d barycenter);
            void resample_boundary(vector<int>& nodes, Vector2d barycenter, int target_no_lines);
            void compute_migration_vectors(phys::FEACase* source, phys::FEACase* target, int pair_idx);
            void initialize();
            void interpolate_nodes(FEACase* active_case, FEACase* source, float fraction, int pair_idx);
            void interpolate_cells(FEACase* source, FEACase* target, FEACase* active_case, float fraction);
            void interpolate(float fraction);
            tuple<int, int> get_nearest_neighbor(Vector2d node, vector<int>* neighbors = 0);
            int get_direct_unvisited_node_neighbor(int node_idx, vector<int>* all_nodes, vector<int>* visited_nodes);
            vector<int> get_indirect_path_to_unvisited_node_neighbor(int node_idx, vector<int>* all_nodes, vector<int>* visited_nodes, int radius);
            vector<int> get_path_to_unvisited_node_neighbor(int node_idx, vector<int>* all_nodes, vector<int>* visited_nodes);
            vector<int> get_unvisited_neighbor_cells(int cell, vector<int>& neighbors, vector<int>* visited_cells);
            pair<int, int> get_cells_with_given_edge(pair<int, int> edge);
            tuple<int, int> get_bound_cells(int start_cell, pair<int, int> neighbor_cells, pair<int, int> prev_line, pair<int, int> next_line);
            void walk_and_collect_bound_nodes(vector<int>* unordered_boundary, vector<int>& walk, int start_node = -1);
            void walk_and_collect_bound_cells(
                FEACase* active_case, map<string, vector<int>>* bound_cells, vector<int>* visited_nodes, pair<int, int> first_line,
                int first_bound_cell, int neighbor_node, string bound_name
            );
            vector<int> get_additional_bound_cells(vector<int>* origin_cells, vector<int>* cells_to_avoid);
            void update_casepaths(string case_folder);

            vector<phys::FEACase> sources, targets, active_cases;
            vector<map<string, vector<Vector2d>>> migration_vectors;
            vector<int> keep_cells = {};
            vector<int> cutout_cells = {};
            vector<int> inactive_cells = {};
            double max_stress_threshold = INFINITY;
            bool maintain_boundary_connection = true;
            int dim_x, dim_y;
            bool dynamic = false;
            string stress_type = "tensile";
        };

        class FEAResults2D {
        public:
            FEAResults2D() = default;
            FEAResults2D(int dim_x, int dim_y) { x = dim_x, y = dim_y; }
            PairSet data;
            map<int, double> data_map;
            int x, y;
            string type;
            double min = INFINITY;
            double max = 0;
        };

        static void start_external_process(
            string base_folder, vector<FILE*>* pipes = 0, bool wait = true, bool verbose = false
        ) {
            std::string command = base_folder + "/run_elmer.bat";
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

        static void call_elmer(
            string case_folder, phys::FEACaseManager* fea_casemanager,
            vector<FILE*>* pipes = 0, bool wait = true, bool verbose = false
        ) {
            for (auto& fea_case : fea_casemanager->active_cases) {
                IO::write_text_to_file(fea_case.name + ".sif\n1", case_folder + "/ELMERSOLVER_STARTINFO");
                start_external_process(case_folder, pipes, wait, verbose);
            }
        }

        static bool load_2d_physics_data(
            vector<string> filenames, FEAResults2D& results, FEACaseManager* fea_casemanager, int dim_x, int dim_y, Vector2d cell_size,
            Vector3d _offset, char* data_type)
        {
            // For each coordinate, retain the maximum value out of all FEA runs
            // We thereby obtain a superposition of stress distributions.
            for (auto& filename : filenames) {
                FEAResults2D single_run_results;
                load_single_VTK_file(filename, &single_run_results, fea_casemanager, dim_x, dim_y, cell_size, _offset, data_type);
                for (auto& [coord, stress] : single_run_results.data_map) {
                    if (single_run_results.data_map[coord] > results.data_map[coord]) {
                        results.data_map[coord] = single_run_results.data_map[coord];
                    }
                }
                if (single_run_results.max > results.max) results.max = single_run_results.max;
                if (single_run_results.min < results.min) results.min = single_run_results.min;
            }
            help::sort(results.data_map, results.data);

            return true;
        }

        // Return the indices of the filled cell neighbors of the given node coordinates
        static vector<int> get_neighbors(int x, int y, int dim_x, int dim_y, uint* _values) {
            vector<pair<int, int>> offsets = { pair(0,0), pair(-1,0), pair(-1, -1), pair(0, -1) };
            vector<int> cell_neighbors;
            for (auto& offset : offsets) {
                int _x = x + offset.first;
                int _y = y + offset.second;
                if (_x == dim_x || _y == dim_y || _x < 0 || _y < 0) continue;
                int neighbor_coord = _x * dim_y + _y;
                if (_values[neighbor_coord] == 1) cell_neighbors.push_back(neighbor_coord);
            }
            return cell_neighbors;
        }

        static void load_nodewise_results(
            string filename, double* results_nodewise, int dim_x, int dim_y, Vector2d cell_size,
            Vector2d offset, string stress_type
        ) {
            // Read data from file
            vtkUnstructuredGridReader* reader = vtkUnstructuredGridReader::New();
            reader->SetFileName(filename.c_str());
            reader->ReadAllScalarsOn();
            reader->Update();
            vtkUnstructuredGrid* output = reader->GetOutput();

            // Get point data (this object contains the physics data)
            vtkPointData* point_data = output->GetPointData();

            // Obtain stress array(s)
            //vtkDoubleArray* vonmises = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars("Vonmises"));
            //vtkDoubleArray* stress_xx = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars("Stress_xx"));
            vtkDoubleArray* displacements = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars("Displacement"));
            if (displacements->GetSize() == 0) return; // If the array is empty, there is no physics data to load.

            // Overwrite grid values with values from results array (only for nodes with coordinates that lie within
            // the FE mesh)
            vtkPoints* points = output->GetPoints();
            double point[2];
            vector<int> coords;
            Vector2d inv_cell_size = Vector2d(1.0 / cell_size(0), 1.0 / cell_size(1));
            for (int i = 0; i < points->GetNumberOfPoints(); i++) {
                point[0] = points->GetData()->GetTuple(i)[0];
                point[1] = points->GetData()->GetTuple(i)[1];
                Vector2d origin_aligned_coord = Vector2d(point[0], point[1]) - offset;
                Vector2d gridscale_coord = inv_cell_size.cwiseProduct(origin_aligned_coord);
                int coord = (round(gridscale_coord[0]) * (dim_y + 1) + round(gridscale_coord[1]));
                coords.push_back(coord);
                double displacement_x = (double)displacements->GetValue(i);
                double displacement_y = (double)displacements->GetValue(i + 1);
                double magnitude = Vector2d(displacement_x, displacement_y).norm();
                magnitude = max(magnitude, 0.0);
                results_nodewise[coord] = magnitude;
            }
            reader->Delete();
        }

        static void write_results_superposition(vector<string> vtk_paths, int dim_x, int dim_y, Vector2d cell_size, Vector3d _offset, string outfile, string stress_type)
        {
            // Initially, populate results array with zeroes (nodes on the grid which are part of the FE mesh will
            // have their corresponding values in the results array overwritten later)
            double* nodewise_results = new double[(dim_x + 1) * (dim_y + 1)];
            help::populate_with_zeroes(nodewise_results, dim_x + 1, dim_y + 1);
            Vector2d offset = Vector2d(_offset(0), _offset(1));
            for (auto& casepath : vtk_paths) {
                double* single_run_results = new double[(dim_x + 1) * (dim_y + 1)]; // Nodes grid has +1 width along each dim
                help::populate_with_zeroes(single_run_results, dim_x + 1, dim_y + 1);
                phys::load_nodewise_results(casepath, single_run_results, dim_x, dim_y, cell_size, offset, stress_type);
                for (int i = 0; i < (dim_x + 1) * (dim_y + 1); i++) {
                    nodewise_results[i] = max(nodewise_results[i], single_run_results[i]);
                }
                delete[] single_run_results;
            }

            // Read data from file
            vtkUnstructuredGridReader* reader = vtkUnstructuredGridReader::New();
            reader->SetFileName(vtk_paths[0].c_str());
            reader->ReadAllScalarsOn();
            reader->Update();
            vtkUnstructuredGrid* output = reader->GetOutput();

            // Get point data (this object contains the physics data)
            vtkPointData* point_data = output->GetPointData();

            // Obtain Von Mises stress array
            //vtkDoubleArray* results_array = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars(data_type));
            vtkDoubleArray* displacements = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars("Displacement"));
            if (displacements->GetSize() == 0) return; // If the array is empty, there is no physics data to load.

            // Overwrite grid values with values from results array (only for nodes with coordinates that lie within
            // the FE mesh)
            vtkPoints* points = output->GetPoints();
            double point[2];
            vector<int> coords;
            Vector2d inv_cell_size = Vector2d(1.0 / cell_size(0), 1.0 / cell_size(1));
            for (int i = 0; i < points->GetNumberOfPoints(); i++) {
                point[0] = points->GetData()->GetTuple(i)[0];
                point[1] = points->GetData()->GetTuple(i)[1];
                Vector2d origin_aligned_coord = Vector2d(point[0], point[1]) - offset;
                Vector2d gridscale_coord = inv_cell_size.cwiseProduct(origin_aligned_coord);
                int coord = (round(gridscale_coord[0]) * (dim_y + 1) + round(gridscale_coord[1]));
                int x = coord / (dim_y + 1);
                int y = coord % (dim_y + 1);
                coords.push_back(coord);
                displacements->SetValue(i, nodewise_results[coord]);
            }

            // Write to vtk OUTFILE
            vtkUnstructuredGridWriter* writer = vtkUnstructuredGridWriter::New();
            writer->SetFileName(outfile.c_str());
            writer->SetFileTypeToASCII();
            writer->SetInputData(reader->GetOutput());
            writer->Update();

            // Delete writer and reader
            delete[] nodewise_results;
            writer->Delete();

            return;
        }

        static bool load_single_VTK_file(
            string filename, FEAResults2D* results, FEACaseManager* fea_casemanager, int dim_x, int dim_y, Vector2d cell_size,
            Vector3d _offset, char* data_type
        ) {
            // Read data from file
            vtkUnstructuredGridReader* reader = vtkUnstructuredGridReader::New();
            reader->SetFileName(filename.c_str());
            reader->ReadAllScalarsOn();
            reader->Update();
            vtkUnstructuredGrid* output = reader->GetOutput();

            // Initially, populate results array with zeroes (nodes on the grid which are part of the FE mesh will
            // have their corresponding values in the results array overwritten later)
            double* results_nodewise = new double[(dim_x + 1) * (dim_y + 1)]; // Nodes grid has +1 width along each dim
            help::populate_with_zeroes(results_nodewise, dim_x + 1, dim_y + 1);

            // Get point data (this object contains the physics data)
            vtkPointData* point_data = output->GetPointData();

            // Obtain Von Mises stress array
            //vtkDoubleArray* results_array = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars(data_type));
            //vtkDoubleArray* stress_xx = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars("Stress_xx"));
            //vtkDoubleArray* displacements = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars("Stress_yy"));
            vtkDoubleArray* displacements = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars("Stress_yy"));
            if (displacements->GetSize() == 0) return false; // If the array is empty, there is no physics data to load.

            // Overwrite grid values with values from results array (only for nodes with coordinates that lie within
            // the FE mesh)
            vtkPoints* points = output->GetPoints();
            double point[2];
            vector<int> coords;
            Vector2d inv_cell_size = Vector2d(1.0 / cell_size(0), 1.0 / cell_size(1));
            Vector2d offset = Vector2d(_offset(0), _offset(1));
            for (int i = 0; i < points->GetNumberOfPoints(); i+=3) {
                point[0] = points->GetData()->GetTuple(i/3)[0];
                point[1] = points->GetData()->GetTuple(i/3)[1];
                Vector2d origin_aligned_coord = Vector2d(point[0], point[1]) - offset;
                Vector2d gridscale_coord = inv_cell_size.cwiseProduct(origin_aligned_coord);
                int coord = (round(gridscale_coord[0]) * (dim_y + 1) + round(gridscale_coord[1]));
                coords.push_back(coord);
                double displacement_x = (double)displacements->GetValue(i);
                double displacement_y = (double)displacements->GetValue(i+1);
                double magnitude = Vector2d(displacement_x, displacement_y).norm();
                magnitude = max(magnitude, 0.0);
                //cout << "displacement x: " << displacement_x << endl;
                //cout << "magnitude: " << magnitude << endl;
                results_nodewise[coord] = magnitude;
            }

            // Create cellwise results distribution by taking maximum of each group of 4 corners of a cell
            double min_stress = 1e30;
            double max_stress = 0;
            for (int i = 0; i < coords.size(); i++) {
                int node_coord = coords[i];
                int x = node_coord / (dim_y + 1);
                int y = node_coord % (dim_y + 1);
                if (x == dim_x || y == dim_y) continue; // Skip coordinates outside cell domain
                int cell_coord = x * dim_y + y;
                if (help::is_in(&fea_casemanager->inactive_cells, cell_coord)) {
                    // Cells marked as 'inactive' are ignored during solution evaluation.
                    results->data_map.insert(pair(cell_coord, -9999));
                    continue;
                }
                Vector4d neighbors;
                neighbors[0] = results_nodewise[node_coord];
                neighbors[1] = results_nodewise[(x + 1) * (dim_y + 1) + y];
                neighbors[2] = results_nodewise[(x + 1) * (dim_y + 1) + (y + 1)];
                neighbors[3] = results_nodewise[x * (dim_y + 1) + (y + 1)];
                if (neighbors.minCoeff() == 0) continue; // Skip cells with corners that have stress value 0
                double cell_stress = neighbors.maxCoeff();
                if (cell_stress > max_stress) max_stress = cell_stress;
                if (cell_stress < min_stress) min_stress = cell_stress;
                //results->data_map.insert(pair(cell_coord, cell_stress));
                results->data_map.insert(pair(cell_coord, neighbors.mean()));
            }
            results->min = min_stress;
            results->max = max_stress;
            delete[] results_nodewise;
            reader->Delete();

            return true;
        }
    };
}



