#pragma once
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <filesystem>
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
#include "timer.h"


//VTK_MODULE_INIT(vtkRenderingOpenGL2)
//#define VTK_DEBUG_LEAKS true

namespace fs = std::filesystem;
using namespace Eigen;


namespace fessga {

    class phys {
    public:

        class FEACase {
        public:
            FEACase() = default;
            FEACase(string _path, int _dim_x, int _dim_y, double _mechanical_threshold) :
                path(_path), dim_x(_dim_x), dim_y(_dim_y), mechanical_threshold(_mechanical_threshold) {
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
            double mechanical_threshold, max_tensile_strength, max_compressive_strength;
        };
        
        class FEACaseManager {
        public:
            FEACaseManager() = default;
            FEACaseManager(
                vector<phys::FEACase> _sources, vector<phys::FEACase> _targets, double _mechanical_threshold,
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
                vector<phys::FEACase> _active_cases, double _mechanical_threshold, int _dim_x, int _dim_y
            ) {
                active_cases = _active_cases;
                mechanical_threshold = _mechanical_threshold;
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
            vector<int> border_nodes = {};
            double mechanical_threshold = INFINITY;
            double max_displacement;
            double max_tensile_strength;
            double max_compressive_strength;
            bool maintain_boundary_connection = true;
            int dim_x, dim_y;
            bool dynamic = false;
            string mechanical_constraint = "";
            int displacement_measurement_cell = -1;
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
            double max_displacement = 0;
            double max_yield_criterion = 0;
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
            Vector3d _offset, string mechanical_constraint,  vector<int>* border_nodes, bool verbose = false, vector<int>* times = 0)
        {
            // For each coordinate, retain the maximum value out of all FEA runs
            // We thereby obtain a superposition of stress distributions.
            for (auto& filename : filenames) {
                FEAResults2D single_run_results;
                load_single_VTK_file(filename, &single_run_results, fea_casemanager, dim_x, dim_y, cell_size, _offset, mechanical_constraint, border_nodes, times, verbose);
                for (auto& [coord, stress] : single_run_results.data_map) {
                    if (single_run_results.data_map[coord] > results.data_map[coord]) {
                        results.data_map[coord] = single_run_results.data_map[coord];
                    }
                }
                if (single_run_results.max > results.max) results.max = single_run_results.max;
                if (single_run_results.min < results.min) results.min = single_run_results.min;
                if (single_run_results.max_displacement > results.max_displacement) results.max_displacement = single_run_results.max_displacement;
                if (single_run_results.max_yield_criterion > results.max_yield_criterion) results.max_yield_criterion = single_run_results.max_yield_criterion;
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

        // Load node data of a specific cell in a vtk channel
        static void load_nodewise_results(
            vtkUnstructuredGrid* vtk_content, vector<double>* corners, int dim_x, int dim_y, Vector2d cell_size,
            Vector2d offset, string channel, map<int, int>* node_coords_map, vector<int>* corner_coords, double sign = 1
        ) {
            // Get point data (this object contains the physics data)
            vtkPointData* point_data = vtk_content->GetPointData();

            // Obtain results array
            vtkDoubleArray* results_array = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars(channel.c_str()));
            if (results_array->GetSize() == 0) return; // If the array is empty, there is no physics data to load.

            vtkPoints* points = vtk_content->GetPoints();
            for (int i = 0; i < 4; i++) {
                int point_idx = node_coords_map->at(corner_coords->at(i));

                if (channel == "Displacement") {
                    double displacement_x = (double)results_array->GetValue(point_idx * 3);
                    double displacement_y = (double)results_array->GetValue(point_idx * 3 + 1);
                    double magnitude = Vector2d(displacement_x, displacement_y).norm();
                    magnitude = max(magnitude, 0.0);
                    corners->push_back(magnitude);
                }
                else {
                    corners->push_back(sign * (double)results_array->GetValue(point_idx));
                }
            }
        }

        // Load all node data of a vtk channel
        static void load_nodewise_results(
            vtkUnstructuredGrid* vtk_content, double* results_nodewise, int dim_x, int dim_y, Vector2d cell_size,
            Vector2d offset, string channel, vector<int>* border_nodes, double sign, bool clamp, map<int, int>* _node_coords_map = 0
        ) {
            if (_node_coords_map) _node_coords_map->clear();

            // Get point data (this object contains the physics data)
            vtkPointData* point_data = vtk_content->GetPointData();

            // Obtain stress array(s)
            //vtkDoubleArray* vonmises = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars("Vonmises"));
            vtkDoubleArray* results_array = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars(channel.c_str()));
            if (results_array->GetSize() == 0) return; // If the array is empty, there is no physics data to load.

            // Overwrite grid values with values from results array (only for nodes with coordinates that lie within
            // the FE mesh)
            vtkPoints* points = vtk_content->GetPoints();
            double point[2];
            map<int, int> node_coords_map;
            Vector2d inv_cell_size = Vector2d(1.0 / cell_size(0), 1.0 / cell_size(1));
            for (int i = 0; i < points->GetNumberOfPoints(); i++) {
                // Compute node coordinate
                point[0] = points->GetData()->GetTuple(i)[0];
                point[1] = points->GetData()->GetTuple(i)[1];
                Vector2d origin_aligned_coord = Vector2d(point[0], point[1]) - offset;
                Vector2d gridscale_coord = inv_cell_size.cwiseProduct(origin_aligned_coord);
                int coord = (round(gridscale_coord[0]) * (dim_y + 1) + round(gridscale_coord[1]));
                
                // Obtain node value
                if (_node_coords_map != 0) (*_node_coords_map)[coord] = i;
                else node_coords_map[coord] = i;
                if (channel == "Displacement") {
                    double displacement_x = (double)results_array->GetValue(i * 3);
                    double displacement_y = (double)results_array->GetValue(i * 3 + 1);
                    double magnitude = Vector2d(displacement_x, displacement_y).norm();
                    magnitude = max(magnitude, 0.0);
                    results_nodewise[coord] = magnitude;
                }
                else {
                    results_nodewise[coord] = sign * (double)results_array->GetValue(i);
                    if (clamp) results_nodewise[coord] = max(0.0, results_nodewise[coord]);
                }
            }
        }

        static void write_results_superposition(
            vector<string> vtk_paths, int dim_x, int dim_y, Vector2d cell_size, Vector3d _offset, string outfile, string mechanical_constraint, vector<int>* border_nodes,
            double max_tensile_strength, double max_compressive_strength
        ) {
            // Initially, populate results array with zeroes (nodes on the grid which are part of the FE mesh will
            // have their corresponding values in the results array overwritten later)
            double* nodewise_tensile_xx = new double[(dim_x + 1) * (dim_y + 1)];
            double* nodewise_tensile_yy = new double[(dim_x + 1) * (dim_y + 1)];
            double* nodewise_compressive_xx = new double[(dim_x + 1) * (dim_y + 1)];
            double* nodewise_compressive_yy = new double[(dim_x + 1) * (dim_y + 1)];
            double* nodewise_xy = new double[(dim_x + 1) * (dim_y + 1)];
            double* nodewise_displacements = new double[(dim_x + 1) * (dim_y + 1)];
            double* nodewise_principle_compressive_stresses = new double[(dim_x + 1) * (dim_y + 1)];
            double* nodewise_principle_tensile_stresses = new double[(dim_x + 1) * (dim_y + 1)];
            double* nodewise_modified_mohr = new double[(dim_x + 1) * (dim_y + 1)];
            double* nodewise_vonmises = new double[(dim_x + 1) * (dim_y + 1)];
            help::populate_with_zeroes(nodewise_tensile_xx, dim_x + 1, dim_y + 1);
            help::populate_with_zeroes(nodewise_tensile_yy, dim_x + 1, dim_y + 1);
            help::populate_with_zeroes(nodewise_compressive_xx, dim_x + 1, dim_y + 1);
            help::populate_with_zeroes(nodewise_compressive_yy, dim_x + 1, dim_y + 1);
            help::populate_with_zeroes(nodewise_xy, dim_x + 1, dim_y + 1);
            help::populate_with_zeroes(nodewise_displacements, dim_x + 1, dim_y + 1);
            help::populate_with_zeroes(nodewise_principle_compressive_stresses, dim_x + 1, dim_y + 1);
            help::populate_with_zeroes(nodewise_principle_tensile_stresses, dim_x + 1, dim_y + 1);
            help::populate_with_zeroes(nodewise_modified_mohr, dim_x + 1, dim_y + 1);
            help::populate_with_zeroes(nodewise_vonmises, dim_x + 1, dim_y + 1);
            Vector2d offset = Vector2d(_offset(0), _offset(1));
            for (auto& casepath : vtk_paths) {
                // Read data from file
                vtkUnstructuredGridReader* reader = vtkUnstructuredGridReader::New();
                reader->SetFileName(casepath.c_str());
                reader->ReadAllScalarsOn();
                reader->Update();
                vtkUnstructuredGrid* output = reader->GetOutput();

                double* single_run_stress_xx = new double[(dim_x + 1) * (dim_y + 1)]; // Nodes grid has +1 width along each dim
                double* single_run_stress_yy = new double[(dim_x + 1) * (dim_y + 1)]; // Nodes grid has +1 width along each dim
                double* single_run_stress_xy = new double[(dim_x + 1) * (dim_y + 1)]; // Nodes grid has +1 width along each dim
                double* single_run_displacements = new double[(dim_x + 1) * (dim_y + 1)]; // Nodes grid has +1 width along each dim
                double* single_run_principal_compressive_stresses = new double[(dim_x + 1) * (dim_y + 1)]; // Nodes grid has +1 width along each dim
                double* single_run_principal_tensile_stresses = new double[(dim_x + 1) * (dim_y + 1)]; // Nodes grid has +1 width along each dim
                double* single_run_modified_mohr = new double[(dim_x + 1) * (dim_y + 1)]; // Nodes grid has +1 width along each dim
                double* single_run_vonmises = new double[(dim_x + 1) * (dim_y + 1)]; // Nodes grid has +1 width along each dim
                help::populate_with_zeroes(single_run_stress_xx, dim_x + 1, dim_y + 1);
                help::populate_with_zeroes(single_run_stress_yy, dim_x + 1, dim_y + 1);
                help::populate_with_zeroes(single_run_stress_xy, dim_x + 1, (dim_y + 1));
                help::populate_with_zeroes(single_run_displacements, dim_x + 1, (dim_y + 1));
                help::populate_with_zeroes(single_run_principal_compressive_stresses, dim_x + 1, (dim_y + 1));
                help::populate_with_zeroes(single_run_principal_tensile_stresses, dim_x + 1, (dim_y + 1));
                help::populate_with_zeroes(single_run_modified_mohr, dim_x + 1, (dim_y + 1));
                help::populate_with_zeroes(single_run_vonmises, dim_x + 1, (dim_y + 1));
                phys::load_nodewise_results(output, single_run_stress_xx, dim_x, dim_y, cell_size, offset, "Stress_xx", border_nodes, 1, false);
                phys::load_nodewise_results(output, single_run_stress_yy, dim_x, dim_y, cell_size, offset, "Stress_yy", border_nodes, 1, false);
                phys::load_nodewise_results(output, single_run_stress_xy, dim_x, dim_y, cell_size, offset, "Stress_xy", border_nodes, 1, false);
                phys::load_nodewise_results(output, single_run_displacements, dim_x, dim_y, cell_size, offset, "Displacement", border_nodes, 1, false);
                phys::load_nodewise_results(output, single_run_vonmises, dim_x, dim_y, cell_size, offset, "Vonmises", border_nodes, 1, false);
                compute_principal_stresses(single_run_stress_xx, single_run_stress_yy, single_run_stress_xy, single_run_principal_compressive_stresses, dim_x, dim_y, "Compressive");
                compute_principal_stresses(single_run_stress_xx, single_run_stress_yy, single_run_stress_xy, single_run_principal_tensile_stresses, dim_x, dim_y, "Tensile");
                for (int i = 0; i < (dim_x + 1) * (dim_y + 1); i++) {
                    single_run_modified_mohr[i] = get_modified_mohr_safety(
                        single_run_principal_tensile_stresses[i], single_run_principal_compressive_stresses[i], max_tensile_strength, max_compressive_strength
                    );
                    single_run_principal_compressive_stresses[i] = max(0.0, -single_run_principal_compressive_stresses[i]);
                    single_run_principal_tensile_stresses[i] = max(0.0, single_run_principal_tensile_stresses[i]);
                    nodewise_tensile_xx[i] = max(nodewise_tensile_xx[i], single_run_stress_xx[i]);
                    nodewise_tensile_yy[i] = max(nodewise_tensile_yy[i], single_run_stress_yy[i]);
                    nodewise_compressive_xx[i] = min(nodewise_compressive_xx[i], single_run_stress_xx[i]);
                    nodewise_compressive_yy[i] = min(nodewise_compressive_yy[i], single_run_stress_yy[i]);
                    nodewise_displacements[i] = max(nodewise_displacements[i], single_run_displacements[i]);
                    nodewise_principle_compressive_stresses[i] = max(nodewise_principle_compressive_stresses[i], single_run_principal_compressive_stresses[i]);
                    nodewise_principle_tensile_stresses[i] = max(nodewise_principle_tensile_stresses[i], single_run_principal_tensile_stresses[i]);
                    nodewise_modified_mohr[i] = max(nodewise_modified_mohr[i], single_run_modified_mohr[i]);
                    nodewise_vonmises[i] = max(nodewise_vonmises[i], single_run_vonmises[i]);
                }

                if (help::is_in(mechanical_constraint, "Mohr")) {
                    // Compute and write ModifiedMohr results to casefile
                    vtkPoints* points = output->GetPoints();
                    vtkDoubleArray* modified_mohr = dynamic_cast<vtkDoubleArray*>(reader->GetOutput()->GetPointData()->GetScalars("Stress_zz"));
                    double point[2];
                    vector<int> coords;
                    Vector2d inv_cell_size = Vector2d(1.0 / cell_size(0), 1.0 / cell_size(1));
                    for (int i = 0; i < points->GetNumberOfPoints(); i++) {
                        point[0] = points->GetData()->GetTuple(i)[0];
                        point[1] = points->GetData()->GetTuple(i)[1];
                        Vector2d origin_aligned_coord = Vector2d(point[0], point[1]) - offset;
                        Vector2d gridscale_coord = inv_cell_size.cwiseProduct(origin_aligned_coord);
                        int coord = (round(gridscale_coord[0]) * (dim_y + 1) + round(gridscale_coord[1]));
                        modified_mohr->SetValue(i, single_run_modified_mohr[coord]);
                    }
                    vtkUnstructuredGridWriter* writer = vtkUnstructuredGridWriter::New();
                    string current_path_base = fs::path(casepath).parent_path().string();
                    casepath.replace(0, current_path_base.size(), fs::path(outfile).parent_path().string());
                    writer->SetFileName(casepath.c_str());
                    writer->SetFileTypeToASCII();
                    writer->SetInputData(reader->GetOutput());
                    writer->Update();
                    writer->Delete();
                }

                delete[] single_run_stress_xx;
                delete[] single_run_stress_yy;
                delete[] single_run_stress_xy;
                delete[] single_run_displacements;
                delete[] single_run_principal_compressive_stresses;
                delete[] single_run_principal_tensile_stresses;
                delete[] single_run_modified_mohr;
                delete[] single_run_vonmises;
                reader->Delete();
            }

            // Read data from file
            vtkUnstructuredGridReader* reader = vtkUnstructuredGridReader::New();
            reader->SetFileName(vtk_paths[0].c_str());
            reader->ReadAllScalarsOn();
            reader->Update();
            vtkUnstructuredGrid* output = reader->GetOutput();

            // Get point data (this object contains the physics data)
            vtkPointData* point_data = output->GetPointData();

            // Obtain value arrays to write to
            vtkDoubleArray* stress_xx = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars("Stress_xx"));
            vtkDoubleArray* stress_yy = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars("Stress_yy"));
            vtkDoubleArray* displacements = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars("Displacement"));
            vtkDoubleArray* principal_compressive_stresses = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars("Stress_xz"));
            vtkDoubleArray* principal_tensile_stresses = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars("Stress_yz"));
            vtkDoubleArray* modified_mohr = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars("Stress_zz"));
            vtkDoubleArray* vonmises = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars("Vonmises"));
            if (displacements->GetSize() == 0) return; // If the array is empty, there is no physics data to load.

            // Overwrite grid values with values from results array (only for nodes with coordinates that lie within
            // the FE mesh)
            vtkPoints* points = output->GetPoints();
            double point[2];
            Vector2d inv_cell_size = Vector2d(1.0 / cell_size(0), 1.0 / cell_size(1));
            for (int i = 0; i < points->GetNumberOfPoints(); i++) {
                point[0] = points->GetData()->GetTuple(i)[0];
                point[1] = points->GetData()->GetTuple(i)[1];
                Vector2d origin_aligned_coord = Vector2d(point[0], point[1]) - offset;
                Vector2d gridscale_coord = inv_cell_size.cwiseProduct(origin_aligned_coord);
                int coord = (round(gridscale_coord[0]) * (dim_y + 1) + round(gridscale_coord[1]));

                // Accumulate the largest absolute values for stress and displacement
                stress_xx->SetValue(i, nodewise_compressive_xx[coord]);
                stress_yy->SetValue(i, nodewise_compressive_yy[coord]);
                displacements->SetValue(i * 3 + 2, nodewise_displacements[coord]);
                principal_compressive_stresses->SetValue(i, nodewise_principle_compressive_stresses[coord]);
                principal_tensile_stresses->SetValue(i, nodewise_principle_tensile_stresses[coord]);
                if (help::is_in(mechanical_constraint, "Mohr")) {
                    modified_mohr->SetValue(i, nodewise_modified_mohr[coord]);
                }
                vonmises->SetValue(i, nodewise_vonmises[coord]);
            }

            // Write to vtk OUTFILE
            vtkUnstructuredGridWriter* writer = vtkUnstructuredGridWriter::New();
            writer->SetFileName(outfile.c_str());
            writer->SetFileTypeToASCII();
            writer->SetInputData(reader->GetOutput());
            writer->Update();

            // Free used memory
            delete[] nodewise_tensile_xx;
            delete[] nodewise_tensile_yy;
            delete[] nodewise_compressive_xx;
            delete[] nodewise_compressive_yy;
            delete[] nodewise_xy;
            delete[] nodewise_displacements;
            delete[] nodewise_principle_compressive_stresses;
            delete[] nodewise_principle_tensile_stresses;
            delete[] nodewise_modified_mohr;
            delete[] nodewise_vonmises;
            writer->Delete();


            return;
        }

        static inline double get_modified_mohr_safety(double theta_1, double theta_3, double max_tensile_strength, double max_compressive_strength) {
            if (theta_1 < 0 && theta_3 < 0) {
                return abs(theta_3) / max_compressive_strength;
            }
            else if ((theta_1 > 0 && theta_3 > 0) || (theta_1 > 0 && theta_1 > abs(theta_3))) {
                return theta_1 / max_tensile_strength;
            }
            else {
                if (theta_1 < 0 && theta_3 > 0) cout << "theta 3 is larger than theta 1. This is a bug\n";
                return ((max_compressive_strength - max_tensile_strength) * theta_1) / (max_compressive_strength * max_tensile_strength) - (theta_3 / max_compressive_strength);
            }
        }

        static void compute_principal_stresses(double* stress_xx, double* stress_yy, double* stress_xy, double* principal_stresses, int dim_x, int dim_y, string mechanical_constraint) {
            for (int node_coord = 0; node_coord < dim_x * dim_y; node_coord++) {
                double xx = stress_xx[node_coord];
                double yy = stress_yy[node_coord];
                double xy = stress_xy[node_coord];
                double direction_wrt_mohr_circle_center = (double)((mechanical_constraint == "Tensile") ? 1 : -1);
                double princ_stress = (xx + yy) / 2.0 + direction_wrt_mohr_circle_center * sqrt(((xx - yy) / 2.0) * ((xx - yy) / 2.0) + xy * xy);
                principal_stresses[node_coord] = princ_stress;
            }
        }

        static bool load_single_VTK_file(
            string filename, FEAResults2D* results, FEACaseManager* fea_casemanager, int dim_x, int dim_y, Vector2d cell_size,
            Vector3d _offset, string mechanical_constraint, vector<int>* border_nodes, vector<int>* times = 0, bool verbose = false
        ) {

            Timer timer1; timer1.start();
            // Read data from file
            vtkUnstructuredGridReader* reader = vtkUnstructuredGridReader::New();
            reader->SetFileName(filename.c_str());
            reader->ReadAllScalarsOn();
            reader->Update();
            vtkUnstructuredGrid* output = reader->GetOutput();

            // Initially, populate results arrays with zeroes to make sure they're initialized (nodes on the grid which are part of the FE mesh will
            // have their corresponding values in the results array overwritten later)
            Vector2d offset = Vector2d(_offset(0), _offset(1));
            double* results1_nodewise = new double[(dim_x + 1) * (dim_y + 1)]; // Nodes grid has +1 width along each dim
            double* results2_nodewise = 0; // Nodes grid has +1 width along each dim
            double* results3_nodewise = 0; // Nodes grid has +1 width along each dim
            double* principal_stresses_nodewise = new double[(dim_x + 1) * (dim_y + 1)];
            map<int, int> node_coords_map;
            help::populate_with_zeroes(results1_nodewise, dim_x + 1, dim_y + 1);
            string mechanical_constraint1, mechanical_constraint2, mechanical_constraint3 = mechanical_constraint;
            double sign = 1;
            if (help::is_in(mechanical_constraint, "ModifiedMohr")) {
                mechanical_constraint1 = "Stress_xx";
                mechanical_constraint2 = "Stress_yy";
                mechanical_constraint3 = "Stress_xy";
                results2_nodewise = new double[(dim_x + 1) * (dim_y + 1)];
                results3_nodewise = new double[(dim_x + 1) * (dim_y + 1)];
                help::populate_with_zeroes(results2_nodewise, dim_x + 1, dim_y + 1);
                help::populate_with_zeroes(results3_nodewise, dim_x + 1, dim_y + 1);
                help::populate_with_zeroes(principal_stresses_nodewise, dim_x + 1, dim_y + 1);
            }
            timer1.stop();
            if (verbose && times != 0) {
                times->push_back(timer1.elapsedMilliseconds());
            }

            Timer timer2; timer2.start();
            // Obtain and process nodewise results
            if (mechanical_constraint == "Displacement") {
                phys::load_nodewise_results(output, results1_nodewise, dim_x, dim_y, cell_size, offset, mechanical_constraint1, border_nodes, 1, false, &node_coords_map);
            }
            else if (help::is_in(mechanical_constraint, "ModifiedMohr")) {
                phys::load_nodewise_results(output, results1_nodewise, dim_x, dim_y, cell_size, offset, mechanical_constraint1, border_nodes, 1, false, &node_coords_map);
                phys::load_nodewise_results(output, results2_nodewise, dim_x, dim_y, cell_size, offset, mechanical_constraint2, border_nodes, 1, false);
                phys::load_nodewise_results(output, results3_nodewise, dim_x, dim_y, cell_size, offset, mechanical_constraint3, border_nodes, 1, false);
                timer2.stop();
                if (verbose && times != 0) {
                    times->push_back(timer2.elapsedMilliseconds());
                }
                double* theta1_nodewise = new double[(dim_x + 1) * (dim_y + 1)];
                double* theta3_nodewise = new double[(dim_x + 1) * (dim_y + 1)];
                help::populate_with_zeroes(theta1_nodewise, dim_x + 1, dim_y + 1);
                help::populate_with_zeroes(theta3_nodewise, dim_x + 1, dim_y + 1);
                compute_principal_stresses(results1_nodewise, results2_nodewise, results3_nodewise, theta1_nodewise, dim_x, dim_y, "Tensile");
                compute_principal_stresses(results1_nodewise, results2_nodewise, results3_nodewise, theta3_nodewise, dim_x, dim_y, "Compressive");
                for (int i = 0; i < (dim_x + 1) * (dim_y + 1); i++) {
                    principal_stresses_nodewise[i] = get_modified_mohr_safety(
                        theta1_nodewise[i], theta3_nodewise[i], fea_casemanager->max_tensile_strength, fea_casemanager->max_compressive_strength
                    );
                }
                delete[] theta1_nodewise;
                delete[] theta3_nodewise;
            }
            else if (help::is_in(mechanical_constraint, "Vonmises")) {
                phys::load_nodewise_results(output, principal_stresses_nodewise, dim_x, dim_y, cell_size, offset, "Vonmises", border_nodes, 1, false, &node_coords_map);
            }
            else {
                // Compute only compressive- or tensile principal stress (depending on the mechanical constraint used)
                phys::load_nodewise_results(output, results1_nodewise, dim_x, dim_y, cell_size, offset, mechanical_constraint1, border_nodes, 1, false, &node_coords_map);
                phys::load_nodewise_results(output, results2_nodewise, dim_x, dim_y, cell_size, offset, mechanical_constraint2, border_nodes, 1, false, &node_coords_map);
                phys::load_nodewise_results(output, results3_nodewise, dim_x, dim_y, cell_size, offset, mechanical_constraint3, border_nodes, 1, false, &node_coords_map);
                compute_principal_stresses(results1_nodewise, results2_nodewise, results3_nodewise, principal_stresses_nodewise, dim_x, dim_y, mechanical_constraint);
            }

            // Create cellwise results distribution by taking mean of each group of 4 corners of a cell
            double min_mechanical_metric = 1e30;
            double max_mechanical_metric = 0;
            double max_yield_criterion = 0;
            for (auto& [node_coord, _] : node_coords_map) {
                //if (verbose) cout << "node coord " << node_coord << ", point index " << _ << endl;

                // Obtain coordinates of the bottomleft corner (= node) of the cell
                int x = node_coord / (dim_y + 1);
                int y = node_coord % (dim_y + 1);
                if (x == dim_x || y == dim_y) continue; // Skip coordinates outside cell domain
                int cell_coord = x * dim_y + y;

                // Check if the cell should be skipped
                if (help::is_in(&fea_casemanager->inactive_cells, cell_coord)) {
                    // Cells marked as 'inactive' are ignored during solution evaluation.
                    results->data_map.insert(pair(cell_coord, -9999));
                    continue;
                }
                if (mechanical_constraint == "Displacement" && cell_coord != fea_casemanager->displacement_measurement_cell) {
                    continue; // If using mechanical constraint 'Displacement', we only consider the displacement measurement cell and thus ignore all other cells.
                }

                // Compute the cell average from the four values defined on its corner nodes
                vector<int> corner_coords = {
                    node_coord,
                    (x + 1) * (dim_y + 1) + y,
                    (x + 1) * (dim_y + 1) + (y + 1),
                    x * (dim_y + 1) + (y + 1)
                };
                Vector4d cell_corners;
                if (mechanical_constraint == "Displacement") {
                    cell_corners[0] = results1_nodewise[corner_coords[0]];
                    cell_corners[1] = results1_nodewise[corner_coords[1]];
                    cell_corners[2] = results1_nodewise[corner_coords[2]];
                    cell_corners[3] = results1_nodewise[corner_coords[3]];
                }
                else {
                    cell_corners[0] = principal_stresses_nodewise[corner_coords[0]];
                    cell_corners[1] = principal_stresses_nodewise[corner_coords[1]];
                    cell_corners[2] = principal_stresses_nodewise[corner_coords[2]];
                    cell_corners[3] = principal_stresses_nodewise[corner_coords[3]];
                }
                double cell_value = cell_corners.mean();

                // Check if the found cell value is the largest found so far
                if (help::is_in(mechanical_constraint, "Compressive")) cell_value = -cell_value;
                cell_value = max(0.0, cell_value);
                results->data_map.insert(pair(cell_coord, cell_corners.mean()));
                if (mechanical_constraint != "Displacement") {
                    if (cell_value > max_mechanical_metric) max_mechanical_metric = cell_value;
                    if (cell_value < min_mechanical_metric) min_mechanical_metric = cell_value;
                    if (cell_value > max_yield_criterion) max_yield_criterion = cell_value;
                }
                if (help::is_in(mechanical_constraint, "Displacement") && cell_coord == fea_casemanager->displacement_measurement_cell) {
                    double* displacements = new double[(dim_x + 1) * (dim_y + 1)];
                    help::populate_with_zeroes(displacements, (dim_x + 1), (dim_y + 1));
                    phys::load_nodewise_results(output, displacements, dim_x, dim_y, cell_size, offset, "Displacement", border_nodes, 1, false);
                    //phys::load_nodewise_results(output, &displacements, dim_x, dim_y, cell_size, offset, "Displacement", &node_coords_map, &corner_coords);
                    //double max_displacement = help::get_max(&displacements);
                    cell_corners[0] = displacements[corner_coords[0]];
                    cell_corners[1] = displacements[corner_coords[1]];
                    cell_corners[2] = displacements[corner_coords[2]];
                    cell_corners[3] = displacements[corner_coords[3]];
                    delete[] displacements;
                    results->max_displacement = cell_corners.maxCoeff();
                    if (results->max_displacement > fea_casemanager->max_displacement) {
                        double relative_displacement = results->max_displacement / fea_casemanager->max_displacement;
                        //cout << "Displacement triggered (" << max_displacement << "). More severe than stress (" << max_mechanical_metric << ")? " << ((relative_displacement * fea_casemanager->mechanical_threshold > max_mechanical_metric) ? "Yes" : "No") << endl;
                        max_mechanical_metric = max(max_mechanical_metric, relative_displacement * fea_casemanager->mechanical_threshold);
                    }
                }
            }
            results->min = min_mechanical_metric;
            results->max = max_mechanical_metric;
            results->max_yield_criterion = max_yield_criterion;
            delete[] results1_nodewise;
            if (results2_nodewise != nullptr) delete[] results2_nodewise;
            if (results3_nodewise != nullptr) delete[] results3_nodewise;
            if (principal_stresses_nodewise != nullptr) delete[] principal_stresses_nodewise;
            reader->Delete();

            return true;
        }
    };
}



