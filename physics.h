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
#include "timer.h"


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
            double max_stress_threshold, max_tensile_strength, max_compressive_strength;
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
            vector<int> border_nodes = {};
            double max_stress_threshold = INFINITY;
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
            Vector2d offset, string mechanical_constraint, vector<int>* border_nodes, double sign, bool clamp, vector<int>* _coords = 0
        ) {
            if (_coords) _coords->clear();

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
            vtkDoubleArray* results_array = dynamic_cast<vtkDoubleArray*>(point_data->GetScalars(mechanical_constraint.c_str()));
            if (results_array->GetSize() == 0) return; // If the array is empty, there is no physics data to load.

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
                //if (help::is_in(border_nodes, coord)) {
                //    results_nodewise[coord] = 0; // Skip nodes on the border of the shape
                //    continue;
                //}
                if (_coords != 0) _coords->push_back(coord);
                else coords.push_back(coord);
                if (mechanical_constraint == "Displacement") {
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
            reader->Delete();
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
                phys::load_nodewise_results(casepath, single_run_stress_xx, dim_x, dim_y, cell_size, offset, "Stress_xx", border_nodes, 1, false);
                phys::load_nodewise_results(casepath, single_run_stress_yy, dim_x, dim_y, cell_size, offset, "Stress_yy", border_nodes, 1, false);
                phys::load_nodewise_results(casepath, single_run_stress_xy, dim_x, dim_y, cell_size, offset, "Stress_xy", border_nodes, 1, false);
                phys::load_nodewise_results(casepath, single_run_displacements, dim_x, dim_y, cell_size, offset, "Displacement", border_nodes, 1, false);
                phys::load_nodewise_results(casepath, single_run_vonmises, dim_x, dim_y, cell_size, offset, "Vonmises", border_nodes, 1, false);
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
                delete[] single_run_stress_xx;
                delete[] single_run_stress_yy;
                delete[] single_run_stress_xy;
                delete[] single_run_displacements;
                delete[] single_run_principal_compressive_stresses;
                delete[] single_run_principal_tensile_stresses;
                delete[] single_run_modified_mohr;
                delete[] single_run_vonmises;
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

                // Accumulate the largest absolute values for stress and displacement
                if (nodewise_tensile_xx[coord] > -nodewise_compressive_xx[coord]) {
                    stress_xx->SetValue(i, nodewise_tensile_xx[coord]);
                }
                else stress_xx->SetValue(i, nodewise_compressive_xx[coord]);
                if (nodewise_tensile_yy[coord] > -nodewise_compressive_yy[coord]) {
                    stress_yy->SetValue(i, nodewise_tensile_yy[coord]);
                }
                else stress_yy->SetValue(i, nodewise_compressive_yy[coord]);
                displacements->SetValue(i * 3 + 2, nodewise_displacements[coord]);
                principal_compressive_stresses->SetValue(i, nodewise_principle_compressive_stresses[coord]);
                principal_tensile_stresses->SetValue(i, nodewise_principle_tensile_stresses[coord]);
                modified_mohr->SetValue(i, nodewise_modified_mohr[coord]);
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
            Vector3d _offset, string mechanical_constraint, vector<int>* border_nodes, vector<int>* times, bool verbose = false
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
            double* principal_stresses_nodewise = 0;
            vector<int> coords;
            help::populate_with_zeroes(results1_nodewise, dim_x + 1, dim_y + 1);
            string mechanical_constraint1, mechanical_constraint2, mechanical_constraint3 = mechanical_constraint;
            double sign = 1;
            if (!help::is_in(mechanical_constraint, "Displacement")) {
                mechanical_constraint1 = "Stress_xx";
                mechanical_constraint2 = "Stress_yy";
                mechanical_constraint3 = "Stress_xy";
                results2_nodewise = new double[(dim_x + 1) * (dim_y + 1)];
                results3_nodewise = new double[(dim_x + 1) * (dim_y + 1)];
                principal_stresses_nodewise = new double[(dim_x + 1) * (dim_y + 1)];
                help::populate_with_zeroes(results2_nodewise, dim_x + 1, dim_y + 1);
                help::populate_with_zeroes(results3_nodewise, dim_x + 1, dim_y + 1);
                help::populate_with_zeroes(principal_stresses_nodewise, dim_x + 1, dim_y + 1);
            }
            timer1.stop();
            if (verbose) {
                cout << "Setup time: " << timer1.elapsedMilliseconds() << " ms\n";
                times->push_back(timer1.elapsedMilliseconds());
            }

            Timer timer2; timer2.start();
            // Obtain and process nodewise results
            if (mechanical_constraint == "Displacement") {
                phys::load_nodewise_results(filename, results1_nodewise, dim_x, dim_y, cell_size, offset, mechanical_constraint1, border_nodes, 1, false, &coords);
            }
            if (mechanical_constraint == "ModifiedMohr") {
                phys::load_nodewise_results(filename, results1_nodewise, dim_x, dim_y, cell_size, offset, mechanical_constraint1, border_nodes, 1, false, &coords);
                phys::load_nodewise_results(filename, results2_nodewise, dim_x, dim_y, cell_size, offset, mechanical_constraint2, border_nodes, 1, false);
                phys::load_nodewise_results(filename, results3_nodewise, dim_x, dim_y, cell_size, offset, mechanical_constraint3, border_nodes, 1, false);
                timer2.stop();
                if (verbose) {
                    cout << "Read time: " << timer2.elapsedMilliseconds() << " ms\n";
                    times->push_back(timer2.elapsedMilliseconds());
                }
                double* theta1_nodewise = new double[(dim_x + 1) * (dim_y + 1)];
                double* theta3_nodewise = new double[(dim_x + 1) * (dim_y + 1)];
                help::populate_with_zeroes(theta1_nodewise, dim_x + 1, dim_y + 1);
                help::populate_with_zeroes(theta3_nodewise, dim_x + 1, dim_y + 1);
                Timer timer3; timer3.start();
                compute_principal_stresses(results1_nodewise, results2_nodewise, results3_nodewise, theta1_nodewise, dim_x, dim_y, "Tensile");
                compute_principal_stresses(results1_nodewise, results2_nodewise, results3_nodewise, theta3_nodewise, dim_x, dim_y, "Compressive");
                for (int i = 0; i < (dim_x + 1) * (dim_y + 1); i++) {
                    principal_stresses_nodewise[i] = get_modified_mohr_safety(
                        theta1_nodewise[i], theta3_nodewise[i], fea_casemanager->max_tensile_strength, fea_casemanager->max_compressive_strength
                    );
                }
                timer3.stop();
                if (verbose) cout << "Calculation time: " << timer3.elapsedMilliseconds() << " ms\n";
            }
            else {
                phys::load_nodewise_results(filename, results1_nodewise, dim_x, dim_y, cell_size, offset, mechanical_constraint1, border_nodes, 1, false, &coords);
                phys::load_nodewise_results(filename, results2_nodewise, dim_x, dim_y, cell_size, offset, mechanical_constraint2, border_nodes, 1, false, &coords);
                phys::load_nodewise_results(filename, results3_nodewise, dim_x, dim_y, cell_size, offset, mechanical_constraint3, border_nodes, 1, false, &coords);
                compute_principal_stresses(results1_nodewise, results2_nodewise, results3_nodewise, principal_stresses_nodewise, dim_x, dim_y, mechanical_constraint);
            }

            Timer timer4; timer4.start();

            // Create cellwise results distribution by taking mean of each group of 4 corners of a cell
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
                if (mechanical_constraint == "Displacement" && cell_coord != fea_casemanager->displacement_measurement_cell) {
                    continue; // If using mechanical constraint 'Displacement', we only consider the displacement measurement cell and thus ignore all other cells.
                }
                Vector4d neighbors;
                int num_neighbors = 4;
                int node_neighbor1 = (x + 1) * (dim_y + 1) + y;
                int node_neighbor2 = (x + 1) * (dim_y + 1) + (y + 1);
                int node_neighbor3 = x * (dim_y + 1) + (y + 1);
                if (mechanical_constraint == "Displacement") {
                    neighbors[0] = results1_nodewise[node_coord];
                    neighbors[1] = results1_nodewise[node_neighbor1];
                    neighbors[2] = results1_nodewise[node_neighbor2];
                    neighbors[3] = results1_nodewise[node_neighbor3];
                }
                else {
                    neighbors[0] = principal_stresses_nodewise[node_coord];
                    neighbors[1] = principal_stresses_nodewise[node_neighbor1];
                    neighbors[2] = principal_stresses_nodewise[node_neighbor2];
                    neighbors[3] = principal_stresses_nodewise[node_neighbor3];
                }
                double cell_value = neighbors.mean();
                /*float intra_cell_variance = (cell_value - neighbors.mean());
                if (help::is_in(border_nodes, node_coord)) { num_neighbors--; }
                if (help::is_in(border_nodes, (x + 1) * (dim_y + 1) + y)) { num_neighbors--; }
                if (help::is_in(border_nodes, (x + 1) * (dim_y + 1) + (y + 1))) { num_neighbors--; }
                if (help::is_in(border_nodes, x * (dim_y + 1) + (y + 1))) { num_neighbors--; }*/
                
                if (help::is_in(mechanical_constraint, "Compressive")) cell_value = -cell_value;
                cell_value = max(0.0, cell_value);
                results->data_map.insert(pair(cell_coord, neighbors.mean()));
                //cout << "cell stress (loader): " << neighbors.mean() << endl;
                if (mechanical_constraint == "Displacement") {
                    max_stress = cell_value;
                }
                else {
                    if (cell_value > max_stress) max_stress = cell_value;
                    if (cell_value < min_stress) min_stress = cell_value;
                }
            }
            timer4.stop();
            if (verbose) cout << "Cell handling time: " << timer4.elapsedMilliseconds() << " ms\n";
            results->min = min_stress;
            results->max = max_stress;
            delete[] results1_nodewise;
            reader->Delete();

            return true;
        }
    };
}



