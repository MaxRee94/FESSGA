#pragma once
#include "physics.h"


using namespace fessga;


Vector2d phys::FEACase::get_node_coords(int idx) {
	int x = idx / dim_y;
	int y = idx % dim_y;
	return Vector2d((double)x, (double)y);
}

void phys::FEACase::compute_node_barycenters() {
	for (auto& [bound_name, lines] : bound_cond_lines) {
		vector<int> _nodes;
		Vector2d coordinate_sum(0, 0);
		for (auto& line : lines) {
			if (!help::is_in(&_nodes, line.first)) {
				int x = line.first / dim_y;
				int y = line.first % dim_y;
				coordinate_sum += Vector2d((double)x, (double)y);
				_nodes.push_back(line.first);
			}
			if (!help::is_in(&_nodes, line.second)) {
				int x = line.second / dim_y;
				int y = line.second % dim_y;
				coordinate_sum += Vector2d((double)x, (double)y);
				_nodes.push_back(line.second);
			}
		}
		node_barycenters[bound_name] = (1.0 / (double)_nodes.size()) * coordinate_sum;
	}
}

Vector2d phys::FEACase::compute_cell_barycenter(vector<int>* cells) {
	Vector2d summed_coords;
	for (auto& cell : *cells) summed_coords += get_node_coords(cell);
	return summed_coords / cells->size();
}

void phys::FEACase::compute_cell_barycenters() {
	for (auto& [bound_name, bound_cell_map] : bound_cond_cells) {
		Vector2d bound_cell_barycenter = compute_cell_barycenter(&bound_cell_map["bound"]);
		Vector2d cutout_cell_barycenter = compute_cell_barycenter(&bound_cell_map["cutout"]);
		Vector2d keep_cell_barycenter = compute_cell_barycenter(&bound_cell_map["keep"]);
		map<string, Vector2d> cell_map;
		cell_map["bound"] = bound_cell_barycenter;
		cell_map["cutout"] = cutout_cell_barycenter;
		cell_map["keep"] = keep_cell_barycenter;
		cell_barycenters[bound_name] = cell_map;
	}
}

void phys::FEACase::compute_barycenters() {
	compute_node_barycenters();
	compute_cell_barycenters();
}


// Resample boundary by adding copies of nodes, starting from the barycenter outwards
void phys::FEACaseManager::resample_boundary(vector<int>& nodes, Vector2d barycenter, int target_no_nodes) {
	vector<int> _nodes = nodes;
	for (int& node : _nodes) {
		if (nodes.size() == target_no_nodes) break;
		nodes.push_back(node);
	}
	if (nodes.size() < target_no_nodes) resample_boundary(nodes, barycenter, target_no_nodes);
}

Vector2d phys::FEACaseManager::get_node_coords(int idx) {
	int x = idx / dim_y;
	int y = idx % dim_y;
	return Vector2d((double)x, (double)y);
}

// Get the index of the nearest node on the grid relative to the given node coordinates
tuple<int, int> phys::FEACaseManager::get_nearest_neighbor(Vector2d coords, vector<int>* _neighbors) {
	double min_dist = INFINITY;
	int global_idx = -1; int local_idx = -1;

	vector<int> neighbors;
	if (!_neighbors) {
		vector<pair<int, int>> offsets = { pair(0,0), pair(1,0), pair(1,1), pair(1,0) };
		for (auto& offset : offsets) {
			int x = (int)coords[0] + offset.first;
			int y = (int)coords[1] + offset.second;
			neighbors.push_back(x * dim_y + y);
		}
	}

	else neighbors = *_neighbors;

	for (int i = 0; i < neighbors.size(); i++) {
		int _neighbor_idx = neighbors[i];
		Vector2d neighbor_coords = get_node_coords(_neighbor_idx);
		double dist = (neighbor_coords - coords).norm();
		if (dist < min_dist) {
			min_dist = dist;
			global_idx = _neighbor_idx;
			local_idx = i;
		}
	}

	return { global_idx, local_idx };
}

// Get vector to the nearest neighbor of the given node
Vector2d phys::FEACaseManager::get_vector2nn(Vector2d coords, vector<int>& neighbors) {
	auto [nn_global_idx, nn_local_idx] = get_nearest_neighbor(coords, &neighbors);
	Vector2d vector2nn = get_node_coords(nn_global_idx) - coords;
	neighbors.erase(neighbors.begin() + nn_local_idx);

	return vector2nn;
}

void phys::FEACaseManager::order_nodes_by_distance(vector<int>* nodes, Vector2d barycenter) {
	map<int, double> distances;
	map<int, int> copied_vals;
	for (auto& _node : *nodes) {
		Vector2d node = get_node_coords(_node);
		double distance = (barycenter - node).norm();
		if (help::get_value(&distances, _node) != -1) {
			if (help::get_value(&copied_vals, _node) != -1) copied_vals[_node]++;
			else copied_vals[_node] = 1;
		}
		distances[_node] = distance;
	}
	PairSet _ordered_nodes;
	help::sort(distances, _ordered_nodes);
	nodes->clear();
	for (auto& [pnt, _] : _ordered_nodes) {
		nodes->push_back(pnt);
		int copies = help::get_value(&copied_vals, pnt);
		for (int i = 0; i < copies; i++) nodes->push_back(pnt);
	}
}

void phys::FEACaseManager::compute_migration_vectors(phys::FEACase* source, phys::FEACase* target, int pair_idx) {
	for (auto& [bound_name, lines] : source->bound_cond_lines) {
		// Compute barycenter alignment vector
		Vector2d align_to_barycenter = target->node_barycenters[bound_name] - source->node_barycenters[bound_name];

		// Convert lines vectors to vectors of nodes
		vector<int> target_nodes;
		for (auto& line : target->bound_cond_lines[bound_name]) {
			if (!help::is_in(&target_nodes, line.first)) target_nodes.push_back(line.first);
			if (!help::is_in(&target_nodes, line.second)) target_nodes.push_back(line.second);
		}
		vector<int> source_nodes;
		for (auto& line : lines) {
			if (!help::is_in(&source_nodes, line.first)) source_nodes.push_back(line.first);
			if (!help::is_in(&source_nodes, line.second)) source_nodes.push_back(line.second);
		}

		// Order source nodes by distance to barycenter
		order_nodes_by_distance(&source_nodes, source->node_barycenters[bound_name]);

		// Resample if the number of nodes of the boundary condition differs between source and target
		if (target_nodes.size() > source_nodes.size()) resample_boundary(source_nodes, source->node_barycenters[bound_name],
			target_nodes.size());
		else if (source_nodes.size() > target_nodes.size()) {
			order_nodes_by_distance(&target_nodes, source->node_barycenters[bound_name]);
			resample_boundary(target_nodes, target->node_barycenters[bound_name], source_nodes.size());
		}
		order_nodes_by_distance(&source_nodes, source->node_barycenters[bound_name]);

		// Compute vectors to align nodes
		vector<Vector2d> _migration_vectors;
		for (auto& src_pnt : source_nodes) {
			Vector2d source_node = source->get_node_coords(src_pnt);
			source_node += align_to_barycenter; // Align the cluster of source nodes to the barycenter of the target cluster
			Vector2d to_nn = get_vector2nn(source_node, target_nodes);
			_migration_vectors.push_back(align_to_barycenter + to_nn);
		}

		// Store bound nodes and morph vectors
		source->bound_cond_nodes[bound_name] = source_nodes;
		target->bound_cond_nodes[bound_name] = target_nodes;
		migration_vectors[pair_idx][bound_name] = _migration_vectors;
	}
}

int phys::FEACaseManager::get_direct_unvisited_node_neighbor(
	int node_idx, vector<int>* all_nodes, vector<int>* visited_nodes
) {
	vector<pair<int, int>> offsets = { pair(0, 1), pair(0, -1), pair(1, 0), pair(-1, 0) };
	int x = node_idx / (dim_y);
	int y = node_idx % (dim_y);
	for (auto& offset : offsets) {
		int neighbor = (x + offset.first) * (dim_y) + (y + offset.second);
		if (help::is_in(visited_nodes, neighbor)) continue;
		if (help::is_in(all_nodes, neighbor)) return neighbor;
	}
	return -1; // Return -1 if we did not find any unvisited direct neighbor
}

vector<int> phys::FEACaseManager::get_indirect_path_to_unvisited_node_neighbor(
	int node_idx, vector<int>* all_nodes, vector<int>* visited_nodes, int radius
) {
	if (radius > 3) return {};

	// Find a neighbor 
	vector<pair<int, int>> offsets = {
		pair(0, radius), pair(0, -radius), pair(radius, 0), pair(-radius, 0),
		pair(radius, radius), pair(-radius, -radius), pair(radius, radius), pair(-radius, -radius)
	};
	int x = node_idx / (dim_y);
	int y = node_idx % (dim_y);
	pair<int, int> neighbor = pair(-1, -1);
	int neighbor_idx;
	for (auto& offset : offsets) {
		neighbor_idx = (x + offset.first) * (dim_y) + (y + offset.second);
		if (help::is_in(visited_nodes, neighbor_idx)) continue;
		if (help::is_in(all_nodes, neighbor_idx)) {
			neighbor = pair(x + offset.first, y + offset.second);
			break;
		}
	}

	// Create path to neighbor
	if (neighbor.first == -1) return get_indirect_path_to_unvisited_node_neighbor(node_idx, all_nodes, visited_nodes, radius + 1);
	else {
		vector<int> path;
		pair<int, int> vec2neighbor = pair(neighbor.first - x, neighbor.second - y);
		int x_sign = vec2neighbor.first != 0 ? vec2neighbor.first / abs(vec2neighbor.first) : 0;
		int y_sign = vec2neighbor.second != 0 ? vec2neighbor.second / abs(vec2neighbor.second) : 0;
		int x_to_travel = abs(vec2neighbor.first);
		int y_to_travel = abs(vec2neighbor.second);

		// Initialize x_frac to total y distance to travel and vice versa, so that first segment will be dedicated to
		// travelling along the axis for which we have the most distance to cover.
		float x_frac = vec2neighbor.second; float y_frac = vec2neighbor.first;
		int x_travelled = 0; int y_travelled = 0;
		while (x_travelled < x_to_travel || y_travelled < y_to_travel) {
			// Use the next path segment to travel along the axis for which a smaller proportion of distance has been travelled so far.
			if (x_frac <= y_frac) x_travelled += 1;
			else y_travelled += 1;
			int _x = x + x_sign * x_travelled;
			int _y = y + y_sign * y_travelled;
			x_frac = (float)x_travelled / (float)abs(vec2neighbor.first);
			y_frac = (float)y_travelled / (float)abs(vec2neighbor.second);
			int node_idx = _x * dim_y + _y;
			path.push_back(node_idx);
		}
		/*cout << "\n\npath: \n";
		help::print_vector(&path);
		cout << "\n\n";*/
		return path;
	}	
}

vector<int> phys::FEACaseManager::get_path_to_unvisited_node_neighbor(
	int node_idx, vector<int>* all_nodes, vector<int>* visited_nodes
) {
	vector<int> path_to_neighbor;
	int neighbor = get_direct_unvisited_node_neighbor(node_idx, all_nodes, visited_nodes);
	if (neighbor == -1) {
		int search_radius = 1;
		path_to_neighbor = get_indirect_path_to_unvisited_node_neighbor(node_idx, all_nodes, visited_nodes, search_radius);
	}
	else path_to_neighbor = { neighbor };
	return path_to_neighbor;
}

// Order the given list of boundary condition nodes such that they describe a walk from one of the endpoints to the other
void phys::FEACaseManager::walk_and_collect_bound_nodes(vector<int>* unordered_boundary, vector<int>& walk, int start_node) {
	int node;
	if (start_node != -1) node = start_node;
	else node = unordered_boundary->at(0);
	int endpoint;
	vector<int> neighbors = { node };
	while (true) {
		for (auto& neighbor: neighbors) walk.push_back(neighbor);
		neighbors = get_path_to_unvisited_node_neighbor(node, unordered_boundary, &walk);
		if (neighbors.size() == 0) {
			endpoint = node;
			break;
		}
		node = neighbors.back();
	}
	if (start_node == -1) { 
		// If start_node == -1, this means we started our walk from a randomly chosen node.
		// We therefore need to re-start the walk from the found endpoint to obtain a complete walk of the boundary condition's path.
		walk.clear();
		walk_and_collect_bound_nodes(unordered_boundary, walk, endpoint);
	}
}

void phys::FEACaseManager::interpolate_nodes(FEACase* active_case, FEACase* source, float fraction, int pair_idx) {
	active_case->bound_cond_nodes.clear();
	for (auto& [bound_name, bound_cond_nodes] : source->bound_cond_nodes) {
		//cout << "\nbound name: " << bound_name << endl;
		
		// Get interpolated boundary nodes
		vector<int> interpolated_bound_cond_nodes;
		for (int i = 0; i < bound_cond_nodes.size(); i++) {
			Vector2d source_coords = source->get_node_coords(bound_cond_nodes[i]);
			//cout << "source coords: " << source_coords.transpose() << endl;
			Vector2d interpolate_vec = fraction * migration_vectors[pair_idx][bound_name][i];
			//cout << "interpolate coords: " << (source_coords + interpolate_vec).transpose() << endl;
			auto [interpolated_node, _] = get_nearest_neighbor(source_coords + interpolate_vec);
			if (!help::is_in(&interpolated_bound_cond_nodes, interpolated_node)) {
				//cout << "adding interpolated node: " << interpolated_node << endl;
				interpolated_bound_cond_nodes.push_back(interpolated_node);
			}
		}

		// Reconstruct boundary lines from node list
		vector<int> boundary_cond_walk;
		walk_and_collect_bound_nodes(&interpolated_bound_cond_nodes, boundary_cond_walk);
		cout << "walk length: " << boundary_cond_walk.size() << endl;
		vector<pair<int, int>> lines;
		for (int i = 0; i < boundary_cond_walk.size() - 1; i++) {
			pair<int, int> line = pair(boundary_cond_walk[i], boundary_cond_walk[i + 1]);
			lines.push_back(line);
		}
		active_case->bound_cond_lines[bound_name] = lines;
		active_case->bound_cond_nodes[bound_name] = interpolated_bound_cond_nodes;

		// Temp
		cout << "No bound cond nodes: " << interpolated_bound_cond_nodes.size() << endl;
		cout << "No bound cond lines: " << lines.size() << endl;
	}
}

vector<int> phys::FEACaseManager::get_unvisited_neighbor_cells(int cell, vector<int>& neighbors, vector<int>* visited_cells) {
	vector<pair<int, int>> offsets = { pair(0,1), pair(1,0), pair(-1, 0), pair(0, -1), pair(-1,-1), pair(-1, 1), pair(1,1), pair(1,-1) };
	int x = cell / (dim_y - 1); int y = cell % (dim_y - 1);
	for (auto& offset : offsets) {
		int _x = x + offset.first;
		int _y = y + offset.second;
		if (_x == (dim_x - 1) || _y == (dim_y - 1) || _x < 0 || _y < 0) continue;
		int neighbor_coord = _x * (dim_y - 1) + _y;
		if (!help::is_in(visited_cells, neighbor_coord)) neighbors.push_back(neighbor_coord);
	}
	return neighbors;
}

pair<int, int> phys::FEACaseManager::get_cells_with_given_edge(pair<int, int> edge) {
	int node1_x = edge.first / dim_y;
	int node1_y = edge.first % dim_y;
	int node2_x = edge.second / dim_y;
	int node2_y = edge.second % dim_y;
	int cell1_x, cell1_y, cell2_x, cell2_y;
	bool vertical = node1_x == node2_x;
	if (vertical) {
		cell1_x = node1_x - 1;
		cell2_x = node1_x;
		cell1_y = min(node1_y, node2_y);
		cell2_y = cell1_y;
	}
	else {
		cell1_y = node1_y - 1;
		cell2_y = node1_y;
		cell1_x = min(node1_x, node2_x);
		cell2_x = cell1_x;
	}

	int cell1 = cell1_x * (dim_y - 1) + cell1_y;
	int cell2 = cell2_x * (dim_y - 1) + cell2_y;

	return pair(cell1, cell2);
}

// Return cell that is not blocked by the given line with respect to the previous cell
tuple<int, int> phys::FEACaseManager::get_bound_cells(
	int prev_cell, pair<int, int> candidate_cells, pair<int, int> prev_line, pair<int, int> next_line
) {
	// Get prev node coords
	int prev_node1_x = prev_line.first / dim_y;
	int prev_node1_y = prev_line.first % dim_y;
	int prev_node2_x = prev_line.second / dim_y;
	int prev_node2_y = prev_line.second % dim_y;

	// Get next node coords
	int next_node1_x = next_line.first / dim_y;
	int next_node1_y = next_line.first % dim_y;
	int next_node2_x = next_line.second / dim_y;
	int next_node2_y = next_line.second % dim_y;

	// Get cell coords
	int candidate_cell1_x = candidate_cells.first / (dim_y - 1);
	int candidate_cell1_y = candidate_cells.first % (dim_y - 1);
	int candidate_cell2_x = candidate_cells.second / (dim_y - 1);
	int candidate_cell2_y = candidate_cells.second % (dim_y - 1);
	int prev_cell_x = prev_cell / (dim_y - 1);
	int prev_cell_y = prev_cell % (dim_y - 1);

	// Choose between candidate bound cells
	bool prevline_vertical = prev_node1_x == prev_node2_x;
	bool nextline_vertical = next_node1_x == next_node2_x;
	bool choose_candidate1 = true;
	if (prevline_vertical == nextline_vertical) {
		if (prevline_vertical) {
			// Continuous vertical line; next cell should have the same x-coord as previous cell
			if (candidate_cell2_x == prev_cell_x) choose_candidate1 = false;
		}
		else {
			// Continuous horizontal line; next cell should have the same y-coord as previous cell
			if (candidate_cell2_y == prev_cell_y) choose_candidate1 = false;
		}
	}
	else {
		// If a direction change has occurred (vertical to horizontal or vice-versa), the next cell has either a different x-
		// AND y-coordinates wrt the previous cell (since it is a diagonal neighbor), or has the same coordinates (because it
		// is enclosed by the two line segments).
		
		// Check if one of the candidate cells has different x AND y coordinates wrt the previous cell
		if (candidate_cell2_x != prev_cell_x && candidate_cell2_y != prev_cell_y) {
			choose_candidate1 = false;
		}
		else if (candidate_cell2_x == prev_cell_x && candidate_cell2_y == prev_cell_y) {
			choose_candidate1 = false;
		}
	}

	// Derive cell coordinates
	int boundcell_x, boundcell_y;
	int keepcell_x, keepcell_y;
	int cutoutcell_x, cutoutcell_y;
	if (choose_candidate1) {
		boundcell_x = candidate_cell1_x; boundcell_y = candidate_cell1_y;
		cutoutcell_x = candidate_cell2_x; cutoutcell_y = candidate_cell2_y;
	}

	int boundcell = boundcell_x * (dim_y - 1) + boundcell_y;
	int cutoutcell = cutoutcell_x * (dim_y - 1) + cutoutcell_y;

	return { boundcell, cutoutcell };
}

vector<int> phys::FEACaseManager::get_additional_bound_cells(vector<int>* origin_cells, vector<int>* cells_to_avoid) {
	vector<int> additional_cells;
	vector<pair<int, int>> offsets = { pair(0,0), pair(1,0), pair(1,1), pair(1,0) };
	for (auto& bound_cell : *origin_cells) {
		Vector2d coords = get_node_coords(bound_cell);
		for (auto& offset : offsets) {
			int x = (int)coords[0] + offset.first;
			int y = (int)coords[1] + offset.second;
			int candidate = x * (dim_y - 1) + y;
			if (!help::is_in(cells_to_avoid, candidate) && !help::is_in(&additional_cells, candidate) &&
				!help::is_in(origin_cells, candidate))
				additional_cells.push_back(candidate);
		}
	}
	return additional_cells;
}

void phys::FEACaseManager::walk_and_collect_bound_cells(
	FEACase* active_case, map<string, vector<int>>* bound_cells, vector<int>* visited_nodes, pair<int, int> first_line,
	int first_bound_cell, int neighbor_node, string bound_name
) {
	pair<int, int> previous_line = first_line;
	int boundcell = first_bound_cell;
	int node = neighbor_node;
	while (true) {
		vector<int> path_to_neighbor = get_path_to_unvisited_node_neighbor(
			node, &active_case->bound_cond_nodes[bound_name], visited_nodes
		);
		if (path_to_neighbor.size()) break; // If there's no neighbor, we've reached the end of the walk.
		for (int i = 1; i < path_to_neighbor.size(); i++) {
			pair<int, int> next_line = pair(path_to_neighbor[i-1], path_to_neighbor[i]);
			pair<int, int> neighbor_cells = get_cells_with_given_edge(next_line);
			auto [_boundcell, cutoutcell] = get_bound_cells(boundcell, neighbor_cells, previous_line, next_line);
			if (_boundcell == -1) continue;
			if (!help::is_in(&bound_cells->at("bound"), _boundcell)) bound_cells->at("bound").push_back(_boundcell);
			if (!help::is_in(&bound_cells->at("cutout"), _boundcell)) bound_cells->at("cutout").push_back(cutoutcell);
			node = path_to_neighbor.back();
			boundcell = _boundcell;
			previous_line = next_line;
		}
	}
}

void phys::FEACaseManager::interpolate_cells(FEACase* source, FEACase* target, FEACase* active_case, float fraction) {
	for (auto& [bound_name, bound_cells] : source->bound_cond_cells) {
		// Get interpolated boundcell- and node barycenters
		Vector2d cell_barycent_diff = target->cell_barycenters[bound_name]["bound"] - source->cell_barycenters[bound_name]["bound"];
		Vector2d cell_barycent = source->cell_barycenters[bound_name]["bound"] + fraction * cell_barycent_diff;
		Vector2d node_barycent_diff = target->node_barycenters[bound_name] - source->node_barycenters[bound_name];
		Vector2d node_barycent = source->node_barycenters[bound_name] + fraction * node_barycent_diff;

		// Get vector from nodes' barycenter to cells' barycenter
		Vector2d barycent_difference = cell_barycent - node_barycent;

		// Get interpolated node closest to node barycenter
		auto [start_node, _] = get_nearest_neighbor(node_barycent, &active_case->bound_cond_nodes[bound_name]);

		// Get neighbor node
		vector<int> visited_nodes = { start_node };
		vector<int> path_to_neighbor = get_path_to_unvisited_node_neighbor(
			start_node, &active_case->bound_cond_nodes[bound_name], &visited_nodes);
		int neighbor_node = path_to_neighbor.back();

		// Check whether formed line is vertical or horizontal
		Vector2d nearest_coords = get_node_coords(start_node);
		Vector2d neighbor_coords = get_node_coords(neighbor_node);
		bool vertical = nearest_coords[0] == neighbor_coords[0];
		int cell_x, cell_y;

		// If formed line is vertical: x-direction of 'barycenter vector' determines side of line on which boundary cell lies
		if (vertical) {
			cell_y = min(nearest_coords[1], neighbor_coords[1]);
			if (barycent_difference[0] < 0) {
				cell_x = nearest_coords[0] - 1;
			}
			else cell_x = nearest_coords[0];
		}
		// Else: y-direction determines side
		else {
			cell_x = min(nearest_coords[0], neighbor_coords[0]);
			if (barycent_difference[1] < 0) {
				cell_y = nearest_coords[1] - 1;
			}
			else cell_y = nearest_coords[1];
		}
		int first_bound_cell = cell_x * (dim_y - 1) + cell_y;
		cout << "first bound cell: " << cell_x << ", " << cell_y << endl;

		// Starting from found boundary cell, walk along boundary condition lines till endpoint. Meanwhile, collect boundary cells.
		map<string, vector<int>> bound_cells;
		bound_cells["bound"] = { first_bound_cell };
		visited_nodes.push_back(neighbor_node);
		walk_and_collect_bound_cells(
			active_case, &bound_cells, &visited_nodes, pair(start_node, neighbor_node),
			first_bound_cell, neighbor_node, bound_name
		);

		// Re-start walk from first cell, walking in opposite direction
		path_to_neighbor = get_path_to_unvisited_node_neighbor(
			start_node, &active_case->bound_cond_nodes[bound_name], &visited_nodes
		);
		if (path_to_neighbor.size() > 0) {
			walk_and_collect_bound_cells(
				active_case, &bound_cells, &visited_nodes, pair(start_node, path_to_neighbor.back()),
				first_bound_cell, path_to_neighbor.back(), bound_name
			);
		}

		// Get keep cells and additional cutout cells
		vector<int> keep_cells = get_additional_bound_cells(&bound_cells["bound"], &bound_cells["cutout"]);
		bound_cells["keep"] = keep_cells;
		vector<int> extra_cutout_cells = get_additional_bound_cells(&bound_cells["cutout"], &bound_cells["bound"]);

		// Store found boundary cells
		help::append_vector(bound_cells["cutout"], &extra_cutout_cells);
		help::append_vector(cells_to_keep, &bound_cells["keep"]);
		help::append_vector(cells_to_keep, &bound_cells["bound"]);
		help::append_vector(cutout_cells, &bound_cells["cutout"]);		
		active_case->bound_cond_cells[bound_name] = bound_cells;
		
		// Temp
		for (auto& [type, cells] : bound_cells) {
			cout << "\nbound name: " << bound_name << ", type: " << type << endl;
			cout << "cells: " << cells.size() << endl;
		}
	}
}

void phys::FEACaseManager::interpolate(float fraction) {
	// Reset
	cells_to_keep.clear();
	cutout_cells.clear();
	for (int i = 0; i < active_cases.size(); i++) active_cases[i].bound_cond_cells.clear();

	// Recompute
	for (int i = 0; i < active_cases.size(); i++) {
		max_stress_threshold = sources[i].max_stress_threshold +
			fraction * (targets[i].max_stress_threshold - sources[i].max_stress_threshold
		);
		interpolate_nodes(&sources[i], &targets[i], fraction, i);
		interpolate_cells(&sources[i], &targets[i], &active_cases[i], fraction);
	}
}

void phys::FEACaseManager::initialize() {
	active_cases = sources;
	max_stress_threshold = sources[0].max_stress_threshold;
	if (dynamic) {
		for (int i = 0; i < sources.size(); i++) {
			sources[i].compute_barycenters();
			targets[i].compute_barycenters();
			compute_migration_vectors(&sources[i], &targets[i], i);
		}
	}
}

void phys::FEACaseManager::update_casepaths(string case_folder) {
	for (auto& fea_case : active_cases) {
		fea_case.path = case_folder + "/" + fea_case.name + ".sif";
	}
}

