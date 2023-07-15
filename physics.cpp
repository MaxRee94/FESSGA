#pragma once
#include "physics.h"


using namespace fessga;


Vector2d phys::FEACase::get_node_coords(int idx) {
	int x = idx / dim_y;
	int y = idx % dim_y;
	return Vector2d((double)x, (double)y);
}

void phys::FEACase::compute_barycenters() {
	for (auto& [bound_name, lines] : bound_conds) {
		vector<int> _nodes;
		Vector2d coordinate_sum(0, 0);
		for (auto& line : lines) {
			if (!help::is_in(&_nodes, line.first)) {
				int x = line.first / dim_y;
				int y = line.first % dim_y;
				coordinate_sum += Vector2d((double)x, (double)y);
			}
			if (!help::is_in(&_nodes, line.second)) {
				int x = line.second / dim_y;
				int y = line.second % dim_y;
				coordinate_sum += Vector2d((double)x, (double)y);
			}
		}
		barycenters[bound_name] = coordinate_sum / (double)(_nodes.size());
	}
}

// Resample boundary by adding copies of nodes, starting from the barycenter outwards
void phys::FEACaseInterpolator::resample_boundary(vector<int>& nodes, Vector2d barycenter, int target_no_nodes) {
	for (int node : nodes) {
		if (nodes.size() == target_no_nodes) break;
		nodes.push_back(node);
	}
	order_nodes_by_distance(&nodes, barycenter);
	if (nodes.size() < target_no_nodes) resample_boundary(nodes, barycenter, target_no_nodes);
}

// Get the index of the nearest node on the grid relative to the given node coordinates
int phys::FEACaseInterpolator::get_nearest_neighbor(Vector2d coords, vector<int>* _neighbors) {
	double min_dist = INFINITY;
	int neighbor_idx;

	vector<int> neighbors;
	vector<pair<int, int>> offsets = { pair(0,0), pair(1,0), pair(1,1), pair(1,0) };
	if (!_neighbors) {
		for (auto& offset : offsets) {
			int x = (int)coords[0] + offset.first;
			int y = (int)coords[1] + offset.second;
			neighbors.push_back(x * source.dim_y + y);
		}
	}
	else neighbors = *_neighbors;

	for (int trgt_idx = 0; trgt_idx < neighbors.size(); trgt_idx++) {
		Vector2d trgt_pnt = source.get_node_coords(trgt_idx);
		double dist = (trgt_pnt - coords).norm();
		if (dist < min_dist) {
			min_dist = dist;
			neighbor_idx = trgt_idx;
		}
	}

	return neighbor_idx;
}

// Get vector to the nearest neighbor of the given node
Vector2d phys::FEACaseInterpolator::get_vector2nn(Vector2d coords, vector<int>& neighbors) {
	int nn_idx = get_nearest_neighbor(coords, &neighbors);
	Vector2d vector2nn = source.get_node_coords(nn_idx) - coords;
	neighbors.erase(neighbors.begin() + nn_idx);

	return vector2nn;
}

void phys::FEACaseInterpolator::order_nodes_by_distance(vector<int>* nodes, Vector2d barycenter) {
	map<int, double> distances;
	for (auto& _node : *nodes) {
		Vector2d node = source.get_node_coords(_node);
		double distance = (barycenter - node).norm();
		distances[_node] = distance;
	}
	PairSet _ordered_nodes;
	help::sort(distances, _ordered_nodes);
	nodes->clear();
	for (auto& [pnt, _] : _ordered_nodes) nodes->push_back(pnt);
}

void phys::FEACaseInterpolator::compute_morph_vectors() {
	source.compute_barycenters();
	target.compute_barycenters();
	for (auto& [bound_name, lines] : source.bound_conds) {
		// Compute barycenter alignment vector
		Vector2d align_to_barycenter = target.barycenters[bound_name] - target.barycenters[bound_name];

		// Convert lines vectors to vectors of nodes
		vector<int> target_nodes;
		for (auto& line : target.bound_conds[bound_name]) {
			if (!help::is_in(&target_nodes, line.first)) target_nodes.push_back(line.first);
			if (!help::is_in(&target_nodes, line.second)) target_nodes.push_back(line.second);
		}
		vector<int> source_nodes;
		for (auto& line : lines) {
			if (!help::is_in(&source_nodes, line.first)) source_nodes.push_back(line.first);
			if (!help::is_in(&source_nodes, line.second)) source_nodes.push_back(line.second);
		}

		// Order source nodes by distance to barycenter
		order_nodes_by_distance(&source_nodes, source.barycenters[bound_name]);

		// Resample if the number of nodes of the boundary condition differs between source and target
		if (target_nodes.size() > source_nodes.size()) resample_boundary(source_nodes, source.barycenters[bound_name], target_nodes.size());
		else if (source_nodes.size() > target_nodes.size()) {
			order_nodes_by_distance(&target_nodes, source.barycenters[bound_name]);
			resample_boundary(target_nodes, target.barycenters[bound_name], source_nodes.size());
		}

		// Compute vectors to align nodes
		vector<Vector2d> _morph_vectors;
		for (auto& src_pnt : source_nodes) {
			Vector2d source_node = source.get_node_coords(src_pnt);
			source_node += align_to_barycenter; // Align the cluster of nodes to the barycenter of the target cluster
			Vector2d to_nn = get_vector2nn(source_node, target_nodes);
		}

		// Store bound nodes and morph vectors
		source.bound_nodes[bound_name] = source_nodes;
		target.bound_nodes[bound_name] = target_nodes;
		morph_vectors[bound_name] = _morph_vectors;
	}
}

void phys::FEACaseInterpolator::morph_nodes(float fraction) {
	interpolated.bound_nodes.clear();
	for (auto& [bound_name, bound_nodes] : source.bound_nodes) {
		// Get interpolated boundary nodes
		vector<int> interpolated_bound_nodes;
		for (int i = 0; i < bound_nodes.size(); i++) {
			Vector2d source_coords = source.get_node_coords(bound_nodes[i]);
			Vector2d morph_vec = fraction * morph_vectors[bound_name][i];
			int interpolated_node = get_nearest_neighbor(source_coords + morph_vec);
			interpolated_bound_nodes.push_back(interpolated_node);
		}

		// Reconstruct boundary lines from node list
		vector<pair<int, int>> lines;
		int i = 0;
		while (i < interpolated_bound_nodes.size() - 1) {
			pair<int, int> line = pair(interpolated_bound_nodes[i], interpolated_bound_nodes[i + 1]);
			lines.push_back(line);
			i++;
		}
		interpolated.bound_conds[bound_name] = lines;
		interpolated.bound_nodes[bound_name] = interpolated_bound_nodes;
	}
}

void phys::FEACaseInterpolator::interpolate(float fraction) {
	interpolated = source;
	interpolated.max_stress_threshold += (target.max_stress_threshold - source.max_stress_threshold) / 2.0;
	morph_nodes(fraction);
}
