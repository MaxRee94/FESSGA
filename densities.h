#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include "helpers.h"

using namespace Eigen;
using namespace std;


class Densities2D {
public:
    Densities2D(int _dim_x, int _dim_y, Vector2d diagonal) {
        construct_grid(_dim_x, _dim_y);
        cell_size = diagonal.cwiseProduct(
            Vector2d(1.0 / (double)dim_x, 1.0 / (double)dim_y)
        );
    }
    Densities2D(int _dim_x, int _dim_y, Vector3d diagonal) {
        construct_grid(_dim_x, _dim_y);
        Vector2d _diagonal2d = Vector2d(diagonal(0), diagonal(1));
        cell_size = _diagonal2d.cwiseProduct(
            Vector2d(1.0 / (double)dim_x, 1.0 / (double)dim_y)
        );
    }
    void construct_grid(int _dim_x, int _dim_y) {
        dim_x = _dim_x;
        dim_y = _dim_y;
        size = dim_x * dim_y;
        count = 0;
        values = new uint[size];
    }
    int get_idx(int x, int y) {
        return x * dim_y + y;
    }
    pair<int, int> get_coords(int idx) {
        pair<int, int> coords(idx / dim_y, idx % dim_y);
        return coords;
    }
    uint operator[](int idx) {
        return values[idx];
    }
    uint operator()(int x, int y) {
        return values[get_idx(x, y)];
    }
    void fill(int idx) {
        if (values[idx] != 1) {
            values[idx] = 1;
            count++;
        }
    }
    void fill(vector<int> indices) {
        for (auto& idx : indices) fill(idx);
    }
    void del(int idx) {
        if (values[idx] != 0) {
            values[idx] = 0;
            count--;
        }
    }
    void del(vector<int> indices) {
        for (auto& idx : indices) del(idx);
    }
    void setzero() {
        for (int i = 0; i < size; i++) values[i] = 0;
    }
    void setone() {
        for (int i = 0; i < size; i++) values[i] = 1;
    }
    void print() {
        for (int y = dim_y - 1; y > -1; y--) {
            for (int x = 0; x < dim_x; x++) {
                cout << values[get_idx(x, y)];
            }
            cout << endl;
        }
    }
    int count = 0;
    int dim_x = 0;
    int dim_y = 0;
    int size = 0;
    Vector2d cell_size;
private:
    uint* values = 0;
};

