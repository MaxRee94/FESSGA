#pragma once
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "meshing.h"



namespace fessga{
	class img {
	public:

		class Image {
		public:
			Image() = default;
			Image(string _path, int _width, int _height, int _channels, grd::Densities2d& _densities) {
				strcpy(path, (_path).c_str());
				width = _width; height = _height; channels = _channels;
				size = width * height * channels;
				densities = _densities;
				pixels_per_cell = height / densities.dim_y;
				vals = new unsigned char[width * height * channels];
			}
			int width, height, channels, size, pixels_per_cell;
			unsigned char* vals = 0;
			grd::Densities2d densities;
			char path[300];
		};

		static void load_distribution_from_image(grd::Densities2d& densities, msh::SurfaceMesh& mesh, const char* filename);
		static void convert_distribution_to_single_channel_image(grd::Densities2d densities, unsigned char* single_channel, Image* image);
		static void singlechannel_to_rgb(unsigned char* single_channel, Image* image);
		static void write_distribution_to_image(grd::Densities2d densities, string path, int x, int y, bool verbose = false);
	};
}
