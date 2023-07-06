#pragma once
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "meshing.h"

#define STB_IMAGE_IMPLEMENTATION
#include "lib/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "lib/stb_image_write.h"


using namespace std;
using namespace fessga;


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

		static void load_distribution_from_image(grd::Densities2d& densities, msh::SurfaceMesh& mesh, const char* filename) {
			if (!IO::file_exists(string(filename))) throw std::runtime_error("Error: image file " + string(filename) + " does not exist");
			
			// Load image
			int width, height, numChannels;
			unsigned char* image = stbi_load(filename, &width, &height, &numChannels, 0);
			if (!image) {
				std::cerr << "Failed to load image." << std::endl;
			}

			// Get size of design domain
			int pixels_per_cell = width / densities.dim_x; // The width is leading in determining the size of the cells
			densities.dim_y = height / pixels_per_cell;
			float mesh_width = mesh.diagonal.maxCoeff();
			float width_per_cell = mesh_width / densities.dim_x;
			Vector3d mesh_size = Vector3d(mesh_width, width_per_cell * densities.dim_y, 0);
			mesh = msh::SurfaceMesh(mesh_size);
			densities = grd::Densities2d(densities.dim_x, mesh.diagonal, densities.output_folder);

			// Extract red values
			int numPixels = width * height;
			int* redValues = new int[numPixels];
			for (int i = 0; i < numPixels; ++i) {
				int index = i * numChannels;
				redValues[i] = image[index];
			}

			// Extract green values
			int* greenValues = new int[numPixels];
			for (int i = 0; i < numPixels; ++i) {
				int index = i * numChannels + 1;
				greenValues[i] = image[index];
			}
			stbi_image_free(image);

			// Sample green values with intervals to derive which cells should be filled and which should be void
			int x_offset = (densities.dim_x / 2);
			x_offset = 0;
			for (int y = 0; y < densities.dim_y; y++) {
				for (int x = 0; x < densities.dim_x; x++) {
					int img_idx = y * width * pixels_per_cell + (x - x_offset) * pixels_per_cell;
					int pixel_count = 0;
					// Count all pixel values in the cell
					for (int i = 0; i < pixels_per_cell * pixels_per_cell; i++) {
						int _y = i / pixels_per_cell;
						int _x = i % pixels_per_cell;
						int _img_idx = img_idx + _y * width + _x;
						pixel_count += greenValues[_img_idx];
					}
					int cell_value = round((float)pixel_count / (float)(255 * (pixels_per_cell * pixels_per_cell)));
					densities.set((x + 1) * densities.dim_y - y - 1, cell_value);
				}
			}
			densities.update_count();
			densities.do_feasibility_filtering();

			// Sample red values with intervals to derive which cells are cutout cells that are not allowed to become filled
			for (int y = 0; y < densities.dim_y; y++) {
				for (int x = 0; x < densities.dim_x; x++) {
					int img_idx = y * width * pixels_per_cell + (x - x_offset) * pixels_per_cell;
					int pixel_count = 0;
					// Count all pixel values in the cell
					for (int i = 0; i < pixels_per_cell * pixels_per_cell; i++) {
						int _y = i / pixels_per_cell;
						int _x = i % pixels_per_cell;
						int _img_idx = img_idx + _y * width + _x;
						pixel_count += redValues[_img_idx];
					}
					int cell_value = round((float)pixel_count / (float)(255 * (pixels_per_cell * pixels_per_cell)));
					bool is_cutout = cell_value;
					if (is_cutout) densities.fea_case.cutout_cells.push_back((x + 1) * densities.dim_y - y - 1);
				}
			}
		}

		static void convert_distribution_to_single_channel_image(grd::Densities2d densities, unsigned char* single_channel, Image* image) {
			for (int y = 0; y < densities.dim_y; y++) {
				for (int x = 0; x < densities.dim_x; x++) {
					int dens_idx = (x + 1) * densities.dim_y - y - 1;
					//int img_idx = (y) * image->width * pixels_per_cell + (x) * pixels_per_cell;
					//single_channel[img_idx] = densities[dens_idx] * 255;
					for (int i = 0; i < image->pixels_per_cell; i++) {
						for (int j = 0; j < image->pixels_per_cell; j++) {
							if (y == image->height - 1) continue;
							int img_idx = (y * image->width * image->pixels_per_cell + i * image->width) + (x * image->pixels_per_cell - j);
							single_channel[img_idx] = densities[dens_idx] * 255;
						}
					}
				}
			}
		}

		static void singlechannel_to_rgb(unsigned char* single_channel, Image* image) {
			// Convert single-channel image to RGB image by copying each value three times
			for (int i = 0; i < image->size; i += 3) {
				image->vals[i] = single_channel[i / 3];
				image->vals[i + 1] = single_channel[i / 3];
				image->vals[i + 2] = single_channel[i / 3];
			}
		}

		static void write_distribution_to_image(grd::Densities2d densities, string path, int x, int y, bool verbose = false) {
			img::Image image = img::Image(path, 1000, 1000, 3, densities);
			unsigned char* single_channel = new unsigned char[image.width * image.height];
			convert_distribution_to_single_channel_image(densities, single_channel, &image);
			singlechannel_to_rgb(single_channel, &image);
			stbi_write_jpg(image.path, image.width, image.height, image.channels, image.vals, 100);
			if (verbose) cout << "Exported distribution to image. Path: " << image.path << endl;
		}
	};
}
