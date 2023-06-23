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
			Image(char* _path, int _width, int _height, int _channels, mesher::Grid3D grid) {
				path = _path;
				width = _width; height = _height; channels = _channels;
				size = width * height * channels;
				pixels_per_cell = height / grid.y;
				vals = new unsigned char[width * height * channels];
			}
			int width, height, channels, size, pixels_per_cell;
			unsigned char* vals = 0;
			char* path;
		};

		static void load_distribution_from_image(const char* filename, uint* densities, mesher::Grid3D grid) {
			// Load image
			int width, height, numChannels;
			unsigned char* image = stbi_load(filename, &width, &height, &numChannels, 0);
			if (!image) {
				std::cerr << "Failed to load image." << std::endl;
			}

			// Extract red values
			int numPixels = width * height;
			int pixels_per_cell = height / grid.y;
			int* redValues = new int[numPixels];
			for (int i = 0; i < numPixels; ++i) {
				int index = i * numChannels;
				redValues[i] = image[index];
			}
			stbi_image_free(image);

			// Sample red values with intervals
			for (int y = 0; y < grid.y; y++) {
				for (int x = 0; x < grid.x; x++) {
					int img_idx = y * width * pixels_per_cell + x * pixels_per_cell + (int)(0.5 * width * pixels_per_cell);
					densities[(x + 1) * grid.y - y - 1] = redValues[img_idx] / 255;
				}
			}
		}

		static void write_distribution_to_image(uint* densities, mesher::Grid3D grid, Image image) {
			// Convert distribution to single-channel image
			unsigned char* single_channel = new unsigned char[image.width * image.height];
			for (int y = 0; y < grid.y; y++) {
				for (int x = 0; x < grid.x; x++) {
					int dens_idx = (x + 1) * grid.y - y - 1;
					//int img_idx = (y) * image.width * pixels_per_cell + (x) * pixels_per_cell;
					//single_channel[img_idx] = densities[dens_idx] * 255;
					for (int i = 0; i < image.pixels_per_cell; i++) {
						for (int j = 0; j < image.pixels_per_cell; j++) {
							if (y == image.height - 1) continue;
							int img_idx = (y * image.width * image.pixels_per_cell + i * image.width) + (x * image.pixels_per_cell - j);
							single_channel[img_idx] = densities[dens_idx] * 255;
						}
					}
				}
			}

			// Convert single-channel image to RGB image by copying each value three times
			for (int i = 0; i < image.size; i+= 3) {
				image.vals[i] = single_channel[i / 3];
				image.vals[i + 1] = single_channel[i / 3];
				image.vals[i + 2] = single_channel[i / 3];
			}

			// Export image to file
			stbi_write_jpg(image.path, image.width, image.height, image.channels, image.vals, 100);
		}
	};
}
