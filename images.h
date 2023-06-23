#pragma once
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "meshing.h"

#define STB_IMAGE_IMPLEMENTATION
#include "lib/stb_image.h"


using namespace std;
using namespace fessga;


namespace fessga{

	class Image {
	public:
		unsigned char* values;
		int width, height, channels;
		float pixels_per_cell;
	};

	class img {
	public:
		static Image load_image(char* filename, mesher::Grid3D grid) {
			Image image;
			image.values = stbi_load(filename, &image.width, &image.height, &image.channels, 1);
			if (image.values == NULL) {
				printf("Error in loading the image\n");
			}
			printf("Loaded image with a width of %dpx, a height of %dpx and %d channels\n", image.width, image.height, image.channels);

			image.pixels_per_cell = (float)(image.width * image.height) / (float)(grid.x * grid.y);

			cout << "pixels per cell: " << image.pixels_per_cell << endl;
			
			return image;
		}

		//static void create_distribution_from_image(char* filename, uint* densities, mesher::Grid3D grid, Image img) {

		//	Image image;
		//	image.values = stbi_load(filename, &image.width, &image.height, &image.channels, 1);
		//	if (image.values == NULL) {
		//		printf("Error in loading the image\n");
		//	}
		//	printf("Loaded image with a width of %dpx, a height of %dpx and %d channels\n", image.width, image.height, image.channels);

		//	image.pixels_per_cell = (float)(image.width * image.height) / (float)(grid.x * grid.y);

		//	cout << "pixels per cell: " << image.pixels_per_cell << endl;

		//	size_t img_size = img.width * img.height * img.channels;
		//	size_t gray_img_size = img.width * img.height;
		//	unsigned char* gray_img = new unsigned char[gray_img_size];
		//	for (unsigned char* p = img.values; p < img.values + img_size - 999000; p += img.channels, gray_img + 1) {
		//		*gray_img = (uint8_t)((*p + *(p + 1) + *(p + 2)) / 3.0);
		//	}
		//	for (int x = 0; x < grid.x; x++) {
		//		for (int y = 0; y < grid.y; y++) {
		//			densities[x * grid.y + y] = (uint)((float)((x * grid.y) * img.pixels_per_cell + (float)y * img.pixels_per_cell));
		//		}
		//	}
		//	cout << "total no pixels: " << img.width * img.height << endl;
		//}

		static void create_distribution_from_image(const char* filename, uint* densities, mesher::Grid3D grid) {
			int width, height, numChannels;
			unsigned char* image = stbi_load(filename, &width, &height, &numChannels, 0);

			if (!image) {
				std::cerr << "Failed to load image." << std::endl;
			}

			// Assuming the image has 3 or 4 channels (RGB or RGBA)
			int numPixels = width * height;
			std::vector<int> redValues(numPixels);

			for (int i = 0; i < numPixels; ++i) {
				int index = i * numChannels;
				redValues[i] = image[index];
			}

			stbi_image_free(image);

			cout << "num pixels: " << numPixels << endl;
			cout << "num cells: " << grid.size2d << endl;

			for (int y = 0; y < grid.y; y++) {
				for (int x = 0; x < grid.x; x++) {
					densities[(x + 1) * grid.y - y - 1] = redValues[
						y * width * (height / grid.y) + x * (height / grid.y) + 0.5 * width * (height / grid.y)
					] / 255;
				}
			}
		}

		static void load_distribution_from_image(char* filename, uint* densities, mesher::Grid3D grid) {
			Image image = load_image(filename, grid);
			//create_distribution_from_image(filename, densities, grid, image);
			create_distribution_from_image(filename, densities, grid);
			stbi_image_free(image.values);
		}
	};
}
