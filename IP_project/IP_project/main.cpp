/*=========================================================================

File:		IP project - Image segmentation & registration code
Language:	C++ & C
Library:	Standard C++ Library
Author:		Yongjin Jeon
Date:		2019-02-20
Version:	1.0
Mail:		yongjin117@naver.com

=========================================================================*/
#include "../Common/memory.h"
#include "../Common/image3d.h"
#include "../Core/raw_io.h"
#include "../Core/raw_io_exception.h"
#include <memory>
#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <limits>
#include <fstream>
#include <ctime>
#include <omp.h>
#include <cmath>

short** pad_image(short ** image_array);
void store_result(short ** image_array);
short** line_detection(short ** image_array);
void minmax_scaling(short ** image_array);
void laplacian_filtering(short** image_array);
short** thresholding(short ** image_array);
short** _3D_CCA(short ** image_array);
void increase_arr_size(int s, short**arr);
void find_min_max(short * arr, short**minmax);

// Fuction for Changing Endian type (Big endian -> Little endian)
template <class T>
std::unique_ptr<mc::image3d<T>> load_image(const std::string& path, const unsigned int w, const unsigned int h, const unsigned int d)
{
	using namespace mc;

	image3d<T>* pImg = nullptr;

	try {
		raw_io<T> io(path.c_str());
		io.setEndianType(raw_io<T>::EENDIAN_TYPE::BIG);
		pImg = io.read(w, h, d);
	}
	catch (raw_io_exception& e) {
		std::cerr << e.what() << std::endl;
	}

	return std::make_unique<image3d<T>>(*pImg);
}

int main()
{
	// Read moving & target images for registration.
	std::unique_ptr<mc::image3d<short>> img1 = load_image<short>("volume1_512x512x56.raw", 512, 512, 56);
	std::unique_ptr<mc::image3d<short>> img2 = load_image<short>("volume2_512x512x58.raw", 512, 512, 58);
	
	// Get the data from loaded image
	//short** image_array = img1->data();
	short** image_array = img2->data();
	// Laplacian filtering with modified laplacian filters
	//laplacian_filtering(image_array);

	// TODO #1 : segmentation of lung region. (use thresholding, CCA)
	short** thresholded = thresholding(image_array);
	_3D_CCA(thresholded);

	// TODO #1-1 : Initial transformation parameter calculation.

	// TODO #2 : edge extraction for both images.

	// TODO #3 : Distance transformation.

	// TODO #4 : Perform iterative REGISTRATION.

	// TODO #5 : Transform moving (floating) image with estimated transformation parameter & generate subtraction image.

	// TODO #6 : store subtraction image (visual purpose).

	// TODO :
	// perform Intensity-based registration (using similarity measure metric with original intensities).
	// Do #[4, 5, 6] and compare with surface-based method.

	return EXIT_SUCCESS;
}

void laplacian_filtering(short ** image_array) {

	// Zero padding
	short ** padding_image = pad_image(image_array);

	/*
	TODO #0 : Semi-isotropic image generation
	Filter 1 - Horizontal			Filter 2 - +45			Filter 3 - Vertical			Filter 4 - -45
	-1  -1  -1				2  -1  -1				-1  2  -1					-1  -1   2
	2   2   2			   -1   2  -1				-1  2  -1					-1   2  -1
	-1  -1  -1			   -1  -1   2				-1  2  -1					 2  -1  -1
	*/
	short** filtered_image = line_detection(padding_image);

	// Store the filtered image
	store_result(filtered_image);

	// Memory free
	free(padding_image);
	free(filtered_image);
}

void minmax_scaling(short ** image_array) {
	// For min max scaling
	short min = 0;
	short max = 0;
	for (int i = 0; i < 56; i++) {
		for (int j = 0; j < 512 * 512; j++) {
			
			// check the minimum value
			if (image_array[i][j] < min) {
				std::cout << "min value changed from " << min << " to " << image_array[i][j] << std::endl;
				min = image_array[i][j];
			}

			// check the maximum value
			if (image_array[i][j] > max) {
				std::cout << "max value changed from " << max << " to " << image_array[i][j] << std::endl;
				max = image_array[i][j];
			}
		}
	}
	std::cout << "min of image array : " << min << std::endl;
	std::cout << "max of image array : " << max << std::endl;
}

short** pad_image(short ** image_array) {
	short ** pad_array = (short**)malloc(sizeof(short*)*56);
	for (int z = 0; z < 56; z++) {
		int original_idx = 0;

		// allocating memory of padded size
		pad_array[z] = (short*)malloc(sizeof(short) * 514 * 514);
		for (int idx = 0; idx < 514 * 514; idx++) {

			// zero padding for the first row and last row
			if (idx < 514 || idx > 263682) {
				pad_array[z][idx] = 0;
			}

			// zero padding at the left and right edge column
			else if ((idx % 514 == 0) || (idx % 514) == 513)
			{
				pad_array[z][idx] = 0;
			}

			// except the top, bottom, left, right edge set the value with the same order of the original image
			else {
				pad_array[z][idx] = image_array[z][original_idx];
				original_idx++;
			}
		}
	}
	return pad_array;
}


short** line_detection(short ** image_array) {

	// Filter for filtering the horizontal, vertical, 45 degrees, 315 degrees, isotropic properties of the image features
	int _H, _V, _45, _315, la, total, processed_idx;
	short ** processed_array = (short**)malloc(sizeof(short*)*56);
		
	for (int z = 0; z < 56; z++) {
		processed_idx = 0;
		_H, _V, _45, _315, total = 0;
		processed_array[z] = (short*)malloc(sizeof(short) * 512 * 512);

		// Filtering with targetting at the center of the neighborhood region
		for (int idx = 0; idx < 514 * 514 - 1028; idx++) {
			if ((idx % 514 == 512) || (idx % 514 == 513)) continue;
			else {
				_H = image_array[z][idx] * (-1)
					+ image_array[z][idx + 1] * (-1)
					+ image_array[z][idx + 2] * (-1)
					+ image_array[z][idx + 514] * 2
					+ image_array[z][idx + 515] * 2
					+ image_array[z][idx + 516] * 2
					+ image_array[z][idx + (514 * 2)] * (-1)
					+ image_array[z][idx + (514 * 2) + 1] * (-1)
					+ image_array[z][idx + (514 * 2) + 2] * (-1);

				_V = image_array[z][idx] * (-1)
					+ image_array[z][idx + 1] * 2
					+ image_array[z][idx + 2] * (-1)
					+ image_array[z][idx + 514] * (-1)
					+ image_array[z][idx + 515] * 2
					+ image_array[z][idx + 516] * (-1)
					+ image_array[z][idx + (514 * 2)] * (-1)
					+ image_array[z][idx + (514 * 2) + 1] * 2
					+ image_array[z][idx + (514 * 2) + 2] * (-1);

				_45 = image_array[z][idx] * 2
					+ image_array[z][idx + 1] * (-1)
					+ image_array[z][idx + 2] * (-1)
					+ image_array[z][idx + 514] * (-1)
					+ image_array[z][idx + 515] * 2
					+ image_array[z][idx + 516] * (-1)
					+ image_array[z][idx + (514 * 2)] * (-1)
					+ image_array[z][idx + (514 * 2) + 1] * (-1)
					+ image_array[z][idx + (514 * 2) + 2] * 2;

				_315 = image_array[z][idx] * (-1)
					+ image_array[z][idx + 1] * (-1)
					+ image_array[z][idx + 2] * 2
					+ image_array[z][idx + 514] * (-1)
					+ image_array[z][idx + 515] * 2
					+ image_array[z][idx + 516] * (-1)
					+ image_array[z][idx + (514 * 2)] * 2
					+ image_array[z][idx + (514 * 2) + 1] * (-1)
					+ image_array[z][idx + (514 * 2) + 2] * (-1);

				/*
				la = image_array[z][idx] 
					+ image_array[z][idx + 1] 
					+ image_array[z][idx + 2] 
					+ image_array[z][idx + 514] 
					+ image_array[z][idx + 515] * (-8)
					+ image_array[z][idx + 516]
					+ image_array[z][idx + (514 * 2)]
					+ image_array[z][idx + (514 * 2) + 1] 
					+ image_array[z][idx + (514 * 2) + 2];
				*/
				
				// we can set sum of every filter or just use one filter
				total = _V + _H + _45 + _315;
				//total = la;

				// taking only positive values
				if (total < 0) { 
					total = 0;
				}
				processed_array[z][processed_idx] = total;
				processed_idx++;	
			}
		}
	}
	return processed_array;
}

short ** thresholding(short** image_array) {
	short ** processed_array = (short**)malloc(sizeof(short*) * 58);

	for (int z = 0; z < 58; z++) {

		processed_array[z] = (short*)malloc(sizeof(short) * 512 * 512);

		for (int idx = 0; idx < 512 * 512; idx++) {

			if ((image_array[z][idx] < (-1024)) || (image_array[z][idx] > (-400))) {
				processed_array[z][idx] = 0;
			}
			else {
				//processed_array[z][idx] = image_array[z][idx];
				processed_array[z][idx] = 1;
			}
		}
	}
	return processed_array;	
}

short** _3D_CCA(short ** image_array) {
	short ** processed_array = (short**)malloc(sizeof(short*) * 58);
	short * equivalency_list = new short[2];
	int pass13_dir_vec[13] = { 513, 512, 511, 1, 0, -1, -511, -512, -513, 513, 512, 511, 1};
	int idx;
	short maximum_label = 0;
	short * minmax = new short[2];
	equivalency_list[0] = 0;
	
	// for z axis,
	for (int z = 0; z <58; z++) {
		printf("current z point : %d\n", z);
		
		// allocate memory for creating label at current z axis
		processed_array[z] = (short*)malloc(sizeof(short) * 512 * 512);

		// for 2d plane,
		for (idx = 0; idx < 512*512; idx++) {
			short *searched_label_list = new short[13];
			int check_step = 0;
			//printf("current idx : %d \n", idx);

			// if current idx value is 0, go to next idx value
			if (image_array[z][idx] != 1) {
				processed_array[z][idx] = 0;
				continue;
			}

			if (z == 0) {

				// if at the z axis 0, idx is 0, and value is 1, we can search nothing.
				if (idx == 0) {
					// so set all searched data into 0
					for (int dir = 0; dir < 13; dir++) {
						searched_label_list[dir] = 0;
					}
				}

				// if idx is not zero, then.. 
				else {
					for (int dir = 0; dir < 13; dir++) {

						// at the axis 0, we can't search up side data  
						if (dir < 9) {
							// so set all upside idx into 0
							searched_label_list[dir] = 0;
							
						}

						// but we can search current z axis data
						else {

							// if idx is at the left edge,
							if (idx % 512 == 0) {
								// we can't search left part so set value at those idx into 0
								if ((dir == 9) || (dir == 12)) {
									searched_label_list[dir] = 0;
								}

								// but we can search at the (idx-511) and (idx-512) value
								else {
									searched_label_list[dir] = processed_array[z][idx - pass13_dir_vec[dir]];
								}
							}

							// if idx is at the right edge,
							else if (idx % 512 == 511) {
								// we can't search right part so set value at those idx into 0
								if (dir == 11) {
									searched_label_list[dir] = 0;
								}

								// but we can search at the (idx-513), (idx-512) and (idx-1) value
								else {
									searched_label_list[dir] = processed_array[z][idx - pass13_dir_vec[dir]];
								}
							}

							// if idx is not at the left and right edge
							else {

								// if idx is on the top edge,
								if ((idx < 512) && (dir<12)) {
									searched_label_list[dir] = 0;
								}

								// in the other case we can search
								else {
									searched_label_list[dir] = processed_array[z][idx - pass13_dir_vec[dir]];
								}
							}
						}
					}
				}
			}

			else if (z > 0) {

				for (int dir = 0; dir < 13; dir++) {

					// if idx is at the left edge, 
					if (idx % 512 == 0) {

						// we can't search left part
						if (dir % 3 == 0) {
							searched_label_list[dir] = 0;
							check_step++;
						}

						else {

							// if idx is at the top corner
							if ((dir % 9 < 3)&&(idx==0)) {
								searched_label_list[dir] = 0;
								check_step++;
							}

							// if idx is at the bottom corner
							else if ((dir % 9 > 5)&& (idx == (512 * 512 - 512))) {
								searched_label_list[dir] = 0;
								check_step++;
								
							}

							else {
								if (dir < 9) {
									searched_label_list[dir] = processed_array[z - 1][idx - pass13_dir_vec[dir]];
									/////////////////
									if (searched_label_list[dir] < 0) {
										printf("z idx : %d \t idx : %d", z, idx);
										printf("searched label : %d", processed_array[z - 1][idx - pass13_dir_vec[dir]]);
									}
									check_step++;
								}
								else {
									searched_label_list[dir] = processed_array[z][idx - pass13_dir_vec[dir]];
									/////////////////
									if (searched_label_list[dir] < 0) {
										printf("z idx : %d \t idx : %d", z, idx);
										printf("searched label : %d", processed_array[z][idx - pass13_dir_vec[dir]]);
									}
									check_step++;
								}
							}
						}
					}

					// if idx is at the right edge
					else if (idx % 512 == 511) {

						// we can't search right part (dir2, dir5, dir8, dir11)
						if ((dir % 3) == 2) {
							searched_label_list[dir] = 0;
							check_step++;
						}

						else {

							// if idx is at the top corner
							if ((dir % 9 < 2)&&(idx==511)) {
								searched_label_list[dir] = 0;
								check_step++;
							}

							// if idx is at the bottom corner
							else if ((dir % 9 > 5)&& (idx == (512 * 512 - 1))) {
								searched_label_list[dir] = 0;
								check_step++;
							}

							else {
								if (dir < 9) {
									searched_label_list[dir] = processed_array[z - 1][idx - pass13_dir_vec[dir]];
									/////////////////
									if (searched_label_list[dir] < 0) {
										printf("step1\n");
										printf("z idx : %d \t idx : %d", z, idx);
										printf("searched label : %d", processed_array[z][idx - pass13_dir_vec[dir]]);
									}
									check_step++;
								}
								else {
									searched_label_list[dir] = processed_array[z][idx - pass13_dir_vec[dir]];
									/////////////////
									if (searched_label_list[dir] < 0) {
										printf("step2\n");
										printf("z idx : %d \t idx : %d", z, idx);
										printf("searched label : %d", processed_array[z][idx - pass13_dir_vec[dir]]);
									}
									check_step++;
								}
							}
						}
					}

					// if idx is not on the edge
					else {

						// if idx is at the top edge
						if (idx < 512) {

							// we can't search upper part
							if (dir % 9 < 3) {
								searched_label_list[dir] = 0;
								check_step++;
							}

							else {
								if (dir < 9) {
									searched_label_list[dir] = processed_array[z - 1][idx - pass13_dir_vec[dir]];
									/////////////////
									if (searched_label_list[dir] < 0) {
										printf("step3\n");
										printf("z idx : %d \t idx : %d", z, idx);
										printf("searched label : %d", processed_array[z - 1][idx - pass13_dir_vec[dir]]);
									}
									check_step++;
								}
								else {
									searched_label_list[dir] = processed_array[z][idx - pass13_dir_vec[dir]];
									/////////////////
									if (searched_label_list[dir] < 0) {
										printf("step4\n");
										printf("z idx : %d \t idx : %d", z, idx);
										printf("searched label : %d", processed_array[z][idx - pass13_dir_vec[dir]]);
									}
									check_step++;
								}
							}
						}

						// if idx is at the bottom edge
						else if (idx > (512 * 512 - 513)) {

							// we can't search down part
							if (dir % 9 > 5) {
								searched_label_list[dir] = 0;
								check_step++;
							}
							else {
								if (dir < 9) {
									searched_label_list[dir] = processed_array[z - 1][idx - pass13_dir_vec[dir]];
									/////////////////
									if (searched_label_list[dir] < 0) {
										printf("step5\n");
										printf("z idx : %d \t idx : %d\t", z, idx);
										printf("searched label : %d\n", processed_array[z - 1][idx - pass13_dir_vec[dir]]);
									}
									check_step++;
								}
								else {
									searched_label_list[dir] = processed_array[z][idx - pass13_dir_vec[dir]];
									/////////////////
									if (searched_label_list[dir] < 0) {
										printf("step6\n");
										printf("z idx : %d \t idx : %d\t", z, idx);
										printf("searched label : %d\n", processed_array[z][idx - pass13_dir_vec[dir]]);
									}
									check_step++;
								}
							}
						}

						else {
							if (dir < 9) {
								check_step++;
								searched_label_list[dir] = processed_array[z - 1][idx - pass13_dir_vec[dir]];
							}
							else {
								check_step++;
								searched_label_list[dir] = processed_array[z][idx - pass13_dir_vec[dir]];
							}
						}
					}
				}
			}

			/*
			Now, we find every possible search point.
			So, we need to check the minimum label and labeling the current idx pixel
			and we also need to set equivalence list
			*/
			find_min_max(searched_label_list, &minmax);

			// if minimum label is 0, 
			if (minmax[0] == 0) {
				// increase the maximum label by add 1
				++maximum_label;

				// set equivalency list
				equivalency_list[maximum_label] = maximum_label;

				// and increase the equivalency list
				increase_arr_size(maximum_label, &equivalency_list);

				// label current idx into 1
				processed_array[z][idx] = maximum_label;
			}

			// if minimum label is bigger than 0,
			else {

				// label current idx into minimum 
				processed_array[z][idx] = minmax[0];

				for (int dir = 0; dir < 13; dir++) {
					if ((searched_label_list[dir] != 0) 
						&& ((minmax[0] != searched_label_list[dir]) 
							|| (minmax[0] == searched_label_list[dir]))) {
						if (searched_label_list[minmax[0] != minmax[0]]) {
							equivalency_list[searched_label_list[dir]] = equivalency_list[minmax[0]];
						}
						else {
							equivalency_list[searched_label_list[dir]] = minmax[0];
						}
					}
				}

				// check if min label and max label is different and if min and max are not same,
				/*
				if (minmax[0] != minmax[1]) {
					for (int dir = 0; dir < 13; dir++) {
						if ((equivalency_list[minmax[0]] != minmax[0])) {
							if ((searched_label_list[dir] != 0) && ((minmax[0] != searched_label_list[dir]) || (minmax[0] == searched_label_list[dir]))) {
								equivalency_list[searched_label_list[dir]] = equivalency_list[minmax[0]];
							}
							else {
								equivalency_list[searched_label_list[dir]] = minmax[0];
							}
						}
						/*
						if ((searched_label_list[dir] != 0) && (minmax[0] != searched_label_list[dir])) {

							// reset the equivalency list
							if (equivalency_list[minmax[0]] != minmax[0]) {
								equivalency_list[searched_label_list[dir]] = equivalency_list[minmax[0]];
							}
							else {
								equivalency_list[searched_label_list[dir]] = minmax[0];
							}
							//equivalency_list[searched_label_list[dir]] = minmax[0];
							//printf("reset %d -> %d\t", searched_label_list[dir], minmax[0]);
						}
						// if searched label is not 0 and not same with min value
					}
				}
				*/
			}
			free(searched_label_list);
		}
	}
	
		for (int z = 0; z < 58; z++) {
		for (int i = 0; i < 512 * 512; i++) {
			if (processed_array[z][i] != equivalency_list[processed_array[z][i]]) {
				processed_array[z][i] = equivalency_list[processed_array[z][i]];
			}
		}
	}
	//store_result(processed_array);
	/*
	for (int z = 0; z < 58; z++) {
		for (int idx = 0; idx < 512 * 512; idx++) {
			if ((processed_array[z][idx] == 189)||(processed_array[z][idx] == 201)) {
				processed_array[z][idx] = 4;
			}
			else if (processed_array[z][idx] == 17) {
				processed_array[z][idx] = 4;
			}
			else {
				processed_array[z][idx] = 0;
			}	
		}
	}	
	*/

	store_result(processed_array);

	
	
	/*
	short * hist_list = new short[maximum_label+1];
	for (int l = 0; l < maximum_label + 1; l++) {
		hist_list[l] = 0;
	}
	for (int z = 0; z < 56; z++) {
		for (int l = 1; l < 512*512 ; l++) {
			hist_list[processed_array[z][l]] = hist_list[processed_array[z][l]] + 1;
		}
	}
	for (int l = 0; l < maximum_label + 1; l++) {
		if (hist_list[l] != 0) {
			std::cout << "number of " << l << "label => "<< hist_list[l] << std::endl;
		}
	}
	*/
	

	return processed_array;
}

void increase_arr_size(int s, short**arr) {
	//int newsize = *s + 1; 
	short* temp = (short*)malloc(sizeof(short)*(s+2));

	for (int i = 0; i < s+1; i++) {
		temp[i] = (*arr)[i];
	}
	
	free(*arr); 
	*arr = temp;
}

void find_min_max(short * arr, short**minmax) {
	
	(*minmax)[0] = arr[0]; // space for minimum
	(*minmax)[1] = arr[0]; // space for maximum

	for (int i = 1; i < 13; i++) {

		// if current min label is 0,
		if ((*minmax)[0] == 0) {
			// then if non zero label found, reset the minimum label
			if ((*minmax)[0] > -(arr[i])) {
				(*minmax)[0] = arr[i];
			}
		}

		// if initial min label is not 0 and current min label is bigger than searched label, reset the minimum label
		else if (((*minmax)[0] != 0) && (arr[i]>0) && ((*minmax)[0] > arr[i])) {
			(*minmax)[0] = arr[i];
		}

		// find max label value
		if (arr[i] > (*minmax)[1]) {
			(*minmax)[1] = arr[i];
		}
	}
}

void store_result(short ** image_array) {
	FILE *file;
	fopen_s(&file, "./Result/cca3d_data_1.raw", "wb");
	printf("\nwriting files ... \n");
	for (int z = 0; z < 58; z++) {
		fwrite(image_array[z], sizeof(short), 512 * 512, file);
	}
	fclose(file);
}

