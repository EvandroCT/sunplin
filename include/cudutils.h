/*************************************************************************
	
	Copyright (C) 2017	Evandro Taquary, Mateus Freitas
	
	This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
	
*************************************************************************/

#include <sys/time.h>

#ifndef _UTILS_H
#define _UTILS_H

/*
 * This CUDA kernel perform a generic modular copy within GPU global memory. 
 * The data ranging from 'source' address to 'source' + 'source_size' 
 * adress are copied in parallel to a location starting at 'destination' until 
 * 'destination' + 'destination_size' in a recycle (modular) fashion,that is, 
 * the data will be copied "byte by byte" until the entire destination range 
 * is repeatively filled by the source data.  
 * */
__global__ void modcpy(void *destination, void *source, size_t destination_size, size_t source_size);

/*
 * Returns the row's index of the equivalent upper triangular matrix where 'i' 
 * is the position of the element relative to the flattened (vectorized) upper 
 * triangular matrix and M is the number of elements in the set.
 * */
__host__ __device__ int row_index( int i, int M );

/*
 * Returns the column's index of the equivalent upper triangular matrix where 'i' 
 * is the position of the element relative to the flattened (vectorized) upper 
 * triangular matrix and M is the number of elements in the set.
 * */
__host__ __device__ int column_index( int i, int M );

/*
 * Macro defined to catch file handling errors.
 * */ 
#define FERR(file) \
		{ \
			if(!file.good()){ \
				cout << "Something went wrong while handling file! Please try again." << endl; \
				cout << "Error: " << __FILE__ ": " << __LINE__ << ", " << endl; \
				exit(EXIT_FAILURE); \
			} \
		}		

/*
 * Macro defined to catch CUDA runtime errors.
 * */ 
#define CHECK(call) \
		{ \
			const cudaError_t error = call; \
			if (error != cudaSuccess) { \
				cout << "Error: " << __FILE__ ": " << __LINE__ << ", "; \
				cout << "code: "<< error << ", reason: " << cudaGetErrorString(error) << endl; \
				exit(EXIT_FAILURE); \
			} \
		}

/*
 * The macros specified below were conceived to allow a precise 
 * time consuming measurement while minimum overhead is required to do 
 * so.
 * */
 
/*
 * Sets up variables needed to allow time measuring. The variables have 
 * local scope, hence this macro has to be ONCE placed within the aimed 
 * code block even if multiple measurements are to be porformed. The arg
 * 'toggle' is a user-defined parameter and specifies wether time 
 * measurements shall be performed.
 * */
#define SETUP_TIMER(toggle) \
			bool timer = toggle; \
			long long start_time; \
			long long end_time; \
			struct timeval tv;

/*
 * Starts time measurement if 'timer' is set to true. Can be called 
 * several times within the same code block provided SETUP_TIMER were 
 * placed before the very first call.
 * */
#define START_TIMER() \
		if(timer){ \
			gettimeofday(&tv, NULL); \
			start_time = tv.tv_sec * 1e6 + tv.tv_usec; \
		}

/*
 * Stops time measurement. Returns time measured in seconds if 'timer' 
 * is set to true.
 * */	
#define STOP_TIMER(time_spent) \
		if(timer){ \
			gettimeofday(&tv, NULL); \
			end_time = tv.tv_sec * 1e6 + tv.tv_usec; \
			time_spent = ((double)(end_time-start_time))/1e6; \
		}


#endif
