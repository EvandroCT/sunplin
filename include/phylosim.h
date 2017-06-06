/*************************************************************************
	
	Copyright (C) 2017	Evandro Taquary
	
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

#ifndef _PHYLOSIM_H
#define _PHYLOSIM_H

#include <curand_kernel.h>
#include "dtree.h"

//create all necessary seeds to massive GPU randomize
__global__ void setup_kernel(long long seed, curandState_t* devStates, ushort N);

//trees' exapansions
__global__ void insertion(DTree tree, curandState_t* devStates);

//generate the patristic distance matrices to all the replics
__global__ void patrix(DTree tree, float* d_matrix);

#endif
