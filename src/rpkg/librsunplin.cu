/*************************************************************************
	
	Copyright (C) 2017	Evandro Taquary
	
	This program is free software: you can redistribute it and/or modify s
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

#include <R.h>
#include <iostream>
#include <htree.h>
#include <cudutils.h>
#include <phylosim.h>

using namespace std;

extern "C" void expansion(int* nreps, char** nw_fname, char** pt_fname, char** o_fname, int* a, int* b, int* c, int* d, int* seed){	
	
	int gpu=0;
	ushort num_reps = *nreps;
	
	HTree tree(gpu,*nw_fname,*pt_fname);	
	DTree replics = tree.gpuRep(*nreps);
	
	*a = tree.getIdxInsSpc();
	*b = tree.getnInsSpc();
	*c = tree.getIdxInsSpc() + tree.getnInsSpc();
	*d = tree.getnNodes()-1;

	curandState_t *devStates;
	cudaDeviceProp device;
	CHECK(cudaGetDeviceProperties(&device,gpu));
	
	int threads = device.warpSize*16; //threads per block; TODO: FIGURE OUT WHICH MULTIPLE IS THE BEST
	int blocks = (num_reps + (threads-1)) / threads;
	dim3 grid(blocks), block(threads);
	
	CHECK(cudaMalloc((void**)&devStates, sizeof(curandState_t)*num_reps));	
	setup_kernel<<<grid,block>>>(*seed,devStates,num_reps);
	CHECK(cudaDeviceSynchronize());
	insertion<<<grid,block>>>(replics,devStates);	
	CHECK(cudaDeviceSynchronize());
	replics.toNewick(tree.getNames(),*o_fname); 
	CHECK(cudaDeviceReset());	
}

extern "C" void chained_patrices(int* nreps, char** nw_fname, char** pt_fname, char** o_fname, int* a, int* b, int* c, int* d, int* seed){	
	
	int gpu=0;
	ushort num_reps = *nreps;
	
	HTree tree(gpu,*nw_fname,*pt_fname);	
	DTree replics = tree.gpuRep(*nreps);
	
	*a = tree.getIdxInsSpc();
	*b = tree.getnInsSpc();
	*c = tree.getIdxInsSpc() + tree.getnInsSpc();
	*d = tree.getnNodes()-1;

	curandState_t *devStates;
	cudaDeviceProp device;
	CHECK(cudaGetDeviceProperties(&device,gpu));
	
	int threads = device.warpSize*16; //threads per block; TODO: FIGURE OUT WHICH MULTIPLE IS THE BEST
	int blocks = (num_reps + (threads-1)) / threads;
	dim3 grid(blocks), block(threads);
	
	CHECK(cudaMalloc((void**)&devStates, sizeof(curandState_t)*num_reps));	
	setup_kernel<<<grid,block>>>(*seed,devStates,num_reps);
	CHECK(cudaDeviceSynchronize());
	insertion<<<grid,block>>>(replics,devStates);	
	
	ushort nleafs = (replics.getnNodes()+1)/2;
	uint msize = nleafs*(nleafs-1)/2;

	float *d_matrices;
	CHECK(cudaMalloc((void**)&d_matrices, sizeof(float)*msize*num_reps));
	patrices<<<num_reps,256,replics.getnNodes()*(sizeof(ushort)*3)>>>(replics, d_matrices);
	
	float *h_matrices = (float*)malloc(sizeof(float)*msize*num_reps);
	CHECK(cudaMemcpy(h_matrices, d_matrices, sizeof(float)*msize*num_reps, cudaMemcpyDeviceToHost));
	
	replics.toCSV(h_matrices, *o_fname);
	CHECK(cudaDeviceReset());
}
