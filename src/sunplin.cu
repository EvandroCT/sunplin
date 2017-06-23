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

#include <htree.h>
#include <cudutils.h>
#include <phylosim.h>
#include <iostream>

using namespace std;

double time_spent[7];

int main(int argc, char *argv[]){

	if(argc < 2 || argc >4){
		cout << "Usage: " << argv[0] << " #replications [newick_file putlist_file]" << endl;
		exit(EXIT_FAILURE);
	}
	char op;
	int gpu=0;	
	int num_reps = atoi(argv[1]);	
	HTree *tree = argc>2 ? new HTree(gpu,argv[2],argv[3]) : new HTree(gpu);
	
	CHECK(cudaSetDevice(gpu));
	
	SETUP_TIMER(true);
	START_TIMER();
	DTree replics = tree->gpuRep(num_reps);
	STOP_TIMER(time_spent[2]);	

	cout << "\n# of species of original tree: " << tree->getIdxInsSpc() << endl;
	cout << "# of PUTs to be inserted: " << tree->getnInsSpc() << endl;
	cout << "# of species after insertion: " << tree->getIdxInsSpc() + tree->getnInsSpc() << endl;
	cout << "# of nodes (including internals): " << tree->getnNodes()-1 << endl << endl;

	/*
	if(replics.compareTo(tree))
		cout << "Data does match!" << endl;
	else
		cout << "Data doesn't match" << endl;
	*/

	curandState_t *devStates;
	cudaDeviceProp device;
	CHECK(cudaGetDeviceProperties(&device,gpu));
	
	int threads = device.warpSize*16; //threads per block; TODO: FIGURE OUT WHICH MULTIPLE IS THE BEST
	int blocks = (num_reps + (threads-1)) / threads;
	dim3 grid(blocks), block(threads);
	
	START_TIMER();
	CHECK(cudaMalloc((void**)&devStates, sizeof(curandState_t)*num_reps));	
	setup_kernel<<<grid,block>>>(1,devStates,num_reps);
	CHECK(cudaDeviceSynchronize());
	insertion<<<grid,block>>>(replics,devStates);	
	CHECK(cudaDeviceSynchronize());
	STOP_TIMER(time_spent[3]);
	
	START_TIMER();
	ushort nleafs = (replics.getnNodes()+1)/2;
	uint msize = nleafs*(nleafs-1)/2;

	float *d_matrix;
	CHECK(cudaMalloc((void**)&d_matrix, sizeof(float)*msize*num_reps));
	patrices<<<num_reps,256,replics.getnNodes()*(sizeof(ushort)*3)>>>(replics, d_matrix);
	CHECK(cudaDeviceSynchronize());
	STOP_TIMER(time_spent[4]);	
	
	START_TIMER();
	float *h_matrix = (float*)malloc(sizeof(float)*msize*num_reps);
	CHECK(cudaMemcpy(h_matrix, d_matrix, sizeof(float)*msize*num_reps, cudaMemcpyDeviceToHost));
	STOP_TIMER(time_spent[5]);

	cout<<"\ntotal time spent to parse the files: "<<time_spent[0]<<"s";
	cout<<"\ntotal time spent to copy backbone tree to GPU: "<<time_spent[1]<<"s";
	cout<<"\ntotal time spent to replicate trees: "<<time_spent[2]<<"s";
	cout<<"\ntotal time spent to expand trees: "<<time_spent[3]<<"s";
	cout<<"\ntotal time spent to generate patrices: "<<time_spent[4]<<"s";
	cout<<"\ntotal time spent to copy patrices back to host: "<<time_spent[5]<<"s\n\n";

/*
	cout<<"\nDo you want to store the generated trees into a file (it may take a long time)? (Y/N): "; cin >> op;

	if(op=='y' || op=='Y'){
		START_TIMER();
		replics.toNewick(tree->getNames()); 
		STOP_TIMER(time_spent[6]);
		cout<<"\ntotal time spent to create output file (versions.tree): "<<time_spent[6]<<"s\n\n";
	}
*/

	replics.print(tree->getNames());

	CHECK(cudaDeviceReset());	
	exit(EXIT_SUCCESS);	
}
