/*************************************************************************
	
	Copyright (C) 2016	Evandro Taquary
	
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

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "modcpy.h"

using namespace std;

#define CHECK(call) \
{ \
        const cudaError_t error = call; \
        if (error != cudaSuccess) \
        { \
                cout << "Error: " << __FILE__ ": " << __LINE__ << ", "; \
                cout << "code: "<< error << ", reason: " << cudaGetErrorString(error) << endl; \
                exit(EXIT_FAILURE); \
        } \
}

#define FERR(file) \
{ \
	if(!file.good()){ \
		cout << "File opening failure. Please try again." << endl; \
		exit(EXIT_FAILURE); \
	} \
}

class SoaTree {
	
	public:
		int *nodeId;			// identification numbers of the nodes
		int *nodeHeight; 		// heights of the nodes on the tree
		int *qtyBeneath;		// amount of taxa under the nodes (size of the subtree)
		int *nodeParent; 		// nodes' parents 
		int *leftChild; 		// nodes' left children
		int *rightChild; 		// nodes' right children
		float *nodeBranch;		// lengths of the nodes' branches (distance to the parent)
		float *distRoot; 		// distances between nodes and root (sum of the paths' branches)
		float *distSpec; 		// distances between all species (patristic distance)
		
		__host__ void setData(int num_elem, void *base_adr);
		__host__ void* getPtr() const {return (void*) nodeId;}
		__host__ size_t getSize(int num_elem){
			size_t size = (6*sizeof(int) + 3*sizeof(float)) * num_elem;	//minimal amount of bytes needed to represent the soa 
			int r = size%sizeof(int4);
			size += r ? sizeof(int4)-r : 0;	//size of the tree padded to multiple of sizeof(int4) (due to a GPU memory aligment requisite)
			return size;
		}		
		__host__ void soalloc(int num_elem) {
			void *ptr = malloc(getSize(num_elem));
			memset(ptr,0,getSize(num_elem));
			setData(num_elem, ptr);
		}
};

void SoaTree::setData(int num_elem, void* base_adr){
	
	nodeId = (int*) base_adr;
	nodeHeight 	= nodeId 		+ num_elem;
	qtyBeneath	= nodeHeight	+ num_elem;
	nodeParent 	= qtyBeneath	+ num_elem;
	leftChild 	= nodeParent	+ num_elem;
	rightChild 	= leftChild 	+ num_elem;
	nodeBranch 	=(float*)(rightChild+ num_elem);
	distRoot 	= nodeBranch	+ num_elem;
	distSpec 	= distRoot 		+ num_elem;	
}

class DTree{
	
	protected:			
		int qtyNodes;			// quantity of nodes on the tree
		int qtyInsSpc;			// quantity of absent species to be inserted
		int idxInsSpc;			// starting idex for insertion of the species
		int idxInsAnc;			// starting idex for insertion of the ancestors	
	public:
		SoaTree devData;		// struct of arrays to hold the trees' data
		__host__ int getQtyNodes() const {return qtyNodes;}
		__host__ int getQtyInsSpc() const {return qtyInsSpc;}
		__host__ int getIdxInsSpc() const {return idxInsSpc;}
		__host__ int getIdxInsAnc() const {return idxInsAnc;}
		virtual __host__ size_t getSize() const = 0;	// force class to be abstract
};

class HTree: public DTree{
	
	public:
		SoaTree hostData;			// struct of arrays to hold the trees' data 		
		vector<string> nodeName; 	// names of taxa fetched from newick and PUT file
		size_t treeSize;			// size of the tree padded to multiple of sizeof(int4) (due to a GPU memory aligment requisite)
		ifstream newick;			// stream object to manipulate input newick file
		//ifstream puts;
		__host__ void setQtyNodes();
		__host__ void setTreeSize();
		__host__ void parseTree();
	public:
		__host__ HTree() : HTree("wellParser.out"){}
		__host__ HTree(string fname);
		__host__ size_t getSize() const {return treeSize;}
		__host__ void* getDevPtr() const {return devData.getPtr();}
};

void HTree::setQtyNodes(){
	
	newick.seekg(0);
	newick >> qtyNodes;
	FERR(newick);
}

void HTree::setTreeSize(){
	//TODO: throw exception
	treeSize = hostData.getSize(qtyNodes);
}

HTree::HTree(string fname){
	
	newick.open(fname);
	FERR(newick);
	setQtyNodes();
	hostData.soalloc(qtyNodes);
	parseTree();

	void* d_tree;
	setTreeSize();
	CHECK(cudaMalloc(&d_tree, treeSize));
	CHECK(cudaMemcpy(d_tree, hostData.getPtr(), treeSize, cudaMemcpyHostToDevice));
	devData.setData(qtyNodes, d_tree);
}

void HTree::parseTree() {
	
	newick >> idxInsAnc >> idxInsSpc >> qtyInsSpc;	
	string taxon_name;
	for(int i=0; i<qtyNodes; i++) newick >> hostData.nodeId[i];
	for(int i=0; i<qtyNodes; i++){newick >> taxon_name; nodeName.push_back(taxon_name);}	
	for(int i=0; i<qtyNodes; i++) newick >> hostData.nodeBranch[i];	
	for(int i=0; i<qtyNodes; i++) newick >> hostData.nodeHeight[i];	
	for(int i=0; i<qtyNodes; i++) newick >> hostData.qtyBeneath[i];
	for(int i=0; i<qtyNodes; i++) newick >> hostData.nodeParent[i];
	for(int i=0; i<qtyNodes; i++) newick >> hostData.leftChild[i];
	for(int i=0; i<qtyNodes; i++) newick >> hostData.rightChild[i];
	FERR(newick);
}

class PhyloGen{
	
	private:
		int qtyReplics;
		void *d_replics;
	public:
		PhyloGen(int qty) : qtyReplics(qty){}
		void gpuRep(const HTree *tree);		
};

void PhyloGen::gpuRep(const HTree *tree){			
	
	size_t tree_size = tree->getSize();
	size_t rep_size = tree_size * qtyReplics;
	CHECK(cudaMalloc(&d_replics, rep_size));
	
	cudaDeviceProp device;
	CHECK(cudaGetDeviceProperties(&device,0));
	int threads = device.warpSize*16;
	int blocks = (rep_size/sizeof(int4) + (threads-1)) / threads;
	dim3 grid = dim3(blocks);
	dim3 block = dim3(threads);
	modcpy<<<grid, block>>>(d_replics, tree->getDevPtr(), rep_size, tree_size);
	CHECK(cudaDeviceSynchronize());
}

int main(int argc, char *argv[]){		

	if(argc < 2 || argc >3){
		cout << "Usage: " << argv[0] << " #replications [file]" << endl;
		exit(EXIT_FAILURE);
	}	
	PhyloGen simulation(atoi(argv[1]));
	HTree *tree = argc==3 ? new HTree(argv[2]) : new HTree();
	simulation.gpuRep(tree);
		
	exit(EXIT_SUCCESS);
}
