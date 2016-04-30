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
	
	private:
		int *idx;			// index of input newickf trees' nodes; or, the subtrees' roots' indices where new nodes shall be inserted to   
		int *height; 		// heights of the nodes on the tree
		int *nSubt;			// amount of taxa within subtree rooted on idx (size of the subtree)
		int *parent; 		// nodes' parents 
		int *lChild; 		// nodes' left children
		int *rChild; 		// nodes' right children
		float *branch;		// lengths of the nodes' branches (distance to the parent)
		//float *dRoot; 	// distances between nodes and root (sum of the paths' branches)
		//float *dSpec; 	// distances between all species (patristic distance)		
	public:
		__host__ SoaTree() = default;
		__host__ SoaTree(int num_nodes) {soalloc(num_nodes);}
		__host__ void* getPtr() const {return (void*) idx;}
		__host__ void setOffs(int num_nodes, void* base);	//set pointers' offsets starting on base accordingly to data structure and number of nodes
		__host__ void setOffs(int num_nodes) {setOffs(num_nodes,idx);}  //set pointers' offsets starting on idx accordingly to data structure and number of nodes
		__host__ size_t getSize(int num_nodes)
		{
			size_t size = (6*sizeof(int) + 2*sizeof(float)) * num_nodes;	//minimal amount of bytes needed to represent the soa 
			int r = size%sizeof(int4);
			size += r ? sizeof(int4)-r : 0;	//size of the tree padded to multiple of sizeof(int4) (due to a GPU memory aligment requisite)
			return size;
		}
		__host__ void soalloc(int num_nodes)
		{
			void *ptr = malloc(getSize(num_nodes));
			memset(ptr,0,getSize(num_nodes));
			setOffs(num_nodes, ptr);
		}
		/* TODO: THROW OVERFLOW EXCEPTION */
		__host__ __device__ int getIdx	 (int i)   const {return idx[i];}
		__host__ __device__ int getHeight(int i)   const {return height[i];}
		__host__ __device__ int getnSubt (int i)   const {return nSubt[i];}
		__host__ __device__ int getParent(int i)   const {return parent[i];}
		__host__ __device__ int getlChild(int i)   const {return lChild[i];}
		__host__ __device__ int getrChild(int i)   const {return rChild[i];}
		__host__ __device__ float getBranch(int i) const {return branch[i];}
		//__host__ __device__ float getdRoot (int i) const {return dRoot[i];}
		/* TODO: THROW OVERFLOW EXCEPTION */
		__host__ void setIdx	(int val, int i)   {idx[i] 	  = val;}
		__host__ void setHeight (int val, int i)   {height[i] = val;}
		__host__ void setnSubt  (int val, int i)   {nSubt[i]  = val;}
		__host__ void setParent (int val, int i)   {parent[i] = val;}
		__host__ void setlChild (int val, int i)   {lChild[i] = val;}
		__host__ void setrChild (int val, int i)   {rChild[i] = val;}
		__host__ void setBranch (float val, int i) {branch[i] = val;}
		//__host__ void setdRoot  (float val, int i) {dRoot[i]  = val;}
};

void SoaTree::setOffs(int num_nodes, void* base) {
	
	idx 	= (int*) base;
	height 	= idx 		+ num_nodes;
	nSubt	= height	+ num_nodes;
	parent 	= nSubt		+ num_nodes;
	lChild 	= parent	+ num_nodes;
	rChild 	= lChild 	+ num_nodes;
	branch 	=(float*)(rChild+num_nodes);
	//dRoot = branch	+ num_nodes;
	//dSpec here
}

class DTree{
	
	protected:
		int nNodes;				// quantity of nodes on the tree (including inserting species)
		int nInsSpc;			// quantity of absent species to be inserted
		int idxInsSpc;			// starting index for insertion of new species
		int idxInsAnc;			// starting index for insertion of new ancestors	
	public:
		SoaTree devData;		// struct of arrays to hold trees' data		
		__host__ __device__ int getnNodes	() const {return nNodes;}
		__host__ __device__ int getnInsSpc	() const {return nInsSpc;}
		__host__ __device__ int getIdxInsSpc() const {return idxInsSpc;}
		__host__ __device__ int getIdxInsAnc() const {return idxInsAnc;}
		__device__ int getIdx	 (int i) const {return devData.getIdx(i);}
		__device__ int getHeight (int i) const {return devData.getHeight(i);}
		__device__ int getnSubt  (int i) const {return devData.getnSubt(i);}
		__device__ int getParent (int i) const {return devData.getParent(i);}
		__device__ int getlChild (int i) const {return devData.getlChild(i);}
		__device__ int getrChild (int i) const {return devData.getrChild(i);}
		__device__ float getBranch (int i) const {return devData.getBranch(i);}
		//__device__ float getdRoot  (int i) const {return devData.getdRoot(i);}
		virtual __host__ size_t getSize() const = 0;	// force class to be abstract
};

class HTree: public DTree{
	
	private:
		SoaTree hostData;			// struct of arrays to hold the trees' data 		
		vector<string> name; 		// names of taxa fetched from newickf and PUT file
		size_t treeSize;			// size of the tree padded to multiple of sizeof(int4) (due to a GPU memory aligment requisite)
		ifstream newickf;			// stream object to manage input newick file
		//ifstream putf;			// stream object to manage input PUT file
		int devId;					// id of the GPU where lies the tree 
		__host__ void setnNodes();
		__host__ void parseTree();
	public:
		__host__ HTree(string nw_fname = "newick.tree"/*, string pt_fname="put.list" */, int dev_id=0);
		__host__ void* gpuRep(int num_reps) const;
		__host__ size_t getSize() const {return treeSize;}
		__host__ bool valiData(void* d_reps, int num_reps);
		/* TODO: CATCH OVERFLOW EXCEPTION */
		__host__ void setIdx	(int val, int i)   {hostData.setIdx(val,i);}
		__host__ void setHeight (int val, int i)   {hostData.setHeight(val,i);}
		__host__ void setnSubt  (int val, int i)   {hostData.setnSubt(val,i);}
		__host__ void setParent (int val, int i)   {hostData.setParent(val,i);}
		__host__ void setlChild (int val, int i)   {hostData.setlChild(val,i);}
		__host__ void setrChild (int val, int i)   {hostData.setrChild(val,i);}
		__host__ void setBranch (float val, int i) {hostData.setBranch(val,i);}
		//__host__ void setdRoot  (float val, int i) {hostData.setdRoot(val,i);}
		__host__ string getName(int i) {return name[i];}
};

HTree::HTree(string nw_fname/*,string pt_fname*/, int dev_id){
	
	devId = dev_id;
	newickf.open(nw_fname);
	FERR(newickf);
	//putf.open(pt_fname);
	//FERR(putf);
	setnNodes();
	treeSize = hostData.getSize(nNodes);
	hostData.soalloc(nNodes);
	parseTree();
	newickf.close();
	//putf.close();
	void * d_tree;
	CHECK(cudaSetDevice(devId));
	CHECK(cudaMalloc(&d_tree, treeSize));
	CHECK(cudaMemcpy(d_tree, hostData.getPtr(), treeSize, cudaMemcpyHostToDevice));
	devData.setOffs(nNodes, d_tree);
}

void HTree::setnNodes(){
	
	newickf.seekg(0);
	newickf >> nNodes;
	FERR(newickf);
}

void HTree::parseTree() {
	
	newickf >> idxInsAnc >> idxInsSpc >> nInsSpc;	
	string taxon_name;
	int ival;
	float fval;
	for(int i=0; i<nNodes; i++) {newickf >> ival; hostData.setIdx(ival,i);}
	for(int i=0; i<nNodes; i++) {newickf >> taxon_name; name.push_back(taxon_name);}	
	for(int i=0; i<nNodes; i++) {newickf >> fval; hostData.setBranch(fval,i);}
	for(int i=0; i<nNodes; i++) {newickf >> ival; hostData.setHeight(ival,i);}	
	for(int i=0; i<nNodes; i++) {newickf >> ival; hostData.setnSubt(ival,i);}
	for(int i=0; i<nNodes; i++) {newickf >> ival; hostData.setParent(ival,i);}
	for(int i=0; i<nNodes; i++) {newickf >> ival; hostData.setlChild(ival,i);}
	for(int i=0; i<nNodes; i++) {newickf >> ival; hostData.setrChild(ival,i);}
	FERR(newickf);
}

void* HTree::gpuRep(int num_reps) const{
	
	size_t rep_size = treeSize * num_reps;
	void * d_replics;
	CHECK(cudaMalloc(&d_replics, rep_size));
	cudaDeviceProp device;
	CHECK(cudaGetDeviceProperties(&device, devId));
	int threads = device.warpSize*16;	//TODO: FIGURE OUT WHICH MULTIPLE IS THE BEST
	int blocks = (rep_size/sizeof(int4) + (threads-1)) / threads;
	dim3 grid = dim3(blocks);
	dim3 block = dim3(threads);
	modcpy<<<grid, block>>>(d_replics, devData.getPtr(), rep_size, treeSize);
	CHECK(cudaDeviceSynchronize());
	return d_replics;
}

bool HTree::valiData(void* d_replics, int num_reps){	
	
	SoaTree tree;	
	void* h_replics = malloc(num_reps * treeSize);
	CHECK(cudaMemcpy(h_replics, d_replics, num_reps*treeSize, cudaMemcpyDeviceToHost));
	cout.precision(3);
	cout.setf(ios::fixed, ios::floatfield);	
	for(int j=0; j<num_reps; j++){
		tree.setOffs(nNodes, h_replics+treeSize*j);
		for(int i=0; i<nNodes; i++){
			cout << "rep " << j+1 << ": " << tree.getIdx(i) << "\t\t" << tree.getBranch(i) << "\t\t" << tree.getHeight(i) << "\t\t" << tree.getParent(i) << "\t\t" << tree.getlChild(i) << "\t\t" << tree.getrChild(i) << "\t\t" << endl;
			if(	tree.getIdx(i) 	  != hostData.getIdx(i)		||
				tree.getHeight(i) != hostData.getHeight(i)	||
				tree.getBranch(i) != hostData.getBranch(i) 	||
				tree.getParent(i) != hostData.getParent(i)	||
				tree.getlChild(i) != hostData.getlChild(i)	||
				tree.getrChild(i) != hostData.getrChild(i)	)
					return false;
		}
	}
	return true;
}

int main(int argc, char *argv[]){
	
	if(argc < 2 || argc >3){
		cout << "Usage: " << argv[0] << " #replications [file]" << endl;
		exit(EXIT_FAILURE);
	}	
	HTree *tree = argc==3 ? new HTree(argv[2]) : new HTree();	
	void* d_replics = tree->gpuRep(atoi(argv[1]));
	if(tree->valiData(d_replics, atoi(argv[1])))
		cout << "Data does match!" << endl;
	else
		cout << "Data doesn't match" << endl;
		
	exit(EXIT_SUCCESS);
}
