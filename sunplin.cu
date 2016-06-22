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
#include <curand_kernel.h>
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
		int *idx;			// index of input newickf trees' nodes   
		int *height; 		// heights of the nodes on the tree
		int *nSubt;			// amount of taxa within subtree rooted on idx (size of the subtree)
		int *parent; 		// nodes' parents; or, the subtrees' roots' indices where new nodes shall be inserted
		int *lChild; 		// nodes' left children
		int *rChild; 		// nodes' right children
		float *branch;		// lengths of the nodes' branches (distance to the parent)
		float *dRoot; 		// distances between nodes and root (sum of the paths' branches)
		int *inseq;			// vector with the sequence of indices of puts to be inserted
		//float *dSpec; 	// distances between all species (patristic distance)		
	public:
		__host__ SoaTree() = default;		
		__host__ SoaTree(int num_nodes, int num_ins) {soalloc(num_nodes, num_ins);}
		__host__ SoaTree(int num_nodes, int num_ins, void* base) {setOffs(num_nodes,num_ins,base);}
		__host__ __device__ void* getPtr() const {return (void*) idx;}
		__host__ __device__ void setOffs(int num_nodes, int num_ins, void* base);	//set pointers' offsets starting on base accordingly to data structure and number of nodes
		__host__ void setOffs(int num_nodes, int num_ins) {setOffs(num_nodes, num_ins, idx);}  //set pointers' offsets starting on idx accordingly to data structure and number of nodes
		__host__ static size_t getSize(int num_nodes, int num_ins)
		{
			size_t size = (6*sizeof(int) + 2*sizeof(float))*num_nodes + sizeof(int)*num_ins; //minimal amount of bytes needed to represent the tree 
			int r = size%sizeof(int4);
			size += r ? sizeof(int4)-r : 0;	//size of the tree padded to multiple of sizeof(int4) (due to a GPU memory aligment requisite)
			return size;
		}
		__host__ void soalloc(int num_nodes, int num_ins)
		{
			void *ptr = malloc(getSize(num_nodes, num_ins));
			memset(ptr,0,getSize(num_nodes, num_ins));
			setOffs(num_nodes, num_ins, ptr);
		}
		/* TODO: THROW OVERFLOW EXCEPTION */
		__host__ __device__ int getIdx	 	(int i) const {return idx[i];}
		__host__ __device__ int getHeight	(int i) const {return height[i];}
		__host__ __device__ int getnSubt 	(int i) const {return nSubt[i];}
		__host__ __device__ int getParent	(int i) const {return parent[i];}
		__host__ __device__ int getlChild	(int i) const {return lChild[i];}
		__host__ __device__ int getrChild	(int i) const {return rChild[i];}
		__host__ __device__ int getInseq 	(int i) const {return inseq[i];}
		__host__ __device__ float getBranch	(int i) const {return branch[i];}
		__host__ __device__ float getdRoot 	(int i) const {return dRoot[i];}
		/* TODO: THROW OVERFLOW EXCEPTION */
		__host__ __device__ void setTreeIdx(int val, int i)   {idx[i] 	 = val;}
		__host__ __device__ void setHeight (int val, int i)   {height[i] = val;}
		__host__ __device__ void setnSubt  (int val, int i)   {nSubt[i]  = val;}
		__host__ __device__ void setParent (int val, int i)   {parent[i] = val;}
		__host__ __device__ void setlChild (int val, int i)   {lChild[i] = val;}
		__host__ __device__ void setrChild (int val, int i)   {rChild[i] = val;}
		__host__ __device__ void setInseq  (int val, int i)   {inseq[i]  = val;}
		__host__ __device__ void setBranch (float val, int i) {branch[i] = val;}
		__host__ __device__ void setdRoot  (float val, int i) {dRoot[i]  = val;}
};

void SoaTree::setOffs(int num_nodes, int num_ins, void* base) {	
	idx 	= (int*) base;
	height 	= idx 		+ num_nodes;
	nSubt	= height	+ num_nodes;
	parent 	= nSubt		+ num_nodes;
	lChild 	= parent	+ num_nodes;
	rChild 	= lChild 	+ num_nodes;	
	branch 	=(float*)(rChild+num_nodes);
	dRoot 	= branch	+ num_nodes;
	inseq 	=(int*)  (dRoot+num_ins);
}

class HTree;

class DTree{
	protected:
		int nNodes;				// quantity of nodes on the tree(s) (including inserting species)
		int nInsSpc;			// quantity of absent species to be inserted
		int idxInsSpc;			// starting index for insertion of new species
		int idxInsAnc;			// starting index for insertion of new ancestors
		int nTrees;				// quantity of trees holded by devData (default=1)
		size_t treeSize;		// size of one tree padded to multiple of sizeof(int4) (due to a GPU memory aligment requisite)
		SoaTree devData;		// struct of arrays to hold trees' data
	public:
		__host__ bool compareTo(HTree *h_tree);
		__host__ DTree() = default;
		__host__ DTree(int nNodes, int nInsSpc, int idxInsSpc, int idxInsAnc, int nTrees, size_t treeSize, void* ptr):
						nNodes(nNodes),
						nInsSpc(nInsSpc),
						idxInsSpc(idxInsSpc), 
						idxInsAnc(idxInsAnc),
						nTrees(nTrees),
						treeSize(treeSize){devData.setOffs(nNodes,nInsSpc,ptr);}
		__host__ __device__ int getnNodes	() const {return nNodes;}
		__host__ __device__ int getnInsSpc	() const {return nInsSpc;}
		__host__ __device__ int getIdxInsSpc() const {return idxInsSpc;}
		__host__ __device__ int getIdxInsAnc() const {return idxInsAnc;}
		__host__ __device__ int getnTrees() const {return nTrees;}
		__host__ __device__ int getIdx	 (int i) const {return devData.getIdx(i);}
		__host__ __device__ int getHeight (int i) const {return devData.getHeight(i);}
		__host__ __device__ int getnSubt  (int i) const {return devData.getnSubt(i);}
		__host__ __device__ int getParent (int i) const {return devData.getParent(i);}
		__host__ __device__ int getlChild (int i) const {return devData.getlChild(i);}
		__host__ __device__ int getrChild (int i) const {return devData.getrChild(i);}
		__host__ __device__ float getBranch (int i) const {return devData.getBranch(i);}
		__device__ int getInseq (int i)   const {return devData.getInseq(i);}
		
		__device__ void setdRoot  (float val, int i)  {devData.setdRoot(val,i);}
		__device__ void setlChild (int val, int i) 	{devData.setlChild(val,i);}
		__device__ void setrChild (int val, int i)	{devData.setrChild(val,i);}
		__device__ void setParent (int val, int i)	{devData.setParent(val,i);}
		__device__ void setBranch (float val, int i)	{devData.setBranch(val,i);}
		__device__ void setInseq (int val, int i)  {devData.setInseq(val,i);}
		__device__ void setTreeIdx(int i){devData.setOffs(nNodes,nInsSpc,devData.getPtr()+treeSize*i);}
		__device__ float getdRoot  (int i) const {return devData.getdRoot(i);}
		__host__ __device__ size_t getSize() const { return treeSize; };
		__host__ void print();		
};

void DTree::print(){

	void *hm = malloc(treeSize*nTrees);	
	CHECK(cudaMemcpy(hm,devData.getPtr(),treeSize*nTrees,cudaMemcpyDeviceToHost));	
	SoaTree ht;
	int i,j;	
	ofstream output("output.txt");
	output.precision(5);
	output.setf(ios::fixed, ios::floatfield);	
	for(i=0; i<nTrees; i++){		
		ht.setOffs(nNodes, nInsSpc, hm+(treeSize*i));		
		for(j=0; j<nNodes; j++) output << ht.getIdx(j) 		<< "\t"; output << endl;
		for(j=0; j<nNodes; j++) output << ht.getBranch(j) 	<< "\t"; output << endl;
		for(j=0; j<nNodes; j++) output << ht.getParent(j) 	<< "\t"; output << endl;
		for(j=0; j<nNodes; j++) output << ht.getlChild(j) 	<< "\t"; output << endl;
		for(j=0; j<nNodes; j++) output << ht.getrChild(j) 	<< "\t"; output << endl;
		for(j=0; j<nNodes; j++) output << ht.getdRoot(j)	<< "\t"; output << endl << endl;
	}
}


class HTree: public DTree{	
	private:
		SoaTree hostData;			// struct of arrays to hold the trees' data 		
		vector<string> name; 		// names of taxa fetched from newickf and PUT file
		ifstream newickf;			// stream object to manage input newick file
		//ifstream putf;			// stream object to manage input PUT file
		int devId;					// id of the GPU where lies the tree 
		__host__ void setParams();
		__host__ void parseTree();
	public:
		__host__ HTree() = default;
		__host__ HTree(int dev_id=0, string nw_fname = "newick.tree"/*, string pt_fname="put.list" */);
		__host__ DTree& gpuRep(int num_reps) const;
		/* TODO: CATCH OVERFLOW EXCEPTION */
		__host__ void setTreeIdx(int val, int i)	{hostData.setTreeIdx(val,i);}
		__host__ void setHeight (int val, int i)   	{hostData.setHeight(val,i);}
		__host__ void setnSubt  (int val, int i)   	{hostData.setnSubt(val,i);}
		__host__ void setParent (int val, int i)   	{hostData.setParent(val,i);}
		__host__ void setlChild (int val, int i)   	{hostData.setlChild(val,i);}
		__host__ void setrChild (int val, int i)	{hostData.setrChild(val,i);}
		__host__ void setBranch (float val, int i)	{hostData.setBranch(val,i);}
		__host__ void setdRoot  (float val, int i)  {hostData.setdRoot(val,i);}
		__host__ void setInseq (int val, int i)  {hostData.setInseq(val,i);}
		
		__host__  int getIdx	 (int i) const {return hostData.getIdx(i);}
		__host__  int getHeight	 (int i) const {return hostData.getHeight(i);}
		__host__  int getnSubt 	 (int i) const {return hostData.getnSubt(i);}
		__host__  int getParent	 (int i) const {return hostData.getParent(i);}
		__host__  int getlChild	 (int i) const {return hostData.getlChild(i);}
		__host__  int getrChild	 (int i) const {return hostData.getrChild(i);}
		__host__  float getBranch(int i) const {return hostData.getBranch(i);}
		__host__  float getdRoot (int i) const {return hostData.getdRoot(i);}
		__host__ string getName	 (int i) const {return name[i];}		
		__host__ int getInseq    (int i) const {return hostData.getInseq(i);}
};

HTree::HTree(int dev_id, string nw_fname/*,string pt_fname*/){	
	devId = dev_id;
	newickf.open(nw_fname);
	FERR(newickf);
	//putf.open(pt_fname);
	//FERR(putf);
	setParams();
	treeSize = hostData.getSize(nNodes, nInsSpc);
	hostData.soalloc(nNodes, nInsSpc);
	parseTree();	
	newickf.close();
	//putf.close();
	void * d_tree;
	CHECK(cudaSetDevice(devId));
	CHECK(cudaMalloc(&d_tree, treeSize));
	CHECK(cudaMemcpy(d_tree, hostData.getPtr(), treeSize, cudaMemcpyHostToDevice));
	devData.setOffs(nNodes, nInsSpc, d_tree);
	nTrees=1;
}

bool DTree::compareTo(HTree *h_tree){	
	if(treeSize != h_tree->getSize() || idxInsSpc != h_tree->getIdxInsSpc() || idxInsAnc != h_tree->getIdxInsAnc())
		return false;
	SoaTree tree;
	size_t rep_size = treeSize * nTrees;	
	void* h_replics = malloc(rep_size);
	CHECK(cudaMemcpy(h_replics, devData.getPtr(), rep_size, cudaMemcpyDeviceToHost));
	cout.precision(3);
	cout.setf(ios::fixed, ios::floatfield);	
	for(int j=0; j<nTrees; j++){
		tree.setOffs(nNodes, nInsSpc, h_replics+treeSize*j);
		for(int i=0; i<nNodes; i++){
			if(	tree.getIdx(i) 	  != h_tree->getIdx(i)		||
				tree.getHeight(i) != h_tree->getHeight(i)	||
				tree.getBranch(i) != h_tree->getBranch(i) 	||
				tree.getParent(i) != h_tree->getParent(i)	||
				tree.getlChild(i) != h_tree->getlChild(i)	||
				tree.getrChild(i) != h_tree->getrChild(i)	)
					return false;
		}
	}
	return true;
}

void HTree::setParams(){	
	newickf.seekg(0);
	newickf >> nNodes >> idxInsAnc >> idxInsSpc >> nInsSpc;
	FERR(newickf);
}

void HTree::parseTree() {	
	string taxon_name;
	int ival;
	float fval;
	for(int i=0; i<nNodes; i++) {newickf >> ival; hostData.setTreeIdx(ival,i);}
	for(int i=0; i<nNodes; i++) {newickf >> taxon_name; name.push_back(taxon_name);}	
	for(int i=0; i<nNodes; i++) {newickf >> fval; hostData.setBranch(fval,i);}
	for(int i=0; i<nNodes; i++) {newickf >> ival; hostData.setHeight(ival,i);}	
	for(int i=0; i<nNodes; i++) {newickf >> ival; hostData.setnSubt(ival,i);}
	for(int i=0; i<nNodes; i++) {newickf >> ival; hostData.setParent(ival,i);}
	for(int i=0; i<nNodes; i++) {newickf >> ival; hostData.setlChild(ival,i);}
	for(int i=0; i<nNodes; i++) {newickf >> ival; hostData.setrChild(ival,i);}
	for(int i=0; i<nNodes; i++) {newickf >> fval; hostData.setdRoot(fval,i);}
	for(int i=0, ival=idxInsSpc; i<nInsSpc; i++){
		hostData.setInseq(ival+i,i);
	}
	FERR(newickf);
}

DTree& HTree::gpuRep(int num_reps) const{	
	size_t rep_size = treeSize * num_reps;
	void *d_replics;
	CHECK(cudaMalloc(&d_replics, rep_size));
	cudaDeviceProp device;
	CHECK(cudaGetDeviceProperties(&device,devId));
	int threads = device.warpSize*16;	//TODO: FIGURE OUT WHICH MULTIPLE IS THE BEST
	int blocks = (rep_size/sizeof(int4) + (threads-1)) / threads;
	dim3 grid = dim3(blocks);
	dim3 block = dim3(threads);
	modcpy<<<grid, block>>>(d_replics,devData.getPtr(),rep_size,treeSize);
	CHECK(cudaDeviceSynchronize());
	return *new DTree(nNodes,nInsSpc,idxInsSpc,idxInsAnc,num_reps,treeSize,d_replics);
}


__global__ void insertion(DTree tree, curandState_t* devStates){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	tree.setTreeIdx(idx);
    curandState state = devStates[idx];    
    unsigned int i,j,t;
    int taxon, mdcc, ancidx;
    float branching;
    float height = tree.getdRoot(0);	//height of the tree (distance from leaf to root)    
        
    if (tree.getnInsSpc() > 1) {
	for (i=0; i<tree.getnInsSpc()-1; i++) {
		j = i + curand(&state) / (UINT_MAX/(tree.getnInsSpc()-i)+1);
		t = tree.getInseq(j);
		tree.setInseq(tree.getInseq(i),j);
		tree.setInseq(t,i);
		}
    }    

    float sum;
	for(i=0; i<tree.getnInsSpc(); i++){
		t = curand(&state);	//path
		mdcc = tree.getParent(tree.getInseq(i));
		branching = curand_uniform(&state) * (height-tree.getdRoot(mdcc));
		taxon = mdcc;
		sum=0;
		do{
			taxon = t&1 ? tree.getlChild(taxon) : tree.getrChild(taxon);
			sum+= tree.getBranch(taxon);
			t>>=1;
		}while(sum<branching);		
		ancidx = tree.getIdxInsAnc()-(tree.getInseq(i)-tree.getIdxInsSpc());	//calculate corresponding ancestor node
		if(taxon==tree.getlChild(tree.getParent(taxon))){
			tree.setrChild(tree.getInseq(i),ancidx);
			tree.setlChild(taxon,ancidx);
			tree.setlChild(ancidx,tree.getParent(taxon));
		}			
		else{
			tree.setlChild(tree.getInseq(i),ancidx);
			tree.setrChild(taxon,ancidx);
			tree.setrChild(ancidx,tree.getParent(taxon));
		}
		tree.setParent(tree.getParent(taxon),ancidx);								 		//set up new ancestor's parent (same of the sister group)
		tree.setParent(ancidx,tree.getInseq(i));									 		//set up PUT's parent
		tree.setParent(ancidx,taxon);												 		//set up sister's new parent
		tree.setBranch(tree.getBranch(taxon)-(sum-branching),ancidx);				 		//set up new ancestor's branch
		tree.setBranch(sum-branching,taxon);											 	//set up sister's new branch length
		tree.setBranch(height-(tree.getdRoot(mdcc)+branching),tree.getInseq(i));	 		//set up PUT's branch length
		tree.setdRoot (tree.getdRoot(tree.getParent(ancidx))+tree.getBranch(ancidx),ancidx);//set up new ancestor's distance to the root		
	}    
}

__global__ void setup_kernel(long long seed, curandState_t* devStates){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;    
    curand_init(seed, idx, 0, &devStates[idx]);
}

int main(int argc, char *argv[]){
	
	if(argc < 2 || argc >3){
		cout << "Usage: " << argv[0] << " #replications [file]" << endl;
		exit(EXIT_FAILURE);
	}	
	int num_reps = atoi(argv[1]);	
	HTree *tree = argc==3 ? new HTree(0,argv[2]) : new HTree(0);
	DTree replics = tree->gpuRep(num_reps);	
	
	if(replics.compareTo(tree))
		cout << "Data does match!" << endl;
	else
		cout << "Data doesn't match" << endl;	
	
	curandState_t *devStates;	
	CHECK(cudaMalloc((void**)&devStates, sizeof(curandState_t)*num_reps));	
	
	setup_kernel<<<1,num_reps>>>(1,devStates);
	insertion<<<1,num_reps>>>(replics,devStates);
		
	replics.print();	
	
	CHECK(cudaDeviceReset());	
	exit(EXIT_SUCCESS);
}
