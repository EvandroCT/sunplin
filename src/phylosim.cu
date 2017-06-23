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

#include <phylosim.h>
#include <cudutils.h>

//create all necessary seeds to massive GPU randomize
__global__ void setup_kernel(long long seed, curandState_t* devStates, ushort N){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int i;
    for(i=idx;i<N;i+=gridDim.x*blockDim.x)
    	curand_init(seed, i, 0, &devStates[i]);
}

//trees' exapansions
__global__ void insertion(DTree tree, curandState_t* devStates){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	curandState state;
	unsigned int i,j,t;
	int taxon, mdcc;
	int ancidx; 	//the put's parent node created to represent the cladogenesis
	int grandpa;
	unsigned int k;

	float depth;	//depth in which the put will be inserted down the subtree rooted at mdcc
	float height;	//height of the tree (distance from leaf to root)
	
	for(k=idx;k<tree.getnTrees();k+=gridDim.x*blockDim.x){	
		tree.setTreeIdx(k);
	    state = devStates[k];
	    height = tree.getdRoot(0); //height of the tree (distance from leaf to root)

	    if (tree.getnInsSpc() > 1) {
		for (i=0; i<tree.getnInsSpc()-1; i++) {
			j = i + curand(&state) / (UINT_MAX/(tree.getnInsSpc()-i)+1);
			t = tree.getInseq(j);
			tree.setInseq(tree.getInseq(i),j);
			tree.setInseq(t,i);
			}
	    }
	    float sum;
	    ushort put; //current put going to be inserted
		for(i=0; i<tree.getnInsSpc(); i++){
			t = curand(&state);	//path
			put = tree.getInseq(i);
			mdcc = tree.getParent(put);	
			depth = curand_uniform(&state) * (height-tree.getdRoot(mdcc));
			taxon = mdcc;
			sum=0;
			do{		
				t>>=1;
				taxon = t&1 ? tree.getlChild(taxon) : tree.getrChild(taxon);
				sum+= tree.getBranch(taxon);			
			}while(sum<depth);
			//after the loop, taxon is the sister clade
			grandpa = tree.getParent(taxon);
			ancidx = tree.getIdxInsAnc()-(put-tree.getIdxInsSpc());	//calculate corresponding ancestor node		
			if(t&1){	//if came from the left
				tree.setrChild(put,ancidx);		//put become the right child
				tree.setlChild(taxon,ancidx);	//the sister clade continue being at left
				tree.setlChild(ancidx,grandpa);//the put's parent node takes place of the sister's clade side
			}			
			else{	//if came from the right
				tree.setlChild(put,ancidx);		//put become the left child
				tree.setrChild(taxon,ancidx);	//the sister clade continue being at right
				tree.setrChild(ancidx,grandpa);//the put's parent node takes place of the sister's clade side
			}
			tree.setParent(grandpa,ancidx);				//set up new ancestor's parent (same of the sister group)
			tree.setSide(t&1,ancidx);									//set up new ancestor's side (same of the sister group)
			tree.setParent(ancidx,put);									//set up PUT's parent
			tree.setSide(!(t&1),put);									//set up PUT's side (the sister's reverse)
			tree.setParent(ancidx,taxon);								//set up sister's new parent
			tree.setBranch(tree.getBranch(taxon)-(sum-depth),ancidx);	//set up new ancestor's branch
			tree.setBranch(sum-depth,taxon);							//set up sister's new branch length
			tree.setBranch(height-(tree.getdRoot(mdcc)+depth),put);		//set up PUT's branch length
			tree.setdRoot (tree.getdRoot(grandpa)+tree.getBranch(ancidx),ancidx);	//set up new ancestor's distance to the root
		}	
	}
}

//generate the patristic distance matrices to all the replics
__global__ void patrices(DTree tree, float* d_matrix){

		tree.setTreeIdx(blockIdx.x);
		uint idx = threadIdx.x;
		ushort row, col, taxon;
		unsigned long long row_bmp, col_bmp; 
		ushort row_len, col_len;
		ushort N = tree.getnNodes();
		ushort nleafs = (N+1)/2;
		uint msize = nleafs*(nleafs-1)/2;

		extern __shared__ ushort s[];

		ushort *parent = s;
		ushort *lchild = parent+N;
		ushort *rchild = lchild+N;

		uint i;

		//separated loops to favor coalesced access
		for(i=idx;i<N;i+=blockDim.x)
				parent[i] = tree.getParent(i);
		for(i=idx;i<N;i+=blockDim.x)
				lchild[i] = tree.getlChild(i);
		for(i=idx;i<N;i+=blockDim.x)
				rchild[i] = tree.getrChild(i);

		__syncthreads();

		for(i=idx;i<msize;i+=blockDim.x)
		{
			row=row_index(i,nleafs);
			col=column_index(i,nleafs);
			row_bmp=0;
			col_bmp=0;
			row_len=0;
			col_len=0;
			for(taxon=row; parent[taxon]!=NOPARENT; taxon=parent[taxon]){
				row_len++;
				row_bmp<<=1;
				row_bmp|=tree.getSide(taxon);
			}
			for(taxon=col; parent[taxon]!=NOPARENT; taxon=parent[taxon]){
				col_len++;
				col_bmp<<=1;
				col_bmp|=tree.getSide(taxon);
			}
			taxon=tree.getnNodes()-1; 		//start with the root
			if((row_bmp&1)==(col_bmp&1)){	//if the LCA isn't the root				
				do{
					taxon = row_bmp&1 ? lchild[taxon] : rchild[taxon]; // either row_bmp or col_bmp (same)
				 	row_bmp>>=1;
				 	col_bmp>>=1;
				 }while((row_bmp&1)==(col_bmp&1));
			}
			d_matrix[blockIdx.x*msize+i] = 2*(tree.getdRoot(row)-tree.getdRoot(taxon));
	}
}
