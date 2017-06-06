/*************************************************************************
	
	Copyright (C) 2017	Evandro C. Taquary
	
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

/* 
 * Data structure to hold both parent's id and the side of which the 
 * node belongs related its parent.
 */ 
typedef struct {
	ushort side	: 1;
	ushort idx	: 15;
} paren_t;

/*
 * All methods inlined for optimization purpose (there's no definition
 * file for the class yet).
 */
class SoaTree {
	private:
		paren_t	*parent;	// nodes' parents or the subtrees' roots' indices where new nodes shall be inserted (MDCC)
		ushort	*lChild;	// nodes' left children
		ushort	*rChild;	// nodes' right children
		float	*branch;	// lengths of the nodes' branches (distance to the parent)
		float	*dRoot;		// distances between nodes and root (sum of the paths' branches)
		ushort	*inseq;		// vector with the sequence of indices of puts to be inserted
	public:
		__host__ SoaTree() = default;
		__host__ SoaTree(int num_nodes, int num_ins) {soalloc(num_nodes,num_ins);}
		__host__ SoaTree(int num_nodes, void* base) {setOffs(num_nodes,base);}
		__host__ __device__ void* getPtr() const {return (void*) parent;}
		__host__ __device__ void setOffs(int num_nodes, void* base) //set pointers' offsets starting on base accordingly to data structure, # nodes and # insertions
		{
			parent 	= (paren_t*) base;
			lChild 	= (ushort*)	(parent+num_nodes);
			rChild 	= lChild	+ num_nodes;	
			branch 	=(float*)	(rChild+num_nodes);
			dRoot 	= branch	+ num_nodes;
			inseq 	=(ushort*)	(dRoot+num_nodes);
		}		
		__host__ void setOffs(int num_nodes) {setOffs(num_nodes, parent);}  //set pointers' offsets starting on the first array accordingly to data structure, # nodes and # insertions		
		__host__ static size_t getSize(int num_nodes, int num_ins)
		{
			size_t size = (3*sizeof(ushort) + 2*sizeof(float))*num_nodes + sizeof(ushort)*num_ins; //minimal amount of bytes needed to represent the tree 
			int r = size%sizeof(int4);
			size += r ? sizeof(int4)-r : 0;	//size of the tree padded to a multiple of sizeof(int4) (due to a GPU memory aligment requisite)
			return size;
		}
		
		__host__ void soalloc(int num_nodes, int num_ins)
		{
			void *ptr = malloc(getSize(num_nodes, num_ins));
			memset(ptr,0,getSize(num_nodes, num_ins));
			setOffs(num_nodes, ptr);
		}
		__host__ __device__ ushort	getParent	(int i) const {return parent[i].idx;}
		__host__ __device__ ushort	getSide		(int i) const {return parent[i].side;}
		__host__ __device__ ushort	getlChild	(int i) const {return lChild[i];}
		__host__ __device__ ushort	getrChild	(int i) const {return rChild[i];}
		__host__ __device__ ushort 	getInseq	(int i) const {return inseq[i];}
		__host__ __device__ float 	getBranch	(int i) const {return branch[i];}
		__host__ __device__ float	getdRoot	(int i) const {return dRoot[i];}

		__host__ __device__ void setParent	(ushort	val, int i)	{parent[i].idx	= val;}
		__host__ __device__ void setSide	(ushort	val, int i)	{parent[i].side	= val;}
		__host__ __device__ void setlChild	(ushort	val, int i)	{lChild[i]		= val;}
		__host__ __device__ void setrChild	(ushort	val, int i)	{rChild[i]		= val;}
		__host__ __device__ void setBranch	(float	val, int i)	{branch[i]		= val;}
		__host__ __device__ void setdRoot	(float	val, int i)	{dRoot[i]		= val;}
		__host__ __device__ void setInseq	(ushort	val, int i)	{inseq[i]		= val;}
};
