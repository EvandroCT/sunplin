/*************************************************************************
	
	Copyright (C) 2017	Evandro C. Taquary, Thiago Santos
	
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

#ifndef _DTREE_H
#define _DTREE_H

#include "soatree.h"
#include <unordered_map>

#define NOCHILD USHRT_MAX		//16 bits
#define NOPARENT USHRT_MAX/2	//15 bits

using namespace std;

class HTree;

class DTree{

	protected:
		void *base;
		ushort nNodes;		// quantity of nodes in the tree(s) (including inserting species)
		ushort nInsSpc;		// quantity of absent species to be inserted
		ushort idxInsSpc;	// starting index for insertion of new species
		ushort idxInsAnc;	// starting index for insertion of new ancestors
		uint nTrees;		// quantity of trees holded by devData (default=1)
		size_t treeSize;	// size of one tree padded to multiple of sizeof(int4) (due to a GPU memory aligment requisite)
		SoaTree devData;	// struct of arrays to hold trees' data
	public:

		__host__ bool compareTo(HTree *h_tree);
		__host__ DTree() = default;
		__host__ DTree(int nNodes, int nInsSpc, int idxInsSpc, int idxInsAnc, int nTrees, size_t treeSize, void* ptr):
						nNodes(nNodes),
						nInsSpc(nInsSpc),
						idxInsSpc(idxInsSpc), 
						idxInsAnc(idxInsAnc),
						nTrees(nTrees),
						treeSize(treeSize){base=ptr; devData.setOffs(nNodes,ptr);}
		
		/* TODO: THROW "OUT OF BOUNDS" EXCEPTION */
		__host__ __device__ ushort	getnNodes	()		const {return nNodes;}
		__host__ __device__ ushort	getnInsSpc	()		const {return nInsSpc;}
		__host__ __device__ ushort	getIdxInsSpc()		const {return idxInsSpc;}
		__host__ __device__ ushort	getIdxInsAnc()		const {return idxInsAnc;}
		__host__ __device__ uint	getnTrees	()		const {return nTrees;}
		__host__ __device__ size_t	getSize		()		const {return treeSize;};

		__host__ __device__ ushort	getParent	(int i)	const {return devData.getParent(i);}
		__host__ __device__ ushort	getSide		(int i)	const {return devData.getSide(i);}
		__host__ __device__ ushort	getlChild	(int i)	const {return devData.getlChild(i);}
		__host__ __device__ ushort	getrChild	(int i)	const {return devData.getrChild(i);}
		__host__ __device__ float	getBranch	(int i)	const {return devData.getBranch(i);}
		__host__ __device__	float	getdRoot	(int i) const {return devData.getdRoot(i);}
		__host__ __device__ ushort	getInseq	(int i)	const {return devData.getInseq(i);}

		__device__ void	setTreeIdx(int i){devData.setOffs(nNodes,base+treeSize*i);}

		/* TODO: THROW "OUT OF BOUNDS" EXCEPTION */
		__device__ void	setParent	(ushort	val, int i)	{devData.setParent(val,i);}
		__device__ void	setSide		(ushort	val, int i)	{devData.setSide(val,i);}
		__device__ void	setlChild	(ushort	val, int i)	{devData.setlChild(val,i);}
		__device__ void	setrChild	(ushort	val, int i)	{devData.setrChild(val,i);}
		__device__ void	setBranch	(float	val, int i)	{devData.setBranch(val,i);}
		__device__ void	setdRoot	(float	val, int i)	{devData.setdRoot(val,i);}
		__device__ void	setInseq	(ushort	val, int i)	{devData.setInseq(val,i);}
		
		/* copy from GPU all the trees holded by the object and print them on the standard output  */
		__host__ void print(unordered_map<int,string> names);
		__host__ void print(unordered_map<int,string> names, int i);
		__host__ void free();
		__host__ void toNewick(unordered_map<int,string> names);
		__host__ string calculateNewick(unordered_map<int,string> names,  SoaTree ht, int idRaiz);
		__host__ void newickToFile(string newick);

};

#endif
