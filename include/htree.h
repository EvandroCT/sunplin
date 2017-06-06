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

#ifndef _HTREE_H
#define _HTREE_H

#include "dtree.h"
#include <fstream>
#include <vector>

using namespace std;

class HTree: public DTree{	
	private:
		SoaTree hostData;					// struct of arrays to hold the trees' data 		
		unordered_map<int, string> name;	// names of taxa fetched from newickf and PUT file
		ifstream newickf;					// stream object to manage input newick file
		ifstream putf;						// stream object to manage input PUT file
		int devId;							// id of the GPU where lies the tree 
		__host__ void setParams(string &fileLine, vector<string> &filePut);
		__host__ void parseTree(string fileLine, vector<string> filePut);
	public:
		__host__ HTree() = default;
		__host__ HTree(int dev_id=0, string nw_fname = "newick.tree", string pt_fname="put.list");
		__host__ DTree& gpuRep(int num_reps) const;
		
		/* TODO: THROW "OUT OF BOUNDS" EXCEPTION */		
		__host__ void setParent (int 	val, int i)	{hostData.setParent(val,i);}
		__host__ void setSide	(int 	val, int i)	{hostData.setSide(val,i);}
		__host__ void setlChild (int 	val, int i)	{hostData.setlChild(val,i);}
		__host__ void setrChild (int 	val, int i)	{hostData.setrChild(val,i);}
		__host__ void setBranch (float 	val, int i)	{hostData.setBranch(val,i);}
		__host__ void setdRoot  (float 	val, int i) {hostData.setdRoot(val,i);}
		__host__ void setInseq 	(int 	val, int i) {hostData.setInseq(val,i);}
		__host__ void setName 	(string val, int i) {name[i]=val;}
		
		/* TODO: THROW "OUT OF BOUNDS" EXCEPTION */
		__host__ ushort	getParent	(int i) const	{return hostData.getParent(i);}
		__host__ ushort	getSide		(int i) const	{return hostData.getSide(i);}
		__host__ ushort	getlChild	(int i) const	{return hostData.getlChild(i);}
		__host__ ushort	getrChild	(int i) const 	{return hostData.getrChild(i);}
		__host__ float	getBranch	(int i) const 	{return hostData.getBranch(i);}
		__host__ float	getdRoot	(int i) const 	{return hostData.getdRoot(i);}
		__host__ ushort	getInseq	(int i) const 	{return hostData.getInseq(i);}
		__host__ string	getName		(int i)			{return name[i];}
		__host__ unordered_map<int, string> getNames(){return name;}
};

#endif
