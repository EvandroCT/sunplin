/*************************************************************************
	
	Copyright (C) 2016	Evandro Taquary, Thiago Santos
	
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

#include <iostream>
#include <string>
#include <fstream>
#include <curand_kernel.h>
#include "modcpy.h"
#include <regex>
#include <unordered_map>
#include <iomanip>
#include <sys/time.h>

using namespace std;

#define START_TIMER() \
		{ \
			gettimeofday(&tv, NULL); \
			start_time = tv.tv_sec * 1000000 + tv.tv_usec; \
		}
//return time measurement in s
#define STOP_TIMER(time_spent) \
		{ \
			gettimeofday(&tv, NULL); \
			end_time = tv.tv_sec * 1000000 + tv.tv_usec; \
			time_spent = ((double)(end_time-start_time))/1000000; \
		}

#define CHECK(call) \
		{ \
			const cudaError_t error = call; \
			if (error != cudaSuccess) { \
				cout << "Error: " << __FILE__ ": " << __LINE__ << ", "; \
				cout << "code: "<< error << ", reason: " << cudaGetErrorString(error) << endl; \
				exit(EXIT_FAILURE); \
			} \
		}

#define FERR(file) \
		{ \
			if(!file.good()){ \
				cout << "Something went wrong while reading the file! Please try again." << endl; \
				cout << "Error: " << __FILE__ ": " << __LINE__ << ", " << endl; \
				exit(EXIT_FAILURE); \
			} \
		}

typedef struct {
	ushort side	: 1;
	ushort idx	: 15;
} paren_t;

#define NOCHILD USHRT_MAX		//16 bits
#define NOPARENT USHRT_MAX/2	//15 bits

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
		__host__ __device__ void setOffs(int num_nodes, void* base);	//set pointers' offsets starting on base accordingly to data structure, # nodes and # insertions
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

void SoaTree::setOffs(int num_nodes, void* base) {
	parent 	= (paren_t*) base;
	lChild 	= (ushort*)	(parent+num_nodes);
	rChild 	= lChild	+ num_nodes;	
	branch 	=(float*)	(rChild+num_nodes);
	dRoot 	= branch	+ num_nodes;
	inseq 	=(ushort*)	(dRoot+num_nodes);
}

class HTree;

class DTree{

	protected:
		void *base;
		ushort nNodes;		// quantity of nodes on the tree(s) (including inserting species)
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
		
		/* TODO: THROW OVER/UNDERFLOW EXCEPTION */
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

		/* TODO: THROW OVER/UNDERFLOW EXCEPTION */
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
		__host__ void free(){CHECK(cudaFree(devData.getPtr()))}
};

void DTree::print(unordered_map<int,string> names){

	size_t rep_size = treeSize*nTrees;	
	void* h_replics = malloc(rep_size);
	CHECK(cudaMemcpy(h_replics, base, rep_size, cudaMemcpyDeviceToHost));	
	SoaTree ht;
	string aux;
	int i,j;	
	cout.precision(4);
	cout.setf(ios::fixed, ios::floatfield);	
	cout << endl;
	for(i=0; i<nTrees; i++){
		cout<<"tree #"<<i<<endl;
		ht.setOffs(nNodes, h_replics+(treeSize*i));		
		for(j=0; j<nNodes; j++){
			aux = names[j]+"("+to_string(j)+")";
			cout << left << setw (40) << aux;
		} 
		cout << endl;
		for(j=0; j<nNodes; j++) {
			aux = ht.getParent(j)!=NOPARENT ? names[ht.getParent(j)]+"("+to_string(ht.getParent(j))+")" : "-1";
			cout << left << setw (40) << aux;
		} 
		cout << endl;
		for(j=0; j<nNodes; j++) {
			aux = ht.getSide(j)==1 ? "left" : "right";
			cout << left << setw (40) << aux;
		} 
		cout << endl;
		for(j=0; j<nNodes; j++) {
			aux = ht.getlChild(j)!=NOCHILD ? names[ht.getlChild(j)]+"("+to_string(ht.getlChild(j))+")" : "-2";
			cout << left << setw (40) << aux;	
		}
		cout << endl;
		for(j=0; j<nNodes; j++) {
			aux = ht.getrChild(j)!=NOCHILD ? names[ht.getrChild(j)]+"("+to_string(ht.getrChild(j))+")" : "-2";
			cout << left << setw (40) << aux;
		}
		cout << endl;
		for(j=0; j<nNodes; j++) cout << left << setw (40) << ht.getBranch(j); cout << endl;
		for(j=0; j<nNodes; j++) cout << left << setw (40) << ht.getdRoot(j); cout << endl;
		for(j=0; j<nInsSpc; j++)cout << names[ht.getInseq(j)] << "("<< ht.getInseq(j) <<") ";
		cout << endl << endl;
	}
}

void DTree::print(unordered_map<int,string> names, int i){

	size_t rep_size = treeSize*nTrees;
	void* h_replics = malloc(rep_size);
	CHECK(cudaMemcpy(h_replics, devData.getPtr(), rep_size, cudaMemcpyDeviceToHost));
	SoaTree ht;
	string aux;
	int j;
	cout.precision(4);
	cout.setf(ios::fixed, ios::floatfield);
	cout << endl;
	cout<<"tree #"<<i<<endl;
	ht.setOffs(nNodes, h_replics+(treeSize*i));
	for(j=0; j<nNodes; j++){
		aux = names[j]+"("+to_string(j)+")";
		cout << left << setw (40) << aux;
	}
	cout << endl;
	for(j=0; j<nNodes; j++) {
		aux = ht.getParent(j)!=NOPARENT ? names[ht.getParent(j)]+"("+to_string(ht.getParent(j))+")" : "-1";
		cout << left << setw (40) << aux;
	}
	cout << endl;
	for(j=0; j<nNodes; j++) {
		aux = ht.getSide(j)==1 ? "left" : "right";
		cout << left << setw (40) << aux;
	}
	cout << endl;
	for(j=0; j<nNodes; j++) {
		aux = ht.getlChild(j)!=NOCHILD ? names[ht.getlChild(j)]+"("+to_string(ht.getlChild(j))+")" : "-2";
		cout << left << setw (40) << aux;
	}
	cout << endl;
	for(j=0; j<nNodes; j++) {
		aux = ht.getrChild(j)!=NOCHILD ? names[ht.getrChild(j)]+"("+to_string(ht.getrChild(j))+")" : "-2";
		cout << left << setw (40) << aux;
	}
	cout << endl;
	for(j=0; j<nNodes; j++) cout << left << setw (40) << ht.getBranch(j); cout << endl;
	for(j=0; j<nNodes; j++) cout << left << setw (40) << ht.getdRoot(j); cout << endl;
	for(j=0; j<nInsSpc; j++)cout << names[ht.getInseq(j)] << "("<< ht.getInseq(j) <<") ";
	cout << endl << endl;

}

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
		
		/* TODO: THROW OVER/UNDERFLOW EXCEPTION */		
		__host__ void setParent (int 	val, int i)	{hostData.setParent(val,i);}
		__host__ void setSide	(int 	val, int i)	{hostData.setSide(val,i);}
		__host__ void setlChild (int 	val, int i)	{hostData.setlChild(val,i);}
		__host__ void setrChild (int 	val, int i)	{hostData.setrChild(val,i);}
		__host__ void setBranch (float 	val, int i)	{hostData.setBranch(val,i);}
		__host__ void setdRoot  (float 	val, int i) {hostData.setdRoot(val,i);}
		__host__ void setInseq 	(int 	val, int i) {hostData.setInseq(val,i);}
		__host__ void setName 	(string val, int i) {name[i]=val;}
		
		/* TODO: THROW OVER/UNDERFLOW EXCEPTION */
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

HTree::HTree(int dev_id, string nw_fname, string pt_fname){	
	long long start_time, end_time;
	struct timeval tv;

	void * d_tree;
	double time_spent;
	devId = dev_id;
	nTrees=1;
	CHECK(cudaSetDevice(devId));

	START_TIMER();
	newickf.open(nw_fname);
	FERR(newickf);
	putf.open(pt_fname);
	FERR(putf);	
	string fileLine;
	vector<string> filePut;
	setParams(fileLine,filePut);	
	hostData.soalloc(nNodes,nInsSpc);
	treeSize = hostData.getSize(nNodes,nInsSpc);
	parseTree(fileLine,filePut);
	newickf.close();
	putf.close();
	STOP_TIMER(time_spent);
	cout<<"\ntotal time spent to parse the files: "<<time_spent<<"s\n";	

	//make a copy of the tree on device side
	START_TIMER();
	CHECK(cudaMalloc(&d_tree, treeSize));
	CHECK(cudaMemcpy(d_tree, hostData.getPtr(), treeSize, cudaMemcpyHostToDevice));	
	STOP_TIMER(time_spent);
	cout<<"\ntotal time spent to copy backbone tree to GPU: "<<time_spent<<"s\n";
	base=d_tree;
	devData.setOffs(nNodes, d_tree);
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
		tree.setOffs(nNodes, h_replics+treeSize*j);
		for(int i=0; i<nNodes; i++){
			if(	tree.getdRoot(i)	!= h_tree->getdRoot(i)	||
				tree.getBranch(i)	!= h_tree->getBranch(i) ||
				tree.getParent(i)	!= h_tree->getParent(i)	||
				tree.getSide(i)		!= h_tree->getSide(i)	||
				tree.getlChild(i)	!= h_tree->getlChild(i)	||
				tree.getrChild(i)	!= h_tree->getrChild(i)	)
					return false;
		}
	}
	return true;
}

void HTree::setParams(string &fileLine, vector<string> &filePut){	

	int fileLines=0;
	int aParen=0, fParen=0, comma=0;
	int quantElementosFile;
	char c;
	string currElement, aux;
	while (newickf.get(c)) {
		fileLine +=c;
		fileLines++;

	}
	nInsSpc = 0; // inicializar durante a construção

	while (getline (putf,aux)) //enquanto end of file for false continua
    {      
      filePut.push_back(aux);                             
      nInsSpc++;
    }

	quantElementosFile = fileLines; // qnts elementos o arquivo tem	
	// primeira varredura apenas para verificar inconsistencias
	for(int i = 0; i < quantElementosFile; i++){ // faz uma varredura no arquivo
		currElement = fileLine[i];
		if(currElement == "(") aParen++;
		if(currElement == ")") fParen++;
		if(currElement == ",") comma++;
	}
	if(aParen != fParen){

		cout<< "Arquivo inconsistente, parentes não balanceados" <<endl;
		throw;
	}	
	idxInsSpc = aParen +1; // nos folhas
	nNodes = (aParen * 2) + (nInsSpc * 2) +2;
	idxInsAnc = nNodes - aParen -1;	
}


void HTree::parseTree(string fileLine, vector<string> filePut) {

	int posParent = -1;
	string leaf =" ", ancestral =" ", currElement=" ", parent = " "; // salva o atual e o ultimo elemento
	string leftChild=" ", rightChild= " ", comprimeRamoLeft ="", comprimeRamoRight = "";
	int auxiliarNumNos =0, auxiliarGeral =0, auxilarPreencherVetor =0; // usado para fazer as trocas de elementos no vetor
	int indexleftChild =-1, indexrightChild =-1;
	bool alphabeticModeOn = false; 
	
	// regex
	int quantFolhas = idxInsSpc;
	
	smatch m;
  	regex e ("\\([^()]+\\)");
  	regex folhas("\\([A-z0-9_+.#]+|,[A-z0-9_+.#]+"); // achar todas as folhas e separar no vetor
  	regex internos("\\)[A-z0-9_+.#]+|\\)[:;]");

  	// fill empty names
    size_t pos = 0;
    int n_unamed=1;
    while ((pos = fileLine.find("):", pos)) != std::string::npos) {
         fileLine.replace(pos, 2, ")#"+to_string(n_unamed)+":");
         pos += to_string(n_unamed).length()+3;
         n_unamed++;
    }
    // fill root's empty name (if absent)
    for(pos=fileLine.length(); fileLine[pos]!=';'; pos--);
    if(fileLine[--pos]==')')
    	fileLine.replace(pos,2,")#"+to_string(n_unamed)+";");
    // fill new ancestors' names
    for(int i=0; i<getnInsSpc();i++)
    	setName("na#"+to_string(i+1),getIdxInsAnc()-i); //the new ancestors' insertions order is backward oriented
  	cout << "FileLine: " << fileLine << endl << endl;	
	for(int i=0;i<nNodes;i++){		
		setParent(NOPARENT,i);
		setlChild(NOCHILD,i);
		setrChild(NOCHILD,i);
		setBranch(0,i);
		setdRoot(0,i);		
	}
  	// preencher vetor com todas as species
	// usando o regex para pegar todos os quantFolhas	

	string copyNewick = fileLine;
	while (std::regex_search (copyNewick,m,folhas)) {
	    for (int i=0; i<m.size(); ++i) {
	    	auxiliarGeral = m.position(i)+1; // posicão do match (sem o '(' ou ',')
	    	leaf = copyNewick[auxiliarGeral++];	    	
	    	while(copyNewick[auxiliarGeral]!=':')
	    		leaf += copyNewick[auxiliarGeral++];	    	
  		}
		setName(leaf,auxilarPreencherVetor++);
	    copyNewick = m.suffix().str();
  	}
  	// preencher vetor com todas as species
	// usando o regex para pegar todos os nos internos	
	auxilarPreencherVetor = quantFolhas + (nInsSpc * 2) + 1;
	copyNewick = fileLine;
	while (std::regex_search (copyNewick,m,internos)) {
		ancestral = "";
	    for (int i=0; i<m.size(); ++i) {

	    	auxiliarGeral = m.position(i) +1; // posicão do match
		    while(copyNewick[auxiliarGeral]!=':' && copyNewick[auxiliarGeral]!=';') {
		    	ancestral += copyNewick[auxiliarGeral++];	    			
		    }		    
  		} 
  		setName(ancestral,auxilarPreencherVetor);
  		auxilarPreencherVetor++;
	    copyNewick = m.suffix().str();
  	}  	  	 
  	setParent(NOPARENT,nNodes-1); // no raiz não tem um pai
  	
	// logica se da no principio de achar todos os nos folhas pares, em cada loop, dai verificamos o seu devido pai
	// e os "eliminamos" da arvore, criando novos filhos folhas.
	// Para isso, estamos usando a biblioteca Redex, para achar os matchs e fazer o replace em seguida.
	// links: http://www.cplusplus.com/reference/regex/regex_search/
	//		  http://www.cplusplus.com/reference/regex/match_results/position/
	//        http://www.cplusplus.com/reference/regex/regex_replace/
	
	
	//regex logica
	// enquanto tivermos nos para buscar, vamos tirar as folhas
	// sobrara no final apenas o pai raiz
 
	int numTotalNos = nNodes-(2*nInsSpc)-1; 
	while(auxiliarNumNos < numTotalNos -1){	

		leftChild = "";
		rightChild = "";
		comprimeRamoLeft = "";
		comprimeRamoRight = "";
		std::regex_search ( fileLine, m, e );
    	
    	currElement = fileLine[m.position(0)]; // primeiro paranteses dos nos folhas achados    	
    	auxiliarGeral = m.position(0);

    	/* read everything until the ':' charactere is reached */
    	while(fileLine[++auxiliarGeral]!=':')
    		if(fileLine[auxiliarGeral]!=' ') leftChild += fileLine[auxiliarGeral];

   		while(fileLine[++auxiliarGeral]!=',')
    		if(fileLine[auxiliarGeral]!=' ') comprimeRamoLeft += fileLine[auxiliarGeral];

    	while(fileLine[++auxiliarGeral]!=':')
    		if(fileLine[auxiliarGeral]!=' ') rightChild += fileLine[auxiliarGeral];

   		while(fileLine[++auxiliarGeral]!=')')
    		if(fileLine[auxiliarGeral]!=' ') comprimeRamoRight += fileLine[auxiliarGeral];

	    auxiliarGeral++;

    	/* fetch name of the internal node (until ':') or of the root (until ';') */
    	parent="";
    	while(fileLine[auxiliarGeral]!=':' && fileLine[auxiliarGeral]!=';') {
	    	parent += fileLine[auxiliarGeral++];	    			
	    }	 

  		// achar o index entao dos filhos tirados e do pai
    	for(int i=0; i<nNodes; i++){
    		if(name[i] == parent){
    			posParent = i;
    			if( (indexleftChild != -1) and (indexrightChild != -1) ) break; // parar se ja achou indexes
    		}
    		else if(name[i]==rightChild){    			
    			indexrightChild = i;
    			if( (indexleftChild != -1) and (posParent != -1) ) break; 
    		}
    		else if(name[i]==leftChild){    			
    			indexleftChild = i;
    			if( (indexrightChild != -1) and (posParent != -1) ) break;
    		}
    	}
    	// preencher vetores
    	setParent(posParent,indexleftChild);
    	setSide(1,indexleftChild);
    	setParent(posParent,indexrightChild);
    	setSide(0,indexrightChild);
    	setlChild(indexleftChild,posParent);
    	setrChild(indexrightChild,posParent);
    	// comprimento do ramo
    	try{
	    	setBranch(atof(comprimeRamoRight.c_str()),indexrightChild);
	    	setBranch(atof(comprimeRamoLeft.c_str()),indexleftChild);
    	}catch(exception e){

    	}

	  	fileLine = m.prefix().str()+m.suffix().str();

	  	posParent = -1;
	  	// reset variaveis
  		rightChild = "";
  		leftChild = "";
  		comprimeRamoLeft = "";
  		comprimeRamoRight = "";
  		indexrightChild = -1;
  		indexleftChild = -1;
		auxiliarNumNos = auxiliarNumNos + 2; // ou seja, foi retirado 2 filhos
	}
	 // preencher novos put
 	string auxiliarPut[2], auxiliar, put;
  	for (int linePut = 0; linePut < nInsSpc; linePut++)
  	{  	
  		auxiliar = filePut[linePut];
  		put = ""; 
  		auxiliarGeral = 0;
  		alphabeticModeOn = false; 		
	    for (int elemenIndex = 0; elemenIndex < auxiliar.length(); elemenIndex++)
	    {	
	       if (isspace(auxiliar[elemenIndex]) and alphabeticModeOn) 
	       {
	       		auxiliarPut[auxiliarGeral++] = put;	       		
	           	put = "";
	           	alphabeticModeOn = false;
	       }else{
	       		if ( !isspace(auxiliar[elemenIndex]) ){ 
	       			alphabeticModeOn = true;
	        		put += auxiliar[elemenIndex];
	        	}	        	
	       }	 
	    }
	    if(put != ""){
	    	auxiliarPut[auxiliarGeral] = put;
	    }
	    //insert no array especies
	    setName(auxiliarPut[0],quantFolhas+linePut);
	    for (int index = 0; index < nNodes; index++)
	    {
	    	if(name[index] == auxiliarPut[1]){
	    		if(index>=getIdxInsSpc())
	    			setParent(index,quantFolhas+linePut);
	    		else // if the MDCC is a leaf, make its parent become the new MDCC 
	    			setParent(getParent(index),quantFolhas+linePut);
	    		break;
	    	}
	    }	
  	} 
	// Calcular comprimento do ramo ate a raiz
	// usando busca em profundidade
	bool folhaDone = false;
	int visited=0;
	setBranch(0,nNodes-1);	//root has no branch
	setdRoot(0,nNodes-1); 	//root has no distance to himself
	int posRamo = getrChild(nNodes-1);//start with the root's right child;
	while(visited<quantFolhas*2-2){
		// primeiramente, faz uma busca profunda, pela esquerda(mas na vdd tanto faz), e busca um no leaf
		// com isso, sabemos a profundidade de todos os outros folhas, restando então apenas os nos internos
		// essa regra se aplica apenas para arvores filogeneticas
		while(not folhaDone){
			setdRoot(getdRoot(getParent(posRamo))+getBranch(posRamo),posRamo);
			if(getrChild(posRamo) == NOCHILD){ // ou seja, não tem filho(leaf)
				folhaDone = true;
				// temos então o comprimento de todos os folhas da arvore
				// atualizar de todas as folhas então				
				for (int i = 0; i < quantFolhas + nInsSpc; i++)
				{
					setdRoot(getdRoot(posRamo),i);
				}
				visited+= quantFolhas;
				posRamo = getParent(posRamo); // volta entao a posição ramo 1 posição, pois chegou no limite da arvore(leaf)				
				break;
			}
			visited++;
			posRamo = getrChild(posRamo); // proximo filho a direita
		}		
		// fazer a busca em profundidade agr para os nos internos
		// se os dois filhos da raiz, ja tiverem seus comprimentos achados,
		// entao significa q a busca em profundidade foi concluida

		// cheka se elemento atual ainda tem filho 
		if(getrChild(posRamo)!=NOCHILD){
			// se tiver filho da direita e o comp dele ainda n foi calculado
			if(getdRoot(getrChild(posRamo))==0){
				// nova posRamo é entao aquele filho da direita
				visited++;
				posRamo = getrChild(posRamo);
				setdRoot(getdRoot(getParent(posRamo))+getBranch(posRamo),posRamo);
			}
			 // se tiver filho da esquerda e o comp dele ainda n foi calculado
			else if(getdRoot(getlChild(posRamo))==0){
				// nova posRamo é entao aquele filho da direita
				visited++;
				posRamo = getlChild(posRamo);
				setdRoot(getdRoot(getParent(posRamo))+getBranch(posRamo),posRamo);
			}
			// ou seja, aquela sub arvore esta concluida
			else{				
				posRamo = getParent(posRamo); // volta entao a posição ramo 1 posição, pois chegou no limite da arvore(leaf)	
			}	
		}
	}
	//setup insertion sequence
	for(int i=0; i<getnInsSpc();i++)
		setInseq(getIdxInsSpc()+i,i);
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

__global__ void setup_kernel(long long seed, curandState_t* devStates, ushort N){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int i;
    for(i=idx;i<N;i+=gridDim.x*blockDim.x)
    	curand_init(seed, i, 0, &devStates[i]);
}

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

__host__ __device__ int row_index( int i, int M ){ // retorna o indice da linha
	M--;
    float m = M;
    float row = (-2*m - 1 + sqrt( (4*m*(m+1) - 8*(float)i - 7) )) / -2;
    if( row == (float)(int) row ) row -= 1;
    return (int) row;
}

__host__ __device__ int column_index( int i, int M ){ // retorna o indice da coluna
    int row = row_index( i, M);
    M--;
    return 1 + (i - M * row + row*(row+1) / 2);
}

__global__ void patrix(DTree tree, float* d_matrix){

		tree.setTreeIdx(blockIdx.x);
		uint idx = threadIdx.x;
		ushort row, col, taxon;
		unsigned long long row_bmp, col_bmp; 
		ushort row_len, col_len;
		ushort N = tree.getnNodes();
		ushort nleafs = (N-1)/2;
		uint msize = nleafs*(nleafs-1)/2;

		extern __shared__ ushort s[];

		ushort *parent = s;
		ushort *lchild = parent+N;
		ushort *rchild = lchild+N;

		 uint i;
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
			taxon=tree.getnNodes()-1; 	//start with the root
			if((row_bmp&1)==(col_bmp&1)){	//if the LCA isn't the root

				//printf("\nrow=%d, col=%d\n",row,col);

				//printf("\nrow_bmp=%llu, col_bmp=%llu\n",row_bmp, col_bmp);

				do{
					taxon = row_bmp&1 ? lchild[taxon] : rchild[taxon]; // either row_bmp or col_bmp (same)
					//printf("taxon: %d\n",taxon);
				 	row_bmp>>=1;
				 	col_bmp>>=1;
				 }while((row_bmp&1)==(col_bmp&1));
			}
			d_matrix[blockIdx.x*msize+i] = tree.getdRoot(row)+tree.getdRoot(col)-2*tree.getdRoot(taxon);
	}
}

int main(int argc, char *argv[]){	

	if(argc < 2 || argc >4){
		cout << "Usage: " << argv[0] << " #replications [newick putlist]" << endl;
		exit(EXIT_FAILURE);
	}
	
	long long start_time, end_time;
	struct timeval tv;

	int gpu=0;
	double time_spent;
	int num_reps = atoi(argv[1]);	
	HTree *tree = argc>2 ? new HTree(gpu,argv[2],argv[3]) : new HTree(gpu);
	
	CHECK(cudaSetDevice(gpu));
	START_TIMER();
	DTree replics = tree->gpuRep(num_reps);
	STOP_TIMER(time_spent);
	cout<<"\ntotal time spent to replicate trees: "<<time_spent<<"s\n";	

	cout << "nNodes: " << tree->getnNodes() << endl;
	cout << "nInsSpc: " << tree->getnInsSpc() << endl;
	cout << "idxInsSpc: " << tree->getIdxInsSpc() << endl;
	cout << "idxInsAnc: " << tree->getIdxInsAnc() << endl << endl;
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
	STOP_TIMER(time_spent);
	cout<<"\ntotal time spent to expand trees: "<<time_spent<<"s\n";	
/*
	replics.print(tree->getNames(),0);
	replics.print(tree->getNames(),1);
	replics.print(tree->getNames(),2);

	START_TIMER();
	float *d_matrix;
	ushort nleafs = (replics.getnNodes()-1)/2;
	ushort msize = nleafs*(nleafs-1)/2;
	CHECK(cudaMalloc((void**)&d_matrix, sizeof(float)*msize*num_reps));
	patrix<<<num_reps,256,replics.getnNodes()*(sizeof(ushort)*3)>>>(replics, d_matrix);
	CHECK(cudaDeviceSynchronize());
	STOP_TIMER(time_spent);
	cout<<"\ntotal time spent to generate patrixes: "<<time_spent<<"s\n";	

	//replics.free();

	START_TIMER();
	float *h_matrix = (float*)malloc(sizeof(float)*msize*num_reps);
	CHECK(cudaMemcpy(h_matrix, d_matrix, sizeof(float)*msize*num_reps, cudaMemcpyDeviceToHost));
	STOP_TIMER(time_spent);
	cout<<"\ntotal time spent to copy patrixes to CPU: "<<time_spent<<"s\n";	
*/
	CHECK(cudaDeviceReset());
	exit(EXIT_SUCCESS);	
}
