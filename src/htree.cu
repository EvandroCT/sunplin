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

#include <htree.h>
#include <cudutils.h>
#include <iostream>
#include <regex>

using namespace std;

extern double time_spent[7];

HTree::HTree(int dev_id, string nw_fname, string pt_fname){
	void * d_tree;
	devId = dev_id;
	nTrees=1;
	CHECK(cudaSetDevice(devId));
	SETUP_TIMER(true);
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
	STOP_TIMER(time_spent[0]);
	//make a copy of the tree on device side
	START_TIMER();
	CHECK(cudaMalloc(&d_tree, treeSize));
	CHECK(cudaMemcpy(d_tree, hostData.getPtr(), treeSize, cudaMemcpyHostToDevice));	
	STOP_TIMER(time_spent[1]);
	base=d_tree;
	devData.setOffs(nNodes, d_tree);
}

//compare argument tree to all the trees within object
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
    while ((pos = fileLine.find("):", pos)) != string::npos) {
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
  	//cout << "FileLine: " << fileLine << endl << endl;	
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
	while (regex_search (copyNewick,m,folhas)) {
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
	while (regex_search (copyNewick,m,internos)) {
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
		regex_search ( fileLine, m, e );
    	
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

//creates 'num_reps' replics of the tree holded by the object, inside GPU Global memory, and return a reference to them
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
