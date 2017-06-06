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

#include <dtree.h>
#include <cudutils.h>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

void DTree::toNewick(unordered_map<int,string> names){ 

	size_t rep_size = treeSize*nTrees;	
	void* h_replics = malloc(rep_size);
	CHECK(cudaMemcpy(h_replics, devData.getPtr(), rep_size, cudaMemcpyDeviceToHost));	
	SoaTree ht;
	int indexThree;
	string  newickFile="";
	for(indexThree=0; indexThree<nTrees; indexThree++){ // total of threes
		ht.setOffs(nNodes, h_replics+(treeSize*indexThree));// get the pointer set for a tree i
		// Ainda precisa resetar as variaveis str para cada nova interação
		newickFile += "Tree ";
		newickFile += to_string(indexThree +1);
		newickFile += "\n";
		newickFile += calculateNewick(names, ht, nNodes -1);
		newickFile += ";\n\n"; 
		// have save the newick in to a file
	}
	
	newickToFile(newickFile);
}

string DTree::calculateNewick(unordered_map<int,string> names, SoaTree ht, int idRaiz ){ 
	string  str_tmp,  str_float;

	//str_tmp = "";
   	//str_float = ""; // idRaiz = nNodes -1 , get the last element in the vector, that's the root
	if ( ht.getlChild(idRaiz)  == NOCHILD) { // left child of the root // Não tem filhos // nz_f1 = filho da esquerda do no
		if ((idRaiz) < 0 || (idRaiz) > (nNodes-1)) // num of nodes
			printf("ERRO %d\n", (idRaiz));
		else
			str_tmp += names[idRaiz]; 
		str_tmp += ":";
		str_float += to_string(ht.getBranch(idRaiz)); 
		str_tmp += str_float;
		return str_tmp;
	} else { // Has child 
		str_tmp +="(";
		str_tmp += calculateNewick(names, ht, ht.getlChild(idRaiz)); 
		str_tmp += ",";
		str_tmp +=  calculateNewick(names, ht, ht.getrChild(idRaiz)); 
		str_tmp += ")";
		str_tmp += names[idRaiz];

		if(nNodes -1 != idRaiz){ // if the element is not the root
			str_tmp +=":";
			str_float += to_string( ht.getBranch(idRaiz)); 
			str_tmp += str_float;
		}
		
		return str_tmp;
	}
	
}

void DTree::newickToFile(string newick ){ 	
	ofstream ofFile;
	ofFile.open("versions.tree");
	ofFile<<newick;
	ofFile.close();
}

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
void DTree::free(){CHECK(cudaFree(devData.getPtr()))}
