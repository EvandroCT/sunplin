/*********************************************************************
11	
12	 Copyright (C) 2013-2016 by Welton Cardoso
13	
14	 This program is free software; you can redistribute it and/or modify
15	 it under the terms of the GNU General Public License as published by
16	 the Free Software Foundation; either version 2 of the License, or
17	 (at your option) any later version.
18	
19	 This program is distributed in the hope that it will be useful,
20	 but WITHOUT ANY WARRANTY; without even the implied warranty of
21	 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
22	 GNU General Public License for more details.
23	
24	 You should have received a copy of the GNU General Public License
25	 along with this program; if not, write to the Free Software
26	 Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
27	
28	 ********************************************************************/

#include <map>
#include <cmath>
#include <ctime>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <utility>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <ctype.h>
#include <sys/timeb.h>
using namespace std;

const double INF = 1e50;
const int MAXN = 100010;
const int MAXGENERATENAMES = 1000;
const int MAXV = 98;
const  double EPS = 1e-9;
int cmp(double a, double b = 0.0){ if(fabs(a-b) < EPS) return 0; return a > b ? 1 : -1; }
struct timeb ini, fim;

class insert{
	public:
		string name;
		int insertionPoint;
		insert(string name = "", int insertionPoint = 0) : name(name), insertionPoint(insertionPoint) { }
};

class subTree{
	public:
		double sumSubTree;
		int sizeSubTree, rootCount;
		subTree(double sumSubTree = 0, int sizeSubTree = 0, int rootCount = 0) : sumSubTree(sumSubTree), sizeSubTree(sizeSubTree), rootCount(rootCount) { }
};

class infoNode{
	public:
		double valor;
		int id;
		infoNode( double valor = 0., int id = 0 ) : valor(valor), id(id) { }
};

bool comp( const infoNode &b, const infoNode &a ) {
	if( cmp( b.valor, a.valor ) > 0 ) return false;
	if( cmp( b.valor, a.valor ) == 0 ) return b.id < a.id;
	return true;
}

int nodedad[ MAXN ];
int treesize[ MAXN ];
int chain[ MAXN ];
int homepos[ MAXN ];
int up[MAXN];
int down[7000][7000];
int sumAux[MAXN];
int noat[MAXN];
int depth[MAXN];
int nodes;
int sometimesGenerates;
int treeAmount;
int amountNewNames = 0;
int camNode = 0;
int pos, cntchain, vezes = 1;
int chainleader[ MAXN ];
char taxon[7000][100];
char line[MAXN];
char taxonCopia[7000][100];
char version;
double bl[MAXN];
double blCopy[MAXN];
double dist[7000][7000];
double custAnt = 0.;
bool dfstop;

void generateNamesSpecies();
void errorType(string s1, string s2, string s3);
void explore( int x, int dad );
void heavy_light( int x, int dad, int k, int p );
void parsingTree();
void menu();
void clearDisplay();
void newSpeciesFile(char local[100], string nameNewFile);
void insertSpecies(int u, string name);
void convertNewickToNexus(int op, fstream &nexus, string nw);
void calcMatrixDistance();
void removeEdge(int fatherAux, int child);
void printToMatrixDistanceInFile(string nameNewFile);
void printTime();
int lca(int a, int b) ;
int dfs_insert_V1( int u );
int searchNameEspecie(char *str);
double dfs_sum_V2( int u );
string grafoToNewick(int u);


subTree sum[MAXN];
vector < insert > espN; 
vector < infoNode > nodeInfo;
vector < insert > espNames; 
vector < vector<int > > grafo, grafoMatrix;
vector < string > nameRand;
vector < int >  father;


int main(int argc, char* argv[]){
	srand ( time(NULL) );
	string local = "";
	ftime(&ini);
	if( argc <= 2 ) menu();
	else{
		if( !strcmp(argv[1],"distMat") ){
			if ( argc > 5 ) menu();
			else if ( argc > 3 ){
				if( !strcmp(argv[3],"-l") && argc > 4 ) local += argv[4];
				else menu();
			}
			else local += "distMatrix";
			fstream file (argv[2], fstream::in | fstream::out);
			if(!file.good()) errorType("Error! The file ",argv[2]," does not exist\n");
			file >> line;
			file.close();
			parsingTree();
			calcMatrixDistance();
			printToMatrixDistanceInFile(local.c_str());
		}
		else if( !strcmp(argv[1],"insertTree") ){
			generateNamesSpecies();
			if ( argc >= 9 || argc < 6) menu();
			else if ( argc > 6 ){
				if( !strcmp(argv[6], "-l" ) && argc > 7 ) local += argv[7];
				else menu();
			}
			else local += "out";
			if(!(argv[5][0] == '1' || argv[5][0] == '2')) errorType("Error! version invalid\n","","");;
			char *p = strstr(argv[2], ".nex");
			if(p != NULL) errorType("Error! The input file must be in Newick format\n","","");
			sscanf(argv[4], "%d", &sometimesGenerates);
			fstream file (argv[2], fstream::in | fstream::out);
			if(!file.good()) errorType("Error! The file ",argv[2]," does not exist\n");
			file >> line;
			file.close();
			parsingTree();
			version = argv[5][0];
			newSpeciesFile(argv[3], local.c_str());
		}
		else menu();
	}
	ftime(&fim);
	printTime();
	return 0;
}

void printTime(){ cout << "Time : "<< ((double) fim.time + ((double) fim.millitm * 0.001)) - ((double) ini.time + ((double) ini.millitm * 0.001)) << "\n"; }

void clearDisplay(){
	#ifdef _WIN32 
		system("cls");
	#endif
	#ifdef __linux__
		system("clear");
	#endif
	#ifdef DJGPP
		system("cls");
	#endif
}

void menu(){
	clearDisplay();
	cout << " To calculate the distance matrix of a(many) tree(s) type:\n";
	cout << "\t sunplin distMat nameTree.tree(nex) [-l nameFile]\n";
	cout << " Where : \n";
	cout << "\t nameTree.tree : local tree(s) containing the phylogeny(ies)\n";
	cout << "\t -l : if you want to specify the filename where the matrix will be saved\n";
	cout << "\t nameFile : name of the file containing the matrix\n\n";
	cout << " To insert missing species in tree type:\n";
	cout << "\t sunplin insertTree nameTree.tree nameFileSpecies numberGTree <version> [-l nameFile]\n";
	cout << " Where : \n";
	cout << "\t nameTree.tree : local tree containing the phylogeny\n";
	cout << "\t nameFileSpecies : local file containing the species to be inserted\n";
	cout << "\t version : insertion method: (1)over the nodes OR (2)over the branch lengths\n";
	cout << "\t -l : if you want to specify the filename where will the(s) new(s) tree(s) will be saved\n";
	cout << "\t nameFile : filename that will contain(s) the new(s) tree(s)\n\n";
	exit(1);
}

void generateNamesSpecies(){
	char alphabet[ ] = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};
	string nameEspecie = "";
	char numberAux[100];
	for(int i = 0; i < MAXGENERATENAMES; i++){
		sprintf(numberAux, "%d", i);
		for( int j = 0; j < 26; j++ ){
			nameEspecie += alphabet[j];
			nameRand.push_back(nameEspecie + numberAux);
			nameEspecie = "";
		}
	}
}

void convertNewickToNexus(int op, fstream &nexus, string nw){
	switch (op){
		case 0: nexus << "begin trees;\n"; break;
		case 1: nexus << "tree" << treeAmount++ << ":\n" << nw << "\n"; break;
		case 2: nexus << "end;"; nexus.close(); break;
	}
}

void errorType(string s1, string s2, string s3){
	if(s1 != "") printf("%s",s1.c_str());
	if(s2 != "") printf("%s",s2.c_str());
	if(s3 != "") printf("%s",s3.c_str());
	exit(1);
}

void parsingTree(){
	int n = 0, i = 0, nodei = -1, leni = -1, atn = 0, branchlengths = 0, done = 0, q;
	int lbrack, rbrack, comma, quant;
	char ch, tmp[15], iname[1000], taxa[1000];
	comma = quant = nodes = rbrack = lbrack = 0;
	quant = strlen(line);
	for(int i = 0; i < quant; i++){
		ch = line[i];
		if(ch == '(') lbrack++;
		if(ch == ')') rbrack++;
		if(ch == ',') comma++;
	}
	if(lbrack != rbrack) errorType("Unbalanced parenthesis.\n","","");
	nodes = lbrack + comma + 1;
	line[quant-1] = 59;
	line[quant] = 0;
	grafo.clear();
	grafoMatrix.clear();
	grafo.resize(nodes+200);
	grafoMatrix.resize(nodes+200);
	father.resize(nodes+200);
	while(line[n] != ';') n++;
	while (i < n) {
		done = 0;
		if (line[i] == '('){ // "(" 
			nodei++;
			if(leni!=-1){
				noat[leni]++;
				grafo[leni+1].push_back(nodei+1);
				grafoMatrix[leni+1].push_back(nodei+1);
				father[nodei+1] = leni+1;
			}
			up[nodei] = leni;
			strcat(taxon[nodei], ".");
			leni = nodei;
			done = 1;
			i++;
		}
		else if (line[i] == ','){ // ","
			done = 1;
			i++;
		}
		else if (line[i] == ')'){ // ")"
			leni = up[leni];
			atn = up[atn];
			done = 1;
			i++;
		}
		else if ((((line[i] >= 65) && (line[i] <= 90)) || ( (line[i] >= 97) && (line[i] <= 122)) || (line[i] == 45) || (line[i] == 95)) && (line[i-1] == ')')){
			string name = "";
			while ((line[i] != ':') && (line[i] != ',') && (line[i] != ')') && (line[i] != '[') && (line[i] != ';') && (line[i] != ']')) {
				name += line[i];
				i++;
			}
			strcpy(taxon[atn], name.c_str());
			done = 1;
		}
		else if (line[i] == '['){
			while (line[i] == ']') i++;	
			i++;
			done = 1;
		}
		else if (line[i] == ':'){
			string num = "";
			i++;
			while(((line[i] >= 48) && (line[i] <= 57)) || (line[i] == 46) || (line[i] == 'E') || (line[i] == '-')){
				num += line[i];
				i++;
			}
			if (num[0] != 0) branchlengths = 1;
			sscanf((char*)num.c_str(), "%lf", &bl[atn]);
			done = 1;
		}
		while ((line[i] == ' ') || (line[i] == '\t') || (line[i] == '\n') || (line[i] == '\r')) i++;
		if (done == 0){
			string tax = "";
			tax += line[i++];
			while ((line[i] != ',') && (line[i] != ')') && (line[i] != ':') && (line[i] != '[')){
				tax += line[i++];
			}
			strcpy(taxa,tax.c_str());
			nodei++;
			atn = nodei;
			up[nodei] = leni;
			strcpy(taxon[nodei], taxa );
			noat[leni]++;
			grafo[leni+1].push_back(nodei+1);
			grafoMatrix[leni+1].push_back(nodei+1);
			father[nodei+1] = leni+1;
		}
	}
	for (i = 0; i <= nodei; i++)
		if (branchlengths == 0) bl[i] = 1.0;
	
	for (i = 0; i < nodes; i++){ 
		if (strcmp(taxon[i], ".") == 0) strcpy(taxon[i], "");
	}
}

/*                        Responsible for the insertion of missing species                  */

string grafoToNewick(int u){
	int tam = grafo[u].size();
	if(!tam){ return ""; }
	string str = "(";
	char num[100];
	for(int i = 0; i < tam; i++){
		string ret = grafoToNewick(grafo[u][i]);
		str += ret;
		if( ret == "") str += taxon[grafo[u][i]-1];
		sprintf(num, "%0.6lf", bl[grafo[u][i]-1]);
		str += ":";
		str += num;
		if(i + 1 != tam) str += ",";
	}
	str += ")";
	return str;
}

int searchNameEspecie(char *str){
	for(int i = 0; i < nodes; i++) 
		if( !strcmp(str, taxon[i]) ) return (i + 1);
	errorType("Error! The specie ",str," does not belong to tree!");
}

void newSpeciesFile(char local[100], string nameNewFile){
	ifstream in( local , ifstream::in );
	if(!in.good()) errorType("Error opening file in place: ",local,"\n");  
	nameNewFile += ".nex";
	fstream filestr (nameNewFile.c_str(), fstream::out);
	convertNewickToNexus(0,filestr,"");
	vector < vector<int > > grafoBackup ( grafo.begin(), grafo.end());
	vector < int > fatherBackup( father.begin(), father.end() );
	for( int i = 0; i <= nodes; i++ ) blCopy[i] = bl[i]; 
	vector < infoNode > :: iterator it;
	int copyNodes = nodes;
	nodeInfo.resize(nodes + 10);
	char name[1000], insertionPoint[1000];
	if( version == '2' ) dfs_sum_V2(1);
	else dfs_insert_V1(1);
	while(in >> name >> insertionPoint){
		int ptInsercao = searchNameEspecie(insertionPoint);
		espN.push_back(insert(name,ptInsercao));
		espNames.push_back(insert(name,ptInsercao));
	}
	for(int i = 0;  i < sometimesGenerates; i++){
		amountNewNames = 0;
		for(int j = 0; j < nodes+10; j++) bl[j] = blCopy[j]; 
		if( version == '2' ) {
			for(int j = 0; j < espN.size(); j++){
				double calc = nodeInfo[(sum[espN[j].insertionPoint-1].sizeSubTree-1) + espN[j].insertionPoint].valor - nodeInfo[espN[j].insertionPoint].valor;
				calc *= (((double)(rand() % MAXV) + 1.0)/100.0);
				calc += nodeInfo[espN[j].insertionPoint].valor;
				it = upper_bound(nodeInfo.begin() + espN[j].insertionPoint, nodeInfo.begin() + espN[j].insertionPoint + (sum[espN[j].insertionPoint-1].sizeSubTree-1), calc, comp);
				insertSpecies(it->id, espN[j].name);
			}
		}
		else
			for( int j = 0; j < espNames.size(); j++) 
				insertSpecies((rand()%sumAux[espNames[j].insertionPoint]) + espNames[j].insertionPoint, espNames[j].name);

		convertNewickToNexus(1,filestr, grafoToNewick(1));
		grafo.clear();
		father.clear();
		grafo.assign(grafoBackup.begin(), grafoBackup.end());
		father.assign(fatherBackup.begin(), fatherBackup.end());
		nodes = copyNodes;
	}
	convertNewickToNexus(2,filestr, "");
}


void removeEdge(int fatherAux, int child){
	for(int i = 0; i < grafo[fatherAux].size(); i++){
		if(grafo[fatherAux][i] == child){
			grafo[fatherAux].erase(grafo[fatherAux].begin() + i);
			break;
		}
	}
}

void insertSpecies(int u, string name){
	int fatherAux = father[u];
	double R, value;
	if(! grafoMatrix[u].size()){
		removeEdge(fatherAux, u);
		grafo[fatherAux].push_back(nodes + 1);
		father[nodes+1] = fatherAux;
		grafo[nodes + 1 ].push_back(u);
		father[u] = nodes + 1;
		grafo[nodes + 1 ].push_back(nodes + 2 );
		father[nodes+2] = nodes+1;
		R = ((double)(rand() % MAXV) + 1.0)/100.0;
		value = bl[u-1];
		value *= R;
		bl[nodes] = value;
		bl[nodes + 1] = bl[u - 1] = fabs(value - bl[u - 1]);
		string novaEspecie = nameRand[amountNewNames++];
		strcpy(taxon[nodes],novaEspecie.c_str());
		strcpy(taxon[nodes+1],name.c_str());
	}
	else if( grafoMatrix[u].size() == 1){
		grafo[u].push_back(nodes + 1);
		father[nodes+1] = u;
		bl[nodes] = bl[grafo[u][0] - 1];
		strcpy(taxon[nodes],name.c_str());
	}
	else{
		int child = rand() % grafo[u].size();
		double division;
		value = sum[u - 1].sumSubTree / double(sum[u - 1].rootCount); 
		R = ((double)(rand() % MAXV) + 1.0)/100.0;
		division = bl[grafo[u][child] - 1];
		division *= R;
		bl[nodes] = division;
		bl[nodes + 1] = fabs(value - division);
		if( bl[nodes + 1]  < EPS ) bl[nodes + 1] = value;
		bl[grafo[u][child] - 1] = fabs(bl[grafo[u][child] - 1] - division); 
		string novaEspecie = nameRand[amountNewNames++];
		strcpy(taxon[nodes],novaEspecie.c_str());
		strcpy(taxon[nodes+1],name.c_str());
		int fi = grafo[u][child];
		grafo[u].erase(grafo[u].begin() + child);
		grafo[u].push_back(nodes + 1 );
		father[nodes + 1] = u;
		grafo[nodes + 1 ].push_back(nodes + 2 );
		father[nodes + 2] = nodes + 1;
		grafo[nodes + 1 ].push_back(fi);
		father[fi] = nodes + 1;
	}
	nodes += 2;
	if(nodes+2 >= grafo.size()) grafo.resize(grafo.size() + 100); 
	if(nodes+2 >= father.size()) father.resize(father.size() + 100); 
}


int dfs_insert_V1( int u ){
	int sAux = 1, poss = 0;
	double ans = 0.;
	sum[u-1].rootCount = (grafoMatrix[u].size() == 0);
	sum[u-1].sumSubTree = 0.;
	for(int i = 0; i < grafoMatrix[u].size(); i++){
		int at = grafoMatrix[u][i];
		sAux += dfs_insert_V1(at);
		sum[u-1].rootCount += sum[at-1].rootCount;
		sum[u-1].sumSubTree += ((double)sum[at-1].rootCount * bl[at-1]) + sum[at-1].sumSubTree;
	}
	sumAux[u] = sAux;
	return sAux;
}

double dfs_sum_V2( int u){
	double sumAux = bl[u-1];
	camNode++;
	custAnt += bl[u-1];
	if( u != 1 ){ nodeInfo[u].valor = custAnt; nodeInfo[u].id = u ; }
	sum[u-1].sizeSubTree = 1;
	sum[u-1].rootCount = (grafoMatrix[u].size() == 0);
	sum[u-1].sumSubTree = 0.;
	for( int i = 0; i < (int)grafoMatrix[u].size(); i++){
		int at = grafoMatrix[u][i];
		sumAux += bl[at-1];
		dfs_sum_V2(at);
		sum[u-1].rootCount += sum[at-1].rootCount;
		sum[u-1].sizeSubTree += sum[at-1].sizeSubTree;
		sum[u-1].sumSubTree += ((double)sum[at-1].rootCount * bl[at-1]) + sum[at-1].sumSubTree;
	}
	return (sumAux);
}

/*                                      ----------------------                                    */

/*                                     Distance matrix calculation                                */
void calcMatrixDistance(){
	int d = 0;
	for(int i = 0; i <= nodes; i++) nodedad[i] = -1;
	explore( 1, 0 );
	pos = 0, cntchain = 0;
	heavy_light( 1, 0, -1, 0 );
	for (int k = 1; k <= nodes; k++){
		if( !noat[k-1] ){
			for (int w = k+1; w <= nodes; w++) {
				if( !noat[w-1] ){
					d = lca(k,w);
					dist[k-1][w-1] = dist[w-1][k-1] = bl[k-1] + bl[w-1] - 2.0*bl[d-1];
				}
			}
		}
	}
}

void printToMatrixDistanceInFile(string nameNewFile){
	nameNewFile += ".txt";
	fstream out (nameNewFile.c_str(), fstream::out);
	for(int i = 0; i < nodes; i++){
		if(!noat[i])
			out << "\t" << taxon[i];
	}
	out << "\n";
	for(int i = 0; i < nodes; i++){
		if(!noat[i]){
			out << taxon[i];
			for(int j = 0; j < nodes; j++){
				if(!noat[j]) out << "\t" << fixed << setprecision(6) << dist[i][j];
			}
			out << "\n";
		}
	}
	out.close();
}

void explore( int x, int dad ) {
	 if( nodedad[x] != -1 ) return;
	 nodedad[x] = dad;
	 treesize[x] = 1;
	 for( int i = 0; i < ( int )grafo[x].size(); ++i ){
		 if( grafo[x][i] != dad ) {
			bl[grafo[x][i] - 1] += bl[ x - 1];
			 explore( grafo[x][i], x );
			 treesize[x] += treesize[ grafo[x][i] ];
		 }
	}
}

void heavy_light( int x, int dad, int k, int p ) {
	 if( p == 0 ){
		 k = cntchain++;
		 chainleader[k] = x;
	 }
	 chain[x] = k;
	 homepos[x] = pos++;
	 int mx = -1;
	 for( int i = 0; i < ( int )grafo[x].size(); ++i )
		 if(grafo[x][i] != dad && (mx == -1 || treesize[grafo[x][i]] > treesize[ grafo[x][mx]])) mx = i;
	 if( mx != -1 ) heavy_light( grafo[x][mx], x, k, p+1 );
	 for( int i = 0; i < ( int )grafo[x].size(); ++i )
		 if( grafo[x][i] != dad && i != mx ) heavy_light( grafo[x][i], x, -1, 0 );
}

int lca(int a, int b) {
	while (chain[a] != chain[b]) {
		if (treesize[chainleader[chain[a]]] >= treesize[chainleader[chain[b]]]) b = nodedad[chainleader[chain[b]]];
		else a = nodedad[chainleader[chain[a]]];
	}
	if (treesize[a] < treesize[b]) return b;
	return a;
}
