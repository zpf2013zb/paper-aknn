#include "gendata.h"
#include "utility.h"
#include "netshare.h"
#include "egtree.h"
#include "KeywordsGenerator.h"
//#pragma comment(lib,"ws2_32.lib")

#define MAX_KEYWORDSN 64
#define MAXLEVEL 10
#define MAXKWD 64

// to be initialized
char **cur_block;
int *cur_block_offset, *cur_node_maxkey, *cur_node_entries;
char *block;
int  root_block_id, num_written_blocks, top_level;
extern int BlkLen;
extern int i_capacity;

int PtMaxKey = 0;
bool PtMaxUsing = false;	// for solving the bug that was to find
float FACTOR;
bool IS_NODE_SIZE_GIVEN = false;

int num_D;
int num_K;


float STEP_SIZE;
float MAX_SEG_DIST;
float AVG_DEG;
extern vector <TreeNode> EGTree;


struct PartAddr {
	int part;
	int addr;
};
map<int, PartAddr> partID;

BTree* initialize(char *treename)
{
	BTree *bt;

	cur_block = new char*[MAXLEVEL]; //ptr to current node at each level
	cur_block_offset = new int[MAXLEVEL];
	cur_node_maxkey = new int[MAXLEVEL];
	cur_node_entries = new int[MAXLEVEL];
	top_level = 0;

	for (int i = 0; i<MAXLEVEL; i++) cur_block[i] = NULL;
	block = new char[BlkLen];

	// number of cached file entries: 128
	bt = new BTree(treename, BlkLen, 128);
	i_capacity = (BlkLen - sizeof(char) - sizeof(int)) / (sizeof(int) + sizeof(int));
	printf("i_capacity=%d\n", i_capacity);
	root_block_id = bt->root_ptr->block;
	num_written_blocks = 0;
	return  bt;
}

void addentry(BTree *bt, int *top_level, int capacity, int level, int key, int *block_id, int RawAddr)
{
	if (cur_block[level] == NULL)   //new node to be created
	{
		if ((*top_level) < level) //new root
			*top_level = level;
		cur_block[level] = new char[BlkLen];

		char l = (char)level;
		memcpy(cur_block[level], &l, sizeof(char));

		cur_block_offset[level] = sizeof(char) + sizeof(int);
		cur_node_entries[level] = 0;
	}
	cur_node_maxkey[level] = key;
	if ((level == 1) && PtMaxUsing) cur_node_maxkey[level] = PtMaxKey>key ? PtMaxKey : key;//Modified by Qin Xu

																						   //copy key as new current entry and also the pointer to lower node
	memcpy(cur_block[level] + cur_block_offset[level], &key, sizeof(int));
	cur_block_offset[level] += sizeof(int);

	//********* (Xblock_id for raw sequential page !)
	int Xblock_id = (level == 1) ? (RawAddr) : (*block_id);
	memcpy(cur_block[level] + cur_block_offset[level], &Xblock_id, sizeof(int));

	cur_block_offset[level] += sizeof(int);
	cur_node_entries[level]++;

	if (cur_node_entries[level] == capacity)   //node is full
	{
		//copy capacity information
		memcpy(cur_block[level] + sizeof(char), &capacity, sizeof(int));
		//add maxkey of this node to upper level node
		bt->file->append_block(cur_block[level]);
		(*block_id)++;
		bt->num_of_inodes++;
		addentry(bt, top_level, capacity, level + 1,
			cur_node_maxkey[level], block_id, 0);
		delete[] cur_block[level];
		cur_block[level] = NULL;
	}
}

void finalize(BTree* bt)
{
	//flush non-empty blocks
	for (int level = 1; level <= top_level; level++)
	{
		if (cur_block[level] != NULL)
		{
			//copy capacity information
			memcpy(cur_block[level] + sizeof(char), &cur_node_entries[level], sizeof(int));
			//add mbr of this node to upper level node
			if (level == top_level)
			{
				//root
				bt->file->write_block(cur_block[level], root_block_id);
				bt->num_of_inodes++;
				bt->root_ptr = NULL;
				bt->load_root();
				//printf("root written, id=%d\n", root_block_id);
			}
			else
			{
				bt->file->append_block(cur_block[level]);
				num_written_blocks++;
				bt->num_of_inodes++;
				addentry(bt, &top_level, i_capacity, level + 1,
					cur_node_maxkey[level], &num_written_blocks, 0);
			}
			delete[] cur_block[level];
			cur_block[level] = NULL;
		}
	}
	delete[] block;
}
//Add by Qin Xu for sort InerNode purpose
struct ComparInerNode
{
	bool operator () (const InerNode& left, const InerNode& right) const
	{
		return left.dis > right.dis;
	}
};

//------extend for egtree
void makeEPtFiles(FILE *ptFile, char* treefile) {
	// construct the extend point file
	set<__int64> vEdgeKey;
	int sizep = sizeof(int) + sizeof(float) + sizeof(unsigned long long);
	PtMaxUsing=true;							 
	BTree* bt=initialize(treefile);
	printf("making EPtFiles\n");
	//------treeNode start from 0
	int treeNodeID = 0;
	int nodeID = 0;
	int RawAddr = 0, key = 0, size, PtSize;
	int numedge = 0;
	for (; treeNodeID < EGTree.size(); treeNodeID++) {
		if (!EGTree[treeNodeID].isleaf) continue;
		// is leaf node ,sort it in ascend order
		sort(EGTree[treeNodeID].leafnodes.begin(), EGTree[treeNodeID].leafnodes.end(), less<int>());
		for (int i = 0; i<EGTree[treeNodeID].leafnodes.size(); i++) {
			nodeID = EGTree[treeNodeID].leafnodes[i];
			__int64 nid = nodeID;
			partID[nodeID].part = treeNodeID;
			// put the adjacent nodes of this node into adjList
			for (int j = 0; j<AdjList[nodeID].size(); j++) {
				int Nj = AdjList[nodeID][j];
				__int64 nje = Nj;
				// Nk can be smaller or greater than Ni !!!
				__int64 mapkey = getKey(nid, nje);
				if (vEdgeKey.find(mapkey) != vEdgeKey.end()) {// this edge has been visited
					continue;
				}
				vEdgeKey.insert(mapkey);
				numedge++;
				edge* e = EdgeMap[mapkey];			
				PtSize = e->pts.size();
				if (PtSize>0) {
					sort(e->pts.begin(), e->pts.end(), ComparInerNode());
					RawAddr = ftell(ptFile);	// set addr to correct amt.
					fwrite(&(e->Ni), 1, sizeof(int), ptFile);
					fwrite(&(e->Nj), 1, sizeof(int), ptFile);
					fwrite(&(e->dist), 1, sizeof(float), ptFile);
					fwrite(&(PtSize), 1, sizeof(int), ptFile);
					fwrite(&(e->pts[0]),e->pts.size(),sizeof(InerNode),ptFile); //???????????是不是要修改
					e->FirstRow = key;
					PtMaxKey = key + e->pts.size() - 1;	// useful for our special ordering !									//printf("(key,value)=(%d,%d)\n",key,RawAddr);
					
					addentry(bt, &top_level, i_capacity, 1, key, &num_written_blocks, RawAddr);
					key += sizeof(int) * 3 + sizeof(float);
					key += PtSize * sizeof(InerNode); //Modified by Qin Xu
				}
				else {
					e->FirstRow = -1; // also later used by AdjFile
				}

			}
		}
	}
	printf("Build-pt The number of edge is:%d\n", numedge);
	finalize(bt);
	bt->UserField=num_D;
	delete bt;
	PtMaxUsing=false;
}

void makeEAdjListFiles(FILE *alFile) { // construct the extend adjacentList file

	printf("making EAdjListFiles\n");
	//------treeNode start from 0
	int treeNodeID = 0;
	int nodeID = 0;
	int size, addr = 0;
	float distsum = 0;
	fwrite(&NodeNum, 1, sizeof(int), alFile);

	// slotted header info.
	int addrss = sizeof(int);
	int adjFixsize = 2 * sizeof(int) + sizeof(float) + sizeof(unsigned long long);

	for (; treeNodeID < EGTree.size(); treeNodeID++) {
		if (!EGTree[treeNodeID].isleaf) continue;
		// is leaf node ,may be no use
		sort(EGTree[treeNodeID].leafnodes.begin(), EGTree[treeNodeID].leafnodes.end(), less<int>());
		for (int i = 0; i<EGTree[treeNodeID].leafnodes.size(); i++) {
			nodeID = EGTree[treeNodeID].leafnodes[i];
			__int64 nid = nodeID;
			partID[nodeID].addr = addrss;
			size = AdjList[nodeID].size();
			fwrite(&(size), 1, sizeof(int), alFile);
			addrss = addrss + sizeof(int);
			for (int k = 0; k<size; k++)
			{
				int Nk = AdjList[nodeID][k];	// Nk can be smaller or greater than Ni !!!
				__int64 Nke = Nk;
				__int64 keyv = getKey(nid, Nke);
				edge* e = EdgeMap[keyv];				

				fwrite(&Nk, 1, sizeof(int), alFile);
				fwrite(&(e->dist), 1, sizeof(float), alFile);
				fwrite(&(e->sumkwds), 1, sizeof(unsigned long long), alFile);
				// pointer information
				fwrite(&(e->FirstRow), 1, sizeof(int), alFile); // use FirstRow for other purpose ...
				
				distsum += e->dist;
			}
			addrss = addrss + adjFixsize * AdjList[nodeID].size();
		}
	}
	distsum = distsum / 2;
	printf("total edge dist: %f\n", distsum);
	//printf("total keywords num:%d\n",num_K);
}

void BuildEBinaryStorage(string fileprefix) { // construct the extend binary storage
	BlkLen = getBlockLength();
	char tmpFileName[255];

	FILE *ptFile, *edgeFile;
	string name;
	name = fileprefix + "\\pfile";
	remove(name.c_str()); // remove existing file
	ptFile = fopen(name.c_str(), "wb+");
	name.clear();
	name = fileprefix + "\\pbtree";
	//sprintf(tmpFileName,"%s\\pbtree",fileprefix);
	//remove(tmpFileName); // remove existing file
	remove(name.c_str()); // remove existing file
	char* c;
	int len = name.length();
	c = new char[len + 1];
	strcpy(c, name.c_str());
	makeEPtFiles(ptFile, c);

	name.clear();
	name = fileprefix + "\\adjlist";

	remove(name.c_str()); // remove existing file
	edgeFile = fopen(name.c_str(), "wb+");
	makeEAdjListFiles(edgeFile);

	// generate reserve the part and address information
	partAddrSave(fileprefix);

	fclose(ptFile);
	fclose(edgeFile);
}
/*
void partAddrSave(string fileprefix) {

	printf("making partAddrFile\n");
	FILE *paFile;
	char tmpFileName[255];
	string name;
	name = fileprefix + "\\partAddr";
	//sprintf(tmpFileName, "%s\\partAddr", fileprefix);
	remove(name.c_str()); // remove existing file
	paFile = fopen(name.c_str(), "w+");

	fwrite(&NodeNum, 1, sizeof(int), paFile);
	map<int, PartAddr>::iterator my_Itr;
	int num = partID.size();
	int i = 0;
	for (my_Itr = partID.begin(); my_Itr != partID.end(); ++my_Itr) {
		int first, pt, addr;
		first = my_Itr->first;
		pt = my_Itr->second.part;
		addr = my_Itr->second.addr;
		
		fwrite(&first,sizeof(int),1, paFile);
		fwrite(&pt,sizeof(int),1, paFile);
		fwrite(&addr, sizeof(int),1,paFile);
		printf("i:%d, nodeid:%d, part:%d, addr:%d\n", i, my_Itr->first, my_Itr->second.part, my_Itr->second.addr);
		i++;
	}
	fclose(paFile);
}
*/

void partAddrSave(string fileprefix) {

	printf("making partAddrFile\n");
	//FILE *paFile;
	char tmpFileName[255];
	string name;
	name = fileprefix + "\\partAddr";
	//sprintf(tmpFileName, "%s\\partAddr", fileprefix);
	remove(name.c_str()); // remove existing file
	ofstream paFile(name.c_str());
	//ofstream owFile(objwetFile.c_str());
	//fwrite(&NodeNum, 1, sizeof(int), paFile);
	map<int, PartAddr>::iterator my_Itr;
	int num = partID.size();
	int i = 0;
	for (my_Itr = partID.begin(); my_Itr != partID.end(); ++my_Itr) {
		int first, pt, addr;
		first = my_Itr->first;
		pt = my_Itr->second.part;
		addr = my_Itr->second.addr;

		paFile << first << " " << pt << " " << addr << endl;

		//fwrite(&first, sizeof(int), 1, paFile);
		//fwrite(&pt, sizeof(int), 1, paFile);
		//fwrite(&addr, sizeof(int), 1, paFile);
		printf("i:%d, nodeid:%d, part:%d, addr:%d\n", i, my_Itr->first, my_Itr->second.part, my_Itr->second.addr);
		i++;
	}
	//fclose(paFile);
	paFile.close();
}

FastArray<float> xcrd, ycrd;

// future work: normalization of edge weight, reassign edge weight, divide into bands ...
// how about scanning whole file and putting edges with both nodeIDs in range ...
void ReadRealNetwork(std::string prefix_name, int _NodeNum)
{
	int eid, Ni, Nj;
	//int eid;
	float dist, x, y;
	//unsigned int Ni, Nj;
	string roadf;
	roadf = prefix_name + "\\road";	

	FILE* roadp = fopen(roadf.c_str(), "r");
	CheckFile(roadp, roadf.c_str());

	NodeNum = 0; // getKey depends on NodeNum so we have to init. it first
	// vertex num start from 0
	while (!feof(roadp)) {
		fscanf(roadp, "%d %d %d %f\n", &eid, &Ni, &Nj, &dist);
		NodeNum = std::max(std::max(Ni, Nj), NodeNum);
	}
	NodeNum++;
	printf("%d nodes read, ", NodeNum);
	PrintElapsed();
	set<int> edgeid;
	fseek(roadp, 0, SEEK_SET);
	AdjList = new FastArray<int>[NodeNum];
	EdgeNum = 0;
	EdgeMap.clear();
	int loop = 0;
	vector<__int64> vist;
	while (!feof(roadp))
	{	
		fscanf(roadp, "%d %d %d %f\n", &eid, &Ni, &Nj, &dist);
		if (Ni < NodeNum && Nj < NodeNum)  	// ignore edges outside the range
		{
			//printf( "%d %d %d %f\n", id, Ni, Nj, dist);
			__int64  a = Ni;
			__int64  b = Nj;
			__int64 key = getKey(a, b);
			if (find(vist.begin(), vist.end(), key) != vist.end()) {
				continue;
			}
			vist.push_back(key);
			edge* e = new edge;
			e->dist = dist;
			e->Ni = Ni<Nj ? Ni : Nj;
			e->Nj = Ni<Nj ? Nj : Ni;	// enforce the constriant Ni<Nj
			AdjList[Ni].push_back(Nj);
			AdjList[Nj].push_back(Ni);
			
			//if(EdgeMap[key]->FirstRow)
			int vists = vist.size();
			//printf(" %d\n", vists);
			EdgeMap[key] = e;	// should be ok
			EdgeNum++;
		}
	}
	printf("%d edges read, ", EdgeNum);
	PrintElapsed();
	fclose(roadp);
}

int GEN_PAIR_CNT = 0;

//Modified by Qin Xu
//for simplicy purpose

void genPairByAd(int& Ni, int& Nj)
{
	int Ri, Rj;
	bool useNormalGenerator = true;	// use a ... generator
	if (useNormalGenerator)
	{
		do
		{
			Ri = rand() % NodeNum;
		} while (AdjList[Ri].size() == 0);
		Rj = AdjList[Ri][rand() % AdjList[Ri].size()];
	}
	Ni = Ri<Rj ? Ri : Rj;
	Nj = Ri<Rj ? Rj : Ri;
}

void printBinary(unsigned long long n)
{
	for (int i = MAX_KEYWORDS - 1; i >= 0; i--)
		cout << ((n >> i) & 1);
	cout<<endl;
}

void GenOutliers(int NumPoint, int avgKeywords)
{
	std::vector<unsigned long long> keys = KeywordsGenerator::Instance().getKeywords(NumPoint, avgKeywords);
	int Ni, Nj;
	num_K = 0;
	num_D = NumPoint;

	for (int z = 0; z<NumPoint; z++)
	{
		genPairByAd(Ni, Nj);
		edge* e = EdgeMap[getKey(Ni, Nj)];
		// 确保POI不落在vertex上
		int randT = 0;
		while (randT == 0 || randT == RAND_MAX) {
			randT = rand();
		}
		float vP = (rand()*e->dist)/RAND_MAX;
		num_K += std::bitset<64>(keys[z]).count();
		unsigned long long tmpvct = keys[z];
		InerNode tempNode;
		tempNode.pid = z;
		tempNode.dis = vP;
		tempNode.vct = tmpvct;
		e->pts.push_back(tempNode);
		// -add sumkwds-
		if (e->pts.size() > 0) {
			e->sumkwds = e->sumkwds | tmpvct;
		}
		else {
			e->sumkwds = tmpvct;
		}		

#ifdef _TEST
		cout << "(" << Ni << "," << Nj << ")" << " ";
		cout << bitset<64>(tmpvct).to_string() << endl;
#endif

	}
	// spread outliers based on distance
	// work for connected graph, may not work for disconnected graphs
	// organize dist based on length and then ...
}

// For test purpose to check if graph is connected
struct StepEvent
{
	//float dist;
	double dist;
	int node;	// posDist for gendata.cc only
};

// ">" means the acsending order
struct StepComparison
{
	bool operator () (const StepEvent& left, const StepEvent& right) const
	{
		return left.dist > right.dist;
	}
};

typedef	priority_queue<StepEvent, vector<StepEvent>, StepComparison> StepQueue;

void ConnectedGraphCheck()
{
	BitStore isVisited;
	isVisited.assign(NodeNum, false);

	StepQueue sQ;
	StepEvent stepX;
	stepX.node = 0;
	sQ.push(stepX);

	int cnt = 0;
	while (!sQ.empty())  	// beware of infty loops
	{
		StepEvent event = sQ.top();
		sQ.pop();
		if (isVisited[event.node]) continue;
		isVisited[event.node] = true;
		cnt++;
		int Ni = event.node;
		for (int k = 0; k<AdjList[Ni].size(); k++)
		{
			int Nk = AdjList[Ni][k];	// Nk can be smaller or greater than Ni !!!
			if (!isVisited[Nk])
			{
				StepEvent newevent = event;
				newevent.node = Nk;
				sQ.push(newevent);
			}
		}
	}

	if (cnt == NodeNum)
		cout << "Road Network is Connected." << endl;
	else
		cout << "Road Network is not Connected." << endl;
}

void getOutliersFromFile(char* prefix_name)
{
	int Ni, Nj;
	int NumOutliers = 0;
	char tmpf[255];
	sprintf(tmpf, "%smap.txt", prefix_name);
	float dist;
	int num;
	unsigned long long keyword;
	float tmpdist;

	FILE* calf = fopen(tmpf, "r");
	CheckFile(calf, tmpf);

	while (!feof(calf))
	{
		fscanf(calf, "%d %d %f %d\n", &Ni, &Nj, &dist, &num);
		for (int i = 0; i<num; ++i)
		{
			fscanf(calf, "%llu %f", &keyword, &tmpdist);
			edge* e = EdgeMap[getKey(Ni, Nj)];
			InerNode tmpNode;
			tmpNode.dis = tmpdist;
			tmpNode.vct = keyword>1 ? keyword : 1;//keywords must less than MAX and more than zero
			e->pts.push_back(tmpNode);
			NumOutliers++;
		}
	}
	num_D = NumOutliers;

	//cout << "Total Num of Outliers Read From File:" << cnt << endl;
	// spread outliers based on distance
	// work for connected graph, may not work for disconnected graphs
	// organize dist based on length and then ...
}

//void reAssignKeywords(int avg)
//{
//    num_K=0;
//    auto iterEdge=EdgeMap.begin();
//    while(iterEdge!=EdgeMap.end())
//    {
//        edge* e=(*iterEdge).second;
//        auto iterIner=e->pts.begin();
//        while(iterIner!=e->pts.end())
//        {
//
//            int keywordnum=random()%(2*avg-1)+1;
//            num_K+=keywordnum;
//            unsigned long long tmpvct = KeywordsGenerator::Instance().getKeywords(keywordnum);
//            (*iterIner).vct=tmpvct;
//            iterIner++;
//        }
//        iterEdge++;
//    }
//}

int mainGenData(string prxfilename, roadParameter rp)
{
	InitClock();	// side effect: init. seeds for random generators
	//string road = prxfilename + "\\road";
	ReadRealNetwork(prxfilename, 0);
	//ConnectedGraphCheck();

	
	GenOutliers(EdgeNum*rp.avgNPt, rp.avgNKwd);
	printf("Avg keyword # per object:%f\n",float(num_K)/num_D);

	//--------------------extend for egtree-------------------
	mainEgtree(NodeNum, EdgeMap);
	BuildEBinaryStorage(prxfilename);
	egtree_save(prxfilename+"\\egtree");
	

	ofstream fout(prxfilename + "\\information");
	fout << "Number of Vertexes:" << NodeNum << endl;
	fout << "Number of Edges:" << EdgeNum << endl;
	fout << "Total Keyowords Number:" << num_K << endl;
	fout << "Total POIs Number:" << num_D << endl;
	fout << "Avg Keywords numbers per POI:" << float(num_K) / num_D << endl;
	fout.close();
	PrintElapsed();
	//*/
	return 0;
}

