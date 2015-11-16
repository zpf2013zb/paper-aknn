#include "btree.h"
#include "utility.h"
#include "netshare.h"
#include "ConfigType.h"
#include "egtree.h"
#include "KeywordsGenerator.h"
#include <random>
#include<algorithm>
#include <bitset>
#include <vector>
#include<stdio.h>
#define MAX_KEYWORDSN 64

int num_D;
int num_K;

int PoiNum;
int KwdNum;
//int distsum;
float STEP_SIZE;
float MAX_SEG_DIST;
float AVG_DEG;

#define MAXLEVEL 10
#define MAXKWD 64

//------------
struct PartAddr {
	int part;
	int addr;
};
//vector <TreeNode> EGTree;

// to be initialized
char **cur_block;
int *cur_block_offset, *cur_node_maxkey, *cur_node_entries;
char *block;
int i_capacity, root_block_id, num_written_blocks, top_level;
int BlkLen;

int PtMaxKey = 0;
bool PtMaxUsing = false;	// for solving the bug that was to find
float FACTOR;
bool IS_NODE_SIZE_GIVEN = false;

map<int, PartAddr> partID;
//????????
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

void addentry(BTree *bt, int *top_level, int capacity, int level, int key, int *block_id, int RawAddr = 0)
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
			cur_node_maxkey[level], block_id);
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
					cur_node_maxkey[level], &num_written_blocks);
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

// Pt FlatFile Field:
//		Header:	Ni(int), Nj(int), dist(float), size(int)
//		Entry:	vP(float)
void makePtFiles(FILE *ptFile, char* treefile)
{
	/*
	PtMaxUsing = true;
	BTree* bt = initialize(treefile);
	printf("making PtFiles\n");

	int RawAddr = 0, key = 0, size;	// changed
	EdgeMapType::iterator iter = EdgeMap.begin();
	while (iter != EdgeMap.end())
	{
		edge* e = iter->second;
		if (e->pts.size()>0)  	// do not index empty groups
		{
			sort(e->pts.begin(), e->pts.end(), ComparInerNode());

			RawAddr = ftell(ptFile);	// set addr to correct amt.
			size = e->pts.size();
			fwrite(&(e->Ni), 1, sizeof(int), ptFile);
			fwrite(&(e->Nj), 1, sizeof(int), ptFile);
			fwrite(&(e->dist), 1, sizeof(float), ptFile);
			fwrite(&(size), 1, sizeof(int), ptFile);
			// ------------M--modify this storage model
			int bytes = 0;
			for (int i = 0; i < size; i++) {
				InerNode in = e->pts[i];
				fwrite(&(in.poid), 1, sizeof(int), ptFile);
				fwrite(&(in.dis), 1, sizeof(float), ptFile);
				fwrite(in.attr, ATTRIBUTE_DIMENSION, sizeof(float), ptFile);
				fwrite(&(in.nOfK), 1, sizeof(int), ptFile);
				fwrite(&in.kwd.begin(), in.nOfK, sizeof(int), ptFile);
				bytes = bytes + (2 + in.nOfK)*sizeof(int) + (ATTRIBUTE_DIMENSION + 1)*sizeof(float);
			}

			//fwrite(&(e->pts[0]),e->pts.size(),sizeof(InerNode),ptFile); //???????????�ǲ���Ҫ�޸�
			e->FirstRow = key;
			PtMaxKey = key + e->pts.size() - 1;	// useful for our special ordering !

												//printf("(key,value)=(%d,%d)\n",key,RawAddr);
			addentry(bt, &top_level, i_capacity, 1, key, &num_written_blocks, RawAddr);
			// -------------------M--modify key value
			//key+=sizeof(int)*3+sizeof(float);
			//key+=e->pts.size()*sizeof(InerNode); //Modified by Qin Xu
			key += bytes;
		}
		else {
			e->FirstRow = -1;		// also later used by AdjFile
		}
		iter++;
	}
	finalize(bt);
	//bt->UserField=num_D;
	bt->UserField = PoiNum;
	delete bt;
	PtMaxUsing = false;
	*/
}

// Adj FlatFile Field:
//		Header:	size(int)
//		Entry:	Nk(int), eDist(float), PtGrpKey(int), PtSize(int)		changed
void makeAdjListFiles(FILE *alFile)
{
	/*
	printf("making alFiles, dependency on makePtFiles\n");
	int key = 0, size, PtSize;
	fwrite(&NodeNum, 1, sizeof(int), alFile);

	// slotted header info.
	int addr = sizeof(int) + sizeof(int)*NodeNum;
	for (int Ni = 1; Ni <= NodeNum; Ni++)
	{
		fwrite(&addr, 1, sizeof(int), alFile);
		addr += sizeof(int) + AdjList[Ni].size()*(2 * sizeof(int) + sizeof(float));
	}

	float distsum = 0;
	for (int Ni = 1; Ni <= NodeNum; Ni++)
	{
		size = AdjList[Ni].size();
		fwrite(&(size), 1, sizeof(int), alFile);

		for (int k = 0; k<AdjList[Ni].size(); k++)
		{
			int Nk = AdjList[Ni][k];	// Nk can be smaller or greater than Ni !!!
			edge* e = EdgeMap[getKey(Ni, Nk)];
			PtSize = e->pts.size();
			fwrite(&Nk, 1, sizeof(int), alFile);
			fwrite(&(e->dist), 1, sizeof(float), alFile);
			fwrite(&(e->FirstRow), 1, sizeof(int), alFile); // use FirstRow for other purpose ...
															//printf("(Ni,Nj,dataAddr)=(%d,%d,%d)\n",Ni,Nk,e->FirstRow);

			distsum += e->dist;
		}
		key = Ni;
	}
	distsum = distsum / 2;
	printf("total edge dist: %f\n", distsum);
	//printf("total keywords num:%d\n",num_K);
	*/
}

void BuildBinaryStorage(const char* fileprefix)
{
	BlkLen = getBlockLength();
	char tmpFileName[255];

	FILE *ptFile, *edgeFile;
	sprintf(tmpFileName, "%s\\pf.p_d", fileprefix);
	remove(tmpFileName); // remove existing file
	ptFile = fopen(tmpFileName, "w+");
	sprintf(tmpFileName, "%s\\pf.p_bt", fileprefix);
	remove(tmpFileName); // remove existing file
	makePtFiles(ptFile, tmpFileName);

	sprintf(tmpFileName, "%s\\adj.al_d", fileprefix);
	remove(tmpFileName); // remove existing file
	edgeFile = fopen(tmpFileName, "w+");
	makeAdjListFiles(edgeFile);

	fclose(ptFile);
	fclose(edgeFile);
}

//------extend for egtree
void makeEPtFiles(FILE *ptFile, char* treefile) {
	// construct the extend point file
	set<int> vEdgeKey;
	int sizep = sizeof(int) + sizeof(float) + sizeof(unsigned long long);
	PtMaxUsing=true;							 
	BTree* bt=initialize(treefile);
	printf("making EPtFiles\n");
	//------treeNode start from 0
	int treeNodeID = 0;
	int nodeID = 0;
	int RawAddr = 0, key = 0, size, PtSize;
	//vector<TreeNode>::iterator it=EGTree.begin();
	for (; treeNodeID < EGTree.size(); treeNodeID++) {
		if (!EGTree[treeNodeID].isleaf) continue;
		// is leaf node ,sort it in ascend order
		sort(EGTree[treeNodeID].leafnodes.begin(), EGTree[treeNodeID].leafnodes.end());
		for (int i = 0; i<EGTree[treeNodeID].leafnodes.size(); i++) {
			nodeID = EGTree[treeNodeID].leafnodes[i];
			partID[nodeID].part = treeNodeID;
			// put the adjacent nodes of this node into adjList
			for (int j = 0; j<AdjList[nodeID].size(); j++) {
				int Nj = AdjList[nodeID][j];
				// Nk can be smaller or greater than Ni !!!
				int mapkey = getKey(nodeID, Nj);
				edge* e = EdgeMap[mapkey];
				if (vEdgeKey.count(mapkey) != 0) {// this edge has not been visited
					continue;
				}
				vEdgeKey.insert(mapkey);
				PtSize = e->pts.size();
				if (PtSize>0) {
					sort(e->pts.begin(), e->pts.end(), ComparInerNode());
					RawAddr = ftell(ptFile);	// set addr to correct amt.
					fwrite(&(e->Ni), 1, sizeof(int), ptFile);
					fwrite(&(e->Nj), 1, sizeof(int), ptFile);
					fwrite(&(e->dist), 1, sizeof(float), ptFile);
					fwrite(&(PtSize), 1, sizeof(int), ptFile);
					fwrite(&(e->pts[0]),e->pts.size(),sizeof(InerNode),ptFile); //???????????�ǲ���Ҫ�޸�
					e->FirstRow = key;
					PtMaxKey = key + e->pts.size() - 1;	// useful for our special ordering !									//printf("(key,value)=(%d,%d)\n",key,RawAddr);
					
					addentry(bt, &top_level, i_capacity, 1, key, &num_written_blocks, RawAddr);
					key+=sizeof(int)*3+sizeof(float);
					key+=PtSize*sizeof(InerNode); //Modified by Qin Xu
					//key += bytes;
				}
				else {
					e->FirstRow = -1;		// also later used by AdjFile
				}

			}
		}
	}
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
	int key = 0, size, PtSize, addr = 0;
	//int *buf = new int[ NodeNum ];
	float distsum = 0;
	fwrite(&NodeNum, 1, sizeof(int), alFile);

	// slotted header info.
	int addr = sizeof(int);
	//vector<TreeNode>::iterator it=EGTree.begin();
	for (; treeNodeID < EGTree.size(); treeNodeID++) {
		if (!EGTree[treeNodeID].isleaf) continue;
		// is leaf node ,may be no use
		sort(EGTree[treeNodeID].leafnodes.begin(), EGTree[treeNodeID].leafnodes.end());
		for (int i = 0; i<EGTree[treeNodeID].leafnodes.size(); i++) {
			nodeID = EGTree[treeNodeID].leafnodes[i];
			partID[nodeID].addr = addr;
			//partID[nodeID] = treeNodeID;
			// put the adjacent nodes of this node into adjList
			//size=AdjList[Ni].size();
			//fwrite(&(size),1,sizeof(int),alFilDe);
			size = AdjList[nodeID].size();
			fwrite(&(size), 1, sizeof(int), alFile);
			//addr = addr + sizeof(int);
			for (int k = 0; k<AdjList[nodeID].size(); k++)
			{
				int Nk = AdjList[nodeID][k];	// Nk can be smaller or greater than Ni !!!
				edge* e = EdgeMap[getKey(nodeID, Nk)];
				PtSize = e->pts.size();

				fwrite(&Nk, 1, sizeof(int), alFile);
				fwrite(&(e->dist), 1, sizeof(float), alFile);
				fwrite(&(e->sumkwds), 1, sizeof(unsigned long long), alFile);
				// pointer information
				fwrite(&(e->FirstRow), 1, sizeof(int), alFile); // use FirstRow for other purpose ...
																//keyword information
				addr = addr + 3*sizeof(int) + sizeof(float) + sizeof(unsigned long long);
				distsum += e->dist;
			}
			key = treeNodeID;
		}
	}
	distsum = distsum / 2;
	printf("total edge dist: %f\n", distsum);
	//printf("total keywords num:%d\n",num_K);
}

void BuildEBinaryStorage(const char* fileprefix) { // construct the extend binary storage
	BlkLen = getBlockLength();
	char tmpFileName[255];

	FILE *ptFile, *edgeFile;
	sprintf(tmpFileName, "%s\\pf.ep_d", fileprefix);
	remove(tmpFileName); // remove existing file
	ptFile = fopen(tmpFileName, "w+");
	sprintf(tmpFileName,"%s.p_bt",fileprefix);
	remove(tmpFileName); // remove existing file
	makeEPtFiles(ptFile, tmpFileName);

	sprintf(tmpFileName, "%s\\adj.eal_d", fileprefix);
	remove(tmpFileName); // remove existing file
	edgeFile = fopen(tmpFileName, "w+");
	makeEAdjListFiles(edgeFile);
	// reserve the part and address information
	partAddrSave(fileprefix);

	fclose(ptFile);
	fclose(edgeFile);
}

void partAddrSave(const char* fileprefix) {

	printf("making partAddrFile\n");
	FILE *paFile;
	char tmpFileName[255];

	sprintf(tmpFileName, "%s\\part.inf", fileprefix);
	remove(tmpFileName); // remove existing file
	paFile = fopen(tmpFileName, "w+");

	fwrite(&NodeNum, 1, sizeof(int), paFile);
	map<int, PartAddr>::iterator my_Itr;
	for (my_Itr = partID.begin(); my_Itr != partID.end(); ++my_Itr) {
		fwrite(&my_Itr->first, 1, sizeof(int), paFile);
		fwrite(&my_Itr->second.part, 1, sizeof(int), paFile);
		fwrite(&my_Itr->second.addr, 1, sizeof(int), paFile);
	}
	fclose(paFile);
}

FastArray<float> xcrd, ycrd;

// future work: normalization of edge weight, reassign edge weight, divide into bands ...
// how about scanning whole file and putting edges with both nodeIDs in range ...
void ReadRealNetwork(std::string prefix_name, int _NodeNum = 0)
{
	int id, Ni, Nj;
	float dist, x, y;

	// side effect: change nodes
	char edgef[255], nodef[255], poif[255];
	///sprintf(edgef, "%s\\road.txt", prefix_name.c_str());
	////sprintf(poif, "%s\\poi.txt", prefix_name.c_str());
	sprintf(edgef, "%s.cedge", prefix_name.c_str());

	FILE* cedge = fopen(edgef, "r");
	CheckFile(cedge, edgef);
	NodeNum = 0;	// getKey depends on NodeNum so we have to init. it first
	// vertex num start from 0
	while (!feof(cedge)) {
		fscanf(cedge, "%d %d %d %f\n", &id, &Ni, &Nj, &dist);
		NodeNum = std::max(std::max(Ni, Nj), NodeNum);
	}
	NodeNum++;
	printf("%d nodes read, ", NodeNum);
	PrintElapsed();

	fseek(cedge, 0, SEEK_SET);
	AdjList = new FastArray<int>[NodeNum];
	EdgeNum = 0;
	while (!feof(cedge))
	{
		fscanf(cedge, "%d %d %d %f\n", &id, &Ni, &Nj, &dist);
		if (Ni <= NodeNum&&Nj <= NodeNum)  	// ignore edges outside the range
		{
			//printf( "%d %d %d %f\n", id, Ni, Nj, dist);
			edge* e = new edge;
			e->dist = dist;
			//Modified by Qin Xu
			e->Ni = Ni<Nj ? Ni : Nj;
			e->Nj = Ni<Nj ? Nj : Ni;	// enforce the constriant Ni<Nj
			AdjList[Ni].push_back(Nj);
			AdjList[Nj].push_back(Ni);
			EdgeMap[getKey(Ni, Nj)] = e;	// should be ok
			EdgeNum++;
		}
	}
	printf("%d edges read, ", EdgeNum);
	PrintElapsed();
	fclose(cedge);
	// --- MM--remove read poi------
	/*
	//-------------read poi--------------
	FILE* cpoi = fopen(poif, "r");
	CheckFile(cpoi, poif);
	PoiNum = 0;	// to record the poi num
	KwdNum = 0;   // to record the kwd num
	int pid, nid, njd, nOfKwd;
	float distMinV;
	float attribute[ATTRIBUTE_DIMENSION];
	memset(attribute, 0.0, sizeof(float)*ATTRIBUTE_DIMENSION);
	int tempMaxNodeNum = 0; // used to verify whether vertex id of poi larger than road
	while (!feof(cpoi)) {
		// read <object id> <vertex id> <vertex id> <distance to edge vertex,little>
		fscanf(cpoi, "%d %d %d %f", &pid, &nid, &njd, &distMinV);
		// read <attributes information,randomly>
		for (int i = 0; i < ATTRIBUTE_DIMENSION; i++) {
			fscanf(cpoi, " %f", &attribute[i]);
		}
		// read <NofKwd>  <textual description/kwd>
		fscanf(cpoi, " %d", &nOfKwd);

		PoiNum = std::max(PoiNum, pid);
		tempMaxNodeNum = std::max(std::max(nid, njd), tempMaxNodeNum);
		if (tempMaxNodeNum > NodeNum) {
			printf("The vertex of poi out of the scope of road!");
		}

		InerNode tempNode;
		tempNode.poid = pid;
		tempNode.dis = distMinV;
		memcpy(tempNode.attr, attribute, sizeof(float)*ATTRIBUTE_DIMENSION);
		tempNode.nOfK = nOfKwd;
		int keywd;
		int loop = 1;
		KwdNum += nOfKwd;
		while (loop<nOfKwd) {
			fscanf(cpoi, " %d", &keywd);
			tempNode.kwd.insert(keywd);
			loop++;
		}
		fscanf(cpoi, " %d\n", &keywd);
		tempNode.kwd.insert(keywd);
		// push poi to edge
		edge* e = EdgeMap[getKey(nid, njd)];
		e->pts.push_back(tempNode);
	}
	fclose(cpoi);

	//--------------------compute the lower&upper bound and kwds
	EdgeMapType::iterator iter = EdgeMap.begin();
	for (; iter != EdgeMap.end(); iter++) {
		edge *e = iter->second;
		for (int i = 0; i<e->pts.size(); i++) {
			InerNode temp = e->pts[i];
			if (i == 0) { //ֱ�ӳ�ʼ��
				for (int j = 0; j<ATTRIBUTE_DIMENSION; j++) {
					e->attrBound[j][0] = temp.attr[j];
					e->attrBound[j][1] = temp.attr[j];
				}
				copy(temp.kwd.begin(), temp.kwd.end(), e->kwds.begin());
			}
			else {
				for (int j = 0; j<ATTRIBUTE_DIMENSION; j++) {
					if (e->attrBound[j][0] > temp.attr[j]) {
						e->attrBound[j][0] = temp.attr[j];
					}
					if (e->attrBound[j][1] < temp.attr[j]) {
						e->attrBound[j][1] = temp.attr[j];
					}
				}
				set_union(temp.kwd.begin(), temp.kwd.end(), e->kwds.begin(), e->kwds.end(), e->kwds.begin());
			}
		}
	} // end for iter

	*/
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
	//cout<<endl;
}

void GenOutliers(int NumPoint, int avgKeywords)
{
	std::vector<unsigned long long> keys = KeywordsGenerator::Instance().getKeywords(NumPoint, avgKeywords);
	int Ni, Nj;
	num_K = 0;
	PoiNum = NumPoint;
	for (int z = 0; z<NumPoint; z++)
	{
		genPairByAd(Ni, Nj);
		edge* e = EdgeMap[getKey(Ni, Nj)];
		//float vP = drand48()*(e->dist);
		float vP = (rand()*1.0)/RAND_MAX *(e->dist);
		num_K += std::bitset<64>(keys[z]).count();
		unsigned long long tmpvct = keys[z];
		InerNode tempNode;
		tempNode.pid = z + 1;
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
	num_D = NumPoint;
	// spread outliers based on distance
	// work for connected graph, may not work for disconnected graphs
	// organize dist based on length and then ...
}



// For test purpose to check if graph is connected
struct StepEvent
{
	float dist;
	//double dist;
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
	//isVisited.assign(NodeNum + 1, false); // more one than nodeNum
	isVisited.assign(NodeNum, false);
	StepQueue sQ;
	StepEvent stepX;
	stepX.node = 1;
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



int main(int argc, char *argv[])
{
	string configFileName = "config.prop";

	ConfigType cr(configFileName, argc, argv);

	cr.ListConfig();

	InitClock();	// side effect: init. seeds for random generators

	ReadRealNetwork(cr.getMapFileName().c_str(), 0);

	ConnectedGraphCheck();

	GenOutliers(EdgeNum*cr.getParameterOutlierDensity(), cr.getParameterAvgKeywordsNumberOfOutliers());
	printf("Avg keyword # per object:%f\n",float(num_K)/num_D);

	//BuildBinaryStorage(cr.getIndexFileName().c_str());
	//--------------------extend for egtree-------------------
	mainFunction(NodeNum, EdgeMap);
	BuildEBinaryStorage(cr.getIndexFileName().c_str());
	//-------------M-- save in this part then destory the temp memory---------
	char tmpFileName[255];
	sprintf(tmpFileName, "%s\\egindex.eg_inx", cr.getIndexFileName().c_str());
	egtree_save(tmpFileName);
	finalize();

	ofstream fout(cr.getIndexFileName() + "_Information");
	fout << "Number of Vertexes:" << NodeNum << endl;
	fout << "Number of Edges:" << EdgeNum << endl;
	fout << "Total Keyowords Number:" << KwdNum << endl;
	fout << "Total POIs Number:" << PoiNum << endl;
	fout << "Avg Keywords numbers per POI:" << float(KwdNum) / PoiNum << endl;
	fout.close();

	PrintElapsed();

	return 0;
}

