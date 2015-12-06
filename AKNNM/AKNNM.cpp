#include <iostream>
#include <unordered_map>
#include <queue>
#include <bitset>
#include <algorithm>
#include "diskbased.h"
#include "egtree.h"
#include "netshare.h"
#include "gendata.h"
#include <fstream>
#include <sstream>
#include <math.h>
#pragma comment(lib,"ws2_32.lib")

using namespace std;

#define EGTD 1
#define EGETD 2
#define TA 3
#define CE 4

#define INFINITE_MAX 100000000.0
#define MAX_KEYWORDN 64 
// define the query path
string basepath = "..\\dataset";

int nOfPageAccessed = 0;
int nOfEdgeExpended = 0;
int nOfEdgeExpendedQ = 0;
int algorithmId;


int BlkLen;
int i_capacity;

int NodeNum;
int EdgeNum;
FastArray<int>* AdjList;
FastArray<point> PtList;
EdgeMapType EdgeMap;	// key: i0*NodeNum+j0

FILE *PtFile, *AdjFile;
BTree *PtTree;

int poicnt = 0;
Query Q;

// for TKDE2005
#define min(a, b) (((a) < (b))? (a) : (b)  )
#define max(a, b) (((a) > (b))? (a) : (b)  )
#define ValAbs(x) (((x) >  0 )? (x) : -(x) )

//typedef map<int, float> DISTMAP;
typedef map<int, float> DISTMAP;
DISTMAP elem_map;
DISTMAP* DistMaps;	// only use cheap storage for range search
struct DistElem {
	int id;
	float dist;
};
typedef vector<DistElem> DistQueue;




BitStore* astar_Visits;		// tmp. bitmap for A* search

BitStore *raVisited, *saVisited; // top-k
int PrintLimit;

typedef vector<int> 	IntVec;
typedef vector<float> 	FloatVec;

FloatVec M_xcrd, M_ycrd, compdist, bestcompdist;

BitStore onSameEdge;
float** partdists;
float bestdist;	// best value found so far

bool isNNquery = true;
bool isEuclidSpace = true;
bool isWeightUse = false;
char ptgrprtfn[255];

DistQueue dQ_graph;
//extern BTree *PtTree;
//extern struct QueryPoint;
//extern QueryPoint* Q.querypts;
//extern int NodeNum;
int NNnum;
bool isSumDist = true;
//extern int Q.k;
//extern Query Q;

// end of TKDE2005

struct result {
	int poid;
	float distance;
};
result rs;

struct PartAddr {
	int part;
	int addr;
};

map<int, PartAddr> paID;


vector<TreeNode> EGT;
set<int> visitedNode; //record the  id of visitedNode treeNodeID

// for EGTD
struct edgeStatu {
	int vi;
	int isVisited; // hvisited = 0; visitedNode = 1;
	float sum;
	vector<float> qTbdist;
};
map<__int64, edgeStatu> edgevisited; // edgeID, isVisited
vector<__int64> queryedge;

// priority queue
struct pnode {
	bool isfilter;
	int nodeID;
	float dist;
	vector<vector<float>> bdist; // border-query points \ record the distance between Q to the border of this node
};

struct comparator
{
	bool operator () (const pnode& left, const pnode& right) const
	{
		return left.dist > right.dist;
	}
};
typedef	priority_queue<pnode, vector<pnode>, comparator> pnodePQ;
pnodePQ pq;

vector<map<int, vector<float>>> queryPath;

struct pvertex {
	int vertexID;
	float sumDist;
	vector<float> vdist; // vertex-query points \ record the distance between Q to the border of this node
};

struct vcomparator
{
	bool operator () (const pvertex& left, const pvertex& right) const
	{
		return left.sumDist > right.sumDist;
	}
};
typedef	priority_queue<pvertex, vector<pvertex>, vcomparator> pvertexPQ;
// pnodePQ pq;

int getAdjListGrpAddr(int NodeID)  	// using AdjFile
{
	int address = 0;
	if (algorithmId < 3) {
		address = paID[NodeID].addr;
		//return paID[NodeID].addr;
	}
	else {
		int addr = sizeof(int) + NodeID*sizeof(int), GrpAddr;
		char* BlockAddr = getFlatBlock(AdjFile, addr / BlkLen);
		memcpy(Ref(GrpAddr), BlockAddr + (addr%BlkLen), sizeof(int));
		address = GrpAddr;
	}
	return address;
}
/*
int partAddrLoad(const char* filename, map<int, PartAddr> &partID) {

	printf("loading partAddrFile\n");
	FILE *paFile;
	paFile = fopen(filename, "r+");
	CheckFile(paFile, filename);
	int nOfNode;
	fread(&nOfNode, sizeof(int), 1, paFile);
	//int loop = 0;
	for (int i = 0; i<nOfNode; i++) {
		int nodeID;
		PartAddr pa;
		int pt, addr;
		fread(&nodeID, sizeof(int), 1, paFile);
		fread(&pt, sizeof(int), 1, paFile);
		fread(&addr, sizeof(int), 1, paFile);
		pa.addr = addr;
		pa.part = pt;
		partID[nodeID] = pa;
		printf("i:%d, nodeid:%d, part:%d, addr:%d\n", i, nodeID, pa.part, pa.addr);
	}
	fclose(paFile);
	return nOfNode;
}
*/

int partAddrLoad(const char* filename, map<int, PartAddr> &partID) {

	printf("loading partAddrFile\n");
	//FILE *paFile;
	//paFile = fopen(filename, "r+");
	//CheckFile(paFile, filename);
	int nOfNode;
	//fread(&nOfNode, sizeof(int), 1, paFile);
	//int loop = 0;
	ifstream paFile(filename);
	string line;
	int loop = 0;
	while (!paFile.eof())
	{
		getline(paFile, line);
		//cout << loop++ << endl;
		if (line == "")
			continue;
		int nid, part, addr; 
		istringstream iss(line);
		// --M-- read the unrelevant coordinates inf
		iss >> nid >> part >> addr;
		PartAddr pa;
		pa.addr = addr;
		pa.part = part;
		partID[nid] = pa;
		loop++;
		//printf("i:%d, nodeid:%d, part:%d, addr:%d\n", loop, nid, pa.part, pa.addr);
	}
	nOfNode = partID.size();
	paFile.close();
	return nOfNode;
}

void comQueryPath(Query Q, vector<TreeNode> &EGT) {
	// 对于Q中的每个查询点，计算它的查询路径并保存到queryPath中
	int i = 0;
	vector<QueryPoint> ::iterator itr = Q.querypts.begin();
	for (; itr != Q.querypts.end(); itr++) {
		QueryPoint q = *itr;
		int pid = paID[q.Ni].part;
		// 计算叶子节点到查询点的距离
		int size = EGT[pid].borders.size();
		float sum = 0.0;
		int posi = find(EGT[pid].leafnodes.begin(), EGT[pid].leafnodes.end(), q.Ni) - EGT[pid].leafnodes.begin();
		int posj = find(EGT[pid].leafnodes.begin(), EGT[pid].leafnodes.end(), q.Nj) - EGT[pid].leafnodes.begin();
		int pi = 0, pj = 0;
		vector<float> dis;
		map<int, vector<float>> nodePath;
		for (int j = 0; j < EGT[pid].borders.size(); j++) {
			pi = j*EGT[pid].leafnodes.size() + posi;
			pj = j*EGT[pid].leafnodes.size() + posj;
			float tempi = EGT[pid].mind[pi];
			float tempj = EGT[pid].mind[pj];
			tempi += q.dist_Ni;
			tempj += (q.distEdge - q.dist_Ni);
			if (tempi < tempj) {
				dis.push_back(tempi);
			}
			else {
				dis.push_back(tempj);
			}
		}
		nodePath[pid] = dis;
		int preid;
		while (pid != 0) {
			//queryPath
			preid = pid;
			pid = EGT[pid].father;
			vector<float> bdis;
			for (int l = 0; l < EGT[pid].borders.size(); l++) {
				int border = EGT[pid].borders[l];
				posi = find(EGT[pid].union_borders.begin(), EGT[pid].union_borders.end(), border) - EGT[pid].union_borders.begin();
				float min = INFINITE_MAX;
				for (int k = 0; k < EGT[preid].borders.size(); k++) {
					int preborder = EGT[preid].borders[k];
					int dist;
					posj = find(EGT[pid].union_borders.begin(), EGT[pid].union_borders.end(), preborder) - EGT[pid].union_borders.begin();
					int index;
					// be careful
					if (posi > posj) {
						index = ((posi + 1)*posi) / 2 + posj;
						dist = EGT[pid].mind[index];
					}
					else {
						index = ((posj + 1)*posj) / 2 + posi;
						dist = EGT[pid].mind[index];
					}
					if (dist < min) min = dist;
				}
				bdis.push_back(min);
			}
			nodePath[pid] = bdis;
		}

		queryPath.push_back(nodePath);
		i++;
	}
	int f = 0;
}


void initialQuery(string fileprefix) {
	nOfPageAccessed = 0;
	nOfEdgeExpended = 0;
	if (algorithmId > 2) return;
	char tmpFileName[255];
	// load egtree
	string name;
	name.clear();
	name = fileprefix + "\\egtree";
	//sprintf(tmpFileName, "%s\\egtree", fileprefix);
	egtree_load(name.c_str(), EGT);
	// load partAddr
	name.clear();
	name = fileprefix + "\\partAddr";
	//sprintf(tmpFileName, "%s\\partAddr", fileprefix);
	//partAddrLoad(name.c_str(), paID);
	comQueryPath(Q, EGT);
}


//计算q到vertex的距离
vector<float> dijkstra_candidate(QueryPoint s, vector<int> &cands) {
	// init
	int AdjGrpAddr, AdjListSize, NewNodeID, PtGrpKey, PtNumOnEdge;
	float EdgeDist, PtDist;

	set<int> todo;
	todo.clear();
	todo.insert(cands.begin(), cands.end());

	unordered_map<int, float> result;
	result.clear();
	set<int> visitedNode;
	visitedNode.clear();
	unordered_map<int, float> q;
	q.clear();
	int vi = s.Ni;
	int vj = s.Nj;
	q[vi] = s.dist_Ni;
	q[vj] = s.distEdge - s.dist_Ni;

	// start
	int min, minpos, adjnode, weight;
	while (!todo.empty() && !q.empty()) {
		min = -1;
		for (unordered_map<int, float>::iterator it = q.begin(); it != q.end(); it++) {
			if (min == -1) {
				minpos = it->first;
				min = it->second;
			}
			else {
				if (it->second < min) {
					min = it->second;
					minpos = it->first;
				}
			}
		}

		// put min to result, add to visitedNode
		result[minpos] = min;
		visitedNode.insert(minpos);
		q.erase(minpos);

		if (todo.find(minpos) != todo.end()) {
			todo.erase(minpos);
		}

		// expand
		AdjGrpAddr = getAdjListGrpAddr(minpos);
		getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);

		for (int i = 0; i < AdjListSize; i++) {
			getVarE(ADJNODE_A, &adjnode, AdjGrpAddr, i);
			//adjnode = graph[minpos].adjnodes[i];
			if (visitedNode.find(adjnode) != visitedNode.end()) {
				continue;
			}			

			getVarE(DIST_A, &weight, AdjGrpAddr, i);
			//weight = graph[minpos].adjweight[i];

			if (q.find(adjnode) != q.end()) {
				if (min + weight < q[adjnode]) {
					q[adjnode] = min + weight;
				}
			}
			else {
				q[adjnode] = min + weight;
			}

		}
	}

	// output
	vector<float> output;
	for (int i = 0; i < cands.size(); i++) {
		output.push_back(result[cands[i]]);
	}

	// return
	return output;
}

float dijkstra(int nodei, int nodej, float upperbound) {
	// init
	int AdjGrpAddr, AdjListSize, NewNodeID, PtGrpKey, PtNumOnEdge;
	float EdgeDist, PtDist;
	set<int> todo;
	todo.clear();	
	todo.insert(nodej);

	unordered_map<int, float> result;
	result.clear();
	set<int> visitedNode;
	visitedNode.clear();
	unordered_map<int, float> q;
	q.clear();
	int vi = nodei;
	//int vj = nodej;
	q[vi] = 0;
	//q[vj] = edgedist;

	// start
	int min, minpos, adjnode, weight;
	while (!todo.empty() && !q.empty()) {
		min = -1;
		for (unordered_map<int, float>::iterator it = q.begin(); it != q.end(); it++) {
			if (min == -1) {
				minpos = it->first;
				min = it->second;
			}
			else {
				if (it->second < min) {
					min = it->second;
					minpos = it->first;
				}
			}
		}
		if (min > upperbound) break;
		// put min to result, add to visitedNode
		result[minpos] = min;
		visitedNode.insert(minpos);
		q.erase(minpos);

		if (todo.find(minpos) != todo.end()) {
			todo.erase(minpos);
		}

		// expand
		AdjGrpAddr = getAdjListGrpAddr(minpos);
		getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);

		for (int i = 0; i < AdjListSize; i++) {
			getVarE(ADJNODE_A, &adjnode, AdjGrpAddr, i);
			//adjnode = graph[minpos].adjnodes[i];
			if (visitedNode.find(adjnode) != visitedNode.end()) {
				continue;
			}

			getVarE(DIST_A, &weight, AdjGrpAddr, i);
			//weight = graph[minpos].adjweight[i];

			if (q.find(adjnode) != q.end()) {
				if (min + weight < q[adjnode]) {
					q[adjnode] = min + weight;
				}
			}
			else {
				q[adjnode] = min + weight;
			}

		}
	}

	// output
	float dis = INFINITE_MAX;
	if (result.count(nodej) > 0) {
		dis = result[nodej];
	}

	// return
	return dis;
}

void EGTDA( ) {
	int AdjGrpAddr, AdjListSize, NewNodeID, PtGrpKey, PtNumOnEdge;
	unsigned long long sumkwd;
	float EdgeDist, PtDist;

	// 初始化结果
	rs.distance = INFINITE_MAX;
	rs.poid = -1;
	edgevisited.clear();
	// step1 将根节点压入到优先队列中
	pnode pni;
	pni.dist = 0;
	pni.nodeID = 0;
	pq.push(pni);

	// step2 while循环，弹出当前最优的节点并展开
	while (!pq.empty()) {
		// step1 弹出top节点
		pnode node = pq.top();
		pq.pop();
		// if n没被过滤 Case1：叶子节点； Case2：中间节点
		int nid = node.nodeID;
		
		//printf("nodeid is:%d\n", nid);
		if ((node.dist < rs.distance) && (EGT[nid].unionkwds&Q.keywords) == Q.keywords) {
			// Case1: 叶子节点
			if (EGT[nid].isleaf) {
				// 计算其中最优的点并更新rs.distance（内部和外部两种情况，内部用Dijkstra外部用边界距离）
				// 计算每个查询点到内部所有vertex的距离
				vector<pvertex> pv;
				for (int i = 0; i < Q.k; i++) {
					// 两种情况，case1：当前该查询点在该节点中，调用Dijkstra;case2:不在，直接相加
					//printf("The nid is:%d, q num is:%d\n", nid, i);
					map<int, vector<float>>::iterator it = queryPath[i].find(nid);
					vector<float> result;
					if (it != queryPath[i].end()) { // 在里面，调用Dijkstra算法
						// note the corresponding relationship between result and leafnodes
						result = dijkstra_candidate(Q.querypts[i], EGT[nid].leafnodes);

						for (int l = 0; l < EGT[nid].leafnodes.size(); l++) {
							if (i == 0) {
								pvertex pt;
								vector<float> vf(Q.k, 0);
								pt.vdist = vf;

								pt.vertexID = EGT[nid].leafnodes[l];
								//pt.vdist.push_back(result[l]);
								pt.vdist[i] = result[l];
								pt.sumDist = result[l];
								pv.push_back(pt);
							}
							else {
								pv[l].vdist[i] = result[l];								
								pv[l].sumDist += result[l];
							}
						}
					}

					else { // 直接相加
						for (int l = 0; l < EGT[nid].leafnodes.size(); l++) {
							if (i == 0) {
								pvertex pt;
								vector<float> vf(Q.k, 0);
								pt.vdist = vf;
								pt.vertexID = EGT[nid].leafnodes[l];

								int dist;
								float min = INFINITE_MAX;
								for (int m = 0; m < EGT[nid].borders.size(); m++) {
									int pos = m*EGT[nid].leafnodes.size() + l;
									dist = EGT[nid].mind[pos] + node.bdist[m][i];
									if (dist < min) min = dist;
								}
								// 可能有问题，vdist未初始化
								pt.vdist[i] = min; 
								pt.sumDist = min;
								pv.push_back(pt);
							}
							else {
								// 计算到q-border-vertex的距离
								int dist;
								float min = INFINITE_MAX;
								for (int m = 0; m < EGT[nid].borders.size(); m++) {
									int pos = m*EGT[nid].leafnodes.size() + l;
									dist = EGT[nid].mind[pos] + node.bdist[m][i];
									if (dist < min) min = dist;
								}
								// 可能有问题，vdist未初始化
								pv[l].vdist[i] = min; 
								pv[l].sumDist += min;
							}
						}
					}
				}

				// 进一步优化在这里************
				// 处理每一个vertex相邻的边，看上面是否有满足条件的POI 
				for (int t = 0; t < pv.size(); t++) {
					//printf("loop %d to num %d\n", t, pv.size());
					pvertex pt = pv[t];
					int vid = pt.vertexID;
					__int64 vide = vid;
					AdjGrpAddr = getAdjListGrpAddr(vid);
					getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);
					for (int n = 0; n < AdjListSize; n++) {
						getVarE(ADJNODE_A, Ref(NewNodeID), AdjGrpAddr, n);
						getVarE(DIST_A, Ref(EdgeDist), AdjGrpAddr, n);
						getVarE(SUMKWD_A, &sumkwd, AdjGrpAddr, n);
						__int64 NewNodeIDe = NewNodeID;
						__int64 edge = getKey(vide, NewNodeIDe);
						// bool isonedge;
						
						map<__int64, edgeStatu>::iterator it = edgevisited.find(edge);
						edgeStatu es;
						if (it != edgevisited.end()) {
							if (edgevisited[edge].isVisited == 0) { //处理一下
								getVarE(PTKEY_A, Ref(PtGrpKey), AdjGrpAddr, n);
								if (PtGrpKey == -1) {
									cout<<"No POI existed on Edge where Q located."<<endl;
								}
								else {
									nOfEdgeExpended++;
									//getVarE(DIST_A, Ref(EdgeDist), AdjGrpAddr, n);
									getFixedF(SIZE_P, Ref(PtNumOnEdge), PtGrpKey);
									//Notice the order POI visitedNode on edge
									unsigned long long pkwd;
									int pid;
									for (int t = 0; t < PtNumOnEdge; t++) {
										getVarE(PT_DIST, Ref(PtDist), PtGrpKey, t);
										getVarE(PT_KWD, &pkwd, PtGrpKey, t);
										getVarE(PT_P, &pid, PtGrpKey, t);
										if ((pkwd&Q.keywords) != Q.keywords) {
											continue;
										}
										else {
											edgeStatu es = it->second;
											float sum = 0;
											for (int f = 0; f < Q.k; f++) {
												// two case:1,q on this edge; 2: q not on this edge
												//if (find(queryedge.begin(), queryedge.end(), edge) != queryedge.end()) { // on edge
												if (queryedge[f] == edge) { // on edge
													float dis = 0;
													dis = fabsf(PtDist - Q.querypts[f].dist_Ni);
													if (dis < (0.5*EdgeDist)) {
														sum += dis;
													}
													else {
														float dist = dijkstra(vid, NewNodeID, dis);
														dist = (EdgeDist - dis) + dist;
														if (dis > dist) dis = dist;
														sum += dis;
													}													
													continue;
												}
												float disi = 0;
												float disj = 0;

												if (vid < es.vi) {
													disi = PtDist + pt.vdist[f];
													disj = EdgeDist - PtDist + es.qTbdist[f];
												}
												else {
													disi = PtDist + es.qTbdist[f];
													disj = EdgeDist - PtDist + pt.vdist[f];
												}
												if (disi > disj) disi = disj;
												sum = sum + disi;
											}
											// have computed the distance between each points to query:sum
											if (sum < rs.distance) {
												rs.distance = sum;
												rs.poid = pid;
											}
										}
									}
								}

								//edgevisited[edge].isVisited = 1;
							}
							// 两个顶点都访问过该条边了，可以删除了
							edgevisited.erase(it);
						}
						else { // 该条边未访问过
							
							if ((sumkwd&Q.keywords) != Q.keywords) {
								es.isVisited = 1;
								edgevisited[edge] = es;
								continue;
							}

							es.isVisited = 0;
							es.vi = vid;
							es.qTbdist = pt.vdist;
							es.sum = pt.sumDist;
							edgevisited[edge] = es;
						}
						//getVarE(ADJNODE_A, Ref(NewNodeID), AdjGrpAddr, n);
						//getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
					}


				}

			}
			// Case2: 中间节点
			else {
				// 拓展该节点，并计算满足条件(关键字和距离)的孩子节点的距离并压入队列
				// 处理每一个孩子节点
				vector<pnode> childnode;
				childnode.clear();

				for (int i = 0; i < Q.k; i++) { 
					// 按照查询路径来处理，处理每个查询点
					for (int j = 0; j < EGT[nid].children.size(); j++) {
						// two cases. 1: q in the node; 2: not in the node
						int cid = EGT[nid].children[j];
						
						if (i == 0) {
							pnode pn;
							pn.isfilter = true;
							pn.nodeID = cid;
							pn.dist = 0;
							for (int f = 0; f < EGT[cid].borders.size(); f++) {
								vector<float> vf(Q.k, 0);
								pn.bdist.push_back(vf);								
							}
							childnode.push_back(pn);
						}
						if ((EGT[cid].unionkwds&Q.keywords) != Q.keywords) continue;
						childnode[j].isfilter = false;
						map<int, vector<float>>::iterator it = queryPath[i].find(cid);
						if (it != queryPath[i].end()) { // do nothing
							continue;						
						}
						else { // compute the distance
							   // first determine the father node in the layer of cid
							int sid;
							bool flag = false;
							for (int l = 0; l < EGT[nid].children.size(); l++) {
								sid = EGT[nid].children[l];
								map<int, vector<float>>::iterator sit = queryPath[i].find(sid);
								if (sit != queryPath[i].end()) { // return this id
									flag = true;
									break;
								}
							}
							if (!flag) {
								sid = nid;
							}
							// then compute the distance from nid
							float mind = INFINITE_MAX;
							for (int k = 0; k < EGT[cid].borders.size(); k++) {
								// 计算每个查询点到所有border的距离
								int bdk = EGT[cid].borders[k];
								int posi = find(EGT[nid].union_borders.begin(), EGT[nid].union_borders.end(), bdk)- EGT[nid].union_borders.begin();
								float min = INFINITE_MAX;
								for (int t = 0; t < EGT[sid].borders.size(); t++) {
									int bdt = EGT[sid].borders[t];
									int posj = find(EGT[nid].union_borders.begin(), EGT[nid].union_borders.end(), bdt)- EGT[nid].union_borders.begin();
									// be careful
									int index;
									float dist;
									if (posi > posj) {
										index = ((posi + 1)*posi) / 2 + posj;
										dist = EGT[nid].mind[index];
									}
									else {
										index = ((posj + 1)*posj) / 2 + posi;
										dist = EGT[nid].mind[index];
									}
									if (flag) { // 
										map<int, vector<float>>::iterator nodeDist = queryPath[i].find(sid);
										dist = dist + nodeDist->second[t];
									}
									else {
										float distance = node.bdist[t][i];
										dist = dist + distance;
									}
									
									if (dist < min) min = dist;
								}
								if (min < mind) mind = min;
								// if is the first then initiate it
								childnode[j].bdist[k][i] = min;
								/*
								if (i == 0) {
									vector<float> qTb(Q.k,0);
									// wrong!!!!!!!!!
									qTb[i] = min;
									//pn.dist = min;
									pn.bdist.push_back(qTb);
								}
								else {
									// wrong!!!!!!!!!
									//pn.dist = pn.dist + min;
									pn.bdist[k][i] = min;
								}
								*/

							}
							childnode[j].dist = childnode[j].dist + mind;
							//pq.push(pn);
						}
					}
				}
				int size = childnode.size();
				for (int loop = 0; loop < size; loop++) {
					if (childnode[loop].isfilter) continue;
					pq.push(childnode[loop]);
				}
			}
		}

	}

	// step3 处理half
	map<__int64, edgeStatu> ::iterator it = edgevisited.begin();
	for (; it != edgevisited.end(); it++) {
		edgeStatu es = it->second;
		// 可能在边上
		__int64 edgeid = it->first;
		vector<__int64> ::iterator itr = find(queryedge.begin(), queryedge.end(), edgeid);
		if ((itr == queryedge.end()) && (es.sum >= rs.distance)) continue;	

		//__int64 vide, vjde;
		__int64 firste = it->first;
		int vid, vjd;
		breakKey(firste, vid, vjd);
		if (vid != es.vi) {
			int temp = vid;
			vid = vjd;
			vjd = temp;
		}
		// 可以优化在这里***********
		AdjGrpAddr = getAdjListGrpAddr(vid);
		getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);
		for (int n = 0; n < AdjListSize; n++) {
			getVarE(ADJNODE_A, Ref(NewNodeID), AdjGrpAddr, n);
			if (NewNodeID != vjd) continue;
			getVarE(DIST_A, Ref(EdgeDist), AdjGrpAddr, n);
			//记录查询边距离
			getVarE(PTKEY_A, Ref(PtGrpKey), AdjGrpAddr, n);
			if (PtGrpKey == -1) {
				//cout<<"No POI existed on Edge where Q located."<<endl;
			}
			else {
				nOfEdgeExpended++;
				getFixedF(SIZE_P, Ref(PtNumOnEdge), PtGrpKey);
				//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
				//Notice the order POI visitedNode on edge
				unsigned long long pkwd;
				int pid;
				for (int j = 0; j < PtNumOnEdge; j++) {
					//poicnt++;
					getVarE(PT_DIST, Ref(PtDist), PtGrpKey, j);
					getVarE(PT_KWD, &pkwd, PtGrpKey, j);
					getVarE(PT_P, &pid, PtGrpKey, j);
					if ((pkwd&Q.keywords) != Q.keywords) {
						continue;
					}
					else {
						float sum = 0;
						for (int f = 0; f < Q.k; f++) {
							// two case:1,q on this edge; 2: q not on this edge
							if (queryedge[f] == edgeid) { // on edge
								float dis = 0;
								dis = fabsf(PtDist - Q.querypts[f].dist_Ni);
								sum += dis;
								continue;
							}
							float dis;
							if (vid < vjd) {
								dis = PtDist;
								//sum = sum + Q.k*dis;
							}
							else {
								dis = EdgeDist - PtDist;
								//sum = sum + Q.k*dis;
							}
							sum = sum + es.qTbdist[f] + dis;								
						}
						if (sum < rs.distance) {
							rs.distance = sum;
							rs.poid = pid;
						}

					}
				}
			}
		}
		

	}
	printf("The result of TDA distance is: %f\n", rs.distance);
	printf("The result of TDA oid is: %d\n", rs.poid);
	printf("The result of TDA edgeexpand is: %d\n", nOfEdgeExpended);
}

bool iscontained(set<unsigned long long> coverkws, unsigned long long kwd) {
	set<unsigned long long> ::iterator it = coverkws.begin();
	bool iscontained = false;
	for (; it != coverkws.end(); it++) {
		unsigned long long key = *it;
		if ((key&kwd) == kwd) {
			iscontained = true;
			break;
		}
	}
	return iscontained;
}

void EGETDA() {

	int AdjGrpAddr, AdjListSize, NewNodeID, PtGrpKey, PtNumOnEdge;
	unsigned long long sumkwd;
	float EdgeDist, PtDist;

	// 初始化结果
	rs.distance = INFINITE_MAX;
	rs.poid = -1;
	edgevisited.clear();
	// step1 将根节点压入到优先队列中
	pnode pni;
	pni.dist = 0;
	pni.nodeID = 0;
	pq.push(pni);

	// step2 while循环，弹出当前最优的节点并展开
	while (!pq.empty()) {
		// step1 弹出top节点
		pnode node = pq.top();
		pq.pop();
		// if n没被过滤 Case1：叶子节点； Case2：中间节点
		int nid = node.nodeID;
		//printf("node is :%d\n", nid);
		//printf("nodeid is:%d\n", nid);
		if ((node.dist < rs.distance) && (EGT[nid].unionkwds&Q.keywords) == Q.keywords) {
			// Case1: 叶子节点
			if (EGT[nid].isleaf) {
				// 计算其中最优的点并更新rs.distance（内部和外部两种情况，内部用Dijkstra外部用边界距离）
				// 计算每个查询点到内部所有vertex的距离
				vector<pvertex> pv;
				for (int i = 0; i < Q.k; i++) {
					// 两种情况，case1：当前该查询点在该节点中，调用Dijkstra;case2:不在，直接相加
					//printf("The nid is:%d, q num is:%d\n", nid, i);
					map<int, vector<float>>::iterator it = queryPath[i].find(nid);
					vector<float> result;
					if (it != queryPath[i].end()) { // 在里面，调用Dijkstra算法
													// note the corresponding relationship between result and leafnodes
						result = dijkstra_candidate(Q.querypts[i], EGT[nid].leafnodes);

						for (int l = 0; l < EGT[nid].leafnodes.size(); l++) {
							if (i == 0) {
								pvertex pt;
								vector<float> vf(Q.k, 0);
								pt.vdist = vf;

								pt.vertexID = EGT[nid].leafnodes[l];
								//pt.vdist.push_back(result[l]);
								pt.vdist[i] = result[l];
								pt.sumDist = result[l];
								pv.push_back(pt);
							}
							else {
								pv[l].vdist[i] = result[l];
								pv[l].sumDist += result[l];
							}
						}
					}

					else { // 直接相加
						for (int l = 0; l < EGT[nid].leafnodes.size(); l++) {
							if (i == 0) {
								pvertex pt;
								vector<float> vf(Q.k, 0);
								pt.vdist = vf;
								pt.vertexID = EGT[nid].leafnodes[l];

								int dist;
								float min = INFINITE_MAX;
								for (int m = 0; m < EGT[nid].borders.size(); m++) {
									int pos = m*EGT[nid].leafnodes.size() + l;
									dist = EGT[nid].mind[pos] + node.bdist[m][i];
									if (dist < min) min = dist;
								}
								// 可能有问题，vdist未初始化
								pt.vdist[i] = min;
								pt.sumDist = min;
								pv.push_back(pt);
							}
							else {
								// 计算到q-border-vertex的距离
								int dist;
								float min = INFINITE_MAX;
								for (int m = 0; m < EGT[nid].borders.size(); m++) {
									int pos = m*EGT[nid].leafnodes.size() + l;
									dist = EGT[nid].mind[pos] + node.bdist[m][i];
									if (dist < min) min = dist;
								}
								// 可能有问题，vdist未初始化
								pv[l].vdist[i] = min;
								pv[l].sumDist += min;
							}
						}
					}
				}

				// 进一步优化在这里************
				// 处理每一个vertex相邻的边，看上面是否有满足条件的POI 
				for (int t = 0; t < pv.size(); t++) {
					//printf("loop %d to num %d\n", t, pv.size());
					pvertex pt = pv[t];
					int vid = pt.vertexID;
					__int64 vide = vid;
					AdjGrpAddr = getAdjListGrpAddr(vid);
					getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);
					for (int n = 0; n < AdjListSize; n++) {
						getVarE(ADJNODE_A, Ref(NewNodeID), AdjGrpAddr, n);
						getVarE(DIST_A, Ref(EdgeDist), AdjGrpAddr, n);
						getVarE(SUMKWD_A, &sumkwd, AdjGrpAddr, n);
						__int64 NewNodeIDe = NewNodeID;
						__int64 edge = getKey(vide, NewNodeIDe);
						// bool isonedge;

						map<__int64, edgeStatu>::iterator it = edgevisited.find(edge);
						edgeStatu es;
						if (it != edgevisited.end()) {
							if (edgevisited[edge].isVisited == 0) { //处理一下
								getVarE(PTKEY_A, Ref(PtGrpKey), AdjGrpAddr, n);
								if (PtGrpKey == -1) {
									cout << "No POI existed on Edge where Q located." << endl;
								}
								else {
									nOfEdgeExpended++;
									//getVarE(DIST_A, Ref(EdgeDist), AdjGrpAddr, n);
									//************ optimize in this part: 1.two vertex 2.kwds
									edgeStatu es = it->second;
									vector<__int64> ::iterator itq = queryedge.begin();
									if ((itq == queryedge.end()) && ((es.sum >= rs.distance) && (pt.sumDist >= rs.distance))) continue;
									getFixedF(SIZE_P, Ref(PtNumOnEdge), PtGrpKey);
									//Notice the order POI visitedNode on edge
									unsigned long long pkwd;
									int pid;
									for (int t = 0; t < PtNumOnEdge; t++) {
										getVarE(PT_DIST, Ref(PtDist), PtGrpKey, t);
										getVarE(PT_KWD, &pkwd, PtGrpKey, t);
										getVarE(PT_P, &pid, PtGrpKey, t);
										if ((pkwd&Q.keywords) != Q.keywords) {
											continue;
										}
										else {
											
											float sum = 0;
											for (int f = 0; f < Q.k; f++) {
												// two case:1,q on this edge; 2: q not on this edge
												//if (find(queryedge.begin(), queryedge.end(), edge) != queryedge.end()) { // on edge
												if (queryedge[f] == edge) { // on edge
													float dis = 0;
													dis = fabsf(PtDist - Q.querypts[f].dist_Ni);
													if (dis < (0.5*EdgeDist)) {
														sum += dis;
													}
													else {
														float dist = dijkstra(vid, NewNodeID, dis);
														dist = (EdgeDist - dis) + dist;
														if (dis > dist) dis = dist;
														sum += dis;
													}
													continue;
												}
												float disi = 0;
												float disj = 0;

												if (vid < es.vi) {
													disi = PtDist + pt.vdist[f];
													disj = EdgeDist - PtDist + es.qTbdist[f];
												}
												else {
													disi = PtDist + es.qTbdist[f];
													disj = EdgeDist - PtDist + pt.vdist[f];
												}
												if (disi > disj) disi = disj;
												sum = sum + disi;
											}
											// have computed the distance between each points to query:sum
											if (sum < rs.distance) {
												rs.distance = sum;
												rs.poid = pid;
											}
										}
									}
								}

								//edgevisited[edge].isVisited = 1;
							}
							// 两个顶点都访问过该条边了，可以删除了
							edgevisited.erase(it);
						}
						else { // 该条边未访问过

							if ((sumkwd&Q.keywords) != Q.keywords) {
								es.isVisited = 1;
								edgevisited[edge] = es;
								continue;
							}

							es.isVisited = 0;
							es.vi = vid;
							es.qTbdist = pt.vdist;
							es.sum = pt.sumDist;
							edgevisited[edge] = es;
						}
						//getVarE(ADJNODE_A, Ref(NewNodeID), AdjGrpAddr, n);
						//getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
					}

				}

			}
			// Case2: 中间节点
			else {
				// 拓展该节点，并计算满足条件(关键字和距离)的孩子节点的距离并压入队列
				// 处理每一个孩子节点
				vector<pnode> childnode;
				childnode.clear();

				for (int i = 0; i < Q.k; i++) {
					// 按照查询路径来处理，处理每个查询点
					for (int j = 0; j < EGT[nid].children.size(); j++) {
						// two cases. 1: q in the node; 2: not in the node
						
						int cid = EGT[nid].children[j];
						
						if (i == 0) {
							pnode pn;
							pn.isfilter = true;
							pn.nodeID = cid;
							pn.dist = 0;
							for (int f = 0; f < EGT[cid].borders.size(); f++) {
								vector<float> vf(Q.k, 0);
								pn.bdist.push_back(vf);
							}
							childnode.push_back(pn);
						}

						if ((EGT[cid].unionkwds&Q.keywords) != Q.keywords || !iscontained(EGT[cid].coverkwds, Q.keywords)) continue;
						childnode[j].isfilter = false;

						map<int, vector<float>>::iterator it = queryPath[i].find(cid);
						if (it != queryPath[i].end()) { // do nothing
							continue;
						}
						else { // compute the distance
							   // first determine the father node in the layer of cid
							int sid;
							bool flag = false;
							for (int l = 0; l < EGT[nid].children.size(); l++) {
								sid = EGT[nid].children[l];
								map<int, vector<float>>::iterator sit = queryPath[i].find(sid);
								if (sit != queryPath[i].end()) { // return this id
									flag = true;
									break;
								}
							}
							if (!flag) {
								sid = nid;
							}
							// then compute the distance from nid
							float mind = INFINITE_MAX;
							for (int k = 0; k < EGT[cid].borders.size(); k++) {
								// 计算每个查询点到所有border的距离
								int bdk = EGT[cid].borders[k];
								int posi = find(EGT[nid].union_borders.begin(), EGT[nid].union_borders.end(), bdk) - EGT[nid].union_borders.begin();
								float min = INFINITE_MAX;
								for (int t = 0; t < EGT[sid].borders.size(); t++) {
									int bdt = EGT[sid].borders[t];
									int posj = find(EGT[nid].union_borders.begin(), EGT[nid].union_borders.end(), bdt) - EGT[nid].union_borders.begin();
									// be careful
									int index;
									float dist;
									if (posi > posj) {
										index = ((posi + 1)*posi) / 2 + posj;
										dist = EGT[nid].mind[index];
									}
									else {
										index = ((posj + 1)*posj) / 2 + posi;
										dist = EGT[nid].mind[index];
									}
									if (flag) { // 
										map<int, vector<float>>::iterator nodeDist = queryPath[i].find(sid);
										dist = dist + nodeDist->second[t];
									}
									else {
										float distance = node.bdist[t][i];
										dist = dist + distance;
									}

									if (dist < min) min = dist;
								}
								if (min < mind) mind = min;
								// if is the first then initiate it
								childnode[j].bdist[k][i] = min;
								/*
								if (i == 0) {
								vector<float> qTb(Q.k,0);
								// wrong!!!!!!!!!
								qTb[i] = min;
								//pn.dist = min;
								pn.bdist.push_back(qTb);
								}
								else {
								// wrong!!!!!!!!!
								//pn.dist = pn.dist + min;
								pn.bdist[k][i] = min;
								}
								*/

							}
							childnode[j].dist = childnode[j].dist + mind;
							//pq.push(pn);
						}
					}
				}
				int size = childnode.size();
				for (int loop = 0; loop < size; loop++) {
					if (childnode[loop].isfilter) continue;
					pq.push(childnode[loop]);
				}
			}
		}

	}

	// step3 处理half
	map<__int64, edgeStatu> ::iterator it = edgevisited.begin();
	for (; it != edgevisited.end(); it++) {
		edgeStatu es = it->second;
		// 可能在边上
		__int64 edgeid = it->first;
		vector<__int64> ::iterator itr = find(queryedge.begin(), queryedge.end(), edgeid);
		if ((itr == queryedge.end()) && (es.sum >= rs.distance)) continue;

		//__int64 vide, vjde;
		__int64 firste = it->first;
		int vid, vjd;
		breakKey(firste, vid, vjd);
		if (vid != es.vi) {
			int temp = vid;
			vid = vjd;
			vjd = temp;
		}
		// 可以优化在这里***********
		AdjGrpAddr = getAdjListGrpAddr(vid);
		getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);
		for (int n = 0; n < AdjListSize; n++) {
			getVarE(ADJNODE_A, Ref(NewNodeID), AdjGrpAddr, n);
			if (NewNodeID != vjd) continue;
			getVarE(DIST_A, Ref(EdgeDist), AdjGrpAddr, n);
			//记录查询边距离
			getVarE(PTKEY_A, Ref(PtGrpKey), AdjGrpAddr, n);
			if (PtGrpKey == -1) {
				//cout<<"No POI existed on Edge where Q located."<<endl;
			}
			else {
				nOfEdgeExpended++;
				getFixedF(SIZE_P, Ref(PtNumOnEdge), PtGrpKey);
				//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
				//Notice the order POI visitedNode on edge
				unsigned long long pkwd;
				int pid;
				for (int j = 0; j < PtNumOnEdge; j++) {
					//poicnt++;
					getVarE(PT_DIST, Ref(PtDist), PtGrpKey, j);
					getVarE(PT_KWD, &pkwd, PtGrpKey, j);
					getVarE(PT_P, &pid, PtGrpKey, j);
					if ((pkwd&Q.keywords) != Q.keywords) {
						continue;
					}
					else {
						float sum = 0;
						for (int f = 0; f < Q.k; f++) {
							// two case:1,q on this edge; 2: q not on this edge
							if (queryedge[f] == edgeid) { // on edge
								float dis = 0;
								dis = fabsf(PtDist - Q.querypts[f].dist_Ni);
								sum += dis;
								continue;
							}
							float dis;
							if (vid < vjd) {
								dis = PtDist;
								//sum = sum + Q.k*dis;
							}
							else {
								dis = EdgeDist - PtDist;
								//sum = sum + Q.k*dis;
							}
							sum = sum + es.qTbdist[f] + dis;
						}
						if (sum < rs.distance) {
							rs.distance = sum;
							rs.poid = pid;
						}

					}
				}
			}
		}


	}
	printf("The result of ETDA distance is: %f\n", rs.distance);
	printf("The result of ETDA oid is: %d\n", rs.poid);
	printf("The result of ETDA edgeexpended is: %d\n", nOfEdgeExpended);
}

int k = 7;
int qk = 3;
int np = 7;
int nOfTest = 1;
struct querycand {
	int nodei;
	int nodej;
	float length;
	bool visited;
};
map<int, querycand>  qc;

void geneQuery(string prxfilename) {
	string name;
	name.clear();

	name = prxfilename + "\\partAddr";
	//sprintf(tmpFileName, "%s\\partAddr", fileprefix);
	NodeNum = partAddrLoad(name.c_str(), paID);
	//NodeNum = paID.size();
	// read info from java transform file, 
	// format: edgeid nodei nodej length

	string inFile = prxfilename + "\\candquery";
	string outFile = prxfilename + "\\queryfile";

	ifstream iFile(inFile.c_str());
	ofstream oFile(outFile.c_str());
	int loop = 0;
	while (!iFile.eof())
	{
		string line;
		getline(iFile, line);
		if (line == "")
			continue;
		istringstream iss(line);
		int edgeid, nodei, nodej;
		float length;
		iss >> edgeid >> nodei >> nodej >> length;
		if (paID[nodei].part != paID[nodej].part) continue;
		querycand qu;
		qu.nodei = nodei;
		qu.nodej = nodej;
		qu.length = length;
		qu.visited = false;
		qc[loop] = qu;
		loop++;
	}
	int size = qc.size();
	// generate Q.keywords
	int l = 0;
	while (l < nOfTest) {
		bitset<MAX_KEYWORDN> keywords_set;
		keywords_set.reset();
		while (keywords_set.count() != qk) {
			int position = rand() % MAX_KEYWORDN;
			keywords_set.set(position);
		}
		unsigned long long kwd = keywords_set.to_ullong();
		oFile << kwd << " " << k;
		// generate Q.querypts
		for (int i = 0; i < k; ) {
			int ni, nj;
			float dis;
			int position = rand() % size;
			querycand qu = qc[position];
			if (qu.visited) continue;
			i++;
			qu.visited = true;
			int randnum = 0;
			while (randnum == 0 || randnum == RAND_MAX) {
				randnum = rand();
			}
			dis = (randnum*qu.length) / RAND_MAX;
			oFile << "," << qu.nodei << " " << qu.nodej << " " << qu.length << " " << dis;
		}
		l++;
	}
	iFile.close();
	oFile.close();
}

//*********** start of TKDE2005
// compute the Euclidean distance
float getDist(float x1, float y1, float x2, float y2) {
	return float(sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1)));
}
// compute the network distance 
// vj: the poi point on the edge
// onSameEdge: whether this edge of VJ is same with the query point locateded in
// eKdist: the length of this edge

inline float getNetworkDist(float vJ, BitStore& onSameEdge, float eKdist) {
	float dist = 0;
	for (int sp = 0; sp<Q.k; sp++) {
		float distL = partdists[sp][0] + vJ;
		float distR = partdists[sp][1] + eKdist - vJ;
		float distT = min(distL, distR);
		if (onSameEdge[sp]) {
			float distQ = ValAbs(Q.querypts[sp].dist_Ni - vJ);
			if (distQ<distT) distT = distQ;
		}

		compdist[sp] = distT;	// for log (log the original one) !
								//if (isWeightUse) distT = distT*weights[sp];

		if (isSumDist)
			dist += distT;
		else
			dist = max(distT, dist);
	}
	return dist;
}

inline void InitDQ(DistQueue& dQ) {
	dQ.empty();
	elem_map.clear();
	//	while (!dQ.empty()) dQ.pop();	// clear kNN-dist heap
}


// the DQ is the storage of the possible results
inline void processDQ(DistQueue& dQ, float& topdist, float newdist, int id) {
	if (elem_map.count(id)>0) {
		if (elem_map[id] <= newdist) return;	// topdist unchanged
		// remove itself
		elem_map[id] = newdist;
		for (int i = 0; i<dQ.size(); i++)
			if (dQ[i].id == id) {
				dQ.erase(dQ.begin() + i);
				break;
			}

		//if (dQ.size()>=NNnum) topdist=dQ[NNnum-1].dist;	// still quite cheap !
	}

	DistElem new_elem;
	new_elem.id = id;		new_elem.dist = newdist;
	if (dQ.size()<NNnum) {		// bestdist still infinity
		elem_map[id] = newdist;	// unsorted here !
		dQ.push_back(new_elem);

		// must be sorted in order for the case below to be valid !
		int pos = dQ.size() - 1;
		while (pos>0 && dQ[pos - 1].dist >= newdist) pos--;
		for (int i = dQ.size() - 1; i>pos; i--) dQ[i] = dQ[i - 1];
		dQ[pos] = new_elem;
	}
	else {
		if (newdist<dQ[NNnum - 1].dist) {
			int old_id = dQ[NNnum - 1].id;
			elem_map.erase(old_id);
			elem_map[id] = newdist;

			int pos = NNnum - 1;
			while (pos>0 && dQ[pos - 1].dist >= newdist) pos--;
			for (int i = NNnum - 1; i>pos; i--) dQ[i] = dQ[i - 1];
			dQ[pos] = new_elem;
		}
	}
	if (dQ.size() >= NNnum) topdist = dQ[NNnum - 1].dist;	// still quite cheap !													
}

// utilize the vj to update the distqueue
inline void UpdateOnePt(float vJ, bool& hasChanged, BitStore& onSameEdge, float eKdist, int ptid) {
	float dist = getNetworkDist(vJ, onSameEdge, eKdist);
	if (dist<bestdist) {	// for S-dist
		hasChanged = true;
		processDQ(dQ_graph, bestdist, dist, ptid);
		bestcompdist.assign(compdist.begin(), compdist.end());
	}
}

// 0 large 1 small // get a possible optimal position for in this edge
inline float getOptPos(int NodeID, int NewNodeID, float eKdist, float& bestPos) {
	float vJ, bestVal, tmpPos, tmpVal;
	int sid = 0, eid = 1;
	if (NewNodeID<NodeID) { sid = 1; eid = 0; }

	for (int sp = 0; sp<Q.k; sp++) {
		partdists[sp][0] = partdists[sp][1] = MAX_DIST;
		if (DistMaps[sp].count(NodeID)>0)
			partdists[sp][sid] = DistMaps[sp][NodeID];
		if (DistMaps[sp].count(NewNodeID)>0)
			partdists[sp][eid] = DistMaps[sp][NewNodeID];
		int qNi = Q.querypts[sp].Ni, qNj = Q.querypts[sp].Nj;
		onSameEdge[sp] = (NodeID == qNi&&NewNodeID == qNj) || (NodeID == qNj&&NewNodeID == qNi);
	}

	// for isSumDist only
	bestPos = 0;
	bestVal = getNetworkDist(bestPos, onSameEdge, eKdist);

	tmpPos = eKdist;
	tmpVal = getNetworkDist(tmpPos, onSameEdge, eKdist);
	if (tmpVal<bestVal) { bestPos = tmpPos; bestVal = tmpVal; }

	int NUM_STEPS = 100;
	if (!isWeightUse) {
		if (isSumDist) {
			for (int sp = 0; sp<Q.k; sp++) {
				if (onSameEdge[sp]) {
					tmpPos = Q.querypts[sp].dist_Ni;
					tmpVal = getNetworkDist(tmpPos, onSameEdge, eKdist);
					if (tmpVal<bestVal) { bestPos = tmpPos; bestVal = tmpVal; }
				}
			}
		}
		else {
			// can find a loose lb for ...
			//bestVal-=eKdist;
			for (int j = 0; j<NUM_STEPS; j++) {	// scan
				tmpPos = j*eKdist / (NUM_STEPS - 1);
				tmpVal = getNetworkDist(tmpPos, onSameEdge, eKdist);
				if (tmpVal<bestVal) { bestPos = tmpPos; bestVal = tmpVal; }
			}
		}
	}
	else {
		for (int j = 0; j<NUM_STEPS; j++) {	// scan
			tmpPos = j*eKdist / (NUM_STEPS - 1);
			tmpVal = getNetworkDist(tmpPos, onSameEdge, eKdist);
			if (tmpVal<bestVal) { bestPos = tmpPos; bestVal = tmpVal; }
		}
	}
	return bestVal;
}


bool UpdatePoints(int NodeID, int NewNodeID, int PtGrpAddr, float eKdist) {
	float vJ;
	bool hasChanged = false;
	unsigned long long sumkwd;
	int sid = 0, eid = 1;
	if (NewNodeID<NodeID) { sid = 1; eid = 0; }

	for (int sp = 0; sp<Q.k; sp++) {
		partdists[sp][0] = partdists[sp][1] = MAX_DIST;
		if (DistMaps[sp].count(NodeID)>0)
			partdists[sp][sid] = DistMaps[sp][NodeID];
		if (DistMaps[sp].count(NewNodeID)>0)
			partdists[sp][eid] = DistMaps[sp][NewNodeID];
		int qNi = Q.querypts[sp].Ni, qNj = Q.querypts[sp].Nj;
		onSameEdge[sp] = (NodeID == qNi&&NewNodeID == qNj) || (NodeID == qNj&&NewNodeID == qNi);
	}

	if (isNNquery) {
		int PtGrpSize, dummy;
		//PtGrpAddr = pointQuery(PtTree, PtGrpKey, dummy);
		nOfEdgeExpended++;
		getFixedF(SIZE_P, Ref(PtGrpSize), PtGrpAddr);
		for (int j = 0; j<PtGrpSize; j++) {	// scan
			int pid = 0;
			getVarE(PT_P, Ref(pid), PtGrpAddr, j);
			getVarE(PT_DIST, Ref(vJ), PtGrpAddr, j);
			getVarE(PT_KWD, Ref(sumkwd), PtGrpAddr, j);
			if ((sumkwd&Q.keywords) != Q.keywords) continue;
			UpdateOnePt(vJ, hasChanged, onSameEdge, eKdist, pid);
		}
	}
	else {
		vJ = 0;
		getOptPos(NodeID, NewNodeID, eKdist, vJ);
		UpdateOnePt(vJ, hasChanged, onSameEdge, eKdist, -1);
	}
	//if (bestdist<MAX_DIST) {PrintElapsed();	exit(0);}
	/*if (hasChanged)
	printf("c: %d %d %f\n",NodeID,NewNodeID,eKdist);*/
	return hasChanged;
}

// top-k
// get the possible minimum distance to this edge and return the value
inline float getMinOptPos(int NodeID, int NewNodeID, float eKdist) {
	float dist = 0, mincpdist, ldist, rdist;
	for (int sp = 0; sp<Q.k; sp++) {
		mincpdist = MAX_DIST;
		ldist = rdist = MAX_DIST;
		if (DistMaps[sp].count(NodeID)>0) ldist = DistMaps[sp][NodeID];
		if (DistMaps[sp].count(NewNodeID)>0) rdist = DistMaps[sp][NewNodeID];
		if (ldist<MAX_DIST) {
			if (rdist<MAX_DIST)
				mincpdist = min(ldist, rdist);
			else
				mincpdist = max(0, ldist - eKdist);
		}
		else {
			if (rdist<MAX_DIST)
				mincpdist = max(0, rdist - eKdist);
			else
				mincpdist = 0;
		}
		if (onSameEdge[sp]) mincpdist = 0;
		mincpdist = max(0, mincpdist);

		//if (isWeightUse) mincpdist = mincpdist*weights[sp];
		if (isSumDist)
			dist += mincpdist;
		else
			dist = max(mincpdist, dist);
	}
	return dist;
}


struct RgnType {
	float partial_dist;
	int crash;
};

typedef map<int, RgnType> REGIONMAP;

// memory requirement: r bitmap and r temp diststore
void ConcurrentExpansion(StepQueue& sQ) {
	int MaxHeapSize = 0, AggPopSize = 0;
	REGIONMAP epsNbrList;
	BitStore epsBits;
	int AdjListSize, NewNodeID, AdjGrpAddr, PtGrpKey;
	float eKdist, xVal, yVal;

	// initialize variables
	//if (!isEuclidSpace) return;
	epsBits.assign(NodeNum, false);
	int maxNbrSize = 0;
	while (!sQ.empty()) {
		if (MaxHeapSize<sQ.size()) MaxHeapSize = sQ.size();
		StepEvent event = sQ.top();
		sQ.pop();	AggPopSize++;

		int NodeID = event.node;
		int id = event.ClusID;

		if (DistMaps[id].count(NodeID)>0) continue;	// already found
		DistMaps[id][NodeID] = event.dist;

		float topbound = event.dist;	// change here !
										//if (isWeightUse) topbound = topbound*weights[id];
		if (topbound >= bestdist) continue;	// for both S-dist and M-dist: bound almost correct 

		AdjGrpAddr = getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);	// read # entries
															//getFixedF(XCRD_A, Ref(xVal), AdjGrpAddr);
															//getFixedF(YCRD_A, Ref(yVal), AdjGrpAddr);

		maxNbrSize = max(maxNbrSize, epsNbrList.size());
		float epsilon = (isSumDist) ? (bestdist / Q.k) : (bestdist);
		if (epsilon>topbound) { 	// any bug ?
			if (epsBits[NodeID] == false) {
				epsBits[NodeID] = true;
				epsNbrList[NodeID].crash = 0;
				epsNbrList[NodeID].partial_dist = 0;
			}
		}

		if (epsNbrList.count(NodeID)>0) {
			RgnType& rgt = epsNbrList[NodeID];
			rgt.crash += 1;	// update crash
			bool willErase = (rgt.crash >= Q.k);

			float curcpdist = event.dist;
			//if (isWeightUse) curcpdist = curcpdist*weights[id];
			if (isSumDist) {	// update cur dist and check if need prune
				rgt.partial_dist += curcpdist;
				//				if (!willErase&&!isWeightUse) 
				//					willErase=((Q.k-rgt.crash)*event.dist+rgt.partial_dist>=bestdist);
			}
			else {	// update cur dist and check if need prune
				rgt.partial_dist = max(curcpdist, rgt.partial_dist);
				//				if (!willErase&&!isWeightUse) 
				//					willErase=(rgt.partial_dist>=bestdist);	// current lb
			}

			// not useful at all
			if (!willErase) { // correct if maxedgedist used ??
				float curXdist = rgt.partial_dist;
				for (int s = 0; s<Q.k; s++)
					if (DistMaps[s].count(NodeID) == 0) {
						// float tmpdist = getDist(xVal, yVal, M_xcrd[s], M_ycrd[s]);
						float tmpdist = 0;
						tmpdist = max(tmpdist, event.dist);
						//tmpdist=event.dist;	// i.e., Euclidean bounds not used

						//if (isWeightUse) tmpdist = tmpdist*weights[s];
						if (isSumDist)
							curXdist += tmpdist;
						else
							curXdist = max(tmpdist, curXdist);
					}
				if (curXdist >= bestdist) willErase = true;
			}
			if (willErase) {
				// print info.
				/*for (int s=0;s<Q.k;s++)
				if (DistMaps[s].count(NodeID)==0)
				if (getDist(xVal,yVal,M_xcrd[s],M_ycrd[s])>event.dist)
				printf("%d %f %f\n",NodeID,event.dist,
				getDist(xVal,yVal,M_xcrd[s],M_ycrd[s]));*/
				epsNbrList.erase(NodeID);
			}
		}

		if (bestdist<MAX_DIST) {	// at least one found
			if (epsNbrList.size() == 0) break;	// all pruned => exit
		}

		for (int z = 0; z<AdjListSize; z++) {
			getVarE(ADJNODE_A, Ref(NewNodeID), AdjGrpAddr, z);
			getVarE(DIST_A, Ref(eKdist), AdjGrpAddr, z);
			// getVarE(DIST_A, Ref(eKdist), AdjGrpAddr, z);
			StepEvent newevent = event;	// copy ...
			newevent.node = NewNodeID;
			newevent.ClusID = event.ClusID;
			newevent.dist = event.dist + eKdist;
			if (DistMaps[id].count(NewNodeID) == 0) sQ.push(newevent);	// propagation

																		// if error, add cond. to (... && not on same edge) [sp]
			bool needUpdate = true;
			for (int sp = 0; sp < Q.k; sp++) {
				if (DistMaps[sp].count(NodeID) == 0 || DistMaps[sp].count(NewNodeID) == 0) {
					needUpdate = false;
					break;
				}
				//needUpdate = false;
			}

			if (needUpdate&&isNNquery) {
				float tmpposition;
				float bestVal = getOptPos(NodeID, NewNodeID, eKdist, tmpposition);
				if (bestVal >= bestdist) needUpdate = false;
			}
			if (!needUpdate) continue;

			getVarE(PTKEY_A, Ref(PtGrpKey), AdjGrpAddr, z);
			if (PtGrpKey >= 0) {	// valid PtGrpKey  (-1 for invalid key)
				bool hasChanged = UpdatePoints(NodeID, NewNodeID, PtGrpKey, eKdist);
			}
		}
	}
	printf("Heap: max_size %d, pop_size %d\n", MaxHeapSize, AggPopSize);
	printf("bestdist: %f\n", bestdist);
	if (!sQ.empty()) printf("topdist: %f\n", sQ.top().dist);
	int totalmapsize = 0;
	for (int sp = 0; sp<Q.k; sp++) {
		//printf("%d) td_sz: %d\n",sp,DistMaps[sp].size());
		totalmapsize += DistMaps[sp].size();
	}
	//printf("totalmapsize: %d, maxNbrSz: %d\n", totalmapsize, maxNbrSize);
	rs.distance = INFINITE_MAX;
	rs.poid = -1;

	if (dQ_graph.begin() != dQ_graph.end()) {
		rs.distance = dQ_graph.begin()->dist;
		rs.poid = dQ_graph.begin()->id;
	}
	
	printf("The result of CE distance is: %f\n", rs.distance);
	printf("The result of CE oid is: %d\n", rs.poid);
	printf("The result of CE edgeexpand is: %d\n", nOfEdgeExpended);
}
// compute the distance between current node to the dest with the A* 
float A_star(StepQueue& aQ, BitStore& isVisited, DISTMAP& curdistmap, int dest, int& PopSize, int& VisitedSize) {
	// gdist: from src to cur ; hdist: from cur to dist
	float x_dest, y_dest, x_cur, y_cur;
	int TmpAdjGrpAddr = getAdjListGrpAddr(dest);
	//getFixedF(XCRD_A, Ref(x_dest), TmpAdjGrpAddr);
	//getFixedF(YCRD_A, Ref(y_dest), TmpAdjGrpAddr);

	while (!aQ.empty()) {
		StepEvent event = aQ.top();
		aQ.pop();	PopSize++;

		int NodeID = event.node;
		if (isVisited[NodeID]) continue;
		isVisited[NodeID] = true;		VisitedSize++;
		curdistmap[NodeID] = event.gdist;	// only gdist means the real dist. from src. !!!

		float eKdist;
		int AdjListSize, NewNodeID, AdjGrpAddr;
		AdjGrpAddr = getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);	// read # entries
		for (int z = 0; z<AdjListSize; z++) {
			getVarE(ADJNODE_A, Ref(NewNodeID), AdjGrpAddr, z);
			getVarE(DIST_A, Ref(eKdist), AdjGrpAddr, z);

			StepEvent newevent = event;	// copy ...
			newevent.node = NewNodeID;
			newevent.gdist += eKdist;
			newevent.hdist = 0;

			TmpAdjGrpAddr = getAdjListGrpAddr(NewNodeID);
			//getFixedF(XCRD_A, Ref(x_cur), TmpAdjGrpAddr);
			//getFixedF(YCRD_A, Ref(y_cur), TmpAdjGrpAddr);
			//newevent.xc = x_cur;		
			//newevent.yc = y_cur;

			//if (isEuclidSpace) newevent.hdist = getDist(x_cur, y_cur, x_dest, y_dest);
			newevent.dist = newevent.gdist + newevent.hdist;

			// pathmax equation for non-monotonic heuristic ?
			if (newevent.dist<event.dist) newevent.dist = event.dist;
			if (isVisited[NewNodeID] == false) aQ.push(newevent);	// propagation		
		}
		if (NodeID == dest) return event.dist;
	}
	return MAX_DIST;
}

struct ValuePair {
	float value;
	int id;
};

FastArray<ValuePair> covernodeset;

void Repair_astar(StepQueue& aQ, BitStore& isVisited, int dest) {
	int tsize = 0;
	StepEvent* tmpAry = new StepEvent[aQ.size()];

	float x_dest, y_dest, x_cur, y_cur;
	int TmpAdjGrpAddr = getAdjListGrpAddr(dest);
	//getFixedF(XCRD_A, Ref(x_dest), TmpAdjGrpAddr);
	//getFixedF(YCRD_A, Ref(y_dest), TmpAdjGrpAddr);

	while (!aQ.empty()) {
		StepEvent event = aQ.top();
		aQ.pop();
		if (isVisited[event.node]) continue;

		// assume xc yc fields initialized by the caller !
		event.hdist = 0;	// for Dijkstra
							//if (isEuclidSpace) event.hdist = getDist(x_dest, y_dest, event.xc, event.yc);
		event.dist = event.gdist + event.hdist;
		tmpAry[tsize] = event;
		tsize++;
	}
	for (int i = 0; i<tsize; i++) aQ.push(tmpAry[i]);
	delete[] tmpAry;
}

void TA_EW(StepQueue* raQ, StepQueue& sQ) {
	float eKdist;
	int AdjListSize, NewNodeID, AdjGrpAddr, PtGrpKey;
	int PopSize = 0, VisitedSize = 0, NodeID;
	float cur_x, cur_y, nextdist;
	int curround, totround = 0;	// # of sorted access
	unsigned long long sumkwd;
	float lb_dist = 0;
	while (!sQ.empty()) {
		StepEvent event = sQ.top();
		sQ.pop();
		NodeID = event.node;
		curround = event.ClusID;
		if (saVisited[curround][NodeID]) continue;
		saVisited[curround][NodeID] = true;	// !!!

		AdjGrpAddr = getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);	// read # entries
		for (int z = 0; z<AdjListSize; z++) {
			getVarE(ADJNODE_A, Ref(NewNodeID), AdjGrpAddr, z);
			getVarE(DIST_A, Ref(eKdist), AdjGrpAddr, z);
			// *** added part
			//getFixedF(XCRD_A, Ref(cur_x), AdjGrpAddr);
			//getFixedF(YCRD_A, Ref(cur_y), AdjGrpAddr);
			StepEvent newevent = event;	// copy ...
			newevent.node = NewNodeID;
			newevent.dist = event.dist + eKdist;
			// *** added part
			//newevent.xc = cur_x;
			//newevent.yc = cur_y;

			if (saVisited[curround][NewNodeID] == false) sQ.push(newevent);	// propagation		
		}
		nextdist = event.dist;

		// get next of curround 
		DistMaps[curround][NodeID] = nextdist;	// assume valid result
		lb_dist = nextdist;

		float mindist = 0, curnetdist = 0;
		// *** added part
		//cur_x = event.xc;
		//cur_y = event.yc;

		for (int s = 0; s<Q.k; s++) {
			//float cpdist = getDist(cur_x, cur_y, M_xcrd[s], M_xcrd[s]);
			float cpdist = 0;
			if (DistMaps[s].count(NodeID)>0) cpdist = DistMaps[s][NodeID];

			//if (isWeightUse) cpdist = cpdist*weights[s];
			if (isSumDist) {
				mindist += lb_dist;
				curnetdist += cpdist;
			}
			else {
				mindist = max(mindist, lb_dist);
				curnetdist = max(curnetdist, cpdist);
			}
		}
		/*
		if (totround%PrintLimit == PrintLimit - 1) {
			printf("%d: %f %f\n", totround, mindist, bestdist);
			PrintElapsed();
		}
		*/
		totround++;
		if (mindist >= bestdist) break;

		bool mustPrune = true;
		for (int z = 0; z<AdjListSize; z++) {
			getVarE(ADJNODE_A, Ref(NewNodeID), AdjGrpAddr, z);
			getVarE(DIST_A, Ref(eKdist), AdjGrpAddr, z);
			getVarE(PTKEY_A, Ref(PtGrpKey), AdjGrpAddr, z);
			getVarE(SUMKWD_A, Ref(sumkwd), AdjGrpAddr, z);
			//if ((sumkwd&Q.keywords) != Q.keywords) mustPrune = true;
			// L0^- filter
			if (isSumDist) {
				if (curnetdist<bestdist + Q.k*eKdist) mustPrune = false;
				break;
			}
			else {
				if (curnetdist<bestdist + eKdist) mustPrune = false;
				break;
			}
		}
		if (mustPrune) continue;

		AdjGrpAddr = getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);	// read # entries
		for (int s = 0; s<Q.k; s++) {
			DISTMAP& curDistMap = DistMaps[s];
			if (curDistMap.count(NodeID) == 0) {
				Repair_astar(raQ[s], raVisited[s], NodeID);
				A_star(raQ[s], raVisited[s], curDistMap, NodeID, PopSize, VisitedSize);
			}
		}
		for (int z = 0; z<AdjListSize; z++) {
			getVarE(ADJNODE_A, Ref(NewNodeID), AdjGrpAddr, z);
			getVarE(DIST_A, Ref(eKdist), AdjGrpAddr, z);
			getVarE(PTKEY_A, Ref(PtGrpKey), AdjGrpAddr, z);
			getVarE(SUMKWD_A, Ref(sumkwd), AdjGrpAddr, z);
			if ((sumkwd&Q.keywords) != Q.keywords) continue;
			// L0 filter
			if (isSumDist) {
				if (curnetdist >= bestdist + Q.k*eKdist) continue;		// cannot have better distance
			}
			else {
				if (curnetdist >= bestdist + eKdist) continue;		// cannot have better distance
			}

			// L1 filter
			if (getMinOptPos(NodeID, NewNodeID, eKdist) >= bestdist) continue;
			if (PtGrpKey >= 0) {	// valid PtGrpKey  (-1 for invalid key)
				for (int s = 0; s<Q.k; s++) {
					DISTMAP& curDistMap = DistMaps[s];
					if (curDistMap.count(NewNodeID) == 0) {
						Repair_astar(raQ[s], raVisited[s], NewNodeID);
						A_star(raQ[s], raVisited[s], curDistMap, NewNodeID, PopSize, VisitedSize);
					}
				}
				bool needUpdate = true;
				if (needUpdate) {
					float tmppossition;
					float bestVal = getOptPos(NodeID, NewNodeID, eKdist, tmppossition);
					if (bestVal >= bestdist) needUpdate = false;
				}
				if (needUpdate) UpdatePoints(NodeID, NewNodeID, PtGrpKey, eKdist);
			}
		}
	}
	//printf("bestdist: %f\n", bestdist);
	//printf("totround: %d, pop size: %d, VisitedSize: %d\n", totround, PopSize, VisitedSize);
	
	rs.distance = INFINITE_MAX;
	rs.poid = -1;

	if (dQ_graph.begin() != dQ_graph.end()) {
		rs.distance = dQ_graph.begin()->dist;
		rs.poid = dQ_graph.begin()->id;
	}
	printf("The result of TA distance is: %f\n", rs.distance);
	printf("The result of TA oid is: %d\n", rs.poid);
	printf("The result of TA edgeexpand is: %d\n", nOfEdgeExpended);
}


inline void printSetting() {
	if (isNNquery) printf(",point"); else printf(",center");
	if (isSumDist) printf(",sum"); else printf(",max");
}


// Only find the necessary dist, not necessarily all pairs dist
void FindSoln(int curAlg) {	// "node" reused as the next nodeID instead !!!
	StepQueue sQ;
	StepEvent stepL, stepR;
	StepQueue *raQ, *saQ;
	raQ = new StepQueue[Q.k];
	saQ = new StepQueue[Q.k];

	StepQueue *mtQ;
	mtQ = new StepQueue[Q.k];

	// initialization
	bestdist = MAX_DIST;
	InitDQ(dQ_graph);	// clear kNN-dist heap

	onSameEdge.assign(Q.k, false);
	partdists = new float*[Q.k];
	for (int i = 0; i<Q.k; i++) partdists[i] = new float[2];
	compdist.assign(Q.k, MAX_DIST);
	bestcompdist.assign(Q.k, MAX_DIST);

	float vP, eDist, eKdist;
	int AdjListSize, Ni, Nj, NewNodeID, TmpAdjGrpAddr;

	while (!sQ.empty()) sQ.pop();
	for (int s = 0; s<Q.k; s++) {
		vP = Q.querypts[s].dist_Ni;
		Ni = Q.querypts[s].Ni;
		Nj = Q.querypts[s].Nj;

		// for greedy search and rapid elimination
		TmpAdjGrpAddr = getAdjListGrpAddr(Ni);
		//getFixedF(XCRD_A, Ref(x1), TmpAdjGrpAddr);
		//getFixedF(YCRD_A, Ref(y1), TmpAdjGrpAddr);

		// lookup eDist from adj. list. of Ni;
		getFixedF(SIZE_A, Ref(AdjListSize), TmpAdjGrpAddr);	// read # entries
		for (int z = 0; z<AdjListSize; z++) {
			getVarE(ADJNODE_A, Ref(NewNodeID), TmpAdjGrpAddr, z);
			getVarE(DIST_A, Ref(eKdist), TmpAdjGrpAddr, z);
			if (NewNodeID == Nj) eDist = eKdist;
		}

		stepL.ClusID = stepR.ClusID = s;
		stepL.dist = vP;			stepL.node = Ni;
		stepR.dist = eDist - vP;	stepR.node = Nj;

		stepL.isSortAcc = true;	stepR.isSortAcc = true;	// ??? 4/3/2004 19:00
		sQ.push(stepL);			sQ.push(stepR);

		stepL.gdist = vP;			stepL.hdist = 0;
		stepR.gdist = eDist - vP;	stepR.hdist = 0;
		stepL.isSortAcc = false;	stepR.isSortAcc = false;	// ??? 4/3/2004 19:00
		mtQ[s].push(stepL);		mtQ[s].push(stepR);
	}


	for (int sp = 0; sp<Q.k; sp++)
		DistMaps[sp].clear(); // clear first (safe before use)


	saVisited = new BitStore[Q.k];
	raVisited = new BitStore[Q.k];
	for (int i = 0; i<Q.k; i++) {
		saVisited[i].assign(NodeNum, false);
		raVisited[i].assign(NodeNum, false);
	}

	if (curAlg == CE) {
		printf("\n[CE");
		printSetting();
		printf("]\n");
		ConcurrentExpansion(sQ);
		//RR_Expand(mtQ);
	}

	// later, change the ... for cache
	if (curAlg == TA) {
		printf("\n[TA");
		printSetting();
		printf("]\n");
		astar_Visits = new BitStore[Q.k];
		for (int i = 0; i<Q.k; i++)
			astar_Visits[i].assign(NodeNum, false);
		TA_EW(raQ, sQ);
	}
}


void queryAlgorithm(string fileprefix) {
	
	clock_t start, finish;
	if (algorithmId == EGTD) {
		//initialQuery(fileprefix);
		OpenDiskComm(basepath, 200);
		initialQuery(fileprefix);
		start = clock();
		EGTDA();
		finish = clock();
		cout << "Time spend of EGTDA is:" << (finish - start) << endl << endl;
		CloseDiskComm();
	}
	else if (algorithmId == EGETD) {
		//initialQuery(fileprefix);
		OpenDiskComm(basepath, 200);
		initialQuery(fileprefix);
		start = clock();
		EGETDA();
		finish = clock();
		//cout << "Time spend of EGETDA is:" << (finish - start) / CLOCKS_PER_SEC << endl;
		cout << "Time spend of EGETDA is:" << (finish - start) << endl << endl;
		CloseDiskComm();
	}
	else if (algorithmId == TA) {
		OpenDiskComm(basepath, 200);
		initialQuery(fileprefix);

		start = clock();
		if (isNNquery) NNnum = 1;	// force to find 1 center only !
		DistMaps = new DISTMAP[Q.k];
		RefreshCache();
		FindSoln(TA);
		finish = clock();
		//cout << "Time spend of EGETDA is:" << (finish - start) / CLOCKS_PER_SEC << endl;
		cout << "Time spend of TA is:" << (finish - start) << endl <<endl;
		CloseDiskComm();
	}
	else {
		OpenDiskComm(basepath, 200);
		initialQuery(fileprefix);

		start = clock();
		if (isNNquery) NNnum = 1;	// force to find 1 center only !
		DistMaps = new DISTMAP[Q.k];
		RefreshCache();
		FindSoln(CE);
		finish = clock();
		//cout << "Time spend of EGETDA is:" << (finish - start) / CLOCKS_PER_SEC << endl;
		cout << "Time spend of CE is:" << (finish - start) << endl << endl;
		CloseDiskComm();
	}
}

//************* end of TKDE2005
int main(int argc, char *argv[]) {

	string dataFilename = basepath;
	
	// set the query:read info from query file, 
	// format: keywords k, nodei nodej lenght dis, ...
	string inFile = dataFilename + "\\queryfile";

	// generate index
	/*
	roadParameter rp;
	rp.avgNKwd = 7;
	rp.avgNPt = 7;
	cout << 1 << endl;
	mainGenData(dataFilename, rp);
	printf("The end of generate data.\n");
	//*/

	/*
	int size; int id;
	float dis;
	unsigned long long ull;
	getFixedF(SIZE_A, Ref(size), 1148);
	printf("Size is:%d\n", size);

	//ADJNODE_A, DIST_A, SUMKWD_A, PTKEY_A
	getVarE(SUMKWD_A, Ref(ull), 1148, 0);
	unsigned long long ull2 = ull; 
	//printf("ull is:%I64d\n", ull);

	getVarE(ADJNODE_A, Ref(id), 1148, 0);
	printf("node id is:%d\n", id);
	*/
	// perform the query	
	///*
	geneQuery(dataFilename);
	
	ifstream iFile(inFile.c_str());
	while (!iFile.eof())
	{
		string line;
		getline(iFile, line);
		if (line == "")
			continue;
		istringstream iss(line);
		
		char c;
		iss >> Q.keywords >> Q.k;
		Q.querypts.clear();
		queryedge.clear();

		for (int i = 0; i < Q.k; i++) {
			QueryPoint q;
			iss >> c >> q.Ni >> q.Nj >> q.distEdge >> q.dist_Ni;
			__int64 ni, nj;
			ni = q.Ni;
			nj = q.Nj;

			__int64 key = getKey(ni, nj);
			queryedge.push_back(key);
			Q.querypts.push_back(q);
		}

		// start query
		for (int i = 1; i < 5; i++) {
			algorithmId = i;
			queryAlgorithm(dataFilename);
		}
		//queryAlgorithm(dataFilename);
	}
	
	
	//printf("The result distance is: %f\n", rs.distance);
	//printf("The result oid is: %d\n", rs.poid);
	//*/
	
	return 0;
}