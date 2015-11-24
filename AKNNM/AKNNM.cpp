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
#pragma comment(lib,"ws2_32.lib")

using namespace std;

#define EGTD 1
#define EGETD 2
#define INFINITE_MAX 100000000.0
#define MAX_KEYWORDN 64 
// define the query path
string basepath = "..\\dataset";

int nOfPageAccessed = 0;
int nOfEdgeExpended = 0;
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
	int addr = paID[NodeID].addr;
	return paID[NodeID].addr;
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
		printf("i:%d, nodeid:%d, part:%d, addr:%d\n", loop, nid, pa.part, pa.addr);
	}
	nOfNode = partID.size();
	paFile.close();
	return nOfNode;
}

void comQueryPath(Query Q, vector<TreeNode> &EGT) {
	// ����Q�е�ÿ����ѯ�㣬�������Ĳ�ѯ·�������浽queryPath��
	int i = 0;
	vector<QueryPoint> ::iterator itr = Q.querypts.begin();
	for (; itr != Q.querypts.end(); itr++) {
		QueryPoint q = *itr;
		int pid = paID[q.Ni].part;
		// ����Ҷ�ӽڵ㵽��ѯ��ľ���
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


//����q��vertex�ľ���
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

void EGTDA( ) {

	int AdjGrpAddr, AdjListSize, NewNodeID, PtGrpKey, PtNumOnEdge;
	unsigned long long sumkwd;
	float EdgeDist, PtDist;

	// ��ʼ�����
	rs.distance = INFINITE_MAX;
	rs.poid = -1;
	edgevisited.clear();
	// step1 �����ڵ�ѹ�뵽���ȶ�����
	pnode pni;
	pni.dist = 0;
	pni.nodeID = 0;
	pq.push(pni);

	// step2 whileѭ����������ǰ���ŵĽڵ㲢չ��
	while (!pq.empty()) {
		// step1 ����top�ڵ�
		pnode node = pq.top();
		pq.pop();
		// if nû������ Case1��Ҷ�ӽڵ㣻 Case2���м�ڵ�
		int nid = node.nodeID;
		
		printf("nodeid is:%d\n", nid);
		if ((node.dist < rs.distance) && (EGT[nid].unionkwds&Q.keywords) == Q.keywords) {
			// Case1: Ҷ�ӽڵ�
			if (EGT[nid].isleaf) {
				// �����������ŵĵ㲢����rs.distance���ڲ����ⲿ����������ڲ���Dijkstra�ⲿ�ñ߽���룩
				// ����ÿ����ѯ�㵽�ڲ�����vertex�ľ���
				vector<pvertex> pv;
				for (int i = 0; i < Q.k; i++) {
					// ���������case1����ǰ�ò�ѯ���ڸýڵ��У�����Dijkstra;case2:���ڣ�ֱ�����
					printf("The nid is:%d, q num is:%d\n", nid, i);
					map<int, vector<float>>::iterator it = queryPath[i].find(nid);
					vector<float> result;
					if (it != queryPath[i].end()) { // �����棬����Dijkstra�㷨
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

					else { // ֱ�����
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
								// ���������⣬vdistδ��ʼ��
								pt.vdist[i] = min; 
								pt.sumDist = min;
								pv.push_back(pt);
							}
							else {
								// ���㵽q-border-vertex�ľ���
								int dist;
								float min = INFINITE_MAX;
								for (int m = 0; m < EGT[nid].borders.size(); m++) {
									int pos = m*EGT[nid].leafnodes.size() + l;
									dist = EGT[nid].mind[pos] + node.bdist[m][i];
									if (dist < min) min = dist;
								}
								// ���������⣬vdistδ��ʼ��
								pv[l].vdist[i] = min; 
								pv[l].sumDist += min;
							}
						}
					}
				}

				// ��һ���Ż�������************
				// ����ÿһ��vertex���ڵıߣ��������Ƿ�������������POI 
				for (int t = 0; t < pv.size(); t++) {
					printf("loop %d to num %d\n", t, pv.size());
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
							if (edgevisited[edge].isVisited == 0) { //����һ��
								getVarE(PTKEY_A, Ref(PtGrpKey), AdjGrpAddr, n);
								if (PtGrpKey == -1) {
									cout<<"No POI existed on Edge where Q located."<<endl;
								}
								else {
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
													sum += dis;
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
							// �������㶼���ʹ��������ˣ�����ɾ����
							edgevisited.erase(it);
						}
						else { // ������δ���ʹ�
							
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
			// Case2: �м�ڵ�
			else {
				// ��չ�ýڵ㣬��������������(�ؼ��ֺ;���)�ĺ��ӽڵ�ľ��벢ѹ�����
				// ����ÿһ�����ӽڵ�
				vector<pnode> childnode;
				childnode.clear();

				for (int i = 0; i < Q.k; i++) { 
					// ���ղ�ѯ·������������ÿ����ѯ��
					for (int j = 0; j < EGT[nid].children.size(); j++) {
						// two cases. 1: q in the node; 2: not in the node
						int cid = EGT[nid].children[j];
											
						if (i == 0) {
							pnode pn;
							pn.nodeID = cid;
							pn.dist = 0;
							for (int f = 0; f < EGT[cid].borders.size(); f++) {
								vector<float> vf(Q.k, 0);
								pn.bdist.push_back(vf);
							}
							childnode.push_back(pn);
						}
			
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
								// ����ÿ����ѯ�㵽����border�ľ���
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
					pq.push(childnode[loop]);
				}
			}
		}

	}

	// step3 ����half
	map<__int64, edgeStatu> ::iterator it = edgevisited.begin();
	for (; it != edgevisited.end(); it++) {
		edgeStatu es = it->second;
		// �����ڱ���
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
		// �����Ż�������***********
		AdjGrpAddr = getAdjListGrpAddr(vid);
		getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);
		for (int n = 0; n < AdjListSize; n++) {
			getVarE(ADJNODE_A, Ref(NewNodeID), AdjGrpAddr, n);
			if (NewNodeID != vjd) continue;
			getVarE(DIST_A, Ref(EdgeDist), AdjGrpAddr, n);
			//��¼��ѯ�߾���
			getVarE(PTKEY_A, Ref(PtGrpKey), AdjGrpAddr, n);
			if (PtGrpKey == -1) {
				//cout<<"No POI existed on Edge where Q located."<<endl;
			}
			else {
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
}


void EGETDA() {
	// �Ľ�1������cover Keyword�� 2�����ö˵��ֵ�����˱��ϵ�ֵ�������ڱ���һ���ڵ��ڲ���ʱ���й���
}

void queryAlgorithm(string fileprefix) {
	initialQuery(fileprefix);
	clock_t start, finish;

	start = clock();
	EGTDA();
	finish = clock();
	cout << "EGTDA:" << (finish - start) / CLOCKS_PER_SEC << endl;
	/*
	start = clock();
	EGETDA();
	finish = clock();
	cout << "EGETDA:" << (finish - start) / CLOCKS_PER_SEC << endl;
	*/
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
	//*/

	// generate query
	//geneQuery(dataFilename);
	///*
	// perform the query
	geneQuery(dataFilename);
	OpenDiskComm(dataFilename, 200);
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
	///*
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
		queryAlgorithm(dataFilename);
	}
	
	CloseDiskComm();
	//*/
	printf("The result distance is: %f\n", rs.distance);
	printf("The result oid is: %d\n", rs.poid);
	return 0;
}