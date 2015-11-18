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


// define the query path
string basepath = "..\\dataset\\";

int nOfPageAccessed = 0;
int nOfEdgeExpended = 0;
int algorithmId;


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
map<int, edgeStatu> edgevisited; // edgeID, isVisited
vector<int> queryedge;

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
	return paID[NodeID].addr;
}

void partAddrLoad(const char* filename, map<int, PartAddr> &partID) {

	printf("loading partAddrFile\n");
	FILE *paFile;
	paFile = fopen(filename, "r+");
	CheckFile(paFile, filename);
	int nOfNode, nodeID;
	fread(&nOfNode, sizeof(int), 1, paFile);
	for (int i = 0; i<nOfNode; i++) {
		PartAddr pa;
		fread(&nodeID, sizeof(int), 1, paFile);
		fread(&pa.part, sizeof(int), 1, paFile);
		fread(&pa.addr, sizeof(int), 1, paFile);
		partID[nodeID] = pa;
	}
	fclose(paFile);
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

		queryPath[i] = nodePath;
		i++;
	}
}


void initialQuery(string fileprefix) {
	nOfPageAccessed = 0;
	nOfEdgeExpended = 0;
	char tmpFileName[255];
	// load egtree
	sprintf(tmpFileName, "%s\\egtree", fileprefix);
	egtree_load(tmpFileName, EGT);
	// load partAddr
	sprintf(tmpFileName, "%s\\partAddr", fileprefix);
	partAddrLoad(tmpFileName, paID);
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

void EGTDA( ) {

	int AdjGrpAddr, AdjListSize, NewNodeID, PtGrpKey, PtNumOnEdge;
	unsigned long long sumkwd;
	float EdgeDist, PtDist;

	// 初始化结果
	rs.distance = INFINITE_MAX;
	rs.poid = -1;

	// step1 将根节点压入到优先队列中
	pnode pn;
	pn.dist = 0;
	pn.nodeID = 0;
	pq.push(pn);

	// step2 while循环，弹出当前最优的节点并展开
	while (!pq.empty()) {
		// step1 弹出top节点
		pnode node = pq.top();
		pq.pop();
		// if n没被过滤 Case1：叶子节点； Case2：中间节点
		int nid = node.nodeID;
		if ((node.dist < rs.distance) && (EGT[nid].unionkwds&Q.keywords) == Q.keywords) {
			// Case1: 叶子节点
			if (EGT[nid].isleaf) {
				// 计算其中最优的点并更新rs.distance（内部和外部两种情况，内部用Dijkstra外部用边界距离）
				// 计算每个查询点到内部所有vertex的距离
				vector<pvertex> pv;
				for (int i = 0; i < Q.k; i++) {
					// 两种情况，case1：当前该查询点在该节点中，调用Dijkstra;case2:不在，直接相加
					map<int, vector<float>>::iterator it = queryPath[i].find(nid);
					vector<float> result;
					if (it != queryPath[i].end()) { // 调用Dijkstra算法
						// note the corresponding relationship between result and leafnodes
						result = dijkstra_candidate(Q.querypts[i], EGT[nid].leafnodes);

						for (int l = 0; l < EGT[nid].leafnodes.size(); l++) {
							if (i == 0) {
								pvertex pt;
								pt.vertexID = EGT[nid].leafnodes[l];

								pt.vdist[i] = result[l];
								pt.sumDist = result[l];
								pv[l] = pt;
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
								pv[l] = pt;
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
					pvertex pt = pv[t];
					int vid = pt.vertexID;

					AdjGrpAddr = getAdjListGrpAddr(vid);
					getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);
					for (int n = 0; n < AdjListSize; n++) {
						getVarE(ADJNODE_A, Ref(NewNodeID), AdjGrpAddr, n);
						int edge = getKey(vid, NewNodeID);
						// bool isonedge;
						
						map<int, edgeStatu>::iterator it = edgevisited.find(edge);
						edgeStatu es;
						if (it != edgevisited.end()) {
							if (edgevisited[edge].isVisited == 0) { //处理一下
								getVarE(PTKEY_A, Ref(PtGrpKey), AdjGrpAddr, n);
								if (PtGrpKey == -1) {
									//cout<<"No POI existed on Edge where Q located."<<endl;
								}
								else {
									getVarE(DIST_A, Ref(EdgeDist), AdjGrpAddr, n);
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
												if (find(queryedge.begin(), queryedge.end(), edge) != queryedge.end()) { // on edge
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
							getVarE(SUMKWD_A, &sumkwd, AdjGrpAddr, n);
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
				for (int i = 0; i < Q.k; i++) { 
					// 按照查询路径来处理，处理每个查询点
					for (int j = 0; j < EGT[nid].children.size(); j++) {
						// two cases. 1: q in the node; 2: not in the node
						int cid = EGT[nid].children[j];
						pnode pn;
						pn.nodeID = cid;
						pn.dist = 0;

						map<int, vector<float>>::iterator it = queryPath[i].find(cid);
						if (it != queryPath[i].end()) { // do nothing
							// wrong!!!!!!!!! vector不能这么赋值
							for (int t = 0; t < EGT[cid].borders.size(); t++) {
								pn.bdist[t].push_back(0);
							}
							//pn.bdist[j] = 0;
						}
						else { // compute the distance
							   // first determine the father node in the layer of cid
							int sid;
							for (int l = 0; l < EGT[nid].children.size(); l++) {
								sid = EGT[nid].children[l];
								map<int, vector<float>>::iterator sit = queryPath[i].find(sid);
								if (sit != queryPath[i].end()) { // return this id
									break;
								}
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
									map<int, vector<float>>::iterator nodeDist = queryPath[i].find(sid);
									dist = dist + nodeDist->second[t];
									if (dist < min) min = dist;
								}
								if (min < mind) mind = min;
								// if is the first then initiate it
								if (i == 0) {
									vector<float> qTb;
									// wrong!!!!!!!!!
									qTb.push_back(min);
									pn.bdist.push_back(qTb);
								}
								else {
									// wrong!!!!!!!!!
									pn.bdist[k].push_back(min);
								}

							}
							pn.dist = pn.dist + mind;
							pq.push(pn);
						}
					}
				}

			}
		}

	}

	// step3 处理half
	map<int, edgeStatu> ::iterator it = edgevisited.begin();
	for (; it != edgevisited.end(); it++) {
		edgeStatu es = it->second;
		if (es.sum < rs.distance) {
			int vid, vjd;
			int first = it->first;
			breakKey(first, vid, vjd);
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

							float sum = es.sum;
							if (vid < vjd) {
								float dis = PtDist;
								sum = sum + Q.k*dis;
							}
							else {
								float dis = EdgeDist - PtDist;
								sum = sum + Q.k*dis;
							}

							if (sum < rs.distance) {
								rs.distance = sum;
								rs.poid = pid;
							}
							else {
								// 后面的距离更大
								break;
							}
						}
					}
				}
			}
		}

	}
}


void EGETDA() {
	// 改进1：利用cover Keyword； 2：利用端点的值来过滤边上的值，可以在遍历一个节点内部点时进行过滤
}

void queryAlgorithm(string fileprefix) {
	initialQuery(fileprefix);
	clock_t start, finish;

	start = clock();
	EGTDA();
	finish = clock();
	cout << "EGTDA:" << (finish - start) / CLOCKS_PER_SEC << endl;

	start = clock();
	EGETDA();
	finish = clock();
	cout << "EGETDA:" << (finish - start) / CLOCKS_PER_SEC << endl;
}

int k = 30;
int qk = 6;
int np = 7;
int nOfTest = 20;
struct querycand {
	int nodei;
	int nodej;
	float length;
	bool visited;
};
map<int, querycand>  qc;

void geneQuery(string prxfilename) {
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
		bitset<MAX_KEYWORDS> keywords_set;
		keywords_set.reset();
		while (keywords_set.count() != qk) {
			int position = rand() % MAX_KEYWORDS;
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
	/*
	// set the query:read info from query file, 
	// format: keywords k, nodei nodej lenght dis, ...
	string inFile = dataFilename + "\\queryfile";

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
			int key = getKey(q.Ni, q.Nj);
			queryedge[i] = key;
			Q.querypts.push_back(q);
		}

		// start query
		queryAlgorithm(dataFilename);
	}
	*/
	roadParameter rp;
	rp.avgNKwd = 7;
	rp.avgNPt = 7;
	mainGenData(dataFilename, rp);

	return 0;
}