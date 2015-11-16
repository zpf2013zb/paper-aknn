#include <iostream>
#include <unordered_map>
#include <queue>
#include <bitset>
#include <algorithm>
#include "diskbased.h"
#include "ConfigType.h"
#include "KeywordsGenerator.h"
#include "egtree.h"
#include <fstream>

using namespace std;
// all paremeters
#define CandidateSet vector<POI>
#define ResultSet vector<int>
#define EA 1
#define EGBU 2
#define EGTD 3
#define INFINITE_MAX 100000000.0
#define visited 4
#define hVisited 5
#define unVisited 6
//#define halfVisited vector<int>
int nOfDomiTest;
int nOfEdgeExpd;
int bottomToUp=1;
float queryEdgeDist;
CandidateSet cS;
ResultSet rS;

// for dijkstra
//halfVQueue hvQ;
map<edgePair, edgeState, eSComparison> visitedState ;
map<int,float> distTQ;
dVQueue dvq;

// for EGBU
struct PartAddr {
	int part;
	int addr;
};
map<int,PartAddr> paID;
vector<TreeNode> EGT;
set<int> Visited; //record the  id of visited treeNodeID

struct TreeNodeC
{
    bool operator () (const TreeNode& left, const TreeNode& right) const
    {
		return left.minDistTQ > right.minDistTQ;
    }
};

typedef	priority_queue<TreeNode,vector<TreeNode>,TreeNodeC> TreeNodeQ;
TreeNodeQ tnq;
vector<int> pathID;


// to verify whether R contains by L
bool LcontianRIKwd(const set<int> L,const set<int> R)
{
    //High bit set to zero
    //R&=tmp;
	
	if(L.size() < R.size()) return false;
	set<int>::iterator it;
	for(it = R.begin(); it != R.end(); it++) {
		if (find(L.begin(), L.end(), (*it)) != L.end()) {
			
		} else {
			return false;
		}
	}
    return  true;
	
	//set<int> temp;
	//set_intersection(L.begin(), L.end(), R.begin(), R.end(), temp.begin());
}
// smaller is better
bool isBeDominate(POI poi, CandidateSet cS, const QueryPoint &Q) {
	for(int i=0; i<cS.size(); ) {
		POI temp = cS[i];
		nOfDomiTest++;
		bool pDt = true;
		bool tDp = true;
		int pEt = 0;
		for(int j=0; j<ATTRIBUTE_DIMENSION; j++) {
			if(Q.subSpace[j]==0) {
			} else {
				if(poi.attr[j]<temp.attr[j]) {
					tDp = false;
				} else if(poi.attr[j] == temp.attr[j]) {
					pEt++;
				} else {
					pDt = false;
				}
			}
		}
		// poi dominate temp
		if(pDt&&pEt!=Q.subSpace.count()) {
			cS.erase(cS.begin()+i);
		}
		if(tDp&&pEt!=Q.subSpace.count()) {
			return true;
		}
		if(pEt == Q.subSpace.count()) {
			i++;
		}
		//nOfDomiTest++;
	}
	return false;
}

bool isBeBoundDominate(float attrBound[ATTRIBUTE_DIMENSION][2], CandidateSet cS, const QueryPoint &Q) {
	for(int i=0; i<cS.size(); ) {
		POI temp = cS[i];
		nOfDomiTest++;
		bool bDt = true;
		bool tDb = true;
		int pEt = 0;
		for(int j=0; j<ATTRIBUTE_DIMENSION; j++) {
			if(Q.subSpace[j]==0) {
			} else {
				if(attrBound[j][0]<temp.attr[j]) {
					tDb = false;
				} else if(attrBound[j][0] == temp.attr[j]) {
					pEt++;
				} else {
					bDt = false;
				}
			}
		}
		// poi dominate temp
		if(tDb&&pEt!=Q.subSpace.count()) {
			return true;
		}
		if(pEt == Q.subSpace.count()) {
			i++;
		}
		//nOfDomiTest++;
	}
	return false;
}


void initialQuery(int algorithmId, FILE *ptAddr) {
	nOfDomiTest = 0;
	nOfEdgeExpd = 0;
	int nOfNode;
	int nodeID;
	if(algorithmId == EA) { // EA

	} else if (algorithmId == EGBU){ // EGBU
		// load partAddr
		fread( &nOfNode, sizeof(int), 1, ptAddr );
		for(int i=0; i<nOfNode; i++) {
			PartAddr pa;
			fread( &nodeID, sizeof(int), 1, ptAddr );
			fread( &pa.part, sizeof(int), 1, ptAddr );
			fread( &pa.addr, sizeof(int), 1, ptAddr );
			paID[nodeID] = pa;
		}
		// load egtree
		egtree_load(EGT);

	} else { // EGTD
		// load partAddr
		fread( &nOfNode, sizeof(int), 1, ptAddr );
		for(int i=0; i<nOfNode; i++) {
			PartAddr pa;
			fread( &nodeID, sizeof(int), 1, ptAddr );
			fread( &pa.part, sizeof(int), 1, ptAddr );
			fread( &pa.addr, sizeof(int), 1, ptAddr );
			paID[nodeID] = pa;
		}
		// load egtree
		egtree_load(EGT);
	}
}

void printResult() {
	for(int i=0; i<rS.size(); i++) {
		cout<<rS[i]<<" ";
	}
}

void EABasedAlgorithm(const QueryPoint &Q) {
	DStepQueue sQ;
    int poicnt=0;
    int roadnodecnt=0;

    //Data
    int Ni=Q.Ni,Nj=Q.Nj;//Ni must less than Nj
    int AdjGrpAddr,AdjListSize,NewNodeID,PtGrpKey,PtNumOnEdge;
    //unsigned long long keywords;
	set<int> tempKwds;
	int *tempKwd;
    float EdgeDist,PtDist;
    //Precess sub edge of Q.Ni,Q.Nj divided by Q
    AdjGrpAddr=getAdjListGrpAddr(Ni);
    getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
    for (int i=0; i<AdjListSize; i++){
        getVarE(ADJNODE_A,Ref(NewNodeID),AdjGrpAddr,i);
        //getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
        if(NewNodeID == Nj){
            getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
            getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);

            nOfEdgeExpd++;

            //cout<<"("<<Ni<<","<<Nj<<") EdgeDist:"<<EdgeDist;
            if(PtGrpKey==-1){
                //cout<<"No POI existed on Edge where Q located."<<endl;
            }
            else {
                getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
                //cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
                //Notice the order POI visited on edge
                for(int j=0; j<PtNumOnEdge; j++){
                    poicnt++;
                    getVarE(PT_DIST,Ref(PtDist),PtGrpKey,j);

                    getVarE(PT_VCT,tempKwd,PtGrpKey,j);
					int num = sizeof(tempKwd)/sizeof(int);
					for(int loop=0; loop<sizeof(tempKwd); loop+sizeof(int)){
						int temp;
						memcpy(&temp,tempKwd+loop,sizeof(int));
						tempKwds.insert(temp);
					}

					if(LcontianRIKwd(tempKwds,Q.kwd)){
                        POI tmp;
						tmp.dist_toquery = fabs(PtDist-Q.dist_Ni);
						getVarE(PT_ATTRIBUTE,tmp.attr,PtGrpKey,j);
						getVarE(PT_P,&tmp.poid,PtGrpKey,j);
						
						if(tmp.dist_toquery <= Q.distCnst) {
							if(isBeDominate(tmp, cS, Q)) {
							} else {
								cS.push_back(tmp);
							}
						}             
                    }// endif
                }
                
            }// endelse
        }
    }
	// record the visited state of edge
    edgePair epi,epj;
	edgeState esi,esj;
	epi.Ni = 0;
	epi.Nj = Q.Ni;
	esi.vState = visited;
	esi.iDisToQuery = 0;
	esi.jDisToQuery = Q.dist_Ni;

	epj.Ni = 0;
	epj.Nj = Q.Nj;
	esj.iDisToQuery = 0;
	esj.jDisToQuery = EdgeDist-Q.dist_Ni;

	visitedState[epi] = esi;
	visitedState[epj] = esj;


	// record the visited node information
	dijkVisit dvi,dvj;
	dvi.N = Q.Ni;
	dvi.disTQ = Q.dist_Ni;

	dvj.N = Q.Nj;
	dvj.disTQ = EdgeDist-Q.dist_Ni;
	
	dvq.push(dvi);
	dvq.push(dvj);

	while(!dvq.empty()) {
		dijkVisit tmpDV = dvq.top();
		dvq.pop();
		if(tmpDV.disTQ <= Q.distCnst) {
			AdjGrpAddr=getAdjListGrpAddr(tmpDV.N);
			getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
			for (int i=0; i<AdjListSize; i++){
				getVarE(ADJNODE_A,Ref(NewNodeID),AdjGrpAddr,i);
				getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
				getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);
				edgePair ep;
				ep.Ni = NewNodeID<tmpDV.N?NewNodeID:tmpDV.N;
				ep.Nj = NewNodeID>tmpDV.N?NewNodeID:tmpDV.N;
				map<edgePair, edgeState, eSComparison>::iterator iter;
				iter = visitedState.find(ep);
				if(iter!=visitedState.end()) {
					if(visitedState[ep].vState == hVisited) {// �������˵�POI���޸�״̬
						visitedState[ep].vState == visited;
						//---��Ӿ���visitedState
						if(tmpDV.N < NewNodeID) {
							visitedState[ep].iDisToQuery = tmpDV.disTQ;
						} else {
							visitedState[ep].jDisToQuery = tmpDV.disTQ;
						}
						//�����˼���POIs
						if(PtGrpKey==-1){
						//cout<<"No POI existed on Edge where Q located."<<endl;
						} else {

							getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
							//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
							//Notice the order POI visited on edge
							for(int j=0; j<PtNumOnEdge; j++){
								poicnt++;
								getVarE(PT_DIST,Ref(PtDist),PtGrpKey,j);

								getVarE(PT_VCT,tempKwd,PtGrpKey,j);
								int num = sizeof(tempKwd)/sizeof(int);
								for(int loop=0; loop<sizeof(tempKwd); loop+sizeof(int)){
									int temp;
									memcpy(&temp,tempKwd+loop,sizeof(int));
									tempKwds.insert(temp);
								}

								if(LcontianRIKwd(tempKwds,Q.kwd)){
									POI tmp;
									//tmp.dist_toquery = fabs(PtDist-Q.dist_Ni);
									float d1,d2;
									d1 = visitedState[ep].iDisToQuery + PtDist;
									d2 = visitedState[ep].jDisToQuery + EdgeDist - PtDist;
									tmp.dist_toquery = d1<d2 ? d1:d2;
									
									getVarE(PT_ATTRIBUTE,tmp.attr,PtGrpKey,j);
									getVarE(PT_P,&tmp.poid,PtGrpKey,j);
						
									if(tmp.dist_toquery <= Q.distCnst) {
										if(isBeDominate(tmp, cS, Q)) {
										} else {
											cS.push_back(tmp);
										}
									}             
								}// endif
							}
						}
					} else if(visitedState[ep].vState == visited) { //�޸ľ���
						//---��Ӿ���visitedState
						if(tmpDV.N < NewNodeID) {
							visitedState[ep].iDisToQuery = tmpDV.disTQ;
						} else {
							visitedState[ep].jDisToQuery = tmpDV.disTQ;
						}
					} else {

					}
				} else {
					edgePair ep;
					edgeState es;
					if(tmpDV.N < NewNodeID) {
						ep.Ni = tmpDV.N;
						ep.Nj = NewNodeID;
						es.iDisToQuery = tmpDV.disTQ;
						es.jDisToQuery = INFINITE_MAX;
					} else {
						ep.Ni = NewNodeID;
						ep.Nj = tmpDV.N;
						es.iDisToQuery = INFINITE_MAX;
						es.jDisToQuery = tmpDV.disTQ;
					}
					if(distTQ.find(NewNodeID)!=distTQ.end()){
						if(distTQ[NewNodeID]>(tmpDV.disTQ+EdgeDist)) {
							distTQ[NewNodeID] = tmpDV.disTQ+EdgeDist;
						}
					} else {
						distTQ[NewNodeID] = tmpDV.disTQ+EdgeDist;
					}
					if(distTQ[NewNodeID]<=Q.distCnst) {
						//�������߼���
						es.vState = visited;
						visitedState[ep] = es;

						if(PtGrpKey==-1){
						//cout<<"No POI existed on Edge where Q located."<<endl;
						} else {
							getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
							//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
							//Notice the order POI visited on edge
							for(int j=0; j<PtNumOnEdge; j++){
								poicnt++;
								getVarE(PT_DIST,Ref(PtDist),PtGrpKey,j);

								getVarE(PT_VCT,tempKwd,PtGrpKey,j);
								int num = sizeof(tempKwd)/sizeof(int);
								for(int loop=0; loop<sizeof(tempKwd); loop+sizeof(int)){
									int temp;
									memcpy(&temp,tempKwd+loop,sizeof(int));
									tempKwds.insert(temp);
								}

								if(LcontianRIKwd(tempKwds,Q.kwd)){
									POI tmp; //��ʱ����ܲ�����ʵ����,Ҳ����Ӱ����
									if(tmpDV.N<NewNodeID) {
										tmp.dist_toquery = tmpDV.disTQ+PtDist;
									} else {
										tmp.dist_toquery = tmpDV.disTQ+EdgeDist-PtDist;
									}
									
								
									getVarE(PT_ATTRIBUTE,tmp.attr,PtGrpKey,j);
									getVarE(PT_P,&tmp.poid,PtGrpKey,j);
						
									if(tmp.dist_toquery <= Q.distCnst) {
										if(isBeDominate(tmp, cS, Q)) {
										} else {
											cS.push_back(tmp);
										}
									}             
								}// endif
							}
						}

					} else {// ���ֻ�в����ڱ��ϣ�hVisited
						es.vState = hVisited;
						visitedState[ep] = es;
					}
				}

				//getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
				//getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);
				if(PtGrpKey==-1){
						//cout<<"No POI existed on Edge where Q located."<<endl;
				} else {

					getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
					//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
					//Notice the order POI visited on edge
					for(int j=0; j<PtNumOnEdge; j++){
						poicnt++;
						getVarE(PT_DIST,Ref(PtDist),PtGrpKey,j);

						getVarE(PT_VCT,tempKwd,PtGrpKey,j);
						int num = sizeof(tempKwd)/sizeof(int);
						for(int loop=0; loop<sizeof(tempKwd); loop+sizeof(int)){
							int temp;
							memcpy(&temp,tempKwd+loop,sizeof(int));
							tempKwds.insert(temp);
						}

						if(LcontianRIKwd(tempKwds,Q.kwd)){
							POI tmp;
							tmp.dist_toquery = fabs(PtDist-Q.dist_Ni);
							getVarE(PT_ATTRIBUTE,tmp.attr,PtGrpKey,j);
							getVarE(PT_P,&tmp.poid,PtGrpKey,j);
						
							if(tmp.dist_toquery <= Q.distCnst) {
								if(isBeDominate(tmp, cS, Q)) {
								} else {
									cS.push_back(tmp);
								}
							}             
						}// endif
					}
				}

			
			}// end for
		}// endif
	}// end while
    
	// ����hVisited

	map<edgePair, edgeState, eSComparison>::iterator iterhV;
	for(iterhV=visitedState.begin(); iterhV!=visitedState.end(); ++iterhV) {
		if(iterhV->second.vState == hVisited) {//����ֻ�о��뵱ǰ�ڵ��������Լ����POI
			int nodei,nodej;
			if(iterhV->second.iDisToQuery == INFINITE_MAX) {
				nodei = iterhV->first.Nj;
				nodej = iterhV->first.Ni;
			} else {
				nodei = iterhV->first.Nj;
				nodej = iterhV->first.Ni;
			}
			//��ȡ���ϵ�POI
			AdjGrpAddr=getAdjListGrpAddr(nodei);
			getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
			for (int i=0; i<AdjListSize; i++){
				getVarE(ADJNODE_A,Ref(NewNodeID),AdjGrpAddr,i);
				//getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
				if(NewNodeID == nodej){
					getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
					getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);

					nOfEdgeExpd++;

					//cout<<"("<<Ni<<","<<Nj<<") EdgeDist:"<<EdgeDist;
					if(PtGrpKey==-1){
						//cout<<"No POI existed on Edge where Q located."<<endl;
					}
					else {
						getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
						//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
						//Notice the order POI visited on edge
						for(int j=0; j<PtNumOnEdge; j++){
							poicnt++;
							getVarE(PT_DIST,Ref(PtDist),PtGrpKey,j);

							getVarE(PT_VCT,tempKwd,PtGrpKey,j);
							int num = sizeof(tempKwd)/sizeof(int);
							for(int loop=0; loop<sizeof(tempKwd); loop+sizeof(int)){
								int temp;
								memcpy(&temp,tempKwd+loop,sizeof(int));
								tempKwds.push_back(temp);
							}

							if(LcontianRIKwd(tempKwds,Q.kwd)){
								POI tmp;
								if(nodei<nodej){
									tmp.dist_toquery = iterhV->second.iDisToQuery+PtDist;
								} else {
									tmp.dist_toquery = iterhV->second.jDisToQuery+EdgeDist-PtDist;
								}
								getVarE(PT_ATTRIBUTE,tmp.attr,PtGrpKey,j);
								getVarE(PT_P,&tmp.poid,PtGrpKey,j);
						
								if(tmp.dist_toquery <= Q.distCnst) {
									if(isBeDominate(tmp, cS, Q)) {
									} else {
										cS.push_back(tmp);
									}
								}             
							}// endif
						}
                
					}// endelse
				}
			}
		}
	}
	

	for(int i=0; i<cS.size(); i++) {
		rS.push_back(cS[i].poid);
		//cout<<cS[i].poid<<" ";
	}
	//cout<<endl;
    //cout<<"POI Read           #:"<<poicnt<<endl;  
    //cout<<"Road Dominate Tested   #:"<<nOfDomiTest<<endl;
    //cout<<"Road Edge Expanded #:"<<nOfEdgeExpd<<endl;
    //base2numedgeexpand+=edgeexpanded;
    //base2numnodeexpand+=roadnodecnt;
}
//����s��cans�ľ���
vector<float> dijkstra_candidate( QueryPoint &Q, int tid, vector<TreeNode> &EGT ) {
	//bottomToUp �������Ϲ���distTQ������·��
	if(bottomToUp == 0) {
		bottomToUp = 1;
		int pid = paID[Q.Ni].part; //�����ѯ��λ��treenode �ڲ�
		pathID.push_back(pid);
		vector<int>::iterator it;
		//�������˵㵽�߽�ľ��룬ȥ��С
		EGT[pid].minDistTQ = INFINITE_MAX;
		for(int i=0; i<EGT[pid].borders.size(); i++) {
			distTQ[i] = INFINITE_MAX;
		}
		int posi,posj,size=EGT[pid].leafnodes.size();
		it = find(EGT[pid].leafnodes.begin(),EGT[pid].leafnodes.end(),Q.Ni);
		posi = (*it);
		it = find(EGT[pid].leafnodes.begin(),EGT[pid].leafnodes.end(),Q.Nj);
		posj = (*it);
		for(int i=0; i<EGT[pid].borders.size(); i++) {
			float distI = Q.dist_Ni+EGT[pid].mind[i*size+posi];
			float distJ = queryEdgeDist-Q.dist_Ni+EGT[pid].mind[i*size+posj];
			if(distI<distJ) {
				EGT[pid].distTQ[i]=distI;
				if(distI<EGT[pid].minDistTQ) EGT[pid].minDistTQ=distI;
			} else {
				EGT[pid].distTQ[i]=distJ;
				if(distJ<EGT[pid].minDistTQ) EGT[pid].minDistTQ=distJ;
			}
		}
		int preID;
		while(pid!=0&&EGT[pid].minDistTQ<Q.distCnst) {
			preID = pid;
			pid=EGT[pid].father;
			pathID.push_back(pid);
			//���㸸�ڵ�ľ���
			//vector<int>::iterator it;
			//�������˵㵽�߽�ľ��룬ȥ��С
			EGT[pid].minDistTQ = INFINITE_MAX;
			for(int i=0; i<EGT[pid].borders.size(); i++) {
				EGT[pid].distTQ[i] = INFINITE_MAX;
			}
			size=EGT[pid].union_borders.size();
			int posb,posu;
			//���ڸò��ÿ��border���������һ�㵽�ò����̾���
			for(int i=0; i<EGT[pid].borders; i++) {
				for(int j=0; j<EGT[preID].borders; j++) {
					//��λ����node��union_border��λ��		
					it = find(EGT[pid].union_borders.begin(),EGT[pid].union_borders.end(),EGT[preID].borders[j]);
					posb = (*it);
					it = find(EGT[pid].union_borders.begin(),EGT[pid].union_borders.end(),EGT[pid].borders[i]);
					posu = (*it);
					float mindPos;
					if(posb<posu) {
						mindPos = (posu*(posu+1))/2+posb;
					} else {
						mindPos = (posb*(posb+1))/2+posu;
					}
					if((EGT[preID].distTQ[j]+EGT[pid].mind[mindPos])<EGT[pid].distTQ[i]) {
						EGT[pid].distTQ[i] = EGT[preID].distTQ[j]+EGT[pid].mind[mindPos];
					}
				}//endfor
				if(EGT[pid].distTQ[i]<EGT[pid].minDistTQ) EGT[pid].minDistTQ = EGT[pid].distTQ[i];
			}//endfor

		}//endwhile
	}//endif
	//����С��ͬ���ڵ�,�������·��
	int commonID = tid;
	vector<int> ctPathID;
	ctPathID.push_back(commonID);
	while(true) {
		if(find(pathID.begin(),pathID.end(),commonID)==pathID.end()) {
			commonID = EGT[commonID].father;
			ctPathID.push_back(commonID);
		}else {
			break;
		}
	}
	//Ѱ��pathID��commonID����һ��
	int pathid = paID[Q.Ni].part;
	while(EGT[pathid].father!=commonID) {
		pathid = EGT[pathid].father;
	}
	//��ʼ��תȥ���������ߣ���ʱ����Ҫ���ϵ��´�����
	//���������ұ�������Ľڵ�
	int upMostID = ctPathID.size()-1;
	if(EGT[upMostID].distTQ.size()!=EGT[upMostID].borders.size()) {
		EGT[upMostID].minDistTQ = INFINITE_MAX;
		for(int i=0; i<EGT[upMostID].borders.size(); i++) {
			EGT[upMostID].distTQ[i] = INFINITE_MAX;
		}
		int posutm,pospth,sizeC = EGT[commonID].union_borders.size();
		for(int i=0; i<EGT[upMostID].borders.size(); i++) {
			for(int j=0; j<EGT[pathid].borders.size(); j++) {
				//��λ
				vector<int>::iterator itr;
				itr = find(EGT[commonID].union_borders.begin(),EGT[commonID].union_borders.end(),EGT[upMostID].borders[i]);
				posutm = (*it);
				itr = find(EGT[commonID].union_borders.begin(),EGT[commonID].union_borders.end(),EGT[pathid].borders[j]);
				pospth = (*it);
				float mindPos;
				if(posutm<pospth) {
					mindPos = (pospth*(pospth+1))/2+posutm;
				} else {
					mindPos = (posutm*(posutm+1))/2+pospth;
				}
				if((EGT[pathid].distTQ[j]+EGT[upMostID].mind[mindPos])<EGT[upMostID].distTQ[i]) {
					EGT[upMostID].distTQ[i] = EGT[pathid].distTQ[j]+EGT[upMostID].mind[mindPos];
				}
			}
			if(EGT[upMostID].distTQ[i]<EGT[upMostID].minDistTQ) EGT[upMostID].minDistTQ = EGT[upMostID].distTQ[i];
		}
	}
	//�ӵ����ڶ�����ʼ
	int uBID,fID;
	for(int k=ctPathID.size()-2; k>=0; k--) {
		uBID = ctPathID[k];
		fID = EGT[uBID].father;
		//preID = pid;
		//pid=EGT[pid].father;
		//pathID.push_back(pid);
		//���㸸�ڵ�ľ���
		//vector<int>::iterator it;
		//�������˵㵽�߽�ľ��룬ȥ��С
		if(EGT[uBID].distTQ.size()==EGT[uBID].borders.size()) continue;
		EGT[uBID].minDistTQ = INFINITE_MAX;
		for(int i=0; i<EGT[uBID].borders.size(); i++) {
			EGT[uBID].distTQ[i] = INFINITE_MAX;
		}
		int sizef=EGT[fID].union_borders.size();
		int posdc,posdf;
		//���ڸò��ÿ��border���������һ�㵽�ò����̾���
		for(int i=0; i<EGT[uBID].borders; i++) {
			for(int j=0; j<EGT[fID].borders; j++) {
				//��λ����node��union_border��λ��		
				it = find(EGT[fID].union_borders.begin(),EGT[fID].union_borders.end(),EGT[uBID].borders[i]);
				posdc = (*it);
				it = find(EGT[fID].union_borders.begin(),EGT[fID].union_borders.end(),EGT[fID].borders[j]);
				posdf = (*it);
				float mindPos;
				if(posdc<posdf) {
					mindPos = (posdf*(posdf+1))/2+posdc;
				} else {
					mindPos = (posdc*(posdc+1))/2+posdf;
				}
				if((EGT[fID].distTQ[j]+EGT[fID].mind[mindPos])<EGT[uBID].distTQ[i]) {
					EGT[uBID].distTQ[i] = EGT[fID].distTQ[j]+EGT[fID].mind[mindPos];
				}
			}
		}//endfor

	}//endfor

	//�����·��


}

void EGBUAlgorithm(const QueryPoint &Q) {
	//process leadnode, two situation:1, not in same part 2, in the same part
	int Ni=Q.Ni,Nj=Q.Nj;//Ni must less than Nj
    int AdjGrpAddr,AdjListSize,NewNodeID,PtGrpKey,PtNumOnEdge;
    //unsigned long long keywords;
	vector<int> tempKwds;
	set<int> visNodes;
	int *tempKwd;
	set<int> edgeSumKwds;
	float edgeSumAttr[ATTRIBUTE_DIMENSION][2];
    float EdgeDist,PtDist;
	// locate the treeNodeID
	int trNodeIDi,trNodeIDj;
	int currTNID;
	float minDist;
	dijkVisit dvi,dvj;

	AdjGrpAddr=getAdjListGrpAddr(Ni);
    getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
    AdjGrpAddr=getAdjListGrpAddr(Ni);
    getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
    for (int i=0; i<AdjListSize; i++){
        getVarE(ADJNODE_A,Ref(NewNodeID),AdjGrpAddr,i);
        //getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
        if(NewNodeID == Nj){
            getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
			//��¼��ѯ�߾���
			queryEdgeDist = EdgeDist;
            getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);
			getVarE(SUMKWD_A,&edgeSumKwds.begin(),AdjGrpAddr,i);
			getVarE(SUMATTR_A,edgeSumAttr,AdjGrpAddr,i);

			edgePair epi,epj;
			edgeState esi,esj;
			epi.Ni = 0;
			epi.Nj = Q.Ni;
			esi.vState = visited;
			esi.iDisToQuery = 0;
			esi.jDisToQuery = Q.dist_Ni;

			epj.Ni = 0;
			epj.Nj = Q.Nj;
			esj.iDisToQuery = 0;
			esj.jDisToQuery = EdgeDist-Q.dist_Ni;

			visitedState[epi] = esi;
			visitedState[epj] = esj;

			if(!LcontianRIKwd(edgeSumKwds,Q.kwd)){
				
				break;
			}
			if(isBeBoundDominate(edgeSumAttr, cS, Q)){
				
				break;
			}
            nOfEdgeExpd++;

            //cout<<"("<<Ni<<","<<Nj<<") EdgeDist:"<<EdgeDist;
            if(PtGrpKey==-1){
                //cout<<"No POI existed on Edge where Q located."<<endl;
            }
            else {
                getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
                //cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
                //Notice the order POI visited on edge
                for(int j=0; j<PtNumOnEdge; j++){
                    poicnt++;
                    getVarE(PT_DIST,Ref(PtDist),PtGrpKey,j);
                    getVarE(PT_VCT,tempKwd,PtGrpKey,j);
					
					int num = sizeof(tempKwd)/sizeof(int);
					for(int loop=0; loop<sizeof(tempKwd); loop+sizeof(int)){
						int temp;
						memcpy(&temp,tempKwd+loop,sizeof(int));
						tempKwds.insert(temp);
					}

					if(LcontianRIKwd(tempKwds,Q.kwd)){
                        POI tmp;
						tmp.dist_toquery = fabs(PtDist-Q.dist_Ni);
						getVarE(PT_ATTRIBUTE,tmp.attr,PtGrpKey,j);
						getVarE(PT_P,&tmp.poid,PtGrpKey,j);
						
						if(tmp.dist_toquery <= Q.distCnst) {
							if(isBeDominate(tmp, cS, Q)) {
							} else {
								cS.push_back(tmp);
							}
						}             
                    }// endif
                }
                
            }// endelse
        }
    }

	trNodeIDi = paID[Ni].part;
	trNodeIDj = paID[Nj].part;
	if(trNodeIDi == trNodeIDj) { //������һ�����ֿ�����
		currTNID = trNodeIDi;
		dvi.N = Q.Ni;
		dvi.disTQ = Q.dist_Ni;

		dvj.N = Q.Nj;
		dvj.disTQ = EdgeDist-Q.dist_Ni;

		dvq.push(dvi);
		dvq.push(dvj);
	} /*else { //��ѯ���ڱ߽������
		if(trNodeIDi < trNodeIDj) {
			currTNID = trNodeIDi;
			dvi.N = Q.Ni;
			dvi.disTQ = Q.dist_Ni;
			dvq.push(dvi);
		} else {
			currTNID = trNodeIDj;
			dvj.N = Q.Nj;
			dvj.disTQ = EdgeDist-Q.dist_Ni;
			dvq.push(dvj);
		}
	}*/
	//------------&&&&&&&&&&&&&&�ڲ��ڵ㴦��
	while(!dvq.empty()) {
		dijkVisit tmpDV = dvq.top();
		dvq.pop();
		if(paID[tmpDV.N].part!=currTNID) continue;
		if(tmpDV.disTQ <= Q.distCnst) {
			AdjGrpAddr=getAdjListGrpAddr(tmpDV.N);
			getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
			for (int i=0; i<AdjListSize; i++){
				getVarE(ADJNODE_A,Ref(NewNodeID),AdjGrpAddr,i);
				getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
				getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);
				edgePair ep;
				ep.Ni = NewNodeID<tmpDV.N?NewNodeID:tmpDV.N;
				ep.Nj = NewNodeID>tmpDV.N?NewNodeID:tmpDV.N;
				map<edgePair, edgeState, eSComparison>::iterator iter;
				iter = visitedState.find(ep);
				if(iter!=visitedState.end()) {
					if(visitedState[ep].vState == hVisited) {// �������˵�POI���޸�״̬
						visitedState[ep].vState == visited;
						//---��Ӿ���visitedState
						if(tmpDV.N < NewNodeID) {
							visitedState[ep].iDisToQuery = tmpDV.disTQ;
						} else {
							visitedState[ep].jDisToQuery = tmpDV.disTQ;
						}
						//�����˼���POIs
						if(PtGrpKey==-1){
						//cout<<"No POI existed on Edge where Q located."<<endl;
						} else {

							getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
							//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
							//Notice the order POI visited on edge
							for(int j=0; j<PtNumOnEdge; j++){
								poicnt++;
								getVarE(PT_DIST,Ref(PtDist),PtGrpKey,j);

								getVarE(PT_VCT,tempKwd,PtGrpKey,j);
								int num = sizeof(tempKwd)/sizeof(int);
								for(int loop=0; loop<sizeof(tempKwd); loop+sizeof(int)){
									int temp;
									memcpy(&temp,tempKwd+loop,sizeof(int));
									tempKwds.insert(temp);
								}

								if(LcontianRIKwd(tempKwds,Q.kwd)){
									POI tmp;
									//tmp.dist_toquery = fabs(PtDist-Q.dist_Ni);
									float d1,d2;
									d1 = visitedState[ep].iDisToQuery + PtDist;
									d2 = visitedState[ep].jDisToQuery + EdgeDist - PtDist;
									tmp.dist_toquery = d1<d2 ? d1:d2;
									
									getVarE(PT_ATTRIBUTE,tmp.attr,PtGrpKey,j);
									getVarE(PT_P,&tmp.poid,PtGrpKey,j);
						
									if(tmp.dist_toquery <= Q.distCnst) {
										if(isBeDominate(tmp, cS, Q)) {
										} else {
											cS.push_back(tmp);
										}
									}             
								}// endif
							}
						}
					} else if(visitedState[ep].vState == visited) { //�޸ľ���
						//---��Ӿ���visitedState
						if(tmpDV.N < NewNodeID) {
							visitedState[ep].iDisToQuery = tmpDV.disTQ;
						} else {
							visitedState[ep].jDisToQuery = tmpDV.disTQ;
						}
					} else {

					}
				} else {
					edgePair ep;
					edgeState es;
					if(tmpDV.N < NewNodeID) {
						ep.Ni = tmpDV.N;
						ep.Nj = NewNodeID;
						es.iDisToQuery = tmpDV.disTQ;
						es.jDisToQuery = INFINITE_MAX;
					} else {
						ep.Ni = NewNodeID;
						ep.Nj = tmpDV.N;
						es.iDisToQuery = INFINITE_MAX;
						es.jDisToQuery = tmpDV.disTQ;
					}
					if(distTQ.find(NewNodeID)!=distTQ.end()){
						if(distTQ[NewNodeID]>(tmpDV.disTQ+EdgeDist)) {
							distTQ[NewNodeID] = tmpDV.disTQ+EdgeDist;
						}
					} else {
						distTQ[NewNodeID] = tmpDV.disTQ+EdgeDist;
					}
					if(distTQ[NewNodeID]<=Q.distCnst) {
						//�������߼���
						es.vState = visited;
						visitedState[ep] = es;



						if(PtGrpKey==-1){
						//cout<<"No POI existed on Edge where Q located."<<endl;
						} else {
							getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
							//��ùؼ�����Ϣ��������Ϣ��ֱ�ӹ���
							if(!LcontianRIKwd(edgeSumKwds,Q.kwd)) {
								//�����߲���
								continue;
							}
							if(isBeBoundDominate(edgeSumAttr, cS, Q)) {
								//�����߲���
								continue;
							}
							//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
							//Notice the order POI visited on edge
							for(int j=0; j<PtNumOnEdge; j++){
								poicnt++;
								getVarE(PT_DIST,Ref(PtDist),PtGrpKey,j);

								getVarE(PT_VCT,tempKwd,PtGrpKey,j);
								int num = sizeof(tempKwd)/sizeof(int);
								for(int loop=0; loop<sizeof(tempKwd); loop+sizeof(int)){
									int temp;
									memcpy(&temp,tempKwd+loop,sizeof(int));
									tempKwds.insert(temp);
								}

								if(LcontianRIKwd(tempKwds,Q.kwd)){
									POI tmp; //��ʱ����ܲ�����ʵ����,Ҳ����Ӱ����
									if(tmpDV.N<NewNodeID) {
										tmp.dist_toquery = tmpDV.disTQ+PtDist;
									} else {
										tmp.dist_toquery = tmpDV.disTQ+EdgeDist-PtDist;
									}
									
								
									getVarE(PT_ATTRIBUTE,tmp.attr,PtGrpKey,j);
									getVarE(PT_P,&tmp.poid,PtGrpKey,j);
						
									if(tmp.dist_toquery <= Q.distCnst) {
										if(isBeDominate(tmp, cS, Q)) {
										} else {
											cS.push_back(tmp);
										}
									}             
								}// endif
							}
						}

					} else {// ���ֻ�в����ڱ��ϣ�hVisited
						if(!LcontianRIKwd(edgeSumKwds,Q.kwd)) {
							//�����߲���
							es.vState = visited;
							visitedState[ep] = es;
							continue;
						}
						if(isBeBoundDominate(edgeSumAttr, cS, Q)) {
							//�����߲���
							es.vState = visited;
							visitedState[ep] = es;
							continue;
						}
						es.vState = hVisited;
						visitedState[ep] = es;
					}
				}

				//getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
				//getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);
				if(PtGrpKey==-1){
						//cout<<"No POI existed on Edge where Q located."<<endl;
				} else {

					getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
					//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
					//Notice the order POI visited on edge
					for(int j=0; j<PtNumOnEdge; j++){
						poicnt++;
						getVarE(PT_DIST,Ref(PtDist),PtGrpKey,j);

						getVarE(PT_VCT,tempKwd,PtGrpKey,j);
						int num = sizeof(tempKwd)/sizeof(int);
						for(int loop=0; loop<sizeof(tempKwd); loop+sizeof(int)){
							int temp;
							memcpy(&temp,tempKwd+loop,sizeof(int));
							tempKwds.insert(temp);
						}

						if(LcontianRIKwd(tempKwds,Q.kwd)){
							POI tmp;
							tmp.dist_toquery = fabs(PtDist-Q.dist_Ni);
							getVarE(PT_ATTRIBUTE,tmp.attr,PtGrpKey,j);
							getVarE(PT_P,&tmp.poid,PtGrpKey,j);
						
							if(tmp.dist_toquery <= Q.distCnst) {
								if(isBeDominate(tmp, cS, Q)) {
								} else {
									cS.push_back(tmp);
								}
							}             
						}// endif
					}
				}

			
			}// end for
		}// endif
	}// end while






	//upward
	//-------------------top treeNodeID is 0
	visitedSet.insert(currTNID);
	TreeNode tn;
	while(!tnq.empty()||currTNID != 0) { //
		if(tnq.empty()) {
			currTNID = EGT[currTNID].father;
			// if currTNID is filter by summary information and distance 
			if(!LcontianRIKwd(EGT[currTNID].union_kwd,Q.kwd)) {
				if(visitedSet.find(currTNID)==visitedSet.end()) visitedSet.insert(currTNID);
				continue;
			}
			if(isBeBoundDominate(EGT[currTNID].attrBound, cS, Q)) {
				if(visitedSet.find(currTNID)==visitedSet.end()) visitedSet.insert(currTNID);
				continue;
			}
			if(EGT[currTNID].minDistTQ > Q.distCnst) {
				if(visitedSet.find(currTNID)==visitedSet.end()) visitedSet.insert(currTNID);
				continue;
			}
			// add unvisited sub node into tnq and compute mindist
			for(int i=0; i<EGT[currTNID].children.size();) {
				int cid = EGT[currTNID].children[i];
				if(visitedSet.find(cid)==visitedSet.end()) {
					//if is not filter by summary information
					if(!LcontianRIKwd(EGT[cid].union_kwd,Q.kwd)) {
						if(visitedSet.find(cid)==visitedSet.end()) visitedSet.insert(cid);
						continue;
					}
					if(isBeBoundDominate(EGT[cid].attrBound, cS, Q)) {
						if(visitedSet.find(cid)==visitedSet.end()) visitedSet.insert(cid);
						continue;
					}
					if(EGT[cid].minDistTQ > Q.distCnst) {
						if(visitedSet.find(cid)==visitedSet.end()) visitedSet.insert(cid);
						continue;
					}
					//compute the minDistTQ,refDistTQ&minDistTQ
					vector<float> distTQ;
					distTQ = dijkstra_candidate( Q, EGT[cid].borders, EGT );
					float minDist = INFINITE_MAX;
					for(int k=0; k<distTQ.size(); k++) {
						if(distTQ[k] < minDist) minDist = distTQ[k];
						EGT[cid].refDistTQ.push_back(Q.distCnst-distTQ[k]);
					}
					EGT[cid].minDistTQ = minDist;
					//----------------&&&&&&&����Dist��
					tnq.push(EGT[cid]);
				}
			}

		}
	
		tn = tnq.top();
		tnq.pop();
		if(tn.isleaf) { //�����Ҷ�ӽڵ㣬����
			// utilize the mind to retrieve this node meet distance constraint
			//---------&&&&&&&&&&&&�����������Լ����node���ϣ���
			vector<float> tempDist;
			for(int k=0; k<tn.leafnodes.size(); k++) {
				tempDist.push_back(INFINITE_MAX);
			}
			for(int k=0; k<tn.refDistTQ.size(); k++) {
				//float mindist = INFINITE_MAX;
				if(tn.refDistTQ[k]<0) continue;
				int pos,size = tn.leafnodes.size();
				float tempD;
				for(int t=0; t<tn.leafnodes.size(); t++) {
					pos = k*size+t;
					tempD = Q.distCnst-tn.refDistTQ[k]+tn.mind[pos];
					if(tempD<=Q.distCnst) {
						if(tempD<tempDist[t]) tempDist[t] = tempD;
					}
				}
			}
			//����dvq
			for(int k=0; k<tempDist.size(); k++) {
				if(tempDist[k]>Q.distCnst) continue;
				//Ϊÿ�����������ĵ�չ��
				/*dijkVisit dv;
				dv.N = tn.leafnodes[k];
				dv.disTQ = tempDis[k];
				dvq.push(dv);*/
				AdjGrpAddr=getAdjListGrpAddr(tn.leafnodes[k]);
				getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
				for (int i=0; i<AdjListSize; i++){
					getVarE(ADJNODE_A,Ref(NewNodeID),AdjGrpAddr,i);
					getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
					getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);
					getVarE(SUMKWD_A,&edgeSumKwds.begin(),AdjGrpAddr,i);
					getVarE(SUMATTR_A,edgeSumAttr,AdjGrpAddr,i);
					edgePair ep;
					ep.Ni = NewNodeID<tn.leafnodes[k]?NewNodeID:tn.leafnodes[k];
					ep.Nj = NewNodeID>tn.leafnodes[k]?NewNodeID:tn.leafnodes[k];
					map<edgePair, edgeState, eSComparison>::iterator iter;
					iter = visitedState.find(ep);
					if(iter!=visitedState.end()) {
						if(visitedState[ep].vState == hVisited) {// �������˵�POI���޸�״̬
							visitedState[ep].vState == visited;
							//---��Ӿ���visitedState
							if(tn.leafnodes[k] < NewNodeID) {
								visitedState[ep].iDisToQuery = tempDis[k];
							} else {
								visitedState[ep].jDisToQuery = tempDis[k];
							}
							//�����˼���POIs
							if(PtGrpKey==-1){
							//cout<<"No POI existed on Edge where Q located."<<endl;
							} else {

								getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
								//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
								//Notice the order POI visited on edge
								for(int j=0; j<PtNumOnEdge; j++){
									poicnt++;
									getVarE(PT_DIST,Ref(PtDist),PtGrpKey,j);

									getVarE(PT_VCT,tempKwd,PtGrpKey,j);
									int num = sizeof(tempKwd)/sizeof(int);
									for(int loop=0; loop<sizeof(tempKwd); loop+sizeof(int)){
										int temp;
										memcpy(&temp,tempKwd+loop,sizeof(int));
										tempKwds.insert(temp);
									}

									if(LcontianRIKwd(tempKwds,Q.kwd)){
										POI tmp;
										//tmp.dist_toquery = fabs(PtDist-Q.dist_Ni);
										float d1,d2;
										d1 = visitedState[ep].iDisToQuery + PtDist;
										d2 = visitedState[ep].jDisToQuery + EdgeDist - PtDist;
										tmp.dist_toquery = d1<d2 ? d1:d2;
									
										getVarE(PT_ATTRIBUTE,tmp.attr,PtGrpKey,j);
										getVarE(PT_P,&tmp.poid,PtGrpKey,j);
						
										if(tmp.dist_toquery <= Q.distCnst) {
											if(isBeDominate(tmp, cS, Q)) {
											} else {
												cS.push_back(tmp);
											}
										}             
									}// endif
								}
							}
						} else if(visitedState[ep].vState == visited) { //�޸ľ���
							//---��Ӿ���visitedState
							if(tmpDV.N < NewNodeID) {
								visitedState[ep].iDisToQuery = tempDis[k];
							} else {
								visitedState[ep].jDisToQuery = tempDis[k];
							}
						} else {
							
					
						}
					} else {
						//�ж��Ƿ�Ҫ���뵽ep��			
							edgeState es;
							if(tn.leafnodes[k] < NewNodeID) {
								ep.Ni = tn.leafnodes[k];
								ep.Nj = NewNodeID;
								es.iDisToQuery = tempDis[k];
								es.jDisToQuery = INFINITE_MAX;
							} else {
								ep.Ni = NewNodeID;
								ep.Nj = tn.leafnodes[k];
								es.iDisToQuery = INFINITE_MAX;
								es.jDisToQuery = tempDis[k];
							}
							if(distTQ.find(NewNodeID)!=distTQ.end()){
								if(distTQ[NewNodeID]>(tempDis[k]+EdgeDist)) {
									distTQ[NewNodeID] = tempDis[k]+EdgeDist;
								}
							} else {
								distTQ[NewNodeID] = tempDis[k]+EdgeDist;
							}
							if(distTQ[NewNodeID]<=Q.distCnst) {
								//�������߼���
								es.vState = visited;
								visitedState[ep] = es;

								if(PtGrpKey==-1){
								//cout<<"No POI existed on Edge where Q located."<<endl;
								} else {
									//��ùؼ�����Ϣ��������Ϣ��ֱ�ӹ���
									if(!LcontianRIKwd(edgeSumKwds,Q.kwd)) {
										//�����߲���
										continue;
									}
									if(isBeBoundDominate(edgeSumAttr, cS, Q)) {
										//�����߲���
										continue;
									}

									//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
									//Notice the order POI visited on edge
									for(int j=0; j<PtNumOnEdge; j++){
										poicnt++;
										getVarE(PT_DIST,Ref(PtDist),PtGrpKey,j);

										getVarE(PT_VCT,tempKwd,PtGrpKey,j);
										int num = sizeof(tempKwd)/sizeof(int);
										for(int loop=0; loop<sizeof(tempKwd); loop+sizeof(int)){
											int temp;
											memcpy(&temp,tempKwd+loop,sizeof(int));
											tempKwds.insert(temp);
										}

										if(LcontianRIKwd(tempKwds,Q.kwd)){
											POI tmp; //��ʱ����ܲ�����ʵ����,Ҳ����Ӱ����
											if(tmpDV.N<NewNodeID) {
												tmp.dist_toquery = tmpDV.disTQ+PtDist;
											} else {
												tmp.dist_toquery = tmpDV.disTQ+EdgeDist-PtDist;
											}
									
								
											getVarE(PT_ATTRIBUTE,tmp.attr,PtGrpKey,j);
											getVarE(PT_P,&tmp.poid,PtGrpKey,j);
						
											if(tmp.dist_toquery <= Q.distCnst) {
												if(isBeDominate(tmp, cS, Q)) {
												} else {
													cS.push_back(tmp);
												}
											}             
										}// endif
									}
								}

							} else {// ���ֻ�в����ڱ��ϣ�hVisited
								//if���ؼ��ֵȹ��˵�����ֱ������visited
								//��ùؼ�����Ϣ��������Ϣ��ֱ�ӹ���
								if(!LcontianRIKwd(edgeSumKwds,Q.kwd)) {
									//�����߲���
									es.vState = visited;
									visitedState[ep] = es;
									continue;
								}
								if(isBeBoundDominate(edgeSumAttr, cS, Q)) {
									//�����߲���
									es.vState = visited;
									visitedState[ep] = es;
									continue;
								}
								es.vState = hVisited;
								visitedState[ep] = es;
							}
						}


					}//endif
				}//endfor

		} else { //����Ҷ�ӽڵ㴦��
			for(int i=0; i<tn.children.size();) {
				int cid = tn.children[i];
				if(visitedSet.find(cid)==visitedSet.end()) {
					//if is not filter by summary information
					if(!LcontianRIKwd(EGT[cid].union_kwd,Q.kwd)) {
						if(visitedSet.find(cid)==visitedSet.end()) visitedSet.insert(cid);
						continue;
					}
					if(isBeBoundDominate(EGT[cid].attrBound, cS, Q)) {
						if(visitedSet.find(cid)==visitedSet.end()) visitedSet.insert(cid);
						continue;
					}
					if(EGT[cid].minDistTQ > Q.distCnst) {
						if(visitedSet.find(cid)==visitedSet.end()) visitedSet.insert(cid);
						continue;
					}
					//compute the minDistTQ
					vector<float> distTQ;
					distTQ = dijkstra_candidate( Q, cid, EGT );
					float minDist = INFINITE_MAX;
					for(int k=0; k<distTQ.size(); k++) {
						if(distTQ[k] < minDist) minDist = distTQ[k];
						EGT[cid].refDistTQ.push_back(Q.distCnst-distTQ[k]);
					}
					EGT[cid].minDistTQ = minDist;
					//----------------&&&&&&&
					tnq.push(EGT[cid]);
				}
			}
		}
	}
	// handel hvisited
	// ����hVisited

	map<edgePair, edgeState, eSComparison>::iterator iterhV;
	for(iterhV=visitedState.begin(); iterhV!=visitedState.end(); ++iterhV) {
		if(iterhV->second.vState == hVisited) {//����ֻ�о��뵱ǰ�ڵ��������Լ����POI
			int nodei,nodej;
			if(iterhV->second.iDisToQuery == INFINITE_MAX) {
				nodei = iterhV->first.Nj;
				nodej = iterhV->first.Ni;
			} else {
				nodei = iterhV->first.Nj;
				nodej = iterhV->first.Ni;
			}
			//��ȡ���ϵ�POI
			AdjGrpAddr=getAdjListGrpAddr(nodei);
			getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
			for (int i=0; i<AdjListSize; i++){
				getVarE(ADJNODE_A,Ref(NewNodeID),AdjGrpAddr,i);
				//getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
				if(NewNodeID == nodej){
					getVarE(DIST_A,Ref(EdgeDist),AdjGrpAddr,i);
					getVarE(PTKEY_A,Ref(PtGrpKey),AdjGrpAddr,i);

					nOfEdgeExpd++;

					//cout<<"("<<Ni<<","<<Nj<<") EdgeDist:"<<EdgeDist;
					if(PtGrpKey==-1){
						//cout<<"No POI existed on Edge where Q located."<<endl;
					}
					else {
						getFixedF(SIZE_P,Ref(PtNumOnEdge),PtGrpKey);
						//cout<<" PtNum on Edge where Q located:"<<PtNumOnEdge<<endl;
						//Notice the order POI visited on edge
						for(int j=0; j<PtNumOnEdge; j++){
							poicnt++;
							getVarE(PT_DIST,Ref(PtDist),PtGrpKey,j);

							getVarE(PT_VCT,tempKwd,PtGrpKey,j);
							int num = sizeof(tempKwd)/sizeof(int);
							for(int loop=0; loop<sizeof(tempKwd); loop+sizeof(int)){
								int temp;
								memcpy(&temp,tempKwd+loop,sizeof(int));
								tempKwds.push_back(temp);
							}

							if(LcontianRIKwd(tempKwds,Q.kwd)){
								POI tmp;
								if(nodei<nodej){
									tmp.dist_toquery = iterhV->second.iDisToQuery+PtDist;
								} else {
									tmp.dist_toquery = iterhV->second.jDisToQuery+EdgeDist-PtDist;
								}
								getVarE(PT_ATTRIBUTE,tmp.attr,PtGrpKey,j);
								getVarE(PT_P,&tmp.poid,PtGrpKey,j);
						
								if(tmp.dist_toquery <= Q.distCnst) {
									if(isBeDominate(tmp, cS, Q)) {
									} else {
										cS.push_back(tmp);
									}
								}             
							}// endif
						}
                
					}// endelse
				}
			}
		}
	}

	for(int i=0; i<cS.size(); i++) {
		rS.push_back(cS[i].poid);
		//cout<<cS[i].poid<<" ";
	}
}

void EGTDAlgorithm(const QueryPoint &Q) {

}

void queryAlgorithm(int algorithmId, const QueryPoint &Q) {
	initialQuery(algorithmId);
	if(algorithmId == EA) { // EA
		EABasedAlgorithm(Q);
	} else if (algorithmId == EGBU){ // EGBU
		EGBUAlgorithm(Q);
	} else { // EGTD
		EGTDAlgorithm(Q);
	}
}

int main() {

	return 0;

}