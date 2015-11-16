#ifndef __NET_SHARED
#define __NET_SHARED

#include <queue>
#include <bitset>
#include <vector>
#include <set>
using namespace std;
// handle
#define MAX_DIST 9999999
#define ATTRIBUTE_DIMENSION 6
int MAX_KEYWORDS = 64;

int NodeNum;
int EdgeNum;
// define the query point

struct QueryPoint
{
	int Ni, Nj;
	float dist_Ni;
	float distEdge;
};

struct Query {
	int k;
	unsigned long long keywords;
	vector<QueryPoint> querypts;
};

ostream& operator<<(ostream& os, const Query& Q)
{
	os << bitset<64>(Q.keywords).to_string() << " "<< Q.k;
	set<QueryPoint> ::iterator iter = Q.querypts.begin();
	for (; iter != Q.querypts.end(); iter++) {
		os << "," << iter->Ni << " " << iter->Nj << " " << iter->dist_Ni;
	}
	return os;
}

struct POI
{
	int Ni, Nj;
	int pre;//Refine use only
	unsigned long long keywords;
	float dist_Ni;
	float dist_Nj;
	float dist_toquery;

};


struct POIComparison
{
	bool operator () (const POI& left, const POI& right) const
	{
		return left.dist_toquery> right.dist_toquery;
	}
};


ostream& operator<<(ostream& os, const POI& poi)
{
	os << "(" << poi.Ni << "," << poi.Nj << ") ";
	os << bitset<64>(poi.keywords).to_string();
	os << " distNi:" << poi.dist_Ni;
	os << " distQ:" << poi.dist_toquery;
	os << " preNode:" << poi.pre << endl;

	return os;
}

struct DStepEvent
{
    double dist;
    int node;
};

struct DStepComparison
{
    bool operator () (const DStepEvent& left, const DStepEvent& right) const
    {
        return left.dist > right.dist;
    }
};

typedef	priority_queue<DStepEvent,vector<DStepEvent>,DStepComparison> DStepQueue;

// for Dijkstra
struct edgePair
{
    int Ni;
    int Nj;
};

struct edgeState 
{
	int vState;
	float iDisToQuery;
	float jDisToQuery;
};

struct eSComparison
{
    bool operator () (const edgePair& left, const edgePair& right) const
    {
		//return left.disToQuery > right.disToQuery;
		if(left.Ni == right.Ni && left.Nj == right.Nj) return true;
		return false;
    }
};

struct dijkVisit 
{
	int N;
	float disTQ;
};
//?????�ǲ��ǰ�����С˳������
struct dVComparison
{
    bool operator () (const dijkVisit& left, const dijkVisit& right) const
    {
		//return left.disToQuery > right.disToQuery;
		return left.disTQ > right.disTQ;
    }
};
typedef	priority_queue<dijkVisit ,vector<dijkVisit>,dVComparison> dVQueue;


struct point
{
    int Ni,Nj,pos;
};

//Modified by Qin Xu
//Denote POI on RoadEdge

// record the information of POIcc
struct InerNode
{
	int pid;
	float dis;
	unsigned long long vct;//Vector of keywords denoted by 64-bit
};

struct edge
{
	int FirstRow;
	int Ni, Nj;
	float dist;
	unsigned long long sumkwds;
	FastArray<InerNode> pts;
	//float dLB,dUB;		// for gendata use only
};

typedef map<int,edge*> EdgeMapType; // map node to edges

// build AdjList on fly
FastArray<int>* AdjList;
FastArray<point> PtList;
EdgeMapType EdgeMap;	// key: i0*NodeNum+j0

// get the key of edge
inline int getKey(int i,int j)
{
    int i0=i<j?i:j,j0=j<i?i:j;	// be careful of the order in other places
    return (i0*NodeNum+j0); // map this edge to unique key
}
// break the index of key of edge
inline void breakKey(int key,int& Ni,int& Nj)
{
    Ni=key/NodeNum;
    Nj=key%NodeNum;
}

void printEdgeKeys()
{
    int Ni,Nj;
    printf("EdgeKeys:\n");
    EdgeMapType::iterator p=EdgeMap.begin();
    while (p!=EdgeMap.end())
    {
        edge* e=p->second;
        breakKey(p->first,Ni,Nj);
        printf("%d %d %d\n",Ni,Nj,(e==NULL));
        p++;
    }
}

#endif // __NET_SHARED

