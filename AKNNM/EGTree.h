#ifndef _EGTREE_H_
#define _EGTREE_H_
//#pragma once
#include<stdio.h>
#include<vector>
#include<stdlib.h>
#include<memory.h>
#include<unordered_map>
#include<map>
#include<set>
#include<deque>
#include<stack>
#include<algorithm>
#include "metis.h"
#include <random>
#include <bitset>
#include "netshare.h"
#include "utility.h"
#include "KeywordsGenerator.h"
//#include<sys/time.h>
using namespace std;

/********************************PreDefinition************************************/
// MACRO for timing----Linux
/*
struct timeval tv;
long long ts, te;
#define TIME_TICK_START gettimeofday( &tv, NULL ); ts = tv.tv_sec * 100000 + tv.tv_usec / 10;
#define TIME_TICK_END gettimeofday( &tv, NULL ); te = tv.tv_sec * 100000 + tv.tv_usec / 10;
#define TIME_TICK_PRINT(T) printf("%s RESULT: %lld (0.01MS)\r\n", (#T), te - ts );
*/
// test in windows

//long long ts, te;
//#define TIME_TICK_START clock_t ts=clock();
//#define TIME_TICK_END clock_t te=clock();
//#define TIME_TICK_PRINT(T) printf("%s RESULT: %lld (0.01MS)\r\n", (#T), te - ts );

// offset 
#define _FILE_OFFSET_BITS 64
// set all edge weight to 1(unweighted graph)
#define ADJWEIGHT_SET_TO_ALL_ONE true
// we assume edge weight is integer, thus (input edge) * WEIGHT_INFLATE_FACTOR = (our edge weight)
#define WEIGHT_INFLATE_FACTOR 100000
// egtree fanout
#define PARTITION_PART 4
// egtree leaf node capacity = tau(in paper)
#define LEAF_CAP 32
#define ATTRIBUTE_DIMENSION 6
#define SKY_PARTITION 5
#define edDis 0.5
#define covThre 10
#define splitBlock 5

/********************************DataStructure************************************/
//-----basic infromation of gtree----------------------
typedef struct {
	double x, y;
	vector<int> adjnodes;
	vector<int> adjweight;
	bool isborder;
	vector<int> egtreepath; // this is used to do sub-graph locating
}Node; 


typedef struct {
	//--------------basic structure of gtree
	vector<int> borders; // the borders of this node
	vector<int> children;
	bool isleaf;
	vector<int> leafnodes;
	int father;
	// ----- min dis -----
	vector<int> union_borders; // for non leaf node, merging borders of children
	vector<float> mind; // min dis, row by row of union_borders 
	// --------------extend structure of egtree for Top-down
	unsigned long long unionkwds; // identical to keyword bitmap
	set<unsigned long long> coverkwds; // the cluster of keywords to cover all points
}TreeNode;

int noe; // number of edges
int nLeafNode;
vector<Node> Nodes;
extern vector<TreeNode> EGTree;

// init status struct
typedef struct {
	int tnid; // tree node id
	set<int> nset; // node set
}Status;

// use for metis
// idx_t = int64_t / real_t = double
idx_t nvtxs; // |vertices|
idx_t ncon; // number of weight per vertex 
idx_t* xadj; // array of adjacency of indices
idx_t* adjncy; // array of adjacency nodes
idx_t* vwgt; // array of weight of nodes
idx_t* adjwgt; // array of weight of edges in adjncy
idx_t nparts; // number of parts to partition
idx_t objval; // edge cut for partitioning solution
idx_t* part; // array of partition vector
idx_t options[METIS_NOPTIONS]; // option array

/********************************Function************************************/
//-------------basic function of gtree
void options_setting(); // METIS setting options
void init_input(int nOfNode, EdgeMapType EdgeMap); // load nodes and edges 
void data_transform_init(set<int> &nset); // transform original data format to that suitable for METIS
void init(int nOfNode, EdgeMapType EdgeMap); // combining two steps above
void finalize(); // free space
unordered_map<int, int> graph_partition(set<int> &nset); // graph partition
void build(EdgeMapType EdgeMap); // egtree construction
void egtree_save(string filename); // dump gtree index to file
void egtree_load(string filename, vector<TreeNode>& EGTree);// load gtree index from file
vector<int> dijkstra_candidate(int s, vector<int> &cands, vector<Node> &graph); // dijkstra search, used for single-source shortest path search WITHIN one gtree leaf node!
void hierarchy_shortest_path_calculation(); // calculate the distance matrix
int mainEgtree(int nOfNode, EdgeMapType EdgeMap); // main function
bool sortBySize(const unsigned long long& left, const unsigned long long& right); //sort by InterNode kwd size
void handleCovKwd(int tn, set<unsigned long long> nodecoverkwds);
#endif