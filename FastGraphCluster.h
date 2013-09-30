#pragma once
#include <utility>
#include <string>
#include <numeric>
#include <fstream>
#include <iostream>
#include <vector>
#include <set>
#include <list>
#include <map>
#include <cstdlib>
#include <time.h>      
#include "emi.h"
using namespace std;

class FiboNode;    //forward declaration
class DegreeArray;
class Pairmatch; 
class Cluster;
class EdgeInfo;
class MisPair;
class FibonacciHeap; 

// node used in heuristic vertex adding
class HeapNode{
public:
	float wsum;	// how many weights are connected to current module
	int index;	// index of the vertex
	int esum;
	HeapNode( float wsum = 0 , int index = 0 , int esum = 0):wsum(wsum),index(index),esum(esum){}
	// the order here are totally reversed. we will get a max heap instead
	bool operator >(const HeapNode &right);
	bool operator >=(const HeapNode &right);
	bool operator <(const HeapNode &right);
	bool operator <=(const HeapNode &right);
};

class FastGraphCluster
{
 public:
  FastGraphCluster(float density, int lowersize,float lowerincrease, int m_nVertex);
  ~FastGraphCluster(void);
  friend class DebugFunc;
  void updateInput(list<Pairmatch * > &active_matches);
  void fastClusterCore(int seedn, float n_overhead, float freq_th, float window_size, float window_size_nfold, ofstream& fout1);
  map <int, EdgeInfo* > *m_neighbor;
  map <int, Cluster* > result_clst;
  //// IBD pairs that are missed by beagle in current windows, key pairs are ordered
  map <pair<int, int>, MisPair* > misPairs;
  
  int clst_topindex;
  int continuous_empty_wins; // max continuous_empty_win allowed 
  int numOfLeftPairs(); // currect window number of IBD pairs that are not used 
  void printAllClst(ofstream& fout, float window_size);
  void dissolve();
  void dissolve(vector< pair <int, int > > &delEdge, vector< pair <int, int > > &addEdge);
  
 protected:
	FiboNode **m_pHeapNode;	// we will use Fibonacci heap
	FiboNode **m_pHeapNode2; // use to reserve memory
	FiboNode **m_pHeapNode3; // use for randomization ???
	int m_nVertex,m_nLowerSize;	// lower size of the modules;
	DegreeArray *seedArray;	// extract min + decrease key
	float m_dLowDen,m_dLowerIncrease;	// lower density of the modules
	float *neighborWeightCnt;	// total sum of the weight for  neighbor
	string intToString(int a);	
	map<int, int > clstID;	

 private:
	int buildCore(int index, set<int> &result, FibonacciHeap &heap, set<int> &changed);
	void extendCore(set<int> &surround,set<int> &core_id, set<int> &node_id, int dn);
	float getWeight(float edgeweight,float vertexweight, float seedW=NULL);
	void deleteClst(int cid, Cluster * i_clst=NULL);
	void initMisFlag();  // set flag = 0
	void refine_bound(vector <float>::iterator &ri, float &pair_start, float &pair_end, float raw_start, float raw_end, vector <float> &freqs);
	void refine_bound(vector <float>::iterator &ri, float &pair_start, float &pair_end, float window_size, float raw_start, float raw_end, vector <float> &freqs, float freq_th, float window_size_nfold, bool verbose = NULL);
	void printMissing(float window_size, float minLen, float freq_th); // print current pairs in current window
};
