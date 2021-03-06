#pragma once
#include <map>
#include <string>
#include <sstream>
#include "FastGraphCluster.h"
#include "emi.h"
using namespace std;

class FastGraphCluster;
class Pairmatch;

struct ComparePair
{
  bool operator() (const pair<int, int> a, const pair<int, int> b) const 
  {return a.second >= b.second;}
};

class DebugFunc
{
 public:
  DebugFunc();
  set<float> pos_to_print;
  vector<int> iter_to_print;
  // static,  call func without an obj
  static void printNeighbor(int index, FastGraphCluster &f);
  static int numOfLeftPairs(FastGraphCluster &f); 
  static void mapToVec(map<int, int> &nodeCount, set< pair<int, int>, ComparePair > &sortCount);
  static void printVec(set<int> &m);
  static void printVec(vector<float> &m);
  static int countTransition(list<Pairmatch * > matches, float WINDOW_SIZE, float &p01, float &p10); // not by reference
};

