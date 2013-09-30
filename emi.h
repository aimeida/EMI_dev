////#define DEBUG
#pragma once
#include <sys/time.h>
#include <math.h>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <list>
#include "FastGraphCluster.h"
#include "DebugFunc.h"
#include "Misc.h"
using namespace std;

class Pairmatch;

class Cluster
{
 public:
 Cluster(set <int> &v, size_t p_start=0, size_t p_end=0):nodes(v),p_start(p_start),p_end(p_end) {};
  set <int> nodes;
  size_t p_start, p_end;
};

class EdgeInfo
{
 public:
 EdgeInfo(float weight, size_t p0, size_t p1):weight(weight),p0(p0),p1(p1) {};
  float weight;
  size_t p0, p1;
};
class MisPair
{
 public:
 MisPair(float p_start, float p_end, int flag, float freq):p_start(p_start),p_end(p_end), flag(flag){ freqs.push_back(freq); };
  float p_start, p_end;
  int flag;
  vector<float> freqs;
};

/*
  class MisPair
  {
 public:
 MisPair(float p_start, float p_end, int flag):p_start(p_start),p_end(p_end), flag(flag){};
  float p_start, p_end;
  int flag;
  };
*/
 
long myclock();

