#define DEBUG
#define DOG
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
using namespace std;

class Pairmatch{
 public:
  Pairmatch(int i1, int i2, size_t p_start, size_t p_end, float weight, float pcm_start, float pcm_end);
  int i1, i2;
  size_t p_start, p_end;
  float weight, pcm_start, pcm_end;
};

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
 MisPair(float p_start, float p_end, int flag):p_start(p_start),p_end(p_end), flag(flag){};
  float p_start, p_end;
  int flag;
};

/* no polymorphism ??
class MisPair
{
 public:
 MisPair(size_t p_start, size_t p_end, int flag):p_start(p_start),p_end(p_end), flag(flag){};
  size_t p_start, p_end;
  int flag;
};
*/

long myclock();

