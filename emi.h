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
 Cluster(set <int> &v):nodes(v) {};
  set <int> nodes;
};

class EdgeInfo
{
 public:
 EdgeInfo(float weight, size_t p0, size_t p1):weight(weight),p0(p0),p1(p1) {};
  float weight;
  size_t p0, p1;
};

//// record IBD pairs that are missed by beagle in current windows 
//class MissingPair
//{
// public:
// MissingIBD(size_t p0, size_t p1, float length, int flag):p0(p0), p1(p1), length(length), flag(flag) {};
//  size_t p0, p1;
//  float length;
//  int flag;
//}

/*
struct EdgeInfo
{
  EdgeInfo(float weight, size_t p0, size_t p1){weight=weight; p0=p0; p1=p1;}
  float weight;
  size_t p0, p1;
};
*/

long myclock();

