#pragma once
#include <list>
#include <fstream>
#include <cstdlib>
#include <map>
#include <vector>
#include <sstream>
using namespace std;

class Pairmatch{
 public:
  Pairmatch(int i1, int i2, float weight); // allowed 
  Pairmatch(int i1, int i2, float weight, float pcm_start, float pcm_end);
  int i1, i2;
  //size_t p_start, p_end;
  float weight, pcm_start, pcm_end;
};

class EdgeInfo
{
 public:
 EdgeInfo(float weight):weight(weight) {};
  //EdgeInfo(float weight, size_t p0, size_t p1):weight(weight),p0(p0),p1(p1) {};
  float weight;
  //size_t p0, p1;
};


class CmdOpt{
 public:
  string len_type, winsize_type;
  float window_size, min_weight, dist2Weight_a, dist2Weight_b, window_size_nfold, emi_weight;
  int iter_count, continuous_empty_wins;
};

float dist2Weight(float dist, float a, float b, float min_weight);
string intToString(int a);
bool read_fam_file(string fam_file, int &m_nVertex, map<string,int> &vertexNameMap, vector<string> &vertexName);
bool read_beagle_input(string input_file, list<Pairmatch * > &matches, map<string,int> &vertexNameMap, CmdOpt &cmdopt);
bool read_emi_input(string input_file, list<Pairmatch * > &matches, float pos_end, float sw_w);
