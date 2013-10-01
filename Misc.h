#pragma once
#include <list>
#include <fstream>
#include <map>
#include <vector>
#include <sstream>
using namespace std;

class Pairmatch{
 public:
  Pairmatch(int i1, int i2, size_t p_start, size_t p_end, float weight, float pcm_start, float pcm_end);
  int i1, i2;
  size_t p_start, p_end;
  float weight, pcm_start, pcm_end;
};

class CmdOpt{
 public:
  string len_type, winsize_type;
  float window_size, min_weight, dist2Weight_a, dist2Weight_b;
  int iter_count;
};

float dist2Weight(float dist, float a, float b, float min_weight);
string intToString(int a);
int read_fam_file(ifstream &fam_list, map<string,int> &vertexNameMap, vector<string> &vertexName);
void read_beagle_input(ifstream &input_seg, list<Pairmatch * > &matches, map<string,int> &vertexNameMap, CmdOpt cmdopt);
