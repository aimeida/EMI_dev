// read pairwise data 
#include "Misc.h"

Pairmatch::Pairmatch(int i1, int i2, float weight):i1(i1),i2(i2),weight(weight),pcm_start(-1), pcm_end(-1){} // -1 for unknow position

Pairmatch::Pairmatch(int i1, int i2, float weight, float pcm_start, float pcm_end):i1(i1),i2(i2),weight(weight),pcm_start(pcm_start), pcm_end(pcm_end){}

float dist2Weight(float dist, float a, float b, float min_weight){
  float t = a * dist + b;
  if ( t < min_weight)
    return min_weight;
  //cerr << "dist2w" << (t < 1 ? t : 1.0) << endl;
  return t < 1 ? t : 1.0;
}

string intToString(int a){
  stringstream ss;
  ss << a ;
  return ss.str();
}

bool read_fam_file(string fam_file, int &m_nVertex, map<string,int> &vertexNameMap, vector<string> &vertexName){
  ifstream fam_list( fam_file.c_str());
  if (!fam_list) 
    return false;
  stringstream ss;
  string line, iid, field1, field2, field3, field4;
  map<int, int> clstMap;
  m_nVertex = 0;
  while ( getline( fam_list , line ) ) {
    ss.clear(); ss.str( line );
    ss >> field1 >> field2;
    iid = field1+" "+field2+".0";
    clstMap[m_nVertex] = -1;
    vertexNameMap[iid] = m_nVertex++;
    vertexName.push_back(iid);
    
    iid = field1+" "+field2+".1";
    clstMap[m_nVertex] = -1;
    vertexNameMap[iid] = m_nVertex++;
    vertexName.push_back(iid);
  }
  fam_list.close();
  return true;
}

bool compare_pairs(Pairmatch *pa, Pairmatch *pb){
  if ( pa->pcm_start < pb->pcm_start) return true;
  else if ( pa->pcm_start > pb->pcm_start) return false;
  else if ( pa->pcm_end < pb->pcm_end) return true;
  else if ( pa->pcm_end > pb->pcm_end) return false;
  else return pa < pb;
}

bool read_beagle_input(string input_file, list<Pairmatch * > &matches, map<string,int> &vertexNameMap, CmdOpt &cmdopt){
  ifstream input_seg( input_file.c_str() );
  if (!input_seg) return false;
  stringstream ss;
  Pairmatch *cur_match;
  string line, field1, field2, field3, field4;
  float logOR, cm_start, cm_end, sw_w;
  size_t pos1, pos2;
  int ia, ib;
  while( getline( input_seg , line ) )
    {
      ss.clear(); ss.str( line );
      ss >> field1 >> field2 >> field3 >> field4 >> pos1 >> pos2 >> logOR >> cm_start >> cm_end;

      ia = vertexNameMap[field1+" "+field2];
      ib = vertexNameMap[field3+" "+field4];      
      if (ia > ib) swap(ia, ib); // such that it's ordered

      if (cmdopt.len_type == "7th" ) {sw_w = logOR; }
      else if (cmdopt.len_type == "cM" ) {sw_w = cm_end - cm_start; }
      else if (cmdopt.len_type == "bp" ) {sw_w = pos2 - pos1;}

      if (cmdopt.winsize_type == "cM"){  
	if ( cm_end - cm_start  < cmdopt.window_size) continue;		
	cur_match = new Pairmatch(ia, ib, dist2Weight(sw_w, cmdopt.dist2Weight_a, cmdopt.dist2Weight_b, cmdopt.min_weight), cm_start, cm_end);      
      }
      else if (cmdopt.winsize_type == "bp") {
	if ( pos2 - pos1 < cmdopt.window_size ) continue;
	cur_match = new Pairmatch(ia, ib, dist2Weight(sw_w, cmdopt.dist2Weight_a, cmdopt.dist2Weight_b, cmdopt.min_weight), (float)pos1, (float) pos2);
	
      }
      matches.push_back( cur_match);
    }     
  input_seg.close();
  return true;
}

bool read_emi_input_byWin(string input_file, list<Pairmatch * > &matches, float pos_end, float sw_w){
  ifstream input_seg( input_file.c_str() );
  if (!input_seg) return false;
  stringstream ss;
  Pairmatch *cur_match;
  int ia, ib;
  float p0, p1;
  string line;
  while( getline( input_seg , line ) ){
    ss.clear(); ss.str( line );
    ss >> ia >> ib >> p0 >> p1;
    if (ia > ib) swap(ia, ib); // such that it's ordered
    if ( p0 < pos_end - 0.001 && p1 >= pos_end) {
      cur_match = new Pairmatch(ia, ib, sw_w);
      matches.push_back(cur_match);
    }
  }
  input_seg.close();
  return true;
}

bool read_emi_input(string input_file, list<Pairmatch * > &matches, float sw_w){
  ifstream input_seg( input_file.c_str() );
  if (!input_seg) return false;
  stringstream ss;
  Pairmatch *cur_match;
  int ia, ib;
  float p0, p1;
  string line;
  while( getline( input_seg , line ) ){
    ss.clear(); ss.str( line );
    ss >> ia >> ib >> p0 >> p1;
    if (ia > ib) swap(ia, ib); // such that it's ordered
    cur_match = new Pairmatch(ia, ib, sw_w, p0, p1);
    matches.push_back(cur_match);
  }
  input_seg.close();
  return true;
}

