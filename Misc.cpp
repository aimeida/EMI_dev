// read pairwise data 
#include "Misc.h"

Pairmatch::Pairmatch(int i1, int i2, size_t p_start, size_t p_end, float weight, float pcm_start, float pcm_end):i1(i1),i2(i2),p_start(p_start),p_end(p_end),weight(weight),pcm_start(pcm_start), pcm_end(pcm_end){}

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

int read_fam_file(ifstream &fam_list, map<string,int> &vertexNameMap, vector<string> &vertexName){
  stringstream ss;
  string line, iid, field1, field2, field3, field4;
  int m_nVertex=0;  
  map<int, int> clstMap;
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
  return m_nVertex;
}

void read_beagle_input(ifstream &input_seg, list<Pairmatch * > &matches, map<string,int> &vertexNameMap, CmdOpt cmdopt){
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
      
      if (cmdopt.len_type == "7th" ) {sw_w = logOR; }
      else if (cmdopt.len_type == "cM" ) {sw_w = cm_end - cm_start; }
      else if (cmdopt.len_type == "bp" ) {sw_w = pos2 - pos1;}

      if (cmdopt.winsize_type == "cM"){  
	if ( cm_end - cm_start  < cmdopt.window_size) continue;		
	cur_match = new Pairmatch(ia, ib, pos1, pos2, dist2Weight(sw_w, cmdopt.dist2Weight_a, cmdopt.dist2Weight_b, cmdopt.min_weight), cm_start, cm_end);      
      }
      else if (cmdopt.winsize_type == "bp") {
	if ( pos2 - pos1 < cmdopt.window_size ) continue;
	cur_match = new Pairmatch(ia, ib, pos1, pos2, dist2Weight(sw_w, cmdopt.dist2Weight_a, cmdopt.dist2Weight_b, cmdopt.min_weight), (float)pos1, (float) pos2);
	
      }
      matches.push_back( cur_match);
    }     
  input_seg.close();
}

