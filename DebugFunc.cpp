// small functions for debug
#include "DebugFunc.h"

int DebugFunc::countTransition(list<Pairmatch * > matches, float WINDOW_SIZE, float &p01, float &p10) 
{
  list<Pairmatch * > active_matches; 
  list< Pairmatch * >::iterator pm_iter, am;  
  float cur_pos = 0,cur_pos_start = 0;
  float min_end;
  int nwin=0;
  size_t n01=0, n10=0, n_all=0;  // Num(non->IBD), Num(IBD->non)
  vector< pair <int, int > > delEdge;
  vector< pair <int, int > > addEdge;
  //cerr << "start0 " << matches.size() << endl;
  pm_iter = matches.begin();
  while (pm_iter != matches.end())
    {
      while ( cur_pos < (*pm_iter)->pcm_start + WINDOW_SIZE ) cur_pos += WINDOW_SIZE;
      //cerr << cur_pos << "\t" << active_matches.size() << "\t"; 
      n_all += active_matches.size();
      for (am = active_matches.begin(); am != active_matches.end(); )
	{
	  if ( (*am)->pcm_end < cur_pos )
	    {
	      delEdge.push_back(make_pair((*am)->i1,(*am)->i2));
	      active_matches.erase( am++ );
	    } else am++;
	}
      while ( pm_iter!= matches.end() && (*pm_iter)->pcm_start <= (cur_pos - WINDOW_SIZE) )
	{
	  if ( (*pm_iter)->pcm_end >= cur_pos )
	    {
	      active_matches.push_back( *pm_iter );
	      addEdge.push_back(make_pair((*pm_iter)->i1,(*pm_iter)->i2));
	    } 
	  matches.erase( pm_iter++ );
	}
      //cerr << delEdge.size() << "\t" << addEdge.size() << endl;
      n10 += delEdge.size();
      n01 += addEdge.size();
      delEdge.clear();
      addEdge.clear();
      cur_pos_start = cur_pos;
      nwin += 1;
    }
  while(!active_matches.empty())
    {
      min_end = 0;
      for (am = active_matches.begin(); am != active_matches.end(); am++ )
	if ( am == active_matches.begin() || (*am)->pcm_end < min_end ) min_end = (*am)->pcm_end;
      while ( min_end >= cur_pos ) cur_pos += WINDOW_SIZE;
      for (am = active_matches.begin(); am != active_matches.end(); )
	{
	  if ( (*am)->pcm_end < cur_pos )   
	    {
	      delEdge.push_back(make_pair((*am)->i1,(*am)->i2));
	      active_matches.erase( am++ );
	    } else am++;
	}
      delEdge.clear();
      cur_pos_start = cur_pos;
      nwin += 1;
    }
  //cerr << "end0 " <<  matches.size() << endl;
  //cerr << n_all << "\t" << n10 << "\t" << n01 << endl;
  p01 = n01/(float)n_all;
  p10 = n10/(float)n_all;
  return nwin;
}

void DebugFunc::printVec(set<int> &m)
{
  for (set <int>::iterator j=m.begin(); j!=m.end(); j++)  
    cerr << *j << "\t";
  cerr << endl;
}

void DebugFunc::printVec(vector<float> &m)
{
  for (vector <float>::iterator j=m.begin(); j!=m.end(); j++)  
    cerr << *j << "\t";
  cerr << endl;
}


void DebugFunc::printNeighbor(int index, FastGraphCluster &f)
{
  map<int, EdgeInfo*>::iterator imap;
  for (imap = f.m_neighbor[index].begin(); imap!= f.m_neighbor[index].end(); imap++)
    if ((imap->second)->weight > 0) cout << imap->first << " ";
  cout << endl;
}

int DebugFunc::numOfLeftPairs(FastGraphCluster &f)
{
  size_t n_used=0, clst_size;  //density=1                                                          
  //  cout << "$$$$"  << result_clst.size() <<  endl;                                               
  for (map <int, Cluster* >::iterator ci=f.result_clst.begin(); ci!=f.result_clst.end(); ci++){
    clst_size = (ci->second)->nodes.size();
    n_used += clst_size*(clst_size-1)/2;
    //cout << "clst_size "<< clst_size << endl;                                                     
  }
  return n_used;
}

void DebugFunc::mapToVec(map<int, int> &nodeCount, set< pair<int, int>, ComparePair > &sortCount)
{
  for (map <int, int>::iterator ni=nodeCount.begin();ni!=nodeCount.end();ni++)
    sortCount.insert(make_pair(ni->second, ni->first)); 
}
