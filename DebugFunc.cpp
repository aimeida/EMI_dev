// small functions for debug
#include "DebugFunc.h"

string DebugFunc::intToString(int a)
{
  stringstream ss;
  ss << a ;
  return ss.str();
}

void DebugFunc::printVec(set<int> m)
{
  for (set <int>::iterator j=m.begin(); j!=m.end(); j++)  
    cout << *j << "\t";
  cout << endl;
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
