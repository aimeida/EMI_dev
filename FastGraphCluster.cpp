#include "DegreeArray.h"
#include "FibonacciHeap.h"
#include "FastGraphCluster.h"
#include "DebugFunc.h"
using namespace std;

#define UNEXPLORED -1
#define getmin(a,b) ((a)<(b)?(a):(b))

extern vector<string> vertexName;
extern float cur_pos_start;
extern float cur_pos;

FastGraphCluster::FastGraphCluster(float density,int lowersize,float lowerincrease, int m_nVertex):m_nLowerSize(lowersize),m_dLowDen(density),m_dLowerIncrease(lowerincrease), m_nVertex(m_nVertex), clst_topindex(0)
{
        int i;
	m_neighbor = new map<int, EdgeInfo* >[m_nVertex]; 
	neighborWeightCnt = new float[m_nVertex];
	m_pHeapNode = new FiboNode*[m_nVertex];
	m_pHeapNode2 = new FiboNode*[m_nVertex];
	for (i=0;i<m_nVertex;i++)
	  {
	    m_pHeapNode[i] = new FiboNode(HeapNode(UNEXPLORED,i,0));
	    m_pHeapNode2[i] = m_pHeapNode[i];
	    clstID[i] = -1;
	  }
}

FastGraphCluster::~FastGraphCluster(void)
{
  for (int i=0; i < m_nVertex; i++)
    {
      delete m_pHeapNode[i];
    }
    delete[] m_pHeapNode; // FiboNode **
    delete[] m_pHeapNode2;
    delete[] neighborWeightCnt;
    delete[] m_neighbor;
    for (map <int, Cluster* >::iterator ci=result_clst.begin(); ci!=result_clst.end(); ci++)
      delete ci->second;
}

void FastGraphCluster::deleteClst(int cid, Cluster * i_clst /*=NULL*/)
{
  if (i_clst==NULL) i_clst = result_clst[cid];  // sometimes too lazy to pass the pointer                  
  for (set<int>::iterator i =i_clst->nodes.begin(); i!=i_clst->nodes.end(); i++)
    clstID[*i] = -1;
  delete i_clst;
  result_clst.erase(cid);
}

template <class K, class V> 
void addKey(map <K,V> &m, K key)
{
  typename map <K,V>::iterator it;
  if ((it=m.find(key)) == m.end()) {
    m[key] = 1;
  } else it->second++;
}

void FastGraphCluster::dissolve(vector< pair <int, int > > &delEdge, vector< pair <int, int > > &addEdge)
{
  map<int, int > changed_clst; 
  Cluster* tmp_clst;
  int ci, ai;
  for (vector< pair <int, int > >::iterator i=delEdge.begin();i!=delEdge.end();i++)
    {
      // edges to delete within cluster
      ci=clstID[(*i).first];      
      if ( ci > -1 &&  ci==clstID[(*i).second])
	addKey(changed_clst, ci);
    }  
  if (!addEdge.empty()){
    for (vector< pair <int, int > >::iterator i=addEdge.begin();i!=addEdge.end();i++)
      {
	if ((ci=clstID[(*i).first]) > -1 )
	  addKey(changed_clst, ci);
	if ((ci=clstID[(*i).second]) > -1 )
	  addKey(changed_clst, ci);
      }  
  }
  if (changed_clst.empty()) return;
  for (map<int, int >::iterator cc = changed_clst.begin(); cc != changed_clst.end(); cc++)
    deleteClst(cc->first);
}

// print all Clusters in current window.
void FastGraphCluster::printAllClst(ofstream& fout)
{
  for (map <int, Cluster* >::iterator ci=result_clst.begin(); ci!=result_clst.end(); ci++)
    {   
      Cluster * cl = ci->second;
      fout << ci->first << "\t" << cur_pos_start << "\t" << cur_pos << "\t" << cl->nodes.size() << "\t";
      for (set <int>::iterator j=cl->nodes.begin(); j!=cl->nodes.end(); j++)
	fout << vertexName[*j] << "\t";
      fout << endl;
    }
}

void FastGraphCluster::updateInput(list<Pairmatch * > &active_matches)
{
  int i;
  float *nw = neighborWeightCnt;
  FiboNode *mp2;
  for (i=0;i < m_nVertex;i++)
    {
      mp2=m_pHeapNode2[i];         // use mp2 to reserve the space, alloc/dealloc expensive
      mp2->key.wsum = UNEXPLORED;
      mp2->key.esum = 0;
      if (clstID[i] > -1){  // already clustered
	m_pHeapNode[i] = NULL;
      } else m_pHeapNode[i] = mp2;
      *nw++ = 0;
    }
  list< Pairmatch * >::iterator am;  
  for (am = active_matches.begin(); am != active_matches.end(); am++)
    {
      if (clstID[(*am)->i1] > -1 || clstID[(*am)->i2] > -1) continue; 
      neighborWeightCnt[(*am)->i1] += (*am)->weight;
      neighborWeightCnt[(*am)->i2] += (*am)->weight;
    }
  float maxWeightDegree = 0;
  for (i=0;i<m_nVertex;i++){
    if( neighborWeightCnt[i] > maxWeightDegree ) maxWeightDegree = neighborWeightCnt[i];
  }
  seedArray = new DegreeArray(neighborWeightCnt,m_nVertex,maxWeightDegree);
}

void FastGraphCluster::fastClusterCore()
{
  int i=0, seedn=5;      
  set<int> result, changed;
  while (!seedArray->empty() && seedArray->top >= m_nLowerSize-1) {
    FibonacciHeap heap;	// local expanding heap
    i = seedArray->getMax();
    expandCore(i, result, heap, changed);
    //cout << cur_pos << "\t" << result.size() << endl;
    //cout << cur_pos << "\t" << result.size() << " " << clstID[*result.begin()]<< endl;

    if (heap.m_nNode > 0){ // there are nodes left after the core
      set< set <int> > core_ext; 
      set<int> surround;
      srand(time(NULL));
      for (int dn=0; dn<seedn; dn++){
	extendCore(surround, result, heap.current_node_id, dn); 
	core_ext.insert(surround);
	surround.clear();
      }
//    if(core_ext.size() > 1 ){
//	cout << "variation " << core_ext.size() << endl;
//	for (set< set <int> >::iterator i=core_ext.begin(); i!=core_ext.end(); i++)
//	  DebugFunc::printVec(*i); 
//      }    
    }

    heap.current_node_id.clear();
    result.clear();
    for ( set<int>::iterator iter = changed.begin() ; iter != changed.end() ; iter++ ) {
      i = *iter;
      if ( m_pHeapNode[i] != NULL ){
	m_pHeapNode[i]->key.wsum = UNEXPLORED;
	m_pHeapNode[i]->key.esum = 0;
      }
    }
    changed.clear();

    // fix me, when to delete heap ??
  }
  if (!seedArray->empty())
    delete seedArray;
}

// used to select the seed node pair
float FastGraphCluster::getWeight(float edgeWeight,float vertexWeight, float seedW/*=NULL*/)
{
  int nbin=10;      // 0.8, 0.9, 1.0 have different values
  //return m_nVertex * edgeWeight + vertexWeight; // similar to dash
  return  m_nVertex * (int)(nbin * edgeWeight) + vertexWeight; 
}

// create core clusters with density==1
int FastGraphCluster::expandCore(int index, set<int> &result, FibonacciHeap &heap, set<int> &changed)
{
  int tEdge=0, imin, i,j,w3, nVertex = 0;
  float density,increase,w,oldWeight;
  FiboNode *ptr = m_pHeapNode[index];
  ptr->key.wsum = 0;
  ptr->key.esum = 0;
  heap.insert(ptr);
  int maxindex = -1; // default 0 may cause problem
  float maxWeight = 0;
  map<int, EdgeInfo*>::iterator imap;    
  for (imap = m_neighbor[index].begin(); imap!=m_neighbor[index].end(); imap++)
    {
      j = imap->first;
      if (m_pHeapNode[j] == NULL) continue;  // already in other clusters
      w = getWeight((imap->second)->weight,neighborWeightCnt[j]);
      if (w > maxWeight ) {
	maxindex = j;
	maxWeight = (float)w;  // it was set to be very large 
      }
    }
  if (maxindex < 0) {   
    seedArray->remove(index,neighborWeightCnt[index]);
    return 0;
  }
  
  // for the expanding procedure, definitely select this edge first
  m_neighbor[index][maxindex]->weight += maxWeight;
  
  while (!heap.empty())
    {
      ptr = heap.extractMin();
      imin = ptr->key.index;
      w3 = ptr->key.esum;      
      tEdge += w3;
      nVertex ++;
        
      if (nVertex > 1) {
  	if ( 2.0*tEdge != (nVertex*(nVertex-1))) {  // fully connected cluster
	  nVertex--;
	  tEdge -= w3;
	  break;	
	}
      }
        
      m_pHeapNode[imin] = NULL;
      seedArray->remove(imin,neighborWeightCnt[imin]);
      result.insert(imin);	// write back result
      
      for (imap = m_neighbor[imin].begin(); imap!=m_neighbor[imin].end(); imap++){
	j = imap->first;
	//cout << "size" << m_neighbor[imin].size() << " " << j << endl;
        w = (imap->second)->weight;
        ptr = m_pHeapNode[j];
	if (ptr==NULL) continue;
	//cout << "UNEXP " << ptr->key.index << endl;
	if (ptr->key.wsum == UNEXPLORED){
	  ptr->key.wsum = w;
	  ptr->key.esum = 1;
	  heap.insert(ptr);
	}else {
	  //cout << "ptr "<< ptr->key.index << " " << ptr->key.wsum << " " << ptr->key.esum  << endl;
	  heap.decrease(ptr,HeapNode( ptr->key.wsum+w,j,ptr->key.esum+1));
	}
          
	oldWeight = neighborWeightCnt[j];        
	if ( ( nVertex == 1 ) && (j == maxindex)) {
	  seedArray->decrease(j,oldWeight,oldWeight + maxWeight - w);
	  neighborWeightCnt[j] -= (w-maxWeight);
	}else {
	  seedArray->decrease(j,oldWeight,oldWeight - w);
	  neighborWeightCnt[j] -= w;
	}
	changed.insert(j);	
      }
    }
  m_neighbor[index][maxindex]->weight -= maxWeight;

  if (result.size() >= m_nLowerSize){
    for (set<int>::iterator ni=result.begin();ni!=result.end();ni++)
      clstID[*ni] = clst_topindex;
    result_clst[clst_topindex++] = new Cluster(result);
  }
  return 1;
}

// add some random variance in the weight 
float addVar(float v0)
{
  float r = 1-(float)(rand()%100)*0.002; // [0,99]==>[0.8,1]
  ////return v0 * 0.9; // debug
  if (v0==0) {
    cerr << "can't be 0 " << endl;
    exit(0);
  }
  return v0 * r;
}

// build local heap from initial nodes, with random variance added
// core nodes in in "core_id", and neighbors of core nodes are "node_id"
void FastGraphCluster::extendCore(set<int> &surround, set<int> &core_id, set<int> &node_id, int dn)
{
  FibonacciHeap heapExt; // Extended heap, start with start_id
  FiboNode *ptr; 
  set<int> current_heap;
  for (set<int>::iterator i = node_id.begin();i!=node_id.end();i++){
    ptr = m_pHeapNode[*i];
    //// fixme: use esum for now, change to wsum later if I'm not lazy
    //cout << "key.wsum " << ptr->key.wsum << " " << ptr->key.esum << " " << *i << endl;
    ptr->key.wsum = addVar((float)ptr->key.esum);
    heapExt.insert(ptr);
    current_heap.insert(*i);
  }
  int imin, es, tEdge, nVertex;
  nVertex = core_id.size();
  tEdge = nVertex * (nVertex-1)/2;
  float density, increase; 
  while (!heapExt.empty()){
    ptr = heapExt.extractMin();
    imin = ptr->key.index;
    es = ptr->key.esum; 
    tEdge += es;
    nVertex ++;
    ////cout << "### " << nVertex << " " << tEdge << endl;
    if (nVertex > 1) {
      density = 2.0*tEdge/(nVertex*(nVertex-1));
      increase = es/(density*(nVertex-1)); 
      //// cout << "### " << density << " " << increase << endl;
      // if ( density < m_dLowDen || increase < m_dLowerIncrease){
      if ( density < m_dLowDen ) {
	nVertex--;
        tEdge -= es;
        break;
      }
    }
    surround.insert(imin);  
    for (map<int, EdgeInfo*>::iterator imap = m_neighbor[imin].begin(); imap!=m_neighbor[imin].end(); imap++){
      int j = imap->first;
      // fixme, should i change m_neighbor[imin]? namely, delete all memebers in the core? faster?
      ptr = m_pHeapNode[j];  // only used to keep wsum(UNEXPLORED) and esum(0)
      if ( ptr==NULL) continue;
      if ( surround.find(j)!=surround.end() || node_id.find(j)!=node_id.end()) continue; 
      float w = addVar((imap->second)->weight);
      if ( current_heap.find(j)==current_heap.end()){
	ptr->key.wsum = w;
	ptr->key.esum = 1;
	heapExt.insert(ptr);
	current_heap.insert(j);
      }else{
	heapExt.decrease(ptr,HeapNode( ptr->key.wsum+w,j,ptr->key.esum+1));
      }
    }
  }
  for (set<int>::iterator ih = current_heap.begin(); ih != current_heap.end(); ih++){
    if (node_id.find(*ih)==node_id.end()){  // esum of left-over nodes remain unchanged.
      m_pHeapNode[*ih]->key.wsum = UNEXPLORED;
      m_pHeapNode[*ih]->key.esum = 0;
    }
  }
}

bool HeapNode::operator >(const HeapNode &right)
{
	return wsum<right.wsum;
}

bool HeapNode::operator <(const HeapNode &right)
{
	return wsum>right.wsum;
}

bool HeapNode::operator <=(const HeapNode &right)
{
	return wsum >= right.wsum;
}

bool HeapNode::operator >=(const HeapNode &right)
{
	return wsum <= right.wsum;
}
