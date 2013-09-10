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

void FastGraphCluster::initMisFlag()
{
  map <pair<int, int>, MisPair*>::iterator i;
  for (i=misPairs.begin(); i!=misPairs.end(); i++)
    (i->second)->flag = 0;
}

void FastGraphCluster::printMissing(float minLen)
{
  // don't erase while iterating, may cause problem
  //cerr << "m3\t" << cur_pos << "\t" << misPairs.size() << endl;
  vector <map <pair<int, int>, MisPair* >::iterator> rms;
  for (map <pair<int, int>, MisPair* >::iterator i=misPairs.begin(); i!=misPairs.end(); i++){
    if ((i->second)->flag == 0){ 
      if ((i->second)->p_end - (i->second)->p_start >= minLen){ 
	cout << vertexName[(i->first).first] << "\t" << vertexName[(i->first).second] << "\t";
	cout << (i->second)->p_start << "\t" << (i->second)->p_end << endl; 
      }
      rms.push_back(i);
    }
  }
  vector <map <pair<int, int>, MisPair* >::iterator>::iterator ri;
  for (ri=rms.begin(); ri!=rms.end();ri++){
    delete (*ri)->second;
    misPairs.erase(*ri);
  }
}

template <class K, class V> 
void addKey(map <K,V> &m, K key, int cnt=1)
{
  typename map <K,V>::iterator it;
  if ((it=m.find(key)) == m.end()) {
    m[key] = cnt;
  } else (it->second) += cnt;
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

void FastGraphCluster::fastClusterCore(int seedn, float freq_th, float minLen, ofstream& fout1)
{
  int i=0;
  set<int> result, changed;
  int freq_thn = (int) seedn * freq_th;
  //cerr << "freq_thn " << freq_thn << endl;
  while (!seedArray->empty() && seedArray->top >= m_nLowerSize-1) {
    FibonacciHeap heap;	// local expanding heap
    i = seedArray->getMax();   
    buildCore(i, result, heap, changed); // update clst_topindex, added to result_clst
    // make sure  with\out extension, core size stays the same.
    fout1 << cur_pos << "\t" << result.size() << endl; 

    if (heap.m_nNode > 0){ // there are nodes left after the core
      map< set <int>, int > core_ext;  // count number of times it appears
      set<int> surround;
      srand(time(NULL));
      //srand(0);
      for (int dn=0; dn<seedn; dn++){
	extendCore(surround, result, heap.current_node_id, dn); 
	if (!surround.empty())
	  addKey(core_ext, surround);
	surround.clear();
      }
      
      if(!core_ext.empty()){
    map <pair<int, int>, MisPair* >::iterator mi;
    map<pair<int,int>, int> pairCount;
    set<int>::iterator ai,ci,bi;
    int id1, id2;
	//cerr << "size1 "<< core_ext.size() << " " << (*core_ext.begin()).size()<< endl;
   	// step 1: missing edges between core_ext and core (look at this first, use freq_thn as filter)
    map<int, int> ext_node;
    set<int> ext_node_set;
    for (map< set <int>, int >::iterator i=core_ext.begin(); i!=core_ext.end(); i++){
       for (ai=(i->first).begin();ai!=(i->first).end();ai++)
          addKey(ext_node, *ai, i->second);
    }
    for (map<int, int>::iterator i = ext_node.begin(); i!=ext_node.end(); i++){
        if (i->second < freq_thn) continue;
        ext_node_set.insert(i->first);   // save time for step 2? while adding pairs
        for (set <int>::iterator j = result.begin(); j!=result.end(); j++){
           if (i->first < *j){
                addKey(pairCount, make_pair(i->first, *j), i->second);
            } else {
                addKey(pairCount, make_pair(*j, i->first), i->second);
            }
        }
    }
    ext_node.clear();
          
	// step 2: missing edges within core_ext

	for (map< set <int>, int >::iterator i=core_ext.begin(); i!=core_ext.end(); i++){
	  if ((i->first).size() <= 1) continue;
	  ci = (i->first).end();
	  ci--;
	  for (ai=(i->first).begin();ai!=ci;ai++){
        if (ext_node_set.find(*ai)==ext_node_set.end()) continue;  // maybe speed up, not sure
	    bi = ai;
	    bi++;
        for (;bi!=(i->first).end();bi++){ // keys are sorted, *ai < *bi,  i hope..
          if (ext_node_set.find(*bi)==ext_node_set.end()) continue;  
	      addKey(pairCount, make_pair(*ai, *bi), i->second); // look up expensive
        }
	  }
	}
          
    ext_node_set.clear();
	MisPair * p_pair;
	for (map<pair<int,int>, int>::iterator pi=pairCount.begin();pi!=pairCount.end();pi++){
	  if (pi->second < freq_thn) continue; //must be frequent enough
	  id1=(pi->first).first; 
	  id2=(pi->first).second;
	  mi = misPairs.find(pi->first);
	  // if exist in pre-window, update position
	  if (mi != misPairs.end()) { 
	    (mi->second)->p_end = cur_pos; 
	    (mi->second)->flag = 1;	    
	  } else { // fixme: window as position for now
	    // if not exist in pre-window and not reported, create new
	    if (m_neighbor[id1].find(id2)== m_neighbor[id1].end()) {
	      p_pair = new MisPair(cur_pos_start, cur_pos, 1);
          misPairs[pi->first] = p_pair;
	    }
	  }
	}
	pairCount.clear();
          
          /*
           for (map<pair<int,int>, int>::iterator pi=pairCount.begin();pi!=pairCount.end();pi++)
           cerr << (pi->first).first <<" " << (pi->first).second << " "<< pi->second << endl;
           */
      }    
      // done with update missing edges
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
  
  printMissing(minLen);
  initMisFlag();
}

// used to select the seed node pair
float FastGraphCluster::getWeight(float edgeWeight,float vertexWeight, float seedW/*=NULL*/)
{
  int nbin=10;      // 0.8, 0.9, 1.0 have different values
  //return m_nVertex * edgeWeight + vertexWeight; // similar to dash
  return  m_nVertex * (int)(nbin * edgeWeight) + vertexWeight; 
}

// create core clusters with density==1
int FastGraphCluster::buildCore(int index, set<int> &result, FibonacciHeap &heap, set<int> &changed)
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
        w = (imap->second)->weight;
        ptr = m_pHeapNode[j];
	if (ptr==NULL) continue;
	if (ptr->key.wsum == UNEXPLORED){
	  ptr->key.wsum = w;
	  ptr->key.esum = 1;
	  heap.insert(ptr);
	}else {
	//cerr << "ptr "<< ptr->key.index << " " << ptr->key.wsum << " " << ptr->key.esum  << endl;
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
  //return v0 * 0.9; // debug
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
    //cerr << "key.wsum " << ptr->key.wsum << " " << ptr->key.esum << " " << *i << endl;
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
    if (nVertex > 1) {
      density = 2.0*tEdge/(nVertex*(nVertex-1));
      increase = es/(density*(nVertex-1)); 
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
