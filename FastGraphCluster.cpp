#include "DegreeArray.h"
#include "FibonacciHeap.h"
#include "FastGraphCluster.h"
#include "DebugFunc.h"
using namespace std;

#define DEBUGMODE 
#define UNEXPLORED -1
#define P_CLST_STOP 0.2
#define FREQ_ERROR 0.3 // add more before counting
#define getmax(a,b) ((a)>(b)?(a):(b))
#define getmin(a,b) ((a)<(b)?(a):(b))

extern vector<string> vertexName;
extern int incre_cutoff;

FastGraphCluster::FastGraphCluster(float density,int lowersize,float lowerincrease, int m_nVertex, int cew):m_nLowerSize(lowersize),m_dLowDen(density),m_dLowerIncrease(lowerincrease), m_nVertex(m_nVertex), continuous_empty_wins(cew), clst_topindex(0)
{
        m_neighbor = new map<int, EdgeInfo* >[m_nVertex]; 
	neighborWeightCnt = new float[m_nVertex];
	m_pHeapNode = new FiboNode*[m_nVertex];
	m_pHeapNode2 = new FiboNode*[m_nVertex];
	for (int i=0;i<m_nVertex;i++) {
	  m_pHeapNode[i] = new FiboNode(HeapNode(UNEXPLORED,i,0));
	  m_pHeapNode2[i] = m_pHeapNode[i];
	  clstID[i] = -1;
	}
}

FastGraphCluster::~FastGraphCluster(void)
{
  cerr << "delete called :-) " << endl;
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
    for (map <pair<int, int>, MisPair* >::iterator mi=misPairs.begin(); mi!=misPairs.end();mi++) 
      delete mi->second;
    misPairs.clear();
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
  for (i=misPairs.begin(); i!=misPairs.end(); i++){
      if ((i->second)->flag < 0 ) {// not updated on this window
          (i->second)->freqs.push_back(0);
          (i->second)->p_end = cur_pos;
      }
      (i->second)->flag -= 1;
  }
}

void FastGraphCluster::refine_bound(vector <float>::iterator &ri, float &pair_start, float &pair_end, float raw_start, float raw_end, vector <float> &freqs)
{
    pair_start = raw_start;
    pair_end = raw_end;
    ri = freqs.begin();
}

void FastGraphCluster::refine_bound(vector <float>::iterator &ri, float &pair_start, float &pair_end, float window_size, float raw_start, float raw_end, vector <float> &freqs, float freq_th, float window_size_nfold, bool verbose)
{
    
    // eg. can be a problem for example like this: 1 1 0 0 0 1
    pair_start=-1;
    pair_end=-1;
    ri = freqs.begin();

    int ni=0, nlen = freqs.size();
    float vj, v_min = 0;
    float sum_freq = 0;
    for (vector <float>::iterator i=freqs.begin(); i!=freqs.end(); i++) sum_freq += *i ;
    
    v_min = getmin(*ri,freqs.back());
    while ( ( v_min < freq_th || sum_freq/nlen < freq_th ) && (nlen >= window_size_nfold) ) {
        if ( *ri == v_min) {
            ni++;
            ri++;
        } else {
            freqs.pop_back();
        }
        sum_freq -= v_min;
        nlen--;
        v_min = getmin(*ri,freqs.back());
        //if (verbose) cerr << "hi " <<  sum_freq << " " << v_min << " " << nlen << " " << freq_th << endl;
    }
    

    if (nlen >= window_size_nfold){
        pair_start = raw_start + window_size * ni;
        pair_end = raw_start + window_size * (ni+nlen);
    }
}

void FastGraphCluster::printMissing(float window_size, float window_size_nfold, float freq_th, ofstream& fout1)
{
    float pair_start, pair_end;
    vector<float>::iterator ri;
  // don't erase while iterating, may cause problem
  vector <map <pair<int, int>, MisPair* >::iterator> rms;
  for (map <pair<int, int>, MisPair* >::iterator i=misPairs.begin(); i!=misPairs.end(); i++){
     if ( abs((i->second)->flag) > continuous_empty_wins){
       rms.push_back(i);
         
       ri = (i->second)->freqs.begin();
       refine_bound(ri, pair_start, pair_end, window_size, (i->second)->p_start, (i->second)->p_end, (i->second)->freqs, freq_th, window_size_nfold);
      
//         if (i->first.first==2522 && i->first.second==2546){
//           cerr <<"before\t" << cur_pos << "\t" << (i->second)->p_start << "\t" << (i->second)->p_end << "\t";
//           DebugFunc::printVec(i->second->freqs);
//         }
//         if (i->first.first==2522 && i->first.second==2546){
//             refine_bound(ri, pair_start, pair_end, window_size, (i->second)->p_start, (i->second)->p_end, (i->second)->freqs, freq_th, window_size_nfold, 1);}
//         else{
//             refine_bound(ri, pair_start, pair_end, window_size, (i->second)->p_start, (i->second)->p_end, (i->second)->freqs, freq_th, window_size_nfold);}
    
	   if (pair_end > 0){
	  ///cout << (i->second)->flag << "\t" << vertexName[(i->first).first] << "\t" << vertexName[(i->first).second] << "\t";
	     fout1 << (i->first).first << "\t" << (i->first).second << "\t" << pair_start << "\t" << pair_end << endl;
           
//	  for (;ri!=(i->second)->freqs.end(); ri++)
//	    fout1 << *ri << "\t" ;
//	  fout1 << endl;
	  }
    }
  }
    
  vector <map <pair<int, int>, MisPair* >::iterator>::iterator mi;
  for (mi=rms.begin(); mi!=rms.end();mi++){
    vector<float>().swap(((*mi)->second)->freqs);
    delete (*mi)->second;
    misPairs.erase(*mi);
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

void FastGraphCluster::dissolve()
{
  for (map<int, int>::iterator i = clstID.begin(); i != clstID.end(); i++)
    i->second = -1;
  for (map <int, Cluster* >::iterator i=result_clst.begin();i!=result_clst.end();i++)
    delete i->second;
  result_clst.clear();
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
      if ( ci > -1 &&  ci==clstID[(*i).second]){
	addKey(changed_clst, ci);
      }
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

void FastGraphCluster::updateInput(list<Pairmatch * > &active_matches)
{
  float *nw = neighborWeightCnt;
  FiboNode *mp2;
  for (int i=0;i < m_nVertex;i++) {
    mp2=m_pHeapNode2[i];         // use mp2 to reserve the space, alloc/dealloc expensive
    mp2->key.wsum = UNEXPLORED;
    mp2->key.esum = 0;
    if (clstID[i] > -1){  // already clustered
      m_pHeapNode[i] = NULL;
    } else m_pHeapNode[i] = mp2;
    *nw++ = 0;
  }
  list< Pairmatch * >::iterator am;  
  for (am = active_matches.begin(); am != active_matches.end(); am++) {
    if (clstID[(*am)->i1] > -1 || clstID[(*am)->i2] > -1) continue; 
    neighborWeightCnt[(*am)->i1] += (*am)->weight;
    neighborWeightCnt[(*am)->i2] += (*am)->weight;
  }
  float maxWeightDegree = 0;
  for (int i=0;i<m_nVertex;i++){
    if( neighborWeightCnt[i] > maxWeightDegree ) maxWeightDegree = neighborWeightCnt[i];
  }
  seedArray = new DegreeArray(neighborWeightCnt,m_nVertex,maxWeightDegree);
}

void FastGraphCluster::updateInput(list<Pairmatch * > &active_matches, list<Pairmatch * > &emi_matches){
}

void FastGraphCluster::fastClusterCore(int seedn, float n_overhead, float freq_th, float window_size, float window_size_nfold, ofstream& fout1)
{
  set<int> result, changed;
  int freq_thn = (int) seedn * (freq_th-FREQ_ERROR);  // allow more segments to be included
  size_t n_clst_all = 0, n_clst_ext = 0;
  //cerr << "clst size " <<  result_clst.size() << endl;  // should be 0
  while (!seedArray->empty() && seedArray->top >= m_nLowerSize-1) {
    FibonacciHeap heap;	// local expanding heap
    int i = seedArray->getMax();   
    buildCore(i, result, heap, changed); // update clst_topindex, added to result_clst
    // make sure  with\out extension, core size stays the same.
    // cerr << cur_pos << "\t" << result.size() << endl;
    if (heap.m_nNode > 0){ // there are nodes left after the core
      map< set <int>, int > core_ext;  // count number of times it appears
      set<int> surround;
#ifdef DEBUGMODE
    srand(0);
#else
    srand(time(NULL));
#endif
      for (int dn=0; dn<seedn; dn++){
        
          int core_ext_size = core_ext.size();
          if (dn >= n_overhead && core_ext_size <= P_CLST_STOP * dn) { // 0.2 is arbitary, need more tests
	    if (core_ext_size) {
	      float n_fold = seedn/(float)dn;
	      for (map< set <int>, int >::iterator i=core_ext.begin(); i!=core_ext.end(); i++) i->second *= n_fold;
	    } 
	    break;
          } // at least 50% calculations can be saved
    
          extendCore(surround, result, heap.current_node_id, dn);
          if (!surround.empty())
              addKey(core_ext, surround);
          surround.clear();
      }
      
    // cerr << core_ext.size() << endl;  // make hist plot 
    
    if(!core_ext.empty()){
    // n_clst_ext += 1;
	map <pair<int, int>, MisPair* >::iterator mi;
	map<pair<int,int>, int> pairCount;
	set<int>::iterator ai,ci,bi;
	int id1, id2;
        
	//step 1: missing edges between core_ext and core (look at this first, use freq_thn as filter)
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
    float new_freq;
    for (map<pair<int,int>, int>::iterator pi=pairCount.begin();pi!=pairCount.end();pi++){
      //if (pi->second > seedn) cerr << "hello " << pi->second << endl;
      if ( pi->second < freq_thn ) continue; //must be frequent enough
      id1=(pi->first).first; 
      id2=(pi->first).second;
      mi = misPairs.find(pi->first);
   
      // if exist in pre-window, update position
      if (mi != misPairs.end()) {
         // if (id1 == 347 && id2==1517) cerr << "%%%%%%%%%%%%%%%%%%%%%%%%" << cur_pos << endl;
          (mi->second)->flag = 0;
          new_freq = pi->second/(float)seedn;
          if ((mi->second)->p_end < cur_pos) {
              (mi->second)->p_end = cur_pos;
              (mi->second)->freqs.push_back(new_freq);
          } else if ( new_freq > (mi->second)->freqs.back() ) {      // pair of segments have been found by other clusters
                (mi->second)->freqs.pop_back();
                (mi->second)->freqs.push_back(new_freq);
          }
    //  if ((cur_pos - (mi->second)->p_start)/window_size > (mi->second)->freqs.size())
    //      cerr << "err " << id1 << " " << id2 << " " << (cur_pos - (mi->second)->p_start)/window_size << " " <<(mi->second)->freqs.size() << endl;
          
      } else if (m_neighbor[id1].find(id2)== m_neighbor[id1].end()) {  // if not exist in pre-window and not reported, create new
              p_pair = new MisPair(cur_pos-window_size, cur_pos, 0, (pi->second)/(float)seedn);
	          misPairs[pi->first] = p_pair;
	  }
        
    }
    pairCount.clear();
   }
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
  printMissing(window_size, window_size_nfold, freq_th, fout1);
  initMisFlag();
  //cerr << "####### " << cur_pos << endl;
  //float n_clst_ext_p = (float)(n_clst_ext+1)/(result_clst.size()+1);
  //cerr << cur_pos << "## " << n_clst_ext << " " << result_clst.size() << " " << n_clst_ext_p << endl;
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
      if (incre_cutoff){
	if ( density < m_dLowDen || increase < m_dLowerIncrease) {
          nVertex--;
          tEdge -= es;
          break;
	}
      } else {
	if ( density < m_dLowDen ) { 
	  nVertex--;
	  tEdge -= es;
	  break;
	}
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
