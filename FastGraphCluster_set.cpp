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
	    m_pHeapNode[i] = new FiboNode(HeapNode(UNEXPLORED,i));
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

bool FastGraphCluster::checkClst(Cluster * cl)
{
  int ai = cl->nodes.size(); // old version, -1
  if(ai < m_nLowerSize || ai*(ai-1)/2.0*m_dLowDen > (float) cl->n_edge)
    return false;
  return true;
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
      ci=clstID[(*i).first];
      // edges to delete within cluster
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

// print all Clusters in current window.
void FastGraphCluster::printAllClst(ofstream& fout)
{
  for (map <int, Cluster* >::iterator ci=result_clst.begin(); ci!=result_clst.end(); ci++)
    {   
      Cluster * cl = ci->second;
      if (cl->n_edge < 3) continue; // at least 3 edges 
#ifdef DEBUG
      fout << ci->first << "\t" << cur_pos_start << "\t" << cur_pos <<"\t" << cl->n_edge << "\t" << cl->nodes.size() << "\t";
#else
      fout << "clst" << "\t" << cur_pos_start << "\t" << cur_pos <<"\t";
#endif

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

// return newly created coreClusters in current window
// fix me: old clusters from previous windows are not updated yet
void FastGraphCluster::fastClusterCore()
{
  int i=0, seedn=10;      

  while (!seedArray->empty() && seedArray->top >= m_nLowerSize-1) {
    i = seedArray->getMax();
    FibonacciHeap heap;	// local expanding heap
    set<int> result;
    expandCore(i, result, heap);
////    // number of nodes in heap is always changing !!!
////    if (!result.empty()){
////      for (set <int>::iterator j=result.begin(); j!=result.end(); j++)
////	cout << vertexName[*j] << "\t";
////    }
    
//     cout << result.size() << "\t" << heap.m_nNode << endl;
//     cout << heap.getMin().index << endl;
//     for (set <int>::iterator j=heap.current_node_id.begin(); j!=heap.current_node_id.end(); j++)
//       cout << m_pHeapNode[*j]->key.index << "\t" << m_pHeapNode[*j]->key.esum << "\t" << m_pHeapNode[*j]->key.wsum << endl;
//     DebugFunc::printVec(heap.current_node_id);
    
    if (heap.m_nNode){
      srand(time(NULL));
      for (int dn=0; dn<seedn; dn++){
	extendCore(heap.getMin().index, heap.current_node_id, dn);  // no need to pass result
	
      }
    }
    // only look at the first cluster for now
    // fix me, when to delete heap ??
    break; 
  }
  if (!seedArray->empty())
    delete seedArray;
}

void FastGraphCluster::fastCluster(ofstream& fout)
{
  int i=0;
  while (!seedArray->empty() && seedArray->top >= m_nLowerSize-1)
    {
      i = seedArray->getMax();
      expand(i);
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
int FastGraphCluster::expandCore(int index, set<int> &result, FibonacciHeap &heap)
{
  int tEdge=0, imin, i,j,w3, nVertex = 0;
  float density,increase,w,oldWeight;
  set <int> changed;

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
      w = ptr->key.wsum;   
      w3 = ptr->key.esum;
      
      tEdge += w3;
      nVertex ++;

      if (nVertex > 1) {
	if ( 2.0*tEdge != (nVertex*(nVertex-1)))  // fully connected cluster
	  {
	    nVertex--;
	    tEdge -= w3;
	    break;	
	  }
      }
      
      m_pHeapNode[imin] = NULL;
      seedArray->remove(imin,neighborWeightCnt[imin]);
      result.insert(imin);	// write back result
      
      for (imap = m_neighbor[imin].begin(); imap!=m_neighbor[imin].end(); imap++)
	{
        j = imap->first;
        w = (imap->second)->weight;
        ptr = m_pHeapNode[j];
        if (ptr==NULL) continue;
	if (ptr->key.wsum == UNEXPLORED){
	  ptr->key.wsum = w;
	  ptr->key.esum = 1;
	  heap.insert(ptr);
	}else  heap.decrease(ptr,HeapNode( ptr->key.wsum+w,j,ptr->key.esum+1));
	oldWeight = neighborWeightCnt[j];
        
	if ( ( nVertex == 1 ) && (j == maxindex) )
	  {
	    seedArray->decrease(j,oldWeight,oldWeight + maxWeight - w);
	    neighborWeightCnt[j] -= (w-maxWeight);
	  }else
	  {
	    seedArray->decrease(j,oldWeight,oldWeight - w);
	    neighborWeightCnt[j] -= w;
	  }
	changed.insert(j);
	}
    }
  m_neighbor[index][maxindex]->weight -= maxWeight;
  
  /*
    for ( set<int>::iterator iter = changed.begin() ; iter != changed.end() ; iter++ )
      {
	i = *iter;
	if ( m_pHeapNode[i] != NULL ){ 
	  m_pHeapNode[i]->key.wsum = UNEXPLORED;
	  m_pHeapNode[i]->key.esum = 0;
	}
      }
    changed.clear();
  */
  return 1;
}

// template <class K> 
// struct CompareWsum {
//   bool operator() (const K a, const K b) const 
//   {return a.second >= b.second;}
// };

// add some random variance in the weight 
float addVar(float v0)
{
  float r = 1-(float)(rand()%100)*0.002; // [0,99], yes, each time different run  
  return v0*r;
}

// build local heap from initial nodes, with random seed added 
int FastGraphCluster::extendCore(int start_id, set<int> &node_id, int dn)
{
  FibonacciHeap heapExt; // Extended heap, start with start_id
  float es, rand_ws;
  FiboNode *ptr; 
  set<int> result;
  //set< pair <int, float >, struct CompareWsum <pair <int, float> > > node_id_new;
  for (set<int>::iterator i = node_id.begin();i!=node_id.end();i++){
    ptr = m_pHeapNode[*i];
    es = (float) ptr->key.esum; 
    node_id_new.insert(make_pair(*i, es));
  }
  for (set< pair <int, float> >::iterator i = node_id_new.begin();i!=node_id_new.end();i++){
      rand_ws = addVar((*i).second); 
      //cout << (*i).first << "## " <<  (*i).second << "## " << rand_ws << endl;
      ptr = 
    }
	 
	 
//  heapExt.insert(ptr);
//  while (!heapExt.empty())
//    {
//      ptr = heapExt.extractMin();
//    }
  

  //cout << "generate2 " << r << endl;
}

int FastGraphCluster::expand(int index)
{
        int tEdge=0, imin, i,j,w3, nVertex = 0;
	set<int> result;
	float density, w_density,increase,w,totalWeight = 0,oldWeight;

	FibonacciHeap heap;	// local expanding heap
	set <int> changed;

	FiboNode *ptr = m_pHeapNode[index];
	ptr->key.wsum = 0;
	ptr->key.esum = 0;
	heap.insert(ptr);

	// start second level heuristic value selection
	int maxindex = -1; // default 0 may cause problem
	float maxWeight = 0;

	map<int, EdgeInfo*>::iterator imap;    

	float seed_WeightCnt = neighborWeightCnt[index];
	//cerr << "start " << neighborWeightCnt[index] << " " << m_neighbor[index].size();
	for (imap = m_neighbor[index].begin(); imap!=m_neighbor[index].end(); imap++)
	  {
	    j = imap->first;
	    if (m_pHeapNode[j] == NULL) continue;
	    w = getWeight((imap->second)->weight,neighborWeightCnt[j], seed_WeightCnt);

	    if (w > maxWeight )
	      {
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
	    w = ptr->key.wsum;   // initial value = 0
	    w3 = ptr->key.esum;

	    totalWeight += w;  // initial value = 0
	    tEdge += w3;
	    nVertex ++;

	    // more than first 2 seeds, check density etc. 
	    if (nVertex > 1)
	      {
		//w_density = 2.0*totalWeight/(nVertex*(nVertex-1));
		density = 2.0*tEdge/(nVertex*(nVertex-1));
		increase = w3/(density*(nVertex-1));  // increase of edge
		
		if ( density < m_dLowDen || increase < m_dLowerIncrease )
		  {
		    // about 1/3 of them stopped because increase < 0.5
		    //cout << "break " << density << " "<< increase << " "<< nVertex << endl;
		    nVertex--;
		    totalWeight -= w;
		    tEdge -= w3;
		    break;	// exceeding the lowest density
		  }
	      }

	    m_pHeapNode[imin] = NULL;
	    // remove merged vertex from seed list
	    seedArray->remove(imin,neighborWeightCnt[imin]);
	    result.insert(imin);	// write back result

    for (imap = m_neighbor[imin].begin(); imap!=m_neighbor[imin].end(); imap++)
    {
        j = imap->first;
        w = (imap->second)->weight;
        ptr = m_pHeapNode[j];
        if (ptr==NULL) continue;
	if (ptr->key.wsum == UNEXPLORED)
	  {	// a new expanded vertex
                ptr->key.wsum = w;
                ptr->key.esum = 1;
                heap.insert(ptr);
	  }else  heap.decrease(ptr,HeapNode( ptr->key.wsum+w,j,ptr->key.esum+1));
	// decrease its weight degree which is independent of the cluster
	oldWeight = neighborWeightCnt[j];
        
	if ( ( nVertex == 1 ) && (j == maxindex) )
	  {
	    seedArray->decrease(j,oldWeight,oldWeight + maxWeight - w);
	    neighborWeightCnt[j] -= (w-maxWeight);
	  }else
	  {
	    seedArray->decrease(j,oldWeight,oldWeight - w);
	    neighborWeightCnt[j] -= w;
	  }
		    changed.insert(j);
    }
    
          
    if (nVertex == 1 ) totalWeight -= maxWeight; // deduce back extra weight for seeding pair
      }
    for ( set<int>::iterator iter = changed.begin() ; iter != changed.end() ; iter++ )
      {
	i = *iter;
	if ( m_pHeapNode[i] != NULL ){
	  m_pHeapNode[i]->key.wsum = UNEXPLORED;
	  m_pHeapNode[i]->key.esum = 0;
	}
      }
    changed.clear();

    m_neighbor[index][maxindex]->weight -= maxWeight;
    if (result.size() >= m_nLowerSize)
      {
    	for (set<int>::iterator ni=result.begin();ni!=result.end();ni++)
	  clstID[*ni] = clst_topindex;
	// cout << "new cluster " << clst_topindex <<" " << tEdge << endl;
	result_clst[clst_topindex++] = new Cluster(result, tEdge);
      }
    return 1;
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
