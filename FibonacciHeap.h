/************************************************************************/
/* Implemented according to <<The introduction to algorithm>>
	operator > and < of Key type must be overloaded!
/************************************************************************/
#pragma once
#include "FastGraphCluster.h"
#include <cstdlib>
#include <cmath>
#include <set>
#include <iostream>

// first version for programming convenience
typedef HeapNode KEYTYPE;

class FiboNode{
public:
	FiboNode(KEYTYPE key,FiboNode *parent = NULL,FiboNode *child = NULL,
		FiboNode *left = NULL ,FiboNode *right = NULL,
		size_t degree = 0,bool mark = false);
	~FiboNode(){}

	KEYTYPE key;	// key value
	FiboNode *parent,*child;	//	child points to any one of its child
	FiboNode *left,*right;	//	link this node into the float linked list of Sibling
	size_t degree;	// the number of children in the child list
	bool mark;	// whether lost a child since the last time x was made the child of another node
};

// in order to get an inverted version. simply change the operate >

class FibonacciHeap
{
public:
	FibonacciHeap(void);
	~FibonacciHeap(void);

	void insert(FiboNode *node);	// make the node yourself first,don't forget to clear the node

	KEYTYPE getMin();	// get the minimum key value
	FiboNode* extractMin();	// extract the min heap

	void decrease(FiboNode *node,KEYTYPE newKey);	// decrease the key
//	void remove(FiboNode *x);
	bool empty(){return minheap == NULL;}
	void combine(FibonacciHeap* rheap);	// combine rheap with ourself
	size_t m_nNode;	// the current number of nodes, i need it public for debugging
	// assume the HeapNode has integer ID ==> otherwise this class is general
	set<int> current_node_id; // record node ids in current heap
private:
	FiboNode *minheap;
	bool *mark;
	void consolidate();
	void cut(FiboNode *x,FiboNode *y);
	void cascadingCut(FiboNode *y);
	void fibHeapLink(FiboNode* y,FiboNode* x);
	void dispose(FiboNode *y);	// dispose the subtree rooted at y
};
