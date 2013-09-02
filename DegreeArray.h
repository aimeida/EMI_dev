#pragma once
#include <iostream>
#include <list>
#include <cstdlib>
using std::list;

class DegreeArray
{
public:
	DegreeArray(const float *degreeWeight,size_t V,float maxWeight);
	~DegreeArray(void);
	size_t extractMax();	// return the vertex index
	size_t getMax();
	bool empty(){return m_nVertex ==0;}
	void decrease(size_t inx,float oValue,float nValue);
	void remove(size_t inx,float value);
	size_t m_nVertex,top;

private:
	list<size_t>::iterator *indexmap;
	list<size_t> *m_pDegree;	// round off the vertex weight to integer
};
