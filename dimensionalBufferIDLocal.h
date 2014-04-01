#ifndef __DIMENSIONAL_BUFFER_IDLOCAL__
#define __DIMENSIONAL_BUFFER_IDLOCAL__

#include <vector>
#include <map>
#include <iostream>

#include "dimensionalBuffer.h"
using namespace std;

#define ID_MAX_CHUNK_COUNT 5

class LSbox;

enum E_POSITIONS
{
	E_FIRST_POSITION = 1,
	E_SECOND_POSITION,
	E_THIRD_POSITION,
	E_FOURTH_POSITION,
	E_LAST_POSITION = -1
};

struct IDChunk
{
	LSbox* local_chunks[ID_MAX_CHUNK_COUNT];
	unsigned int total_elements;
	IDChunk() : total_elements(0)
	{}
	IDChunk& operator=(const IDChunk& other)
	{
	    if (this != &other)
	    {
	    	total_elements = other.total_elements;
	    }
	    return *this;
	}
	bool insertAtPosition(E_POSITIONS position, LSbox* element)
	{
		if(position == E_LAST_POSITION)
		{
			if(total_elements >= ID_MAX_CHUNK_COUNT)
			{
				return false;
			}
			else
			{
				local_chunks[total_elements] = element;
				total_elements++;
			}
		}
		else
		{
			int actual_position = position - E_FIRST_POSITION;
			if(actual_position < 0 || actual_position >= ID_MAX_CHUNK_COUNT)
			{
				return false;
			}
			else
			{
				LSbox* element_to_write = element;
				for(int i=actual_position; i<total_elements+1 && i < ID_MAX_CHUNK_COUNT; i++)
				{
					LSbox* swap = local_chunks[i];
					local_chunks[i] = element_to_write;
					element_to_write = swap;
				}
			}
			total_elements++;
			if(total_elements > ID_MAX_CHUNK_COUNT)
				total_elements = ID_MAX_CHUNK_COUNT;
		}
		return true;
	}
	LSbox* getElementAt(int pos)
	{
		if(pos < total_elements)
			return local_chunks[pos];
		else
			return NULL;
	}
	inline void clear()
	{
		total_elements = 0;
	}
};

class DimensionalBufferIDLocal	: public DimensionalBuffer<IDChunk>
{
public:
	void clear()
	{
		for(auto& iterator : this->m_values)
			iterator.clear();
	}
};



#endif //__DIMENSIONAL_BUFFER_IDLOCAL__
