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
	unsigned int total_chunks;
	IDChunk() : total_chunks(0)
	{}
	IDChunk& operator=(const IDChunk& other)
	{
	    if (this != &other)
	    {
	    	total_chunks = other.total_chunks;
	    }
	    return *this;
	}
	bool insertAtPosition(E_POSITIONS position, LSbox* element)
	{
		if(position == E_LAST_POSITION)
		{
			if(total_chunks >= ID_MAX_CHUNK_COUNT)
			{
				return false;
			}
			else
			{
				local_chunks[total_chunks] = element;
				total_chunks++;
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
				for(int i=actual_position; i<total_chunks+1 && i < ID_MAX_CHUNK_COUNT; i++)
				{
					LSbox* swap = local_chunks[i];
					local_chunks[i] = element_to_write;
					element_to_write = swap;
				}
			}
			total_chunks++;
			if(total_chunks > ID_MAX_CHUNK_COUNT)
				total_chunks = ID_MAX_CHUNK_COUNT;
		}
		return true;
	}
	LSbox* getElementAt(int pos)
	{
		if(pos < total_chunks)
			return local_chunks[pos];
		else
			return NULL;
	}
};

class DimensionalBufferIDLocal	: public DimensionalBuffer<IDChunk>
{
public:
	void clear()
	{
		this->clearValues(IDChunk());
	}
};



#endif //__DIMENSIONAL_BUFFER_IDLOCAL__
