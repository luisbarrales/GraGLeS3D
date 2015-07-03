/*
	GraGLeS 2D A grain growth simulation utilizing level set approaches
    Copyright (C) 2015  Christian Miessen, Nikola Velinov

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef		__GRAIN_BOUNDARY__
#define		__GRAIN_BOUNDARY__

#include "dimensionalBufferReal.h"
#include "dimensionalBufferIDLocal.h"
#include "triangle.h"

#include <vector>
#include <map>
using namespace std;

class LSbox;

/*!
 * \class ExplicitGrainBoundary
 * \brief Class that holds all information required to describe a grain boundary.
 *
 * This class is a utility class that is used to encapsulate the explicit artifacts of a
 * grain boundary. It is responsible for generating the explicit boundary from a given
 * distance function, identify the neighboring grains to the given grain, construct all sectors
 * on the boundary as well as compute information for the grain like perimeter and volume.
 */

class ExplicitGrainBoundary
{
public:
	ExplicitGrainBoundary(LSbox* owner);
	~ExplicitGrainBoundary();
	inline void addDirectNeighbor(unsigned int id) {m_directNeighbors.push_back(id);}
	bool generateHull();
private:

	LSbox*					m_owningGrain;
	vector<unsigned int>	m_directNeighbors;
	vector<Triangle>		m_grainHull;

};

#endif		//__GRAIN_BOUNDARY__
