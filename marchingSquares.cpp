#include "marchingSquares.h"


MarchingSquaresAlgorithm::MarchingSquaresAlgorithm(DimensionalBufferVar& distance_buffer, DimensionalBufferIDLocal& id_local, LSbox* current_grain) :

	m_DistanceBuffer(distance_buffer),
	m_IDLocal(id_local),
	m_CurrentGrain(current_grain),
	m_top(-1), m_bottom(-1), m_left(-1), m_right(-1)
{}
bool MarchingSquaresAlgorithm::generateContour(std::vector<SPoint>& contour_output, std::vector<GrainJunction>& junction_output)
{
	int startX = -1, startY = -1;
	//Find starting point on contour
	for(int i=m_DistanceBuffer.getMinY(), done = 0; i<m_DistanceBuffer.getMaxY() && !done; i++)
		for(int j=m_DistanceBuffer.getMinX(); j<m_DistanceBuffer.getMaxX(); j++)
		{
			if(isInside(i,j))
			{
				startX = j;
				startY = i;
				done = 1;
				break;
			}
		}
	if(startX == startY && startX == -1)
	{
		return false;
	}
	
	//Start point found determine four corners of square
	m_left =  	startX - 1;
	m_right = 	startX;
	m_top = 	startY - 1;
	m_bottom = 	startY;
	E_MOVEMENT_DIRECTIONS next_step = NO_DIRECTION;
	E_MOVEMENT_DIRECTIONS previous_step = NO_DIRECTION;
	do
	{
		int stateMask = 0;

		if(isInside(m_top, m_left))
			stateMask |= TOP_LEFT;
		if(isInside(m_top, m_right))
			stateMask |= TOP_RIGHT;
		if(isInside(m_bottom, m_left))
			stateMask |= BOTTOM_LEFT;
		if(isInside(m_bottom, m_right))
			stateMask |= BOTTOM_RIGHT;
		SPoint next(-1,-1,-1);
		switch( stateMask )
		{
			case TOP_LEFT:
				next_step = UP;
				//generate point on m_top side of square
				next = generatePoint(next_step);
				insertPoint(contour_output, next);
				break;
			case TOP_RIGHT:
				next_step = RIGHT;
				//generate point on m_right side of square
				next = generatePoint(next_step);
				insertPoint(contour_output, next);
				break;
			case TOP_LEFT | TOP_RIGHT:
				next_step = RIGHT;
				//generate point on m_right side of square
				next = generatePoint(next_step);
				insertPoint(contour_output, next);
				break;
			case BOTTOM_LEFT:
				next_step = LEFT;
				//generate point on the m_left side of square
				next = generatePoint(next_step);
				insertPoint(contour_output, next);
				break;
			case BOTTOM_LEFT | TOP_LEFT:
				next_step = UP;
				//generate point on m_top side of square
				next = generatePoint(next_step);
				insertPoint(contour_output, next);
				break;
			case BOTTOM_LEFT | TOP_RIGHT:
				if(previous_step == UP)
					next_step = RIGHT;
				else
					next_step = LEFT;
				break;
			case BOTTOM_LEFT | TOP_RIGHT | TOP_LEFT:
				next_step = RIGHT;
				//generate point on m_right side of square
				next = generatePoint(next_step);
				insertPoint(contour_output, next);
				break;
			case BOTTOM_RIGHT:
				next_step = DOWN;
				//generate point on m_bottom side of square
				next = generatePoint(next_step);
				insertPoint(contour_output, next);
				break;
			case BOTTOM_RIGHT | TOP_LEFT:
				if(previous_step == RIGHT)
					next_step = DOWN;
				else
					next_step = UP;
				break;
			case BOTTOM_RIGHT | TOP_RIGHT:
				next_step = DOWN;
				//generate point on m_bottom side of square
				next = generatePoint(next_step);
				insertPoint(contour_output, next);
				break;
			case BOTTOM_RIGHT | TOP_RIGHT | TOP_LEFT:
				next_step = DOWN;
				//generate point on m_bottom side of square
				next = generatePoint(next_step);
				insertPoint(contour_output, next);
				break;
			case BOTTOM_RIGHT | BOTTOM_LEFT:
				next_step = LEFT;
				//generate point on the m_left side of square
				next = generatePoint(next_step);
				insertPoint(contour_output, next);
				break;
			case BOTTOM_RIGHT | BOTTOM_LEFT | TOP_LEFT:
				next_step = UP;
				//generate point on m_top side of square
				next = generatePoint(next_step);
				insertPoint(contour_output, next);
				break;
			case BOTTOM_RIGHT | BOTTOM_LEFT | TOP_RIGHT:
				next_step = LEFT;
				//generate point on the m_left side of square
				next = generatePoint(next_step);
				insertPoint(contour_output, next);
				break;
		}
		//Generate junctions if existing
		generateJunction(junction_output, stateMask);
		switch(next_step)
		{
			case UP:
				m_top--;
				m_bottom--;
				break;
			case DOWN:
				m_top++;
				m_bottom++;
				break;
			case LEFT:
				m_left--;
				m_right--;
				break;
			case RIGHT:
				m_left++;
				m_right++;
				break;
		}
		previous_step = next_step;
	}
	while(m_left != startX || m_bottom != startY);

	if(!contour_output.empty())
	{
		if(contour_output[0].x != contour_output[contour_output.size()-1].x || contour_output[0].y != contour_output[contour_output.size()-1].y)
		{
			contour_output.push_back(contour_output[0]);
		}
	}
	return true;
}
bool MarchingSquaresAlgorithm::isInside(int row, int column)
{
	if(row < m_DistanceBuffer.getMinY() || row >= m_DistanceBuffer.getMaxY() ||
	   column < m_DistanceBuffer.getMinX() || column >= m_DistanceBuffer.getMaxX())
	   return false;
	else
		return m_DistanceBuffer.getValueAt(row, column) >= 0;
}
void MarchingSquaresAlgorithm::insertPoint(std::vector<SPoint>& output, SPoint p)
{
	if ( output.empty() || p.x != output[output.size()-1].x || p.y != output[output.size()-1].y )
	{
		output.push_back(p);
	}
}
SPoint MarchingSquaresAlgorithm::generatePoint(E_MOVEMENT_DIRECTIONS dir)
{
	SPoint result(-1,-1,-1);
	double delta;
	double position;
	switch(dir)
	{
		case LEFT:
			result.x = m_left;
			if(m_top < m_DistanceBuffer.getMinY())
			{
				result.y = m_bottom;
			}
			else
			{
				delta = m_DistanceBuffer.getValueAt(m_bottom, m_left) - m_DistanceBuffer.getValueAt(m_top, m_left);
				position = m_bottom - m_DistanceBuffer.getValueAt(m_bottom,m_left)/delta;
				result.y = position;
			}
			break;
		case RIGHT:
			result.x = m_right;
			if(m_bottom >= m_DistanceBuffer.getMaxY())
			{
				result.y = m_top;
			}
			else
			{
				delta = m_DistanceBuffer.getValueAt(m_top, m_right) - m_DistanceBuffer.getValueAt(m_bottom, m_right);
				position = m_top + m_DistanceBuffer.getValueAt(m_top,m_right)/delta;
				result.y = position;
			}
			break;
		case UP:
			result.y = m_top;
			if(m_right >= m_DistanceBuffer.getMaxX())
			{
				result.x = m_left;
			}
			else
			{
				delta = m_DistanceBuffer.getValueAt(m_top, m_left) - m_DistanceBuffer.getValueAt(m_top, m_right);
				position = m_left + m_DistanceBuffer.getValueAt(m_top,m_left)/delta;
				result.x = position;
			}
			break;
		case DOWN:
			result.y = m_bottom;
			if(m_left < m_DistanceBuffer.getMinX())
			{
				result.x = m_right;
			}
			else
			{
				delta = m_DistanceBuffer.getValueAt(m_bottom, m_right) - m_DistanceBuffer.getValueAt(m_bottom, m_left);
				position = m_right - m_DistanceBuffer.getValueAt(m_bottom, m_right)/delta;
				result.x = position;
			}
			break;
	}
	return result;
}
void	MarchingSquaresAlgorithm::generateJunction(std::vector<GrainJunction>& junctions, int state_mask) const
{
	LSbox* boxes[4] = {NULL, NULL, NULL, NULL};
	int junction_order = 0;
	if(state_mask & BOTTOM_LEFT)
		insertDistinctPointer(m_CurrentGrain, boxes, junction_order);
	else
		insertDistinctPointer(m_IDLocal.getValueAt(m_bottom, m_left).getElementAt(0),
						  boxes, junction_order);
	if(state_mask & TOP_LEFT)
		insertDistinctPointer(m_CurrentGrain, boxes, junction_order);
	else
		insertDistinctPointer(m_IDLocal.getValueAt(m_top, m_left).getElementAt(0),
						  boxes, junction_order);
	if(state_mask & BOTTOM_RIGHT)
		insertDistinctPointer(m_CurrentGrain, boxes, junction_order);
	else
		insertDistinctPointer(m_IDLocal.getValueAt(m_bottom, m_right).getElementAt(0),
						  boxes, junction_order);
	if(state_mask & TOP_RIGHT)
		insertDistinctPointer(m_CurrentGrain, boxes, junction_order);
	else
		insertDistinctPointer(m_IDLocal.getValueAt(m_top, m_right).getElementAt(0),
						  boxes, junction_order);
	if(junction_order > 2)
	{
		junctions.push_back(GrainJunction(boxes, junction_order, SPoint(m_left+0.5, m_top+0.5, 0.0)));
	}
}

void	MarchingSquaresAlgorithm::insertDistinctPointer(LSbox* pointer, LSbox** array, int& elem_count) const
{
	bool found = false;
	for(int i=0; i<elem_count; i++){
		if(array[i] == pointer)
		{
			found = true;
			break;
		}
	}
	if(false == found)
	{
		array[elem_count] = pointer;
		elem_count++;
	}
}
