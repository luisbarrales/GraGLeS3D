#ifndef		__GRAIN_JUNCTION__
#define		__GRAIN_JUNCTION__

#include "spoint.h"
#include "box.h"

/*!
 * \enum E_JUNCTION_TYPE
 * \brief Enumeration that is used to identify the type of junction described by the
 * structure.
 */
enum E_JUNCTION_TYPE
{
	E_TRIPPLE_JUNCTION,
	E_QUADRUPLE_JUNCTION,
	E_INVALID_JUNCTION
};

#define MAXIMUM_JUNCTION_ORDER	4
/*!
 * \struct GrainJunction
 * \brief Structure that encapsulates a grain junction. It carries information about the
 * order of the junction as well as all the grains active at that junction.
 */
struct GrainJunction
{
	E_JUNCTION_TYPE	junction_type;
	SPoint coordinates;
	LSbox* grains[MAXIMUM_JUNCTION_ORDER];
	GrainJunction(LSbox** boxes, int count, SPoint coords) : coordinates(coords)
	{
		if(count == 3)
			junction_type = E_TRIPPLE_JUNCTION;
		else if(count == 4)
			junction_type = E_QUADRUPLE_JUNCTION;
		else
			junction_type = E_INVALID_JUNCTION;
		for(int i=0; i < count && i < MAXIMUM_JUNCTION_ORDER; i++)
		{
			grains[i] = boxes[i];
		}
	}
};

#endif		//__GRAIN_JUNCTION__
