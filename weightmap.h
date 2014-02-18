#ifndef WEIGHTMAP_H
#define WEIGHTMAP_H
#include "ggLS.h"

class grainhdl;
class LSbox;

/*!
 * \class Weightmap
 * \brief Class that encapsulates weightmap for gridpoints.
 * Each entry represents a weight value corresponding to the 3 closest LSBoxes that influence the point.
 */
class Weightmap{
/*!
 * \struct mapkey
 * \brief Internal structure that represents the key with which the actual value is fetched.
 * In essence represents an ordered 3element tuple.
 */
		struct mapkey
		{
			//TODO: Investigate pointer usage and possibly migrate to Box IDs.
			LSbox* first; 	/*!< First or the "smallest address" of an LSBox.*/
			LSbox* second; 	/*!< Second or the "second smallest address" of an LSBox.*/
			
			mapkey(vector<LSbox*> IDs);

			bool operator<(const mapkey & other) const; /*!< Overloaded operator that compares
														 *ordered 3 tuples required for the map
														 *STL container.*/
		};
private:
	grainhdl* 					m_pHandler;	/*!< Pointer to the owning handler.*/
	std::map<mapkey, double> 	m_Weights;	/*!< Actual weight information.*/
	/*!
	 * \brief Method to compute weights.
	 * \param rep The ordered tuple containing the three active LSBoxes.
	 * \param me Unknown usage.
	 * \param ST Unknown usage.
	 */
	double 	computeWeights( mapkey rep, LSbox* me, double* ST);
public:
	friend class LSbox;
	/*!
	 * \brief Basic constructor.
	 * \param owner Specifies the owning \b grainhdl instance.
	 */
	Weightmap(grainhdl* owner);
	/*!
		 * \brief Basic destructor.
	*/
	
	double isTriplePoint(vector<LSbox*> IDs);
	/*!
	 * \brief Method to return a weight if Triple Point exist else returns 0.0.
	 * \param IDs IDs of the neighboring cells.
	 * 
	 * 
	 */
	
	
	~Weightmap();
	/*!
	 * \brief Method to return weights.
	 * \param IDs IDs of the neighboring cells.
	 * \param me Unknown usage.
	 * \param ST Surface Tension Array.
	 */
	double loadWeights(vector<LSbox*> IDs, LSbox* me, double* ST);
	/*!
	 * \brief Debug method to plot the weightmap.
	 * \param length Unknown usage.
	 * \param ID Unknown usage.
	 * \param ST Unknown usage.
	 * \param zeroBox Unknown usage.
	 */
	void plotWeightmap(int length, LSbox*** ID, double* ST, LSbox * zeroBox);
};

#endif
