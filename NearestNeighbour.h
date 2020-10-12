
#ifndef NEIGHBOUR_NEARESTNEIGHBOUR_H_
#define NEIGHBOUR_NEARESTNEIGHBOUR_H_

#include <string>
#include <vtkDistancePolyDataFilter.h>

#include "Neighbour.h"
#include "SimulatedShape.h"
#include "Frontend.h"

using namespace std;

/**
 * The NearestNeighbour class finds the nearest point on the outer mesh for ever point
 * on the inner mesh. The distance between these 2 points is measured as the "thickness'
 * and the largest value is found to determine the greatest thickness of the geometry.
 */
class NearestNeighbour : public Neighbour {

public:

// CONSTRUCTORS & DESTRUCTOR

	/**
	 * Default constructor
	 */
	NearestNeighbour();

	/**
	 * Copy constructor
	 */
	NearestNeighbour(const NearestNeighbour&);

	/**
	 * Equality operator
	 */
	NearestNeighbour& operator = (const NearestNeighbour&);

	/**
	 * Destructor
	 */
	~NearestNeighbour();


// FUNCTION

	/**
	 * Returns greatest thickness value between the inner and outer meshes
	 */
	double findThickest();

private:

	/**
	 * thicknessList
	 * String containing the inner point, outer point, and current thickness
	 * measurement from these points
	 */
	string thicknessList;

// HELPER FUNCTION

	/**
	 * Returns closes point on outer mesh from given inner mesh point
	 * @param point on inner mesh
	 */
	double* findClosestPoint(double [3]);
};

#endif /* NEARESTNEIGHBOUR_H_ */
