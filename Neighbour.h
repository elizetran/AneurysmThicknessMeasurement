/**
 * The Neighbour class is the abstract class which is to be inherited by NormalNeighbour
 * and NearestNweighbour measuring algorithms. All functions and data required by both
 * classes outlined in Neighbour.
 */

#ifndef NEIGHBOUR_H_
#define NEIGHBOUR_H_

#include <string>

#include "SimulatedShape.h"

using namespace std;

/**
 * The neighbour class is the parent class for both the NormalNeighbour and NearestNeighbour
 * classes. It contains information required by both algorithms, such as the mesh and global
 * thickness value.
 */
class Neighbour {

public:

// CONSTRUCTORS & DESTRUCTORS

	/**
	 * Default constructor
	 */
	Neighbour();

	/**
	 * Copy constructor
	 */
	Neighbour(const Neighbour&);

	/**
	 * Equality operator
	 */
	Neighbour& operator = (const Neighbour& src);

	/**
	 * Destructor
	 * To be overridden by children class
	 */
	virtual ~Neighbour();


// SETTERS & GETTERS

	/**
	 * Setter for aorta mesh
	 * @param mesh to be copied
	 */
	void setAorta(SimulatedShape*);

	/**
	 * Setter for inner thickest point
	 * @param point on inner mesh where thickest measurement is made
	 */
	void setInnerThickestPoint(double [3]);

	/**
	 * Setter for outer thickest point
	 * @param point on outer mesh where thickest measurement is made
	 */
	void setOuterThickestPoint(double [3]);

	/**
	 * Setter for thickness
	 * @param value of thickest measurement calculated in whole mesh
	 */
	void setThickness(double);

	/**
	 * Setter for thicknessList
	 * @param string containing all data for thickness measurement and location
	 */
	void setThicknessList(string);

	/**
	 * Returns aorta mesh
	 */
	SimulatedShape *getAorta();

	/**
	 * Returns point on inner mesh where thickest measurement occurs
	 */
	double* getInnerThickestPoint();

	/**
	 * Returns point on outer mesh where thickest measurement occurs
	 */
	double* getOuterThickestPoint();

	/**
	 * Returns maximum thickness calculated throughout the aorta mesh
	 */
	double getThickness();

	/**
	 * Returns tring containing all data for thickness measurement
	 * and location
	 */
	string getThicknessList();


// FUNCTIONS: general measurement functions

	/**
	 * Returns the distance between 2 points
	 * @param point1: starting point for distance measurement
	 * @param point2: ending point for distance measurement
	 */
	double calculateDistance(double [3], double [3]);

	/**
	 * Copies data of one point into another and returns altered point
	 * @param destination point to be changed
	 * @param source point holding data to be copied
	 */
	double* copyPoint(double[3], double[3]);


private:

// MEMBER VARIABLES

	/**
	 * aorta
	 * mesh for inner and outer aorta wall
	 */
	SimulatedShape *aorta;

	/**
	 * thickness
	 * Maximum thickness calculated throughout whole aorta mesh
	 */
	double thickness;

	/**
	 * innerThickestPoint
	 * point on inner mesh where maximum thickness measurement was calculated
	 */
	double *innerThickestPoint;

	/**
	 * outerThickestPoint
	 * point on outer mesh where maximum thickness measurement was calculated
	 */
	double *outerThickestPoint;

	/**
	 * thicknessList
	 * string containing data about thickness measurements and the location of measurement
	 * for all points in mesh
	 */
	string thicknessList;

};

#endif /* NEIGHBOUR_H_ */
