

#include <string>

#include <vtkCleanPolyData.h>
#include <vtkDistancePolyDataFilter.h>
#include <vtkPointLocator.h>

#include "NearestNeighbour.h"
#include "Frontend.h"

using namespace std;

NearestNeighbour::NearestNeighbour()
: Neighbour() {
}

NearestNeighbour::NearestNeighbour(const NearestNeighbour& src)
: Neighbour(src){
}

NearestNeighbour& NearestNeighbour::operator = (const NearestNeighbour& src) {
	Neighbour::operator = (src);
	return *this;
}

NearestNeighbour::~NearestNeighbour() {

}


// member functions

// point must be given from inner mesh
double* NearestNeighbour::findClosestPoint(double point [3]) {
	double *closestPoint = (double*) malloc(3* sizeof(double));
	vtkIdType i;

	vtkSmartPointer <vtkPointLocator> locator = vtkSmartPointer <vtkPointLocator> :: New();


	locator -> SetDataSet(getAorta() -> getOuter());
	locator -> BuildLocator();

	i = locator -> FindClosestPoint(point);
	getAorta() -> getOuter() -> GetPoint(i, closestPoint);

	return closestPoint;
}

double NearestNeighbour::findThickest() {
	SimulatedShape *temp = new SimulatedShape(*getAorta());
	int numOfPoints = temp -> getInner() -> GetNumberOfPoints();
	double innerPoint[3];
	double outerPoint[3];
	string list;

	double currThickness;
	for (vtkIdType i = 0 ; i < numOfPoints; i++) {
		temp -> getInner() -> GetPoint(i, innerPoint);
		copyPoint(outerPoint, findClosestPoint(innerPoint));
		currThickness = calculateDistance(innerPoint, outerPoint);
		if (currThickness > getThickness()) {
			setThickness(currThickness);
			setInnerThickestPoint(innerPoint);
			setOuterThickestPoint(outerPoint);
		}
		list += to_string(innerPoint[0]) + "," + to_string(innerPoint[1]) + "," + to_string(innerPoint[2])
			+ ",," + to_string(outerPoint[0]) + "," + to_string(outerPoint[1]) + "," + to_string(outerPoint[2])
			+ ",," + to_string(currThickness) + "\n";
	}
	setThicknessList(list);
	return getThickness();
}
