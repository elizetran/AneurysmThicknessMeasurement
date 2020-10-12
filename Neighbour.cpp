
#include "Neighbour.h"

Neighbour::Neighbour() {
	aorta = new SimulatedShape();
	thickness = 0;
	innerThickestPoint = new double [3];
	outerThickestPoint = new double [3];
	//originIsInner = false;
}

Neighbour::Neighbour(const Neighbour& src) {
	aorta = new SimulatedShape();
	aorta = src.aorta;
	thickness = src.thickness;
	//originIsInner = src.originIsInner;

	innerThickestPoint = new double [3];
	outerThickestPoint = new double [3];
	copyPoint(innerThickestPoint , src.innerThickestPoint);
	copyPoint(outerThickestPoint , src.outerThickestPoint);
}

Neighbour& Neighbour::operator = (const Neighbour& src) {
	if (this != &src) {
		aorta = src.aorta;
		thickness = src.thickness;
		innerThickestPoint = src.innerThickestPoint;
		outerThickestPoint = src.outerThickestPoint;
		//originIsInner = src.originIsInner;
	}
	return *this;
}

Neighbour::~Neighbour() {
	delete aorta;
	delete [] innerThickestPoint;
	delete [] outerThickestPoint;
}

// setters and getters
void Neighbour::setAorta(SimulatedShape *src) {
	aorta = src;
}

SimulatedShape *Neighbour::getAorta() {
	return aorta;
}

void Neighbour::setThickness(double thickness) {
	this -> thickness = thickness;
}

double Neighbour::getThickness() {
	return thickness;
}

void Neighbour::setInnerThickestPoint(double* point) {
	copyPoint(innerThickestPoint, point);
}

double* Neighbour::getInnerThickestPoint() {
	return innerThickestPoint;
}

void Neighbour::setOuterThickestPoint(double* point) {
	copyPoint(outerThickestPoint, point);
}

double* Neighbour::getOuterThickestPoint() {
	return outerThickestPoint;
}

void Neighbour::setThicknessList(string src) {
	if (thicknessList.compare(src) != 0) {
		thicknessList.clear();
		thicknessList = src;
	}
}

string Neighbour::getThicknessList(){
	return thicknessList;
}

// member variables
double Neighbour::calculateDistance(double point1 [3], double point2 [3]) {
	double distance;
	distance = pow((point1[0] - point2[0]), 2) + pow((point1[1] - point2[1]), 2)
				+ pow((point1[2] - point2[2]), 2);

	return sqrt(distance);
}

double* Neighbour::copyPoint(double location [3], double src[3]) {
	for (int i = 0; i < 3; i++)
		location[i] = src[i];
	return location;
}

