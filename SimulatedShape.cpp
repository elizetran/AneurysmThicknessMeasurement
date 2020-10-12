#include <stdio.h>

#include <vtkAdaptiveSubdivisionFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkAxesActor.h>
#include <vtkCaptionActor2D.h>
#include <vtkCleanPolyData.h>
#include <vtkClipClosedSurface.h>
#include <vtkClipPolyData.h>
#include <vtkBox.h>
#include <vtkCylinder.h>
#include <vtkCylinderSource.h>
#include <vtkLineSource.h>
#include <vtkNamedColors.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkPlane.h>
#include <vtkPlaneCollection.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSphere.h>
#include <vtkSphereSource.h>
#include <vtkTextProperty.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkTubeFilter.h>


#include "SimulatedShape.h"
using namespace std;

SimulatedShape::SimulatedShape() {
	inner = vtkSmartPointer <vtkPolyData>::New();
	outer = vtkSmartPointer <vtkPolyData>::New();
}

SimulatedShape::SimulatedShape(const SimulatedShape &src) {
	outer = vtkSmartPointer <vtkPolyData>::New();
	outer -> DeepCopy(src.outer);
	inner = vtkSmartPointer <vtkPolyData>::New();
	inner -> DeepCopy(src.inner);
}

SimulatedShape& SimulatedShape::operator = (const SimulatedShape &src) {
	if (this != &src) {
		outer -> DeepCopy(src.outer);
		inner -> DeepCopy(src.inner);
	}
	return *this;
}

SimulatedShape::~SimulatedShape() {
	inner -> Delete();
	outer -> Delete();
}

// setters and getters
void SimulatedShape::setInner(vtkSmartPointer <vtkPolyData> src) {
	inner -> DeepCopy(src);
}

void SimulatedShape::setOuter(vtkSmartPointer <vtkPolyData> src) {
	outer -> DeepCopy(src);
}

const vtkSmartPointer <vtkPolyData> SimulatedShape::getOuter() {
	return outer;
}

const vtkSmartPointer <vtkPolyData> SimulatedShape::getInner() {
	return inner;
}

// public functions
void SimulatedShape::build2Cylinders() {
	appendPolyData(innerCylinder(), inner);
	appendPolyData(outerCylinder(), outer);
}

void SimulatedShape::build2CylindersWith2Spheres() {
	vtkSmartPointer <vtkCylinder> cylCut = vtkSmartPointer <vtkCylinder> :: New();
	cylCut -> SetRadius(2);
	appendPolyData(cutOutOfShape(cylCut, outerSphere()), outer);

	addSeparatedCylinders(4, 2, outer);

	cylCut -> SetRadius(1);
	appendPolyData(cutOutOfShape(cylCut, innerSphere()), inner);

	addSeparatedCylinders(3, 1, inner);
}

void SimulatedShape::build2CylindersWithSphere() {
	vtkSmartPointer <vtkCylinder> cylCut = vtkSmartPointer <vtkCylinder> :: New();
	cylCut -> SetRadius(2);
	appendPolyData(cutOutOfShape(cylCut, outerSphere()), outer);

	addSeparatedCylinders(4, 2, outer);
	appendPolyData(innerCylinder(), inner);

// *** not clipping outer & unsure why
// *** apparently vtkClipPolyData might have problems with clipping vtkCylinderSource
// *** http://vtk.1045678.n5.nabble.com/ClipPolyData-fails-to-clip-CylinderSource-td5740148.html#none
// *** https://stackoverflow.com/questions/21549620/vtk-clipping-does-not-always-work-with-cylinders-spheres-cones
//	sphereCut -> SetRadius(4);
//	cutOutOfShape(sphereCut, sphere);
}

void SimulatedShape::buildSkewedSphere() {
	vtkSmartPointer <vtkCylinder> cylCut = vtkSmartPointer <vtkCylinder> :: New();
	cylCut -> SetRadius(2);
	vtkSmartPointer <vtkSphere> sphereCut = vtkSmartPointer <vtkSphere> :: New();
	sphereCut -> SetRadius(4);

	vtkSmartPointer <vtkTransform> transform = vtkSmartPointer <vtkTransform> :: New();
	vtkSmartPointer <vtkTransformPolyDataFilter> transPD = vtkSmartPointer <vtkTransformPolyDataFilter> :: New();
	transform -> Translate(1, 0, 0);
	transPD -> SetInputConnection(outerSphere() -> GetOutputPort());
	transPD -> SetTransform(transform);
	transPD -> Update();

	appendPolyData(cutOutOfShape(cylCut, transPD), outer);
	appendSkewedSeparatedCylinder(4, 1, outerCylinder(), outer);
	appendPolyData(innerCylinder(), inner);
}

void SimulatedShape::build2SkewedSpheres() {
	vtkSmartPointer <vtkCylinder> cylCut = vtkSmartPointer <vtkCylinder> :: New();
	cylCut -> SetRadius(2);

	vtkSmartPointer <vtkTransform> transform = vtkSmartPointer <vtkTransform> :: New();
	vtkSmartPointer <vtkTransformPolyDataFilter> transPD = vtkSmartPointer <vtkTransformPolyDataFilter> :: New();
	transform -> Translate(1, 0, 0);
	transPD -> SetInputConnection(outerSphere() -> GetOutputPort());
	transPD -> SetTransform(transform);
	transPD -> Update();

	appendPolyData(cutOutOfShape(cylCut, transPD), outer);
	appendSkewedSeparatedCylinder(4, 1, outerCylinder(), outer);

	cylCut -> SetRadius(1);
	transPD -> SetInputConnection(innerSphere() -> GetOutputPort());
	appendPolyData(cutOutOfShape(cylCut, transPD), inner);

	appendSkewedSeparatedCylinder(3, 1, innerCylinder(), inner);
}

void SimulatedShape::buildSelfIntersectingTest() {
	vtkSmartPointer <vtkCylinder> cylCut = vtkSmartPointer <vtkCylinder> :: New();
	cylCut -> SetRadius(2);
	vtkSmartPointer <vtkSphere> sphereCut = vtkSmartPointer <vtkSphere> :: New();
	sphereCut -> SetRadius(4);

	vtkSmartPointer <vtkTransform> transform = vtkSmartPointer <vtkTransform> :: New();
	vtkSmartPointer <vtkTransformPolyDataFilter> transPD = vtkSmartPointer <vtkTransformPolyDataFilter> :: New();
	transform -> Translate(1, 0, 0);
	transPD -> SetInputConnection(outerSphere() -> GetOutputPort());
	transPD -> SetTransform(transform);
	transPD -> Update();

	appendPolyData(cutOutOfShape(cylCut, transPD), outer);
	appendSkewedSeparatedCylinder(4, 1, outerCylinder(), outer);

	vtkSmartPointer <vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
	sphere -> SetRadius(1);
	sphere -> SetThetaResolution(100);
	sphere -> SetPhiResolution(100);

	transform -> Translate(2.5, 0, 0);
	transPD -> SetInputConnection(sphere -> GetOutputPort());
	transPD -> Update();
	appendPolyData(transPD, inner);
}

// helper functions

void SimulatedShape::appendPolyData (vtkSmartPointer <vtkPolyDataAlgorithm> input, vtkSmartPointer <vtkPolyData> mesh) {
	vtkSmartPointer <vtkAppendPolyData> appendFilter = vtkSmartPointer <vtkAppendPolyData> :: New();
	appendFilter -> AddInputConnection(input -> GetOutputPort());
	appendFilter -> AddInputData(mesh);
	appendFilter -> Update();

	vtkSmartPointer <vtkTriangleFilter> triangleFilter = vtkSmartPointer <vtkTriangleFilter> :: New();
	triangleFilter -> SetInputConnection(appendFilter -> GetOutputPort());
	triangleFilter -> Update();

	vtkSmartPointer <vtkAdaptiveSubdivisionFilter> subdivideFilter = vtkSmartPointer <vtkAdaptiveSubdivisionFilter> :: New();
	subdivideFilter -> SetInputConnection(triangleFilter -> GetOutputPort());
	subdivideFilter -> SetMaximumEdgeLength(0.5);
	subdivideFilter -> Update();

	vtkSmartPointer <vtkCleanPolyData> cleanFilter = vtkSmartPointer <vtkCleanPolyData> :: New();
	cleanFilter -> SetInputConnection (subdivideFilter -> GetOutputPort());
	cleanFilter -> Update();

	mesh -> ShallowCopy(cleanFilter -> GetOutput());
}

vtkSmartPointer <vtkClipPolyData> SimulatedShape::cutOutOfShape(vtkSmartPointer <vtkImplicitFunction> hole, vtkSmartPointer <vtkPolyDataAlgorithm> shape) {
	vtkSmartPointer <vtkClipPolyData> cut = vtkSmartPointer <vtkClipPolyData> :: New();
	cut -> SetInputConnection(shape -> GetOutputPort());
	cut -> SetClipFunction(hole);
	cut -> Update();

	return cut;
}

vtkSmartPointer <vtkPolyDataAlgorithm> SimulatedShape::innerCylinder() {
	vtkSmartPointer <vtkCylinderSource> innerCyl = vtkSmartPointer <vtkCylinderSource>::New();
	innerCyl -> SetRadius(1);
	innerCyl -> SetHeight(15);
	innerCyl -> SetResolution(100);
	innerCyl -> CappingOff();

//	vtkSmartPointer <vtkTriangleFilter> triangleFilter = vtkSmartPointer <vtkTriangleFilter> :: New();
//	triangleFilter -> SetInputConnection(innerCyl -> GetOutputPort());
//	triangleFilter -> Update();
//
//	vtkSmartPointer <vtkAdaptiveSubdivisionFilter> subdivideFilter = vtkSmartPointer <vtkAdaptiveSubdivisionFilter> :: New();
//	subdivideFilter -> SetInputConnection(triangleFilter -> GetOutputPort());
//	subdivideFilter -> SetMaximumEdgeLength(0.5);
//	//subdivideFilter -> SetMaximumNumberOfTriangles(10000);
//	subdivideFilter -> Update();

	return innerCyl;
}

vtkSmartPointer <vtkPolyDataAlgorithm> SimulatedShape::outerCylinder() {
	vtkSmartPointer <vtkCylinderSource> outerCyl = vtkSmartPointer <vtkCylinderSource>::New();
	outerCyl -> SetRadius(2);
	outerCyl -> SetHeight(15);
	outerCyl -> SetResolution(100);
	outerCyl -> CappingOff();

//	vtkSmartPointer <vtkTriangleFilter> triangleFilter = vtkSmartPointer <vtkTriangleFilter> :: New();
//	triangleFilter -> SetInputConnection(outerCyl -> GetOutputPort());
//	triangleFilter -> Update();
//
//	vtkSmartPointer <vtkAdaptiveSubdivisionFilter> subdivideFilter = vtkSmartPointer <vtkAdaptiveSubdivisionFilter> :: New();
//	subdivideFilter -> SetInputConnection(triangleFilter -> GetOutputPort());
//	subdivideFilter -> SetMaximumEdgeLength(0.5);
////	subdivideFilter -> SetMaximumNumberOfTriangles(10000);
//	subdivideFilter -> Update();

	return outerCyl;
}

vtkSmartPointer <vtkSphereSource> SimulatedShape::outerSphere() {
	vtkSmartPointer <vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
	sphere -> SetRadius(4);
	sphere -> SetThetaResolution(100);
	sphere -> SetPhiResolution(100);
	return sphere;
}

vtkSmartPointer <vtkSphereSource> SimulatedShape::innerSphere() {
	vtkSmartPointer <vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
	sphere -> SetRadius(3);
	sphere -> SetThetaResolution(100);
	sphere -> SetPhiResolution(100);
	return sphere;
}

void SimulatedShape::addSeparatedCylinders(double sphereRad, double cylRad, vtkSmartPointer <vtkPolyData> mesh) {
	vtkSmartPointer <vtkTransform> outerTrans = vtkSmartPointer <vtkTransform> :: New();
	vtkSmartPointer <vtkTransformPolyDataFilter> transPD = vtkSmartPointer <vtkTransformPolyDataFilter> :: New();
	vtkSmartPointer <vtkCylinderSource> buildingCyl = vtkSmartPointer <vtkCylinderSource> :: New();

	double calcHeight = (15 -2*sqrt(pow(sphereRad, 2) - pow(cylRad, 2)))/2;
	double calcTrans = sqrt(pow(sphereRad, 2) - pow(cylRad, 2)) + calcHeight/2;

	buildingCyl -> SetRadius(cylRad);
	buildingCyl -> SetResolution(100);
	buildingCyl -> CappingOff();
	buildingCyl -> SetHeight(calcHeight);

	outerTrans -> Translate(0, calcTrans, 0);
	transPD -> SetTransform(outerTrans);
	transPD -> SetInputConnection(buildingCyl -> GetOutputPort());
	appendPolyData(transPD, mesh);

	outerTrans -> Translate(0, -2*calcTrans, 0);
	transPD -> SetTransform(outerTrans);
	transPD -> SetInputConnection(buildingCyl -> GetOutputPort());
	appendPolyData(transPD, mesh);
}

void SimulatedShape::appendSkewedSeparatedCylinder(double radius, double offset, vtkSmartPointer <vtkPolyDataAlgorithm> shape, vtkSmartPointer <vtkPolyData> mesh) {
	vtkSmartPointer<vtkPlane> plane1 = vtkSmartPointer<vtkPlane>::New();
	plane1->SetOrigin(0, (sqrt(pow(radius, 2) - pow((radius - offset), 2)) + sqrt(pow(radius, 2) - pow(offset, 2)))/2, 0);

	// hard-coded normal
	plane1->SetNormal(-1, (3*sqrt(7) + sqrt(15))/4, 0);

	appendPolyData(cutOutOfShape(plane1, shape), mesh);

	vtkSmartPointer<vtkPlane> plane2 = vtkSmartPointer<vtkPlane>::New();
	plane2->SetOrigin(0, (sqrt(pow(radius, 2) - pow((radius - offset), 2)) + sqrt(pow(radius, 2) - pow(offset, 2)))/-2, 0);
	plane2->SetNormal(1, (sqrt(pow(radius, 2) - pow((radius - offset), 2)) + 3*sqrt(pow(radius, 2) - pow(offset, 2)))/4, 0);

	vtkSmartPointer <vtkClipPolyData> tempCut = vtkSmartPointer <vtkClipPolyData> :: New();
	tempCut = cutOutOfShape(plane2, shape);
	tempCut -> InsideOutOn();
	appendPolyData(tempCut, mesh);
}

void SimulatedShape::trimLumen() {
	vtkSmartPointer <vtkPlane> plane = vtkSmartPointer <vtkPlane> :: New();
	plane -> SetOrigin(0, 0,115);
	plane -> SetNormal(0, 0, 1);

	vtkSmartPointer <vtkClipPolyData> cut = vtkSmartPointer <vtkClipPolyData> :: New();
	cut -> SetInputData(getInner());
	cut -> SetClipFunction(plane);
	cut -> InsideOutOn();
	cut -> Update();

	setInner(cut->GetOutput());
}

