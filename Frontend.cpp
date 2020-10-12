
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>

#include <vtkAdaptiveSubdivisionFilter.h>
#include <vtkAxesActor.h>
#include <vtkCenterOfMass.h>
#include <vtkColor.h>
#include <vtkDistancePolyDataFilter.h>
#include <vtkLineSource.h>
#include <vtkLookupTable.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkScalarBarActor.h>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSTLReader.h>
#include <vtkTextProperty.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkXMLPolyDataWriter.h>

#include "Frontend.h"

using namespace std;

Frontend::Frontend() {
	renderer = vtkSmartPointer <vtkRenderer> :: New();
	renderWindow = vtkSmartPointer <vtkRenderWindow> :: New();
	renderWindowInteractor = vtkSmartPointer <vtkRenderWindowInteractor> :: New();
	widget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();

	double* rgb;
	renderer -> SetBackground(1.0, 1.0, 1.0);
	renderer -> ResetCamera(-7,7,-7,7,-7,7);
	renderWindow-> AddRenderer(renderer);
	renderWindowInteractor -> SetRenderWindow(renderWindow);

	addAxes();
}

Frontend::~Frontend() {
	renderer -> Delete();
	renderWindowInteractor -> Delete();
	renderWindow -> Delete();
	widget -> Delete();
}

void Frontend::addSimulatedShape(SimulatedShape *src) {

	vtkSmartPointer <vtkPolyDataMapper> outerMapper = vtkSmartPointer <vtkPolyDataMapper>::New();
	vtkSmartPointer <vtkPolyDataMapper> innerMapper = vtkSmartPointer <vtkPolyDataMapper>::New();

	outerMapper -> SetInputData (src -> getOuter());
	innerMapper -> SetInputData (src -> getInner());
	outerMapper -> ScalarVisibilityOn();
	innerMapper -> ScalarVisibilityOn();

	vtkSmartPointer <vtkActor> outerActor = vtkSmartPointer <vtkActor>::New();
	vtkSmartPointer <vtkActor> innerActor = vtkSmartPointer <vtkActor>::New();

	outerActor -> SetMapper(outerMapper);
	//outerActor -> GetProperty() -> SetColor(0.0, 1, 0.0);
	outerActor -> GetProperty() -> SetOpacity(0.5);

	innerActor -> SetMapper(innerMapper);
	innerActor -> GetProperty() -> SetColor(0.0, 0.0, 1.0);
	innerActor -> GetProperty() -> SetOpacity(0.5);

	renderer -> AddActor(outerActor);
	renderer -> AddActor(innerActor);
}

void Frontend::display() {
	renderWindow -> Render();
	renderWindowInteractor -> Start();
}

void Frontend::highlightLine(double point1 [3], double point2 [3]) {
	vtkSmartPointer <vtkLineSource> line = vtkSmartPointer <vtkLineSource> :: New();
	line -> SetPoint1(point1);
	line -> SetPoint2(point2);
	line -> Update();

	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(line->GetOutputPort());

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetLineWidth(2);
	actor -> GetProperty() -> SetColor(1, 0, 0);
	actor -> GetProperty() -> SetOpacity(1);

	renderer -> AddActor(actor);
}

// helper functions
void Frontend::addAxes() {
	vtkSmartPointer <vtkAxesActor> axes = vtkSmartPointer <vtkAxesActor> :: New();

	widget -> SetOutlineColor( 0.9300, 0.5700, 0.1300 );
	widget -> SetOrientationMarker( axes );
	widget -> SetInteractor( renderWindowInteractor );
	widget -> SetViewport( 0.0, 0.0, 0.4, 0.4 );
	widget -> SetEnabled( 1 );
	widget -> InteractiveOn();
}

void Frontend::writeThicknessOutput(string input) {
	FILE *file = fopen("thicknessList.csv", "w");

	fprintf(file, ",Origin Point,,,,Target Point,,,Thickness\n");
	fprintf(file, "x,y,z,,x,y,z\n");
	fprintf(file, "%s\n", input.c_str());

	fclose(file);
}

void Frontend::addNormalThicknessDisplay(NormalNeighbour *src) {
	vtkSmartPointer <vtkPolyDataMapper> targetMapper = vtkSmartPointer <vtkPolyDataMapper>::New();
	vtkSmartPointer <vtkPolyData> temp = vtkSmartPointer <vtkPolyData> :: New();

	src -> getAorta() -> getInner() -> GetPointData() -> SetActiveScalars("Distance");
	src -> getAorta() -> getOuter() -> GetPointData() -> SetActiveScalars("Distance");

	if(src -> getOriginIsInner()) {
		targetMapper -> SetInputData (src -> getAorta() -> getOuter());
		temp -> DeepCopy(src -> getAorta() -> getInner());

	} else {
		targetMapper -> SetInputData (src -> getAorta() -> getInner());
		temp -> DeepCopy(src -> getAorta() -> getOuter());
	}
	vtkSmartPointer <vtkActor> targetActor = vtkSmartPointer <vtkActor>::New();
	targetActor -> SetMapper(targetMapper);

	vtkSmartPointer <vtkPolyDataMapper> originMapper = vtkSmartPointer <vtkPolyDataMapper> :: New();
	originMapper -> SetInputData(temp);

	vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
	scalarBar->SetTitle("Thickness");
	scalarBar -> GetTitleTextProperty() -> SetLineSpacing(3);
	scalarBar->SetNumberOfLabels(4);

	vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
	lut -> SetTableRange (
			temp -> GetPointData() -> GetScalars() -> GetRange()[0],
			temp -> GetPointData() -> GetScalars() -> GetRange()[1]);
	lut->Build();

	scalarBar -> SetLookupTable(lut);

	originMapper -> SetLookupTable(lut);
	originMapper -> UseLookupTableScalarRangeOn();

	vtkSmartPointer<vtkActor> originActor = vtkSmartPointer<vtkActor>::New();
	originActor -> SetMapper(originMapper);

	targetActor -> GetProperty() -> SetOpacity(0.5);
	originActor -> GetProperty() -> SetOpacity(0.5);

	renderer -> AddActor(targetActor);
	renderer -> AddActor(originActor);
	renderer -> AddActor2D(scalarBar);
}

//void Frontend::addNormalThicknessDisplay(NormalNeighbour *src) {
//	vtkSmartPointer <vtkPolyDataMapper> outerMapper = vtkSmartPointer <vtkPolyDataMapper>::New();
//	vtkSmartPointer <vtkPolyDataMapper> innerMapper = vtkSmartPointer <vtkPolyDataMapper>::New();
//	vtkSmartPointer <vtkActor> outerActor = vtkSmartPointer <vtkActor>::New();
//	vtkSmartPointer <vtkActor> innerActor = vtkSmartPointer <vtkActor>::New();
//	vtkSmartPointer<vtkLookupTable> outerLut = vtkSmartPointer<vtkLookupTable>::New();
//	vtkSmartPointer<vtkLookupTable> innerLut = vtkSmartPointer<vtkLookupTable>::New();
//	vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
//
//	src -> getAorta() -> getInner() -> GetPointData() -> SetActiveScalars("Distance");
//	src -> getAorta() -> getOuter() -> GetPointData() -> SetActiveScalars("Distance");
//
//	outerMapper -> SetInputData(src -> getAorta() -> getOuter());
//	innerMapper -> SetInputData(src -> getAorta() -> getInner());
//
//	outerLut -> SetTableRange (
//			src -> getAorta() -> getOuter() -> GetPointData() -> GetScalars() -> GetRange()[0],
//			src -> getAorta() -> getOuter() -> GetPointData() -> GetScalars() -> GetRange()[1]);
//	outerLut->Build();
//	innerLut -> SetTableRange (
//			src -> getAorta() -> getInner() -> GetPointData() -> GetScalars() -> GetRange()[0],
//			src -> getAorta() -> getInner() -> GetPointData() -> GetScalars() -> GetRange()[1]);
//	innerLut->Build();
//
//	scalarBar->SetTitle("Thickness");
//	scalarBar->SetNumberOfLabels(4);
//	scalarBar -> SetLookupTable(innerLut);
//
//	outerMapper -> SetLookupTable(outerLut);
//	outerMapper -> UseLookupTableScalarRangeOn();
//
//	innerMapper -> SetLookupTable(innerLut);
//	innerMapper -> UseLookupTableScalarRangeOn();
//
//	outerActor -> SetMapper(outerMapper);
//	innerActor -> SetMapper(innerMapper);
//
//	outerActor -> GetProperty() -> SetOpacity(0.5);
//	innerActor -> GetProperty() -> SetOpacity(0.5);
//
//	renderer -> AddActor(innerActor);
//	renderer -> AddActor(outerActor);
//	renderer -> AddActor2D(scalarBar);
//}

int Frontend::selectShape() {
	int choice = -1;
	cout << "Please select a test geometry:\n" <<
			"1. Cylinders\n" <<
			"2. Spheres in inner and outer mesh\n" <<
			"3. Sphere in outer mesh\n" <<
			"4. Skewed spheres in inner and outer mesh\n" <<
			"5. Skewed sphere in outer mesh\n" <<
			"6. Self-intersecting measurement test mesh" << endl;
	cin >> choice;
	return choice;
}

//int Frontend::selectOrigin() {
//	int choice = -1;
//	cout << "Please select a mesh to display distance on:\n" <<
//			"1. Inner mesh\n" <<
//			"2. Outer mesh" << endl;
//	cin >> choice;
//	return choice;
//}

void Frontend::menuTestGeometry(SimulatedShape *simShape) {
	int choice = 0;
	do {
		choice = selectShape();
		switch(choice) {
		case 1:
			simShape -> build2Cylinders();
			cout << "Building 2 cylinders..." << endl;
			break;
		case 2:
			simShape -> build2CylindersWith2Spheres();
			cout << "Building 2 cylinders with 2 spheres..." << endl;
			break;
		case 3:
			simShape -> build2CylindersWithSphere();
			cout << "Building 2 cylinders with 1 sphere in outer mesh..." << endl;
			break;
		case 4:
			simShape -> build2SkewedSpheres();
			cout << "Building 2 cylinders with 2 skewed spheres..." << endl;
			break;
		case 5:
			simShape -> buildSkewedSphere();
			cout << "Building 2 cylinders with 1 skewed sphere in outer mesh..." << endl;
			break;
		case 6:
			simShape -> buildSelfIntersectingTest();
			cout << "Building shape to test thickness measurement intersecting origin mesh..." << endl;
			break;
		default:
			cout << "Invalid choice! Please try again." << endl;
		}
	} while (choice < 1 || choice > 6);
}

//void Frontend::menuOriginTarget(NormalNeighbour *normN) {
//	int choice = 0;
//	do {
//		choice = selectOrigin();
//		switch(choice) {
//		case 1:
//			normN -> setOriginIsInnerToTrue();
//			cout << "Inner mesh is origin..." << endl;
//			break;
//		case 2:
//			normN -> setOriginIsInnerToFalse();
//			cout << "Outer mesh is origin..." << endl;
//			break;
//		default:
//			cout << "Invalid choice! Please try again." << endl;
//		}
//	} while (choice < 1 || choice > 2);
//}

vtkSmartPointer <vtkPolyDataAlgorithm> Frontend::readFileOFF(char* fileName)
{
	vtkSmartPointer <vtkTriangleFilter> triangles = vtkSmartPointer <vtkTriangleFilter> :: New();

	vtkSmartPointer <vtkPolyData> poly = vtkPolyData::New();
	vtkSmartPointer <vtkPoints> points = vtkPoints::New();
	vtkSmartPointer <vtkCellArray> polys = vtkCellArray::New();
	float p[3];

	FILE *in = fopen(fileName,"r");

	if(!in)
		return triangles;

	int nPts,nPolys;

	fscanf(in,"OFF\n%d %d 0\n",&nPts,&nPolys);


	// read the points from the file float p[3];
	int i,j;


	for(i=0;i<nPts;i++)
	{
		fscanf(in,"%f %f %f\n",&p[0],&p[1],&p[2]);
		points->InsertPoint(i,p);
	}

	for(i=0;i<nPolys;i++)
	{
		int nVertices;
		fscanf(in,"%d ",&nVertices);
		polys->InsertNextCell(nVertices);
		int vertex_id;

		for(j=0; j < nVertices; j++)
		{
			fscanf(in,"%d ",&vertex_id);
			polys->InsertCellPoint(vertex_id);
		}
	}

	poly->SetPoints(points);
	points->Delete();
	poly->SetPolys(polys);
	polys->Delete();

	triangles->SetInputData(poly);
	triangles->Update();

	return triangles;
}

SimulatedShape* Frontend::requestData() {

	SimulatedShape *temp = new SimulatedShape();
	char *wallFile = new char[256];
	char *lumenFile = new char[256];

	cout << "Please input filename for aorta wall:" << endl;
	cin >> wallFile;

	cout << "Please input filename for aorta lumen:" << endl;
	cin >> lumenFile;

	temp -> setOuter(readFileOFF(wallFile) -> GetOutput());
	temp -> setInner(readFileOFF(lumenFile) -> GetOutput());

	centerModel(temp);

	cout << "Files read..." << endl;

//	subdivideGeometry(temp -> getOuter());
//	subdivideGeometry(temp -> getInner());

	return temp;
}

void Frontend::centerModel(SimulatedShape *input) {
	double *centerPoint;

	vtkSmartPointer <vtkCenterOfMass> centerFilter = vtkSmartPointer <vtkCenterOfMass> :: New();
	centerFilter -> SetInputData(input -> getInner());
	centerFilter -> SetUseScalarsAsWeights(false);
	centerFilter -> Update();
	centerPoint = centerFilter -> GetCenter();

	vtkSmartPointer <vtkTransform> transform = vtkSmartPointer <vtkTransform> :: New();
	transform -> Translate(-centerPoint[0], -centerPoint[1], -centerPoint[2]);

	vtkSmartPointer <vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer <vtkTransformPolyDataFilter> :: New();

	transformFilter -> SetInputData(input -> getInner());
	transformFilter -> SetTransform(transform);
	transformFilter -> Update();
	input -> setInner(transformFilter -> GetOutput());

	transformFilter -> SetInputData(input -> getOuter());
	transformFilter -> Update();
	input -> setOuter(transformFilter -> GetOutput());
}

void Frontend::saveResults(NormalNeighbour *normN) {
	vtkSmartPointer <vtkXMLPolyDataWriter> writer = vtkSmartPointer <vtkXMLPolyDataWriter> :: New();
	writer -> SetInputData(normN -> getAorta() -> getOuter());
	writer -> SetFileName("wallThickness.vtp");
	writer -> Update();

	writer -> SetInputData(normN -> getAorta() -> getInner());
	writer -> SetFileName("lumenThickness.vtp");
	writer -> Update();
}

void Frontend::subdivideGeometry(vtkSmartPointer <vtkPolyData> mesh) {
	vtkSmartPointer <vtkAdaptiveSubdivisionFilter> subdivideFilter = vtkSmartPointer <vtkAdaptiveSubdivisionFilter> :: New();
	subdivideFilter -> SetInputData(mesh);
	subdivideFilter -> SetMaximumEdgeLength(0.1);
	subdivideFilter -> Update();

	mesh -> DeepCopy(subdivideFilter -> GetOutput());
}

void Frontend::writeMeshPointData(string input) {
	FILE *file = fopen("pointData.csv", "w");

	fprintf(file, "Inner Pts,Inner Polys,Outer Pts,Outer Polys,Mean Thickness,Median Thickness,Max Thickness, Inner SA\n");
	fprintf(file, "%s\n", input.c_str());

	fclose(file);
}


