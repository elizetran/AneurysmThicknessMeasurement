

#ifndef FRONTEND_H_
#define FRONTEND_H_

#include <string>

#include <vtkDistancePolyDataFilter.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkPointSetAlgorithm.h>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

#include "SimulatedShape.h"
#include "NormalNeighbour.h"
#include "Neighbour.h"

using namespace std;

/**
 * The Frontend class handles output of data to use. This includes writing all thicknesses
 * measured to a csv file, displaying the geometry, and outlining the location of the
 * thickest measurement, etc.
 */
class Frontend {

public:

// CONSTRUCTORS & DESTRUCTOR:

	/**
	 * Default constructor
	 */
	Frontend();

	//Frontend(const Frontend&);

	//Frontend& operator = (const Frontend);

	/**
	 * Destructor
	 */
	~Frontend();

// FUNCTIONS:

	/**
	 * Assigns a shape to be displayed
	 * @param SimulatedShape* pointer to source geometry
	 */
	void addSimulatedShape(SimulatedShape*);

	/**
	 * Displays shape in an interactive window
	 */
	void display();

	/**
	 * Write list of all thicknesses measured at each point to a csv file
	 * @param string string formatted in the form: x_inner,y_inner,z_inner,, x_outer,y_outer,z_outer,,thickness
	 */
	void writeThicknessOutput(string);

	/**
	 * Highlights line defined by 2 points passed as args
	 * @param source point to start finite line from
	 * @param target point to end finite line at
	 */
	void highlightLine(double[3], double [3]);

	/**
	 * Displays the simulated shape with thickness-based, colour-coded inner mesh
	 * @param NormalNeighbour* object with geometry and thickness data
	 */
	void addNormalThicknessDisplay(NormalNeighbour*);

	/**
	 * Presents all options of test geometry to user
	 * @param SimulatedShape* valid pointer at SimulatedShape object to build test geometry into
	 */
	void menuTestGeometry(SimulatedShape*);

	/**
	 * Requests files for wall and lumen of aorta as .off files
	 * Used for NormalNeighbour results only
	 */
	SimulatedShape* requestData();

	/**
	 * Saves outer and inner meshes of normal neighbour as a .vtp files
	 */
	void saveResults(NormalNeighbour *);

	/**
	 * Save number of points, number of polys, and thickness in pointData.csv
	 */
	void writeMeshPointData(string);


private:

// HELPER FUNCTIONS:

	/**
	 * Adds axes to orient geometry as user rotates in display
	 */
	void addAxes();

	/**
	 * Displays menu of available test geometries and return choice
	 */
	int selectShape();

	/**
	 * Reads .off files and return them as PolyData
	 * @param valid directory or filename of .off file
	 */
	vtkSmartPointer <vtkPolyDataAlgorithm> readFileOFF(char*);

	/**
	 * Finds center of inner mesh and move both inner and outer meshes to global origin
	 * @param aorta mesh to be centered
	 */
	void centerModel(SimulatedShape *);

	/**
	 * Subdivides geometry to refine inner and outer meshes
	 * @param polydata to be subdivided
	 */
	void subdivideGeometry(vtkSmartPointer <vtkPolyData>);


//MEMBER VARIABLES
	vtkSmartPointer <vtkRenderer> renderer;
	vtkSmartPointer <vtkRenderWindowInteractor> renderWindowInteractor;
	vtkSmartPointer <vtkRenderWindow> renderWindow;
	vtkSmartPointer <vtkOrientationMarkerWidget> widget;

};

#endif //FRONTEND_H_
