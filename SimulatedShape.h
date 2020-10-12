
#ifndef SIMULATEDSHAPE_H_
#define SIMULATEDSHAPE_H_

#include <vtkAppendPolyData.h>
#include <vtkClipPolyData.h>
#include <vtkCylinderSource.h>
#include <vtkImplicitFunction.h>
#include <vtkPolyData.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTransformPolyDataFilter.h>

/**
 * The SimulatedShape class builds various regular shapes to simulate the shape of an
 * aorta. These shapes will be used to test thickness-measuring algorithm  since these
 * have known thicknesses.
 */
class SimulatedShape {

public:

// CONSTRUCTORS & DESTRUCTORS

	/**
	 * Default constructor
	 */
	SimulatedShape();

	/**
	 * Copy constructor
	 */
	SimulatedShape(const SimulatedShape&);

	/**
	 * Equality operator
	 */
	SimulatedShape& operator = (const SimulatedShape&);

	/**
	 * Destructor
	 */
	~SimulatedShape();


// SETTERS & GETTERS

	/**
	 * Setter for inner mesh
	 */
	void setInner(vtkSmartPointer <vtkPolyData>);

	/**
	 * Setter for outer mesh
	 */
	void setOuter(vtkSmartPointer <vtkPolyData>);

	/**
	 * Getter for inner mesh
	 */
	const vtkSmartPointer <vtkPolyData> getInner();

	/**
	 * Getter for outer mesh
	 */
	const vtkSmartPointer <vtkPolyData> getOuter();


// FUNCTIONS: build aorta-like geometries to test on

	/**
	 * Builds test geometry with inner mesh as a cylinder and outer mesh as a cylinder
	 * Thickness is consistently 1
	 */
	void build2Cylinders();

	/**
	 * Builds test geometry with both meshes as a cylinder with a sphere aligned at its center.
	 * Outer mesh is slightly larger.
	 * Thickness is consistently 1
	 */
	void build2CylindersWith2Spheres();

	/**
	 * Builds test geometry with both meshes as a cylinder and a sphere with a skewed
	 * center in the x-axis.
	 * Outer mesh is slightly larger.
	 * Thickness is consistently 1
	 */
	void build2SkewedSpheres();

	/**
	 * Builds test geometry with inner mesh as cylinder and outer mesh as cylinder
	 * with center-aligned sphere.
	 * Expected largest, normal thickness is 3.
	 */
	void build2CylindersWithSphere();

	/**
	 * Builds test geometry with inner mesh as cylinder and outer mesh as cylinder
	 * with sphere with a skewed center of 1 in the x-axis.
	 * Expected largest, normal thickness is 4.
	 */
	void buildSkewedSphere();

	/**
	 * Builds test geometry with inner mesh as a small sphere with a skewed centre in the
	 * x-axis. Outer mesh is a cylinder with a skewed-centred sphere.
	 */
	void buildSelfIntersectingTest ();

	void trimLumen();
	//temp function

private:

// HELPER FUNCTIONS: pre-built sizes of regular shapes making up geometries

	/**
	 * Builds cylinder with a radius of 1 and height of 15
	 */
	vtkSmartPointer <vtkPolyDataAlgorithm> innerCylinder();

	/**
	 * Builds cylinder with radius of 2 and height of 15
	 */
	vtkSmartPointer <vtkPolyDataAlgorithm> outerCylinder();

	/**
	 * Builds sphere with a radius of 3
	 */
	vtkSmartPointer <vtkSphereSource> outerSphere();

	/**
	 * Builds sphere with a radius of 4
	 */
	vtkSmartPointer <vtkSphereSource> innerSphere();


// HELPER FUNCTIONS: editing and combining of pre-built regular shapes to geometries

	/**
	 * Appends 2 PolyData shapes together.
	 * @param PolyDataAlgorithm pointer to be appended
	 * @param PolyDataAlgorithm pointer to append to
	 */
	void appendPolyData (vtkSmartPointer <vtkPolyDataAlgorithm>, vtkSmartPointer <vtkPolyData>);

	/**
	 * Clips PolyData result from clipping an implicit function out of a PolyData
	 * @param ImplicitFunction defining how to clip
	 * @param PolyDataAlgorithm containing shape to be clipped
	 */
	vtkSmartPointer <vtkClipPolyData> cutOutOfShape(vtkSmartPointer <vtkImplicitFunction>,  vtkSmartPointer <vtkPolyDataAlgorithm>);

	/**
	 * cutOutOfShape replacement using cylinders as PolyData
	 * @param double radius of sphere to be used as clipping template
	 * @param double radius of cylinder to be clipped
	 */
	void addSeparatedCylinders(double, double, vtkSmartPointer <vtkPolyData>);

	/**
	 * cutOutOfShape replacement using cylinders as PolyData for skewed sphere
	 * @param double radius of sphere to be used as clipping template
	 * @param double x-axis skew of cylinder to be used as clipping template
	 * @param PolyData to be clipped
	 * @param PolyData to add clipped cylinders to
	 */
	void appendSkewedSeparatedCylinder (double, double, vtkSmartPointer<vtkPolyDataAlgorithm>, vtkSmartPointer<vtkPolyData>);


// MEMBER VARIABLES
	vtkSmartPointer <vtkPolyData> inner; /**< simulates inner wall of an aorta as a PolyData */

	vtkSmartPointer <vtkPolyData> outer; /**< simulates outer wall of an aorta as a PolyData */
};

#endif /* SIMULATESHAPE_H_ */
