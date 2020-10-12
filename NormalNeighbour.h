// https://vtk.org/Wiki/VTK/Examples/Cxx/PolyData/PolyDataCellNormals
// https://vtk.org/Wiki/VTK/Examples/Cxx/PolyData/PolyDataExtractNormals

#ifndef NEIGHBOUR_NORMALNEIGHBOUR_H_
#define NEIGHBOUR_NORMALNEIGHBOUR_H_

#include "Neighbour.h"
#include "SimulatedShape.h"

#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkPolyDataNormals.h>
#include <vtkSmartPointer.h>

/**
 * Unused class, besides function: findHausdorffThickest()
 */
class NormalNeighbour : public Neighbour {
public:
	NormalNeighbour();
	NormalNeighbour(const NormalNeighbour&);
	NormalNeighbour& operator = (const NormalNeighbour&);
	~NormalNeighbour();

	void setPointNormals(vtkSmartPointer <vtkPolyDataNormals>);
	void setCellNormals(vtkSmartPointer <vtkPolyDataNormals>);
	void setFloatArray(vtkSmartPointer<vtkFloatArray>);
	void setDoubleArray(vtkSmartPointer<vtkDoubleArray>);
	void setOriginIsInnerToTrue();
	void setOriginIsInnerToFalse();

	vtkSmartPointer <vtkPolyDataNormals> getPointNormals();
	vtkSmartPointer <vtkPolyDataNormals> getCellNormals();
	vtkSmartPointer<vtkFloatArray> getFloatArray();
	vtkSmartPointer<vtkDoubleArray> getDoubleArray();
	bool getOriginIsInner();
	double getMeanThickness();
	double getMedianThickness();
	double getInnerSurfaceArea();

	double findThickest(vtkSmartPointer <vtkPolyData>, vtkSmartPointer <vtkPolyData>, double);
	void findCurvatureThickest();

	/**
	 * Uses the Hausdorff class to calculate a thicknes for each point.
	 * Inner mesh is origin.
	 */
	void findHausdorffThickest();

	/**
	 * Reduces the number of points in aorta mesh. Used for point density testing.
	 */
	void decimate(double);

	/**
	 * Increases the number of points in aorta mesh. Used for point density testing.
	 */
	void subdivide(double);

private:
	void findIntersectingCells(double [3], double [3], vtkSmartPointer <vtkPolyData>, double);
	void calculateAllThicknesses(vtkSmartPointer <vtkIdList>, double [3], vtkSmartPointer <vtkPolyData>);
	vtkSmartPointer <vtkPolyDataNormals> findPointNormals(vtkSmartPointer <vtkPolyData>);
	vtkSmartPointer <vtkPolyDataNormals> findCellNormals(vtkSmartPointer <vtkPolyData> );
	bool isCellNormals(vtkPolyData*);
	bool passesThroughLumen(double point1 [3], double point2 [3]);

	/**
	 * Mean, median, and surface area all used for assessing the affect of point density
	 * on thickness calculations.
	 */
	void calculateMeanThickness();
	void calculateMedianThickness();
	void calculateSurfaceArea();

	vtkSmartPointer <vtkPolyDataNormals> pointNormals;
	vtkSmartPointer <vtkPolyDataNormals> cellNormals;
	vtkSmartPointer<vtkFloatArray> floatArray;
	vtkSmartPointer<vtkDoubleArray> doubleArray;
	bool originIsInner;
	double meanThickness;
	double medianThickness;
	double innerSurfaceArea;

};

#endif /* NORMALNEIGHBOUR_H_ */
