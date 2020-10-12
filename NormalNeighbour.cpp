
#include <vtkAdaptiveSubdivisionFilter.h>
#include <vtkCellCenters.h>
#include <vtkCellLocator.h>
#include <vtkCellData.h>
#include <vtkCurvatures.h>
#include <vtkDataArray.h>
#include <vtkDecimatePro.h>
#include <vtkGlyph3D.h>
#include <vtkHausdorffDistancePointSetFilter.h>
#include <vtkLineSource.h>
#include <vtkMassProperties.h>
#include <vtkPCACurvatureEstimation.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkSmartPointer.h>
#include <vtkSortDataArray.h>
#include <vtkTriangle.h>

#include "NormalNeighbour.h"
#include "Hausdorff.h"


NormalNeighbour::NormalNeighbour()
: Neighbour() {
	pointNormals = vtkSmartPointer <vtkPolyDataNormals> :: New();
	cellNormals = vtkSmartPointer <vtkPolyDataNormals> :: New();
	floatArray = vtkSmartPointer<vtkFloatArray> :: New();
	doubleArray = vtkSmartPointer<vtkDoubleArray> :: New();
	originIsInner = true;
	meanThickness = 0;
	medianThickness = 0;
	innerSurfaceArea = 0;
}

NormalNeighbour::NormalNeighbour(const NormalNeighbour &src)
: Neighbour(src) {
	pointNormals = vtkSmartPointer <vtkPolyDataNormals> :: New();
	cellNormals = vtkSmartPointer <vtkPolyDataNormals> :: New();
	floatArray = vtkSmartPointer<vtkFloatArray> :: New();
	doubleArray = vtkSmartPointer<vtkDoubleArray> :: New();

	pointNormals = src.pointNormals;
	cellNormals = src.cellNormals;
	floatArray = src.floatArray;
	doubleArray = src.doubleArray;
	originIsInner = src.originIsInner;
	meanThickness = src.meanThickness;
	medianThickness = src.medianThickness;
	innerSurfaceArea =src.innerSurfaceArea;
}

NormalNeighbour& NormalNeighbour::operator = (const NormalNeighbour& src) {
	if(this != &src) {
		pointNormals -> Delete();
		cellNormals -> Delete();
		floatArray -> Delete();
		doubleArray -> Delete();

		Neighbour::operator = (src);

		pointNormals = src.pointNormals;
		cellNormals = src.cellNormals;
		floatArray = src.floatArray;
		doubleArray = src.doubleArray;
		originIsInner = src.originIsInner;
		meanThickness = src.meanThickness;
		medianThickness = src.medianThickness;
		innerSurfaceArea =src.innerSurfaceArea;
	}
	return *this;
}

NormalNeighbour::~NormalNeighbour() {
	pointNormals -> Delete();
	cellNormals -> Delete();
	floatArray -> Delete();
	doubleArray -> Delete();
}

// setters
void NormalNeighbour::setPointNormals(vtkSmartPointer <vtkPolyDataNormals> src) {
	pointNormals -> Delete();
	pointNormals = src;
}

void NormalNeighbour::setCellNormals(vtkSmartPointer <vtkPolyDataNormals> src) {
	cellNormals -> Delete();
	cellNormals = src;
}

void NormalNeighbour::setFloatArray(vtkSmartPointer<vtkFloatArray> src) {
	floatArray = src;
}

void NormalNeighbour::setDoubleArray(vtkSmartPointer<vtkDoubleArray> src) {
	doubleArray = src;
}

void NormalNeighbour::setOriginIsInnerToTrue() {
	originIsInner = true;
}

void NormalNeighbour::setOriginIsInnerToFalse() {
	originIsInner = false;
}


// getters
vtkSmartPointer <vtkPolyDataNormals> NormalNeighbour::getPointNormals (){
	return pointNormals;
}

vtkSmartPointer <vtkPolyDataNormals> NormalNeighbour::getCellNormals() {
	return cellNormals;
}

vtkSmartPointer<vtkFloatArray> NormalNeighbour::getFloatArray() {
	return floatArray;
}

vtkSmartPointer<vtkDoubleArray> NormalNeighbour::getDoubleArray() {
	return doubleArray;
}

bool NormalNeighbour::getOriginIsInner() {
	return originIsInner;
}

double NormalNeighbour::getMeanThickness() {
	return meanThickness;
}

double NormalNeighbour::getMedianThickness() {
	return medianThickness;
}

double NormalNeighbour::getInnerSurfaceArea() {
	return innerSurfaceArea;
}

// member functions


// helper functions
vtkSmartPointer <vtkPolyDataNormals> NormalNeighbour::findPointNormals(vtkSmartPointer <vtkPolyData> input) {
	vtkSmartPointer <vtkPolyDataNormals> normals = vtkSmartPointer <vtkPolyDataNormals> :: New();
	normals -> SetInputData(input);
	normals -> ComputePointNormalsOn();
	normals -> ComputeCellNormalsOff();
	normals -> Update();

	return normals;
}

vtkSmartPointer <vtkPolyDataNormals> NormalNeighbour::findCellNormals(vtkSmartPointer <vtkPolyData> input) {
	vtkSmartPointer <vtkPolyDataNormals> normals = vtkSmartPointer <vtkPolyDataNormals> :: New();
	normals -> SetInputData(input);
	normals -> ComputePointNormalsOff();
	normals -> ComputeCellNormalsOn();
	normals -> Update();

	return normals;
}

// reference = outer ... current = inner
// change all outer to target
double NormalNeighbour::findThickest(vtkSmartPointer <vtkPolyData> origin, vtkSmartPointer <vtkPolyData> target, double scalarParam) {
	vtkSmartPointer <vtkPolyData> copyInner = vtkSmartPointer <vtkPolyData> :: New();
	copyInner -> DeepCopy(origin);

	vtkSmartPointer <vtkPolyDataNormals> normals = vtkSmartPointer <vtkPolyDataNormals> :: New();
	bool hasCellNormals = isCellNormals(origin);
	if(!hasCellNormals)
		normals = findCellNormals(origin);
	else
		normals -> SetInputData(origin);

	origin -> DeepCopy(normals -> GetOutput());
	//normals -> FlipNormalsOn();
	normals -> Update();

	vtkSmartPointer <vtkCellCenters> originCenters = vtkSmartPointer <vtkCellCenters> :: New();
	originCenters -> SetInputConnection(normals -> GetOutputPort());
	originCenters -> Update();

	vtkSmartPointer <vtkFloatArray> normalArray =
			vtkFloatArray::SafeDownCast(normals -> GetOutput() -> GetCellData() -> GetNormals());

	for (vtkIdType i = 0; i < normalArray -> GetNumberOfTuples(); i++) {
		double centerPoint[3];
		double normArray[3];

		originCenters -> GetOutput() -> GetPoint(i, centerPoint);
		normalArray -> GetTuple(i, normArray);

//		cout << "CenterPoint: ";
//		for(int i = 0; i < 3; i++)
//			cout << centerPoint[i] << ", ";
//		cout << endl;

		findIntersectingCells(centerPoint, normArray, target, scalarParam);
	}

	return getThickness();
}


// arg1: call centersRef -> GetOutput() -> GetPoint(i, centerPoint) and pass centerPoint
// arg2: call normalArray->GetTuple(i,normalArray) and pass normalArray
void NormalNeighbour::findIntersectingCells(double centerPoint[3], double normalArray[3], vtkSmartPointer <vtkPolyData> target, double scalarParam) {
	double normalPoint[3];
	double tolerance = 0.0;

	vtkSmartPointer <vtkIdList> ids = vtkSmartPointer <vtkIdList> :: New();
	vtkSmartPointer <vtkCellLocator> locator = vtkSmartPointer <vtkCellLocator> :: New();

	locator -> SetDataSet(target);
	locator -> BuildLocator();

	for(int i = 0; i < 3; i++)
		normalPoint[i] = centerPoint[i] + normalArray[i]*scalarParam;

	locator -> FindCellsAlongLine(centerPoint, normalPoint, tolerance, ids);

//	cout << ids -> GetNumberOfIds() << endl;
	calculateAllThicknesses(ids, centerPoint, target);
}

void NormalNeighbour::calculateAllThicknesses(vtkSmartPointer <vtkIdList> ids, double centerPoint[3], vtkSmartPointer <vtkPolyData> target) {
	double currentPoint [3];
	string thicknessData;

	vtkSmartPointer <vtkCellCenters> targetCenters = vtkSmartPointer <vtkCellCenters> :: New();
	targetCenters -> SetInputData(target);
	targetCenters -> Update();

	// only enters for loop for points at a specific y level --> these are the only points who have any ids in id list
	// all of the points at this level only have 1 id in the id list
	for(vtkIdType i = 0; i < ids -> GetNumberOfIds(); i++) {
		targetCenters -> GetOutput() -> GetPoint(ids -> GetId(i), currentPoint);

		double distance = calculateDistance(currentPoint, centerPoint);

//		cout << "CURRENT POINT: ";
//		for(int i = 0; i < 3; i++)
//			cout << currentPoint[i] << ", ";
//		cout << endl;
//
//		cout << "DISTNACE: " << distance << endl;

		if (distance > getThickness()) {
			setThickness(distance);
			setInnerThickestPoint(centerPoint);
			setOuterThickestPoint(currentPoint);
		}
		thicknessData = to_string(centerPoint[0]) + "," + to_string(centerPoint[1]) + "," + to_string(centerPoint[2])
					+ ",," + to_string(currentPoint[0]) + "," + to_string(currentPoint[1]) + "," + to_string(currentPoint[2])
					+ ",," + to_string(distance) + "\n";
		setThicknessList(getThicknessList() + thicknessData);
	}
}

bool NormalNeighbour::isCellNormals(vtkPolyData* polydata)
{
  cout << "Looking for cell normals..." << endl;

  // Count points
  vtkIdType numCells = polydata->GetNumberOfCells();
  cout << "There are " << numCells << " cells." << endl;

  // Count triangles
  vtkIdType numPolys = polydata->GetNumberOfPolys();
  cout << "There are " << numPolys << " polys." << endl;

  // Double normals in an array
  vtkDoubleArray* normalDataDouble =
    vtkDoubleArray::SafeDownCast(polydata->GetCellData()->GetArray("Normals"));

  if(normalDataDouble)
    {
    int nc = normalDataDouble->GetNumberOfTuples();
    cout << "There are " << nc
            << " components in normalDataDouble" << endl;
    setDoubleArray(normalDataDouble);
    return true;
    }

  // Double normals in an array
  vtkSmartPointer <vtkFloatArray> normalDataFloat =
    vtkFloatArray::SafeDownCast(polydata->GetCellData()->GetArray("Normals"));

  if(normalDataFloat)
    {
    int nc = normalDataFloat->GetNumberOfTuples();
    cout << "There are " << nc
            << " components in normalDataFloat" << endl;

    setFloatArray(normalDataFloat);
    return true;
    }

  // Point normals
  vtkDoubleArray* normalsDouble =
    vtkDoubleArray::SafeDownCast(polydata->GetCellData()->GetNormals());

  if(normalsDouble)
    {
    cout << "There are " << normalsDouble->GetNumberOfComponents()
              << " components in normalsDouble" << endl;
    return true;
    }

  ////////////////////////////////////////////////////////////////
  // Point normals
  vtkFloatArray* normalsFloat =
    vtkFloatArray::SafeDownCast(polydata->GetCellData()->GetNormals());

  if(normalsFloat)
    {
    cout << "There are " << normalsFloat->GetNumberOfComponents()
              << " components in normalsFloat" << endl;
    return true;
    }

  /////////////////////////////////////////////////////////////////////
  // Generic type point normals
  vtkDataArray* normalsGeneric = polydata->GetCellData()->GetNormals(); //works
  if(normalsGeneric)
    {
    cout << "There are " << normalsGeneric->GetNumberOfTuples()
              << " normals in normalsGeneric" << endl;

    double testDouble[3];
    normalsGeneric->GetTuple(0, testDouble);

    cout << "Double: " << testDouble[0] << " "
              << testDouble[1] << " " << testDouble[2] << endl;

    return true;
    }


  // If the function has not yet quit, there were none of these types of normals
  cout << "Normals not found!" << endl;
  return false;

}

void NormalNeighbour::findCurvatureThickest() {
//	vtkSmartPointer <vtkPolyData> temp = vtkSmartPointer <vtkPolyData> :: New();
//	vtkSmartPointer<vtkCurvatures> curvaturesFilter = vtkSmartPointer<vtkCurvatures>::New();
//	curvaturesFilter->SetInputData(getAorta() -> getOuter());
//
//	curvaturesFilter->SetCurvatureTypeToMinimum();
//	curvaturesFilter -> Update();
//
//	curvaturesFilter -> GetOutput() -> GetPointData() -> SetActiveScalars("Normal");
//	temp = curvaturesFilter -> GetOutput();
//	//vtkSmartPointer <vtkDataArray> minCurve = vtkDataArray :: SafeDownCast(temp -> GetPointData() -> GetScalars());
//
//	cout << "~~~~~~MIN CURVE~~~~~~" << endl;
//	curvaturesFilter -> Print(cout);
//
//	for (vtkIdType i = 0; i < temp -> GetPointData() -> GetScalars() -> GetSize(); i++) {
//		cout << temp -> GetPointData() -> GetScalars() -> GetTuple(i)[0] << ", " << endl;
//	}
//
//	curvaturesFilter->SetCurvatureTypeToMaximum();
//	curvaturesFilter -> Update();
//	temp = curvaturesFilter -> GetOutput();
//
//	curvaturesFilter -> GetOutput() -> GetPointData() -> SetActiveScalars("Maximum_Curvature");
//	temp = curvaturesFilter -> GetOutput();
//	vtkSmartPointer <vtkDataArray> maxCurve = curvaturesFilter -> GetOutput() -> GetPointData() -> GetNormals();
//

	vtkSmartPointer <vtkPolyData> temp = vtkSmartPointer <vtkPolyData> :: New();
	temp -> DeepCopy(getAorta() -> getInner());

	const vtkNew<vtkIdList> vertices_n;
	const vtkNew<vtkIdList> vertices;

	double n_f[3]; // normal of facet (could be stored for later?)
	double n_n[3]; // normal of edge
	double t[3];   // to store the cross product of n_f n_n

	double ore[3]; // origin of e
	double end[3]; // end of e
	double oth[3]; //     third vertex necessary for comp of n

	double vn0[3];
	double vn1[3]; // vertices for computation of neighbour's n
	double vn2[3];

	double e[3]; // edge (oriented)

	int n;

	const int F = temp -> GetNumberOfCells();
	cout << "F: " << F << endl;
	for (int f = 0; f < F; f++)
	{
		getAorta() -> getInner() -> GetCellPoints(f, vertices);
		const int nv = vertices->GetNumberOfIds();

		cout << "nv: " << nv << endl;
		for (int v = 0; v < nv; v++)
		{
			double Af = vtkTriangle::TriangleArea(ore, end, oth);

			// find 3 corners of n: in order!
			temp -> GetCellPoints(n, vertices_n);
			temp -> GetPoint(vertices_n->GetId(0), vn0);
			temp -> GetPoint(vertices_n->GetId(1), vn1);
			temp -> GetPoint(vertices_n->GetId(2), vn2);

			Af += double(vtkTriangle::TriangleArea(vn0, vn1, vn2));
			// compute normal of n
			vtkTriangle::ComputeNormal(vn0, vn1, vn2, n_n);
			// the cosine is n_f * n_n
			const double cs = vtkMath::Dot(n_f, n_n);
			// the sin is (n_f x n_n) * e
			vtkMath::Cross(n_f, n_n, t);
			cout << "t: " << t[0] << ", " << t[1] << ", " << t[2] << endl;
			const double sn = vtkMath::Dot(t, e);
			// signed angle in [-pi,pi]
			if (sn != 0.0 || cs != 0.0)
			{
			  const double angle = atan2(sn, cs);
			  //Hf = length * angle;
			}
			else
			{
			 // Hf = 0.0;
			}
	  }
	}
}

//void NormalNeighbour::findHausdorffThickest() {
////	vtkSmartPointer <vtkHausdorffDistancePointSetFilter> distFilter = vtkSmartPointer <vtkHausdorffDistancePointSetFilter> :: New();
//	Hausdorff *distFilter = new Hausdorff();
//	if(originIsInner) {
//		distFilter -> SetInputData(0, getAorta() -> getInner());
//		distFilter -> SetInputData(1, getAorta() -> getOuter());
//	} else {
//		distFilter -> SetInputData(0, getAorta() -> getOuter());
//		distFilter -> SetInputData(1, getAorta() -> getInner());
//	}
//	distFilter -> SetTargetDistanceMethodToPointToCell();
//	distFilter -> Update();
//
//	cout << "Origin to target distance: " << distFilter -> GetRelativeDistance()[0] << endl;
//	cout << "Target to origin distance: " << distFilter -> GetRelativeDistance()[1] << endl;
//
//	if(originIsInner) {
//		getAorta() -> setInner(distFilter -> GetPolyDataOutput());
//		setThickness(distFilter -> GetRelativeDistance()[1]);
//
//	} else {
//		getAorta() -> setOuter(distFilter -> GetPolyDataOutput());
//		setThickness(distFilter -> GetRelativeDistance()[0]);
//	}
//}

void NormalNeighbour::findHausdorffThickest() {
//	vtkSmartPointer <vtkHausdorffDistancePointSetFilter> distFilter = vtkSmartPointer <vtkHausdorffDistancePointSetFilter> :: New();
	Hausdorff *distFilter = new Hausdorff();

	distFilter -> SetInputData(0, getAorta() -> getInner());
	distFilter -> SetInputData(1, getAorta() -> getOuter());

	distFilter -> SetTargetDistanceMethodToPointToCell();
	distFilter -> Update();

	cout << "Inner to outer distance: " << distFilter -> GetRelativeDistance()[0] << endl;
	cout << "Outer to inner distance: " << distFilter -> GetRelativeDistance()[1] << endl;

	getAorta() -> setInner(distFilter -> GetPolyDataOutput());
//	setThickness(distFilter -> GetHausdorffDistance());

//	vtkSmartPointer <vtkPolyData> outer = vtkSmartPointer <vtkPolyData> :: New();
//	outer -> DeepCopy(distFilter -> GetOutput(1));
//	getAorta() -> setOuter(outer);

	setThickness(distFilter -> GetRelativeDistance()[0]);

//	calculateMeanThickness();
//	calculateMedianThickness();
//	calculateSurfaceArea();
}

void NormalNeighbour::calculateMeanThickness() {
	double sum = 0;

	vtkDoubleArray *thicknessArr = vtkDoubleArray :: SafeDownCast(
			getAorta() -> getInner() -> GetPointData() -> GetArray("Distance"));
	cout << thicknessArr -> GetNumberOfTuples() << endl;

	for(int i = 0; i < thicknessArr -> GetNumberOfTuples(); i++) {
		sum += thicknessArr -> GetTuple(i)[0];
		cout << i << ": " << thicknessArr -> GetTuple(i)[0] << endl;
	}

	meanThickness = sum /  thicknessArr -> GetNumberOfTuples();
	cout << "Mean: " << meanThickness << endl;
}

void NormalNeighbour::calculateMedianThickness() {
	int medianIndex = 0;
	vtkSmartPointer <vtkSortDataArray> sorter = vtkSmartPointer <vtkSortDataArray> :: New();
	vtkSmartPointer <vtkPolyData> temp = vtkSmartPointer <vtkPolyData> :: New();

	temp -> DeepCopy(getAorta() -> getInner());
	vtkDoubleArray *thicknessArr = vtkDoubleArray :: SafeDownCast(
			temp -> GetPointData() -> GetArray("Distance"));

	cout << thicknessArr -> GetNumberOfTuples() << endl;
	sorter -> Sort(thicknessArr);

	medianIndex = thicknessArr -> GetNumberOfTuples() / 2 + 1;
	medianThickness = thicknessArr -> GetTuple(medianIndex)[0];
	cout << "Median: " << medianThickness << endl;
}

void NormalNeighbour::subdivide(double numOfPasses) {
	vtkSmartPointer <vtkAdaptiveSubdivisionFilter> subdivideFilter = vtkSmartPointer <vtkAdaptiveSubdivisionFilter> :: New();

	subdivideFilter -> SetInputData(getAorta() -> getInner());
	subdivideFilter -> SetMaximumNumberOfPasses(numOfPasses);
	subdivideFilter -> Update();
	getAorta() -> setInner(subdivideFilter -> GetOutput());

	subdivideFilter -> SetInputData(getAorta() -> getOuter());
	subdivideFilter -> SetMaximumNumberOfPasses(numOfPasses);
	subdivideFilter -> Update();
	getAorta() -> setOuter(subdivideFilter -> GetOutput());
}

void NormalNeighbour::decimate(double reduction){
	vtkSmartPointer <vtkDecimatePro> decimator = vtkSmartPointer <vtkDecimatePro> :: New();

	decimator -> SetInputData(getAorta() -> getInner());
	decimator -> SetTargetReduction(reduction);
	decimator -> PreserveTopologyOn();
	decimator -> Update();
	getAorta() -> setInner(decimator -> GetOutput());

	decimator -> SetInputData(getAorta() -> getOuter());
	decimator -> SetTargetReduction(reduction);
	decimator -> PreserveTopologyOn();
	decimator -> Update();
	getAorta() -> setOuter(decimator -> GetOutput());
}

void NormalNeighbour::calculateSurfaceArea() {
	vtkSmartPointer <vtkMassProperties> massFilter = vtkSmartPointer <vtkMassProperties> :: New();
	massFilter -> SetInputData(getAorta() -> getInner());
	innerSurfaceArea = massFilter -> GetSurfaceArea();
	cout << "Inner mesh's surface area: " << innerSurfaceArea << endl;
}

//bool NormalNeighbour::passesThroughLumen(double point1 [3], double point2 [3]) {
//	vtkSmartPointer <vtkLineSource> line = vtkSmartPointer <vtkLineSource>::New();
//	line -> SetPoint1(point1);
//	line -> SetPoint2(point2);
//	line -> Update();
//
//	vtkSmartPointer <vtkSelectEnclosedPoints> sep =vtkSmartPointer <vtkSelectEnclosedPoints>::New();
//	sep -> CheckSurfaceOff();
//
//	sep -> SetSurfaceData(getAorta() -> getInner());
//
//	vtkSmartPointer<vtkPolyData> pointsPolyData = vtkSmartPointer<vtkPolyData>::New();
//	pointsPolyData->SetPoints(line -> GetPoints());
//	sep -> SetInputData(pointsPolyData);
//
//	sep -> Update();
//
//	vtkSmartPointer <vtkDataArray> insideArray
//		= vtkDataArray::SafeDownCast(sep->GetOutput()->GetPointData()->GetArray("SelectedPoints"));
//
//	cout << "Number of points inside: " << insideArray -> GetSize() << endl;
//	if(insideArray -> GetSize() == 0)
//		return false;
//	return true;
//}
