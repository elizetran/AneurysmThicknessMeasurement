#include "Hausdorff.h"

#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"

#include "vtkCellLocator.h"
#include "vtkGenericCell.h"
#include "vtkKdTreePointLocator.h"
#include "vtkPointSet.h"
#include "vtkSmartPointer.h"

#include <vtkAdaptiveSubdivisionFilter.h>
#include <vtkLineSource.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkTriangleFilter.h>

vtkStandardNewMacro(Hausdorff);

Hausdorff::Hausdorff()
{
  this->RelativeDistance[0] = 0.0;
  this->RelativeDistance[1] = 0.0;
  this->HausdorffDistance = 0.0;

  this->SetNumberOfInputPorts(2);
  this->SetNumberOfInputConnections(0, 1);
  this->SetNumberOfInputConnections(1, 1);

  this->SetNumberOfOutputPorts(2);

  this->TargetDistanceMethod = POINT_TO_CELL;
}

Hausdorff::~Hausdorff() {}

int Hausdorff::RequestData(vtkInformation* vtkNotUsed(request),
		  vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  // Get the info objects
  vtkInformation* inInfoA = inputVector[0]->GetInformationObject(0);
  vtkInformation* inInfoB = inputVector[1]->GetInformationObject(0);
  vtkInformation* outInfoA = outputVector->GetInformationObject(0);

  if (inInfoA == nullptr || inInfoB == nullptr)
  {
    return 0;
  }

  // Get the input
  vtkPointSet* inputA = vtkPointSet::SafeDownCast(inInfoA->Get(vtkDataObject::DATA_OBJECT()));
  vtkPointSet* inputB = vtkPointSet::SafeDownCast(inInfoB->Get(vtkDataObject::DATA_OBJECT()));
  vtkPointSet* outputA = vtkPointSet::SafeDownCast(outInfoA->Get(vtkDataObject::DATA_OBJECT()));

  if (inputA->GetNumberOfPoints() == 0 || inputB->GetNumberOfPoints() == 0)
  {
    return 0;
  }

  // Re-initialize the distances
  this->RelativeDistance[0] = 0.0;
  this->RelativeDistance[1] = 0.0;
  this->HausdorffDistance = 0.0;

  vtkSmartPointer<vtkKdTreePointLocator> pointLocatorA =
    vtkSmartPointer<vtkKdTreePointLocator>::New();
  vtkSmartPointer<vtkKdTreePointLocator> pointLocatorB =
    vtkSmartPointer<vtkKdTreePointLocator>::New();

  vtkSmartPointer<vtkCellLocator> cellLocatorA = vtkSmartPointer<vtkCellLocator>::New();
  vtkSmartPointer<vtkCellLocator> cellLocatorB = vtkSmartPointer<vtkCellLocator>::New();

  double dist;
  double currentPoint[3];
  double closestPoint[3];
  double slope[3];
  double offCurrP [3];
  double offCloseP [3];
  double buffer = 1.5;
  int numOfPointsA = inputA -> GetNumberOfPoints();
  bool invalidClosestPoint;
  vtkIdType cellId;
  vtkIdType closestPointId;
  vtkSmartPointer <vtkIdList> ids = vtkSmartPointer <vtkIdList> :: New();
  vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
  vtkSmartPointer <vtkPolyData> tempB = vtkSmartPointer <vtkPolyData> :: New();
  int subId;

  vtkSmartPointer<vtkDoubleArray> distanceAToB = vtkSmartPointer<vtkDoubleArray>::New();
  distanceAToB->SetNumberOfComponents(1);
  distanceAToB->SetNumberOfTuples(inputA->GetNumberOfPoints());
  distanceAToB->SetName("Distance");

//  vtkSmartPointer<vtkDoubleArray> distanceBToA = vtkSmartPointer<vtkDoubleArray>::New();
//  distanceBToA->SetNumberOfComponents(1);
//  distanceBToA->SetNumberOfTuples(inputB->GetNumberOfPoints());
//  distanceBToA->SetName("Distance");

  // Find the nearest neighbors to each point and add edges between them,
  // if they do not already exist and they are not self loops

  for (int i = 0; i < numOfPointsA; i++)
  {
	  //cout << i - 1 << "/" << numOfPointsA << ": "<< distanceAToB -> GetValue(i-1) << endl;
	  tempB -> DeepCopy(inputB);

	  cellLocatorA->SetDataSet(inputA);
	  cellLocatorA->BuildLocator();
	  cellLocatorB->SetDataSet(inputB);
	  cellLocatorB->BuildLocator();

	  inputA->GetPoint(i, currentPoint);

	  do {
		invalidClosestPoint = true;

		cellLocatorB->FindClosestPoint(currentPoint, closestPoint, cell, cellId, subId, dist);


		dist = std::sqrt(
				pow(closestPoint[0] - currentPoint[0], 2) +
				pow(closestPoint[1] - currentPoint[1], 2) +
				pow(closestPoint[2] - currentPoint[2], 2));

		// if distance is long enough for check, do intersection check
		if(dist > 2*buffer) {
			for(int j = 0; j < 3; j++ ) {
				offCurrP[j] = currentPoint[j] + (closestPoint[j] - currentPoint[j])* buffer / dist;
				offCloseP[j] = closestPoint[j] - (closestPoint[j] - currentPoint[j])* buffer / dist;
			}

			cellLocatorA -> FindCellsAlongLine(offCurrP, offCloseP, 0.000, ids);

			// if check passes, set distance value and end loop
			if (ids -> GetNumberOfIds() == 0) {
//				cout << i << ": " << count << ": PASS" << endl;
//				cout << "0. Current Point: " << currentPoint[0] << ", " << currentPoint[1] << ", " << currentPoint[2] << endl;
//				cout << "0. Closest Point: " << closestPoint[0] << ", " << closestPoint[1] << ", " << closestPoint[2] << endl;
//				cout << "0. Slope: " << slope[0] << ", " << slope[1] << ", " << slope[2] << endl;
//				cout << dist << ". Off Current Point: " << offCurrP[0] << ", " << offCurrP[1] << ", " << offCurrP[2] << endl;
//				cout << dist << ". Off Closest Point: " << offCloseP[0] << ", " << offCloseP[1] << ", " << offCloseP[2] << endl;


				invalidClosestPoint = false;
				distanceAToB->SetValue(i, dist);

				if (dist > this->RelativeDistance[0])
				{
				     this->RelativeDistance[0] = dist;
				}
				break;

			// if check fails, delete current B point from search scope and re-build locators
			} else {
				tempB -> DeleteCell(cellId);
				tempB -> RemoveDeletedCells();

				if (tempB -> GetNumberOfCells() < 1 || tempB -> GetNumberOfPoints() < 1) {
					invalidClosestPoint = false;
					distanceAToB->SetValue(i, 0);
					break;
				}

				cellLocatorB->SetDataSet(tempB);
				cellLocatorB->BuildLocator();
			 }

		// if distance too short for check, set distance value and end loop
		} else {
			invalidClosestPoint = false;
			distanceAToB->SetValue(i, dist);

			if (dist > this->RelativeDistance[0])
			{
			    this->RelativeDistance[0] = dist;
			}
		}
	  } while(invalidClosestPoint);
  }


//		if (ids -> GetNumberOfIds() == 0) {
//			cout << "Here" << endl;
//			if (count > 1) {
//				cout << "Count: " << count << endl;
//				cout << "(A) " << i << "'s # of ids: "<< ids -> GetNumberOfIds() << endl;
//				cout << "0. Current Point: " << currentPoint[0] << ", " << currentPoint[1] << ", " << currentPoint[2] << endl;
//				cout << "0. Closest Point: " << closestPoint[0] << ", " << closestPoint[1] << ", " << closestPoint[2] << endl;
//				cout << "0. Slope: " << slope[0] << ", " << slope[1] << ", " << slope[2] << endl;
//				cout << "0. Off Current Point: " << offCurrP[0] << ", " << offCurrP[1] << ", " << offCurrP[2] << endl;
//				cout << "0. Off Closest Point: " << offCloseP[0] << ", " << offCloseP[1] << ", " << offCloseP[2] << endl;
//				cout << endl;
//			}
//
//			invalidClosestPoint = false;
//
//			dist = std::sqrt(std::pow(currentPoint[0] - closestPoint[0], 2) +
//			  std::pow(currentPoint[1] - closestPoint[1], 2) +
//			  std::pow(currentPoint[2] - closestPoint[2], 2));
//
//			if (dist > 40) {
//				cout << dist << ". Current Point: " << currentPoint[0] << ", " << currentPoint[1] << ", " << currentPoint[2] << endl;
//				cout << dist << ". Closest Point: " << closestPoint[0] << ", " << closestPoint[1] << ", " << closestPoint[2] << endl;
//
//				cout << dist << ". Off Current Point: " << offCurrP[0] << ", " << offCurrP[1] << ", " << offCurrP[2] << endl;
//				cout << dist << ". Off Closest Point: " << offCloseP[0] << ", " << offCloseP[1] << ", " << offCloseP[2] << endl;
//
//				//dist = -2;
//			}
//
//			distanceAToB->SetValue(i, dist);
//
//			if (dist > this->RelativeDistance[0])
//			{
//		      this->RelativeDistance[0] = dist;
//			}
//		} else {
//			cout << "(A) " << i << "'s # of ids: "<< ids -> GetNumberOfIds() << endl;
//			cout << "Count: " << count << endl;
//			distanceAToB->SetValue(i, -1);
//			if (this->TargetDistanceMethod == POINT_TO_POINT)
//			{
//				tempB -> DeletePoint(closestPointId);
//				tempB -> Modified();
//
//				pointLocatorB->SetDataSet(tempB);
//				pointLocatorB->BuildLocator();
//			}
//			else {
//				tempB -> DeleteCell(cellId);
//				tempB -> RemoveDeletedCells();
//
//				cellLocatorB->SetDataSet(tempB);
//				cellLocatorB->BuildLocator();
//			}
//			if (tempB -> GetNumberOfCells() == 0 || tempB -> GetNumberOfPoints() == 0) {
//				distanceAToA->SetValue(i, dist);
//				invalidClosestPoint = false;
//			}
//
//		}
//	  } while(invalidClosestPoint);
//  }

//  for (int i = 0; i < inputB->GetNumberOfPoints(); i++)
//  {
//	 vtkSmartPointer <vtkPolyData> tempA = vtkSmartPointer <vtkPolyData> :: New();
//	 tempA -> DeepCopy(inputA);
//
//	 if (this->TargetDistanceMethod == POINT_TO_POINT)
//	 {
//		 pointLocatorA->SetDataSet(inputA);
//	 	 pointLocatorA->BuildLocator();
//	 	 pointLocatorB->SetDataSet(inputB);
//	 	 pointLocatorB->BuildLocator();
//	 }
//	 else
//	 {
//	 	 cellLocatorA->SetDataSet(inputA);
//	 	 cellLocatorA->BuildLocator();
//	 	 cellLocatorB->SetDataSet(inputB);
//	 	 cellLocatorB->BuildLocator();
//	 }
//
//    inputB->GetPoint(i, currentPoint);
//
//    do {
//    	ids -> Reset();
//    	invalidClosestPoint = true;
//		if (this->TargetDistanceMethod == POINT_TO_POINT)
//		{
//		  vtkIdType closestPointId = pointLocatorA->FindClosestPoint(currentPoint);
//		  inputA->GetPoint(closestPointId, closestPoint);
//
//		}
//		else
//		{
//		  cellLocatorA->FindClosestPoint(currentPoint, closestPoint, cell, cellId, subId, dist);
//		}
//
//		for(int j = 0; j < 3; j++ ) {
//			slope[j] = closestPoint[j] - currentPoint[j];
//			offCurrP[j] = currentPoint[j] + 0.1 * slope[j];
//			offCloseP[j] = currentPoint[j] + (1 - 0.1) * slope[j];
//		}
//
//		cellLocatorB -> FindCellsAlongLine(offCurrP, offCloseP, 0.0, ids);
//		if (ids -> GetNumberOfIds() == 0) {
//
//			invalidClosestPoint = false;
//
//			dist = std::sqrt(std::pow(currentPoint[0] - closestPoint[0], 2) +
//			  std::pow(currentPoint[1] - closestPoint[1], 2) +
//			  std::pow(currentPoint[2] - closestPoint[2], 2));
//			distanceBToA->SetValue(i, dist);
//
//			if (dist > this->RelativeDistance[1])
//			{
//				if(dist < 1) {
//					cout << "1. Current Point: " << currentPoint[0] << ", " << currentPoint[1] << ", " << currentPoint[2] << endl;
//					cout << "1. Closest Point: " << closestPoint[0] << ", " << closestPoint[1] << ", " << closestPoint[2] << endl;
//
////					cout << "1. Slope: " << slope[0] << ", " << slope[1] << ", " << slope[2] << endl;
////					cout << "1. Off Current Point: " << offCurrP[0] << ", " << offCurrP[1] << ", " << offCurrP[2] << endl;
////					cout << "1. Off Closest Point: " << offCloseP[0] << ", " << offCloseP[1] << ", " << offCloseP[2] << endl;
//				}
//
//		        this->RelativeDistance[1] = dist;
//			}
//		}
//	    else
//	    {
////	    	cout << "(B) " << i << "'s # of ids: "<< ids -> GetNumberOfIds() << endl;
////	    	distanceBToA->SetValue(i, -1);
//			if (this->TargetDistanceMethod == POINT_TO_POINT)
//			{
//				tempA -> DeletePoint(closestPointId);
//				tempA -> Modified();
//
//				pointLocatorA->SetDataSet(tempA);
//				pointLocatorA->BuildLocator();
//			}
//			else {
//				tempA -> DeleteCell(cellId);
//				tempA -> RemoveDeletedCells();
//				tempA -> Modified();
//
//				cellLocatorA->SetDataSet(tempA);
//				cellLocatorA->BuildLocator();
//			}
//
//			if (tempA -> GetNumberOfCells() == 0 || tempA -> GetNumberOfPoints() == 0) {
//				distanceBToA->SetValue(i, dist);
//				invalidClosestPoint = false;
//			}
//	    }
//    } while (invalidClosestPoint);
//  }

  if (this->RelativeDistance[0] >= RelativeDistance[1])
  {
    this->HausdorffDistance = this->RelativeDistance[0];
  }
  else
  {
    this->HausdorffDistance = this->RelativeDistance[1];
  }

  vtkSmartPointer<vtkDoubleArray> relativeDistanceAtoB = vtkSmartPointer<vtkDoubleArray>::New();
  relativeDistanceAtoB->SetNumberOfComponents(1);
  relativeDistanceAtoB->SetName("RelativeDistanceAtoB");
  relativeDistanceAtoB->InsertNextValue(RelativeDistance[0]);

  vtkSmartPointer<vtkDoubleArray> hausdorffDistanceFieldDataA =
    vtkSmartPointer<vtkDoubleArray>::New();
  hausdorffDistanceFieldDataA->SetNumberOfComponents(1);
  hausdorffDistanceFieldDataA->SetName("HausdorffDistance");
  hausdorffDistanceFieldDataA->InsertNextValue(HausdorffDistance);

  vtkSmartPointer<vtkDoubleArray> hausdorffDistanceFieldDataB =
    vtkSmartPointer<vtkDoubleArray>::New();
  hausdorffDistanceFieldDataB->SetNumberOfComponents(1);
  hausdorffDistanceFieldDataB->SetName("HausdorffDistance");
  hausdorffDistanceFieldDataB->InsertNextValue(HausdorffDistance);

  outputA->DeepCopy(inputA);
  outputA->GetPointData()->AddArray(distanceAToB);
  outputA->GetFieldData()->AddArray(relativeDistanceAtoB);
  outputA->GetFieldData()->AddArray(hausdorffDistanceFieldDataA);

  return 1;
}



int Hausdorff::FillInputPortInformation(int port, vtkInformation* info)
{
  // The input should be two vtkPointsSets
  if (port == 0)
  {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  if (port == 1)
  {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  return 0;
}
