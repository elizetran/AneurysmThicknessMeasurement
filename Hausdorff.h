#ifndef Hausdorff_h_
#define Hausdorff_h_

#include <vtkFiltersModelingModule.h>
#include <vtkPointSetAlgorithm.h>
#include <vtkSmartPointer.h>
#include <vtkPointSet.h>

//#include "vtkPythonProgrammableFilter.h"

class VTKFILTERSMODELING_EXPORT Hausdorff : public vtkPointSetAlgorithm /*, public vtkPythonProgrammableFilter*/ {
public:
	Hausdorff();
	~Hausdorff();

	static Hausdorff* New();
	vtkTypeMacro(Hausdorff, vtkPointSetAlgorithm);
	//void PrintSelf(ostream& os, vtkIndent indent) override;

	enum DistanceMethod
	{
	  POINT_TO_POINT,
	  POINT_TO_CELL
	};

	vtkSetMacro(TargetDistanceMethod, int);
	vtkGetMacro(TargetDistanceMethod, int);
	void SetTargetDistanceMethodToPointToPoint() { this->SetTargetDistanceMethod(POINT_TO_POINT); }
	void SetTargetDistanceMethodToPointToCell() { this->SetTargetDistanceMethod(POINT_TO_CELL); }
	const char* GetTargetDistanceMethodAsString();

	double *GetRelativeDistance() { return RelativeDistance; }
	double GetHausdorffDistance() { return HausdorffDistance; }

protected:
	int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
	int FillInputPortInformation(int port, vtkInformation* info) override;

	int TargetDistanceMethod;   //!< point-to-point if 0, point-to-cell if 1
	double RelativeDistance[2]; //!< relative distance between inputs
	double HausdorffDistance;   //!< hausdorff distance (max(relative distance))

private:
	Hausdorff(const Hausdorff&) = delete;
	void operator=(const Hausdorff&) = delete;

};

inline const char* Hausdorff::GetTargetDistanceMethodAsString()
{
  if (this->TargetDistanceMethod == POINT_TO_POINT)
  {
    return "PointToPoint";
  }
  else
  {
    return "PointToCell";
  }
}

#endif
