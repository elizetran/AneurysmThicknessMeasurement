
#include "SimulatedShape.h"
#include "Frontend.h"
#include "NearestNeighbour.h"
#include "NormalNeighbour.h"

int main (void) {
	SimulatedShape *simShape = new SimulatedShape();
	Frontend *fe = new Frontend();

	simShape = fe -> requestData();
	simShape -> trimLumen();
//	fe -> menuTestGeometry(simShape);


///* NEAREST NEIGHBOUR TESTING */
//	NearestNeighbour *nearN = new NearestNeighbour();
//	nearN -> setAorta(simShape);
//	nearN -> findThickest();
//
//	for (int i = 0; i < 3; i++)
//		cout << nearN -> getInnerThickestPoint()[i] << " ";
//	cout << " with thickness of: " << nearN -> getThickness() << endl;
//
//	// writing to file does not run from terminal, only app -> contents -> macos -> geometrymodule.exe
//	fe -> writeThicknessOutput(nearN -> getThicknessList());
//	fe -> highlightLine(nearN -> getOuterThickestPoint(), nearN -> getInnerThickestPoint());
//	fe -> addSimulatedShape(simShape);
//	fe -> writeThicknessOutput(nearN -> getThicknessList());
//
//	cout << "Thickest point: " ;
//	for (int i = 0; i < 3; i++)
//		cout << nearN -> getInnerThickestPoint()[i] << " ";


/* NORMAL NEIGHBOUR TESTING	*/
	string pointData;
	NormalNeighbour *normN = new NormalNeighbour();
	normN -> setAorta(simShape);


	/* in -> out: scalarParam = 4
	   out -> in: scalarParam = -4 */
	//normN -> findThickest(simShape -> getOuter(), simShape -> getInner(), -2);
//	normN -> findCurvatureThickest();

	/* POINT DENSITY TESTING */
//  normN -> subdivide(1);
//	double surfaceArea = 0;
//
//	for (int i = 0; i < 15; i++) {
//		cout << "Run " << i << endl;
//		normN -> findHausdorffThickest();
//
//		if (i == 0) {
//			fe -> addNormalThicknessDisplay(normN);
//			fe -> saveResults(normN);
//			surfaceArea = normN -> getInnerSurfaceArea();
//		}
//
//		pointData += to_string(normN -> getAorta() -> getInner() -> GetNumberOfPoints())
//				+ "," + to_string(normN -> getAorta() -> getInner() -> GetNumberOfPolys())
//				+ "," + to_string(normN -> getAorta() -> getOuter() -> GetNumberOfPoints())
//				+ "," + to_string(normN -> getAorta() -> getOuter() -> GetNumberOfPolys())
//				+ "," + to_string(normN -> getMeanThickness())
//				+ "," + to_string(normN -> getMedianThickness())
//				+ "," + to_string(normN -> getThickness())
//				+ "," + to_string(normN -> getInnerSurfaceArea()) + "\n";
//		normN -> decimate(0.1);
//	}
//
//	pointData += "Surface Area," + to_string(surfaceArea);
//
//	fe -> writeMeshPointData(pointData);

	normN -> findHausdorffThickest();
	fe -> saveResults(normN);
	fe -> addNormalThicknessDisplay(normN);

/* DISPLAY */
	cout << "Thickness: " << normN -> getThickness() << endl;
	fe -> display();

	return EXIT_SUCCESS;
}
