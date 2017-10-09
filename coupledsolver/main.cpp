#include "vtkAutoInit.h" 
VTK_MODULE_INIT(vtkRenderingOpenGL2); // VTK was built with vtkRenderingOpenGL2
VTK_MODULE_INIT(vtkInteractionStyle);


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include "json/json.h"



#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <vtkProperty.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkMath.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLStructuredGridWriter.h>

#include <vtkPolyDataMapper.h>
#include <vtkStructuredGridWriter.h>
#include <vtkStructuredGridGeometryFilter.h>
using namespace std;

class Point
{
public:
	Point(double x, double y, double z);
	Point(array<double, 3> p);
	double& Point::operator[] (int coordi);
	inline double x() const
	{
		return pCoodinates_[0];
	}
	inline double y() const
	{
		return pCoodinates_[1];
	}
	inline double z() const
	{
		return pCoodinates_[2];
	}

private:
	array<double, 3> pCoodinates_ = { 0,0,0 };

};

Point::Point(double x, double y, double z)
	:pCoodinates_(array<double, 3> {x,y,z})
{
}

Point::Point(array<double, 3> p)
	:pCoodinates_(p)
{
}

double& Point::operator[] (int coordi)
{
	if (coordi >= 0 && coordi < 3)
	{
		return pCoodinates_[coordi];
	}
	else
	{
		cout << "index out of the range of point coorinates, should be 0-2\n";
		exit(0);
	}
}

ostream & operator<<(ostream &os, const Point &p)
{
	os << "Point ( " << p.x() << " " << p.y() << " " << p.z() << "\n";
	return os;
}

class Face
{
public:
	Face(array<Point, 4> pts);
	const Point& operator[] (int pindex) const;
private:
	array<Point, 4> points_;
};

Face::Face(array<Point, 4> pts)
	:points_(pts)
{
}

const Point& Face::operator[] (int pindex) const
{
	if (pindex >= 0 && pindex < 4)
	{
		return points_[pindex];
	}
	else
	{
		cout << "index out of the range of face points, should be 0-3\n";
		exit(0);
	}
}

ostream & operator<<(ostream &os, const Face &f)
{
	os << "Face :\n" << f[0] << f[1] << f[2] << f[3] << "\n";
	return os;
}

class Cell
{
public:
	Cell(array<Face, 6> faces);
	const Face& operator[] (int findex) const;
private:
	array<Face, 6> faces_;
};

Cell::Cell(array<Face, 6> faces)
	:faces_(faces)
{
}

const Face& Cell::operator[] (int findex) const
{
	if (findex >= 0 && findex < 6)
	{
		return faces_[findex];
	}
	else
	{
		cout << "index out of the range of cell faces, should be 0-5\n";
		exit(0);
	}
}

ostream & operator<<(ostream &os, const Cell &c)
{
	os << "Cell :\n" << c[0] << c[1] << c[2] << c[3] << c[4] << c[5] << "\n";
	return os;
}

class Mesh
{
public:
	Mesh(string meshfile);
	void render();
	void writeVTK();
	
private:
	int mNcell_;
	int mNface_;
	int mNpoint_;
	int mNinnerface_;

	vector<Cell> mCells_;
	vector<Face> mFaces_;
	vector<Point> mPoints_;

};

Mesh::Mesh(string meshfile)
{
	Json::Value root;
	std::ifstream meshdata(meshfile, std::ifstream::binary);
	meshdata >> root;
	const Json::Value MeshConstants = root["Mesh"];

	mNcell_ = MeshConstants["Ncell"].asInt();
	mNface_ = MeshConstants["Nface"].asInt();
	mNpoint_ = MeshConstants["Npoint"].asInt();
	mNinnerface_ = MeshConstants["Ninnerface"].asInt();

	const Json::Value MeshPoints = root["Points"];
	for (int index = 0; index < MeshPoints.size(); ++index)
	{
		double pCoordx = MeshPoints[index][0].asDouble();
		double pCoordy = MeshPoints[index][1].asDouble();
		double pCoordz = MeshPoints[index][2].asDouble();
		Point pt(pCoordx, pCoordy, pCoordz);
		mPoints_.push_back(pt);
	}

	const Json::Value MeshFaces = root["Faces"];
	for (int index = 0; index < MeshFaces.size(); ++index)
	{
		int fPt0 = MeshFaces[index][0].asInt();
		int fPt1 = MeshFaces[index][1].asInt();
		int fPt2 = MeshFaces[index][2].asInt();
		int fPt3 = MeshFaces[index][3].asInt();
		array<Point, 4> facearray = { mPoints_[fPt0],mPoints_[fPt1],mPoints_[fPt2],mPoints_[fPt3] };
		Face face(facearray);
		mFaces_.push_back(face);
	}

	const Json::Value MeshCells = root["Cells"];
	for (int index = 0; index < MeshCells.size(); ++index)
	{
		int cf0 = MeshCells[index][0].asInt();
		int cf1 = MeshCells[index][1].asInt();
		int cf2 = MeshCells[index][2].asInt();
		int cf3 = MeshCells[index][3].asInt();
		int cf4 = MeshCells[index][4].asInt();
		int cf5 = MeshCells[index][5].asInt();
		array<Face, 6> cellarray = { mFaces_[cf0],mFaces_[cf1],mFaces_[cf2],mFaces_[cf3],mFaces_[cf4],mFaces_[cf5] };
		Cell cell(cellarray);
		mCells_.push_back(cell);
	}


}

void Mesh::render()
{
	// Create a grid
	vtkSmartPointer<vtkStructuredGrid> structuredGrid =
		vtkSmartPointer<vtkStructuredGrid>::New();

	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();

	unsigned int gridSize = mNpoint_;
	unsigned int counter = 0;
	// Create a 5x5 grid of points
	for (unsigned int i = 0; i < gridSize; i++)
	{
		double pCoordx = mPoints_[i][0];
		double pCoordy = mPoints_[i][1];
		double pCoordz = mPoints_[i][2];
		points->InsertNextPoint(pCoordx, pCoordy, pCoordz);
		counter++;
	}

	// Specify the dimensions of the grid
	structuredGrid->SetDimensions(6, 6, 2);

	structuredGrid->SetPoints(points);

	structuredGrid->Modified();

	vtkSmartPointer<vtkStructuredGridWriter> writer =
		vtkSmartPointer<vtkStructuredGridWriter>::New();
	writer->SetFileName("output.vtk");
#if VTK_MAJOR_VERSION <= 5
	writer->SetInput(structuredGrid);
#else
	writer->SetInputData(structuredGrid);
#endif
	writer->Write();

	// Create a mapper and actor
	vtkSmartPointer<vtkDataSetMapper> gridMapper =
		vtkSmartPointer<vtkDataSetMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	gridMapper->SetInputConnection(structuredGrid->GetProducerPort());
#else
	gridMapper->SetInputData(structuredGrid);
#endif

	vtkSmartPointer<vtkActor> gridActor =
		vtkSmartPointer<vtkActor>::New();
	gridActor->SetMapper(gridMapper);
	gridActor->GetProperty()->EdgeVisibilityOn();
	gridActor->GetProperty()->SetEdgeColor(0, 0, 1);

	// Create a renderer, render window, and interactor
	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Add the actor to the scene
	renderer->AddActor(gridActor);
	renderer->SetBackground(.3, .6, .3); // Background color green

										 // Render and interact
	renderWindow->Render();
	renderWindowInteractor->Start();

}

void Mesh::writeVTK()
{

}

class Solver
{
public:
	Solver(Mesh m);

private:
	Mesh solMesh_;
};

Solver::Solver(Mesh m)
	:solMesh_(m)
{}




int main(int argc, char * argv[])
{
	array<double, 3> testparray = { 0.1,0.3,0.5 };
	Point testp11(testparray[0], testparray[1], testparray[2]);
	Point testp1(testparray);
	cout << testp1;
	cout << testp11;
	array<double, 3> testparray2 = { 0.2,0.4,0.6 };
	array<double, 3> testparray3 = { 0.1,0.3,0.6 };
	array<double, 3> testparray4 = { 0.2,0.4,0.5 };
	array<double, 3> testparray5 = { 0.3,0.6,0.9 };
	Point testp2(testparray2);
	Point testp3(testparray3);
	Point testp4(testparray4);
	Point testp5(testparray5);
	array<Point, 4> fpts1 = { testp1,testp2,testp3,testp4 };
	array<Point, 4> fpts2 = { testp2,testp3,testp4,testp5 };
	array<Point, 4> fpts3 = { testp1,testp3,testp4,testp5 };
	array<Point, 4> fpts4 = { testp1,testp2,testp4,testp5 };
	array<Point, 4> fpts5 = { testp1,testp2,testp3,testp5 };
	array<Point, 4> fpts6 = { testp2,testp3,testp4,testp5 };
	Face testf1(fpts1);
	Face testf2(fpts2);
	Face testf3(fpts3);
	Face testf4(fpts4);
	Face testf5(fpts5);
	Face testf6(fpts6);
	cout << testf1;
	array<Face, 6> cfaces = { testf1,testf2,testf3,testf4,testf5,testf6 };
	Cell testc(cfaces);
	cout << testc;

	Mesh meshtest("mesh.json");
	meshtest.render();


	cin.get();

	return 0;
}