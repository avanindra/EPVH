/*
    Copyright (c) 2012, avanindra <email>
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
        * Neither the name of the <organization> nor the
        names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY avanindra <email> ''AS IS'' AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL avanindra <email> BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include "display3droutines.h"
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkRenderer.h>
#include <vtkPolyDataMapper.h>
#include "vtkProperty.h"
#include <vtkIdList.h>
#include <vtkTriangleFilter.h>
#include "vtkVersion.h"

namespace tr
{

Display3DRoutines::Display3DRoutines()
{

}


void Display3DRoutines::displayPolyData( vtkSmartPointer< vtkPolyData > mesh )
{
  vtkSmartPointer< vtkRenderer > renderer = vtkSmartPointer< vtkRenderer >::New();
  
  vtkSmartPointer< vtkPolyDataMapper > meshMapper = vtkSmartPointer< vtkPolyDataMapper >::New();
 
#if VTK_MAJOR_VERSION < 6  
  meshMapper->SetInput( mesh );
#else 
  meshMapper->SetInputData( mesh );
#endif  
  
  vtkSmartPointer< vtkActor > meshActor = vtkSmartPointer< vtkActor >::New();
  
  meshActor->SetMapper( meshMapper );
  
  renderer->AddActor( meshActor );/*

  meshActor->GetProperty()->SetAmbient(1.0);
  meshActor->GetProperty()->SetDiffuse(0.0);
  meshActor->GetProperty()->SetSpecular(0.0);*/
  
  vtkSmartPointer< vtkRenderWindow > window = vtkSmartPointer< vtkRenderWindow >::New();
  vtkSmartPointer< vtkRenderWindowInteractor > interactor = vtkSmartPointer< vtkRenderWindowInteractor >::New();
  vtkSmartPointer< vtkInteractorStyleTrackballCamera > style = vtkSmartPointer< vtkInteractorStyleTrackballCamera >::New();  
  window->AddRenderer( renderer );
  window->SetInteractor( interactor );
  interactor->SetInteractorStyle( style );
  
  window->Render();
  
  interactor->Start();
}


void Display3DRoutines::displayPolygons(vtkSmartPointer< vtkPolyData > mesh)
{

}


//void Display3DRoutines::displayCarvePolyhedron( CarvePolyhedron& polyhedron )
//{
/*
  carve::mesh::MeshSet<3> *carveMesh = polyhedron.getPolyhedron();
  
   
  int numVertices = carveMesh->vertex_storage.size();
  
  vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();
  
  points->Allocate( numVertices );
  
  for( int vv = 0; vv < numVertices; vv++ )
  {
    double p[ 3 ];
    
    p[ 0 ] = carveMesh->vertex_storage[ vv ].v[ 0 ];
    p[ 1 ] = carveMesh->vertex_storage[ vv ].v[ 1 ];
    p[ 2 ] = carveMesh->vertex_storage[ vv ].v[ 2 ];
    
    points->InsertNextPoint( p );
  }
  
  int numMeshes = carveMesh->meshes.size();
  
  int numFaces = 0;
  
  for( int mm = 0; mm < numMeshes; mm++ )
  {
    numFaces += carveMesh->meshes[ mm ]->faces.size();
  }
  
  vtkSmartPointer< vtkPolyData > polydataMesh = vtkSmartPointer< vtkPolyData >::New();
  
  polydataMesh->Allocate( numFaces );
  
    std::map<carve::mesh::MeshSet<3>::vertex_t*, unsigned int > vertexToIndex_map;
    std::vector<carve::mesh::MeshSet<3>::vertex_t>::iterator it = carveMesh->vertex_storage.begin();
    for (int i = 0; it != carveMesh->vertex_storage.end(); ++i, ++it) {
        carve::mesh::MeshSet<3>::vertex_t *vertex = &(*it);
        vertexToIndex_map[vertex] = i;
    }
    
    
  vtkSmartPointer< vtkIdList > poly = vtkSmartPointer< vtkIdList >::New();  
  
  for( carve::mesh::MeshSet<3>::face_iter iter = carveMesh->faceBegin() ; iter != carveMesh->faceEnd(); iter++ )
  {
    carve::mesh::MeshSet<3>::face_t *face = *iter;
    
    int numFacetPoints = face->nVertices();
    
    std::vector< carve::mesh::Vertex<3>* > verts;
    
    face->getVertices(verts);
    
    poly->Reset();
    
    for( int vv = 0; vv < numFacetPoints; vv++ )
    {
      poly->InsertNextId( vertexToIndex_map[ verts[ vv ] ] );
    }
    
    polydataMesh->InsertNextCell( VTK_POLYGON , poly );

  }
  
  polydataMesh->SetPoints( points );
  
      vtkSmartPointer< vtkTriangleFilter > triangulator = vtkSmartPointer< vtkTriangleFilter >::New();

#if VTK_MAJOR_VERSION < 6      
    triangulator->SetInput( polydataMesh );
#else
    triangulator->SetInputData( polydataMesh );
#endif    
    
    triangulator->Update();
    polydataMesh = triangulator->GetOutput();
  
  displayPolyData( polydataMesh ); */
  
//} 


void Display3DRoutines::displayPointSet(std::vector< Eigen::Vector3d >& vertices, std::vector< Eigen::Vector3d >& colors)
{
	vtkSmartPointer< vtkPolyData > dataSet = vtkSmartPointer< vtkPolyData >::New();
	vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();
	vtkSmartPointer< vtkUnsignedCharArray > vColors = vtkSmartPointer< vtkUnsignedCharArray >::New();

	int numVerts = vertices.size();

	points->Allocate(numVerts);
	dataSet->Allocate(numVerts);

	vColors->SetNumberOfComponents(3);
	vColors->SetName("Colors");

	// Add the three colors we have created to the array
	vtkSmartPointer< vtkIdList > vert = vtkSmartPointer< vtkIdList >::New();

	int id = 0;

	for (int vv = 0; vv < numVerts; vv++)
	{
		unsigned char col[] = {(unsigned char) (colors[vv](0) * 255),
        				(unsigned char) (colors[vv](1) * 255),
        				(unsigned char) (colors[vv](2) * 255) };



		vColors->InsertNextTypedTuple(col);

		points->InsertNextPoint(vertices[vv](0), vertices[vv](1), vertices[vv](2));

		vert->Reset();

		vert->InsertNextId(vv);

		dataSet->InsertNextCell(VTK_VERTEX, vert);

	}

	dataSet->SetPoints(points);
	dataSet->GetPointData()->SetScalars(vColors);

	displayPolyData(dataSet);
}

Display3DRoutines::~Display3DRoutines()
{

}

}
