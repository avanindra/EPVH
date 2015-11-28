/*
    Copyright (c) 2011, <copyright holder> <email>
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

    THIS SOFTWARE IS PROVIDED BY <copyright holder> <email> ''AS IS'' AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL <copyright holder> <email> BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include "math3d.h"

#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkIdList.h"
#include "limits"
#include "vtkPoints.h"

namespace tr{

Math3D::Math3D()
{

}

cv::Point3f Math3D::ray_plane_intersection(cv::Vec3f &normal, float origin_distance, cv::Point3f &point, cv::Vec3f &direction)
{
  
  float _coeff1 = normal[ 0 ] * direction[ 0 ] + normal[ 1 ] * direction[ 1 ] 
                  + normal[ 2 ] * direction[ 2 ];
  
  float _coeff2 = normal[ 0 ] * point.x  + 
                  normal[ 1 ] * point.y  + 
		  normal[ 2 ] * point.z  + origin_distance;		  
   
  float _depth = - _coeff2 / _coeff1 ;
  
  cv::Point3f intersection_point;
  
  intersection_point.x = point.x + _depth * direction[0];
  intersection_point.y = point.y + _depth * direction[1];
  intersection_point.z = point.z + _depth * direction[2];
  
  return intersection_point;
  
}

void Math3D::copy_poly_data( vtkPolyData *src , vtkPolyData *dst )
{
	dst->Reset();

	dst->Allocate( src->GetNumberOfCells() );

	dst->SetPoints( src->GetPoints() );

	dst->GetPointData()->SetTCoords( src->GetPointData()->GetTCoords() );

	dst->GetPointData()->SetNormals( src->GetPointData()->GetNormals() );

	int num_cells = src->GetNumberOfCells();

	vtkIdList *triangle = vtkIdList::New();

	for( int cc = 0; cc < num_cells; cc++)
	{
		triangle->Reset();

		src->GetCellPoints( cc , triangle );

		dst->InsertNextCell( VTK_TRIANGLE , triangle );
	}
}



double Math3D::lineIntersection( Point3d &point1, Vector3d &direction1, Point3d &point2, Vector3d &direction2 )
{
  double d1 = 0 , d2 = 0;
  
  double denom = ( direction2.x() * direction1.y() - direction1.x() * direction2.y() );
  
  double numer = ( point2( 1 ) - point1( 1 ) ) * direction2( 0 ) + ( point1( 0 )  - point2( 0 ) ) * direction2( 1 ) ;
  
  d1 = ( numer / denom );

  
   //Vector3d vec1 = ( point2 - point1 ).cross( direction2 );
   //Vector3d vec2 = direction1.cross( direction2 ); 
   //
   //d1 = vec1.norm() / vec2.norm();  
  
  return d1;
}

bool pointDIstancePredicate( const std::pair< int , double > &p1 , const std::pair< int , double > &p2 )
{
	return p1.second < p2.second;
}

void Math3D::denseCenter( std::vector< tr::Point3d > &points ,  tr::Point3d &center , double denseFraction )
{
	if( points.size() == 0 )
	{
	  return;
	}

	tr::Point3d average( 0 , 0 , 0 );	

	int numPts = points.size();

	std::vector< double > distances( numPts , 0 );

	std::vector< std::pair< int , double > > sortedPoints( numPts );

	for( int pp = 0; pp < numPts ; pp++ )
	{
		average += points[ pp ];
	}

	average /= numPts;

	int numIters = 1;

	int numThresholdPoint = numPts * denseFraction;

	if( numThresholdPoint < 20 )
	{
		center = average;

		return;
	}

	for( int iter = 0; iter < numIters ; iter++ )
	{

	   double maxD = 0;

	   for( int pp = 0; pp < numPts ; pp++ )
	   {
		  sortedPoints[ pp ].first = pp;
		  sortedPoints[ pp ].second = ( points[ pp ] - average ).norm();
	   }

	   std::sort( sortedPoints.begin() , sortedPoints.end() , pointDIstancePredicate );

	   average = tr::Point3d::Zero();

	   for( int pp = 0; pp < numThresholdPoint ; pp++ )
	   {
		   average += points[ sortedPoints[ pp ].first ];		
	   }

	   average /= numThresholdPoint;

	   center = average;
	}


}

void Math3D::denseCenter( std::vector< tr::Point3d > &points ,  tr::Point3d &center , double aabb[ 6 ] , double denseFraction )
{
	if( points.size() == 0 )
	{
	  return;
	}

	tr::Point3d average( 0 , 0 , 0 );	

	int numPts = points.size();

	std::vector< double > distances( numPts , 0 );

	std::vector< std::pair< int , double > > sortedPoints( numPts );

	for( int pp = 0; pp < numPts ; pp++ )
	{
		average += points[ pp ];
	}

	average /= numPts;

	int numIters = 1;

	int numThresholdPoint = numPts * denseFraction;

	if( numThresholdPoint < 20 )
	{
		center = average;

		return;
	}

	for( int iter = 0; iter < numIters ; iter++ )
	{

	   double maxD = 0;

	   for( int pp = 0; pp < numPts ; pp++ )
	   {
		  sortedPoints[ pp ].first = pp;
		  sortedPoints[ pp ].second = ( points[ pp ] - average ).norm();
	   }

	   std::sort( sortedPoints.begin() , sortedPoints.end() , pointDIstancePredicate );

	   average = tr::Point3d::Zero();

	   for( int pp = 0; pp < numThresholdPoint ; pp++ )
	   {
		   average += points[ sortedPoints[ pp ].first ];		
	   }

	   average /= numThresholdPoint;

	   center = average;
	}

	aabb[ 0 ] = std::numeric_limits< double >::max();
	aabb[ 1 ] = std::numeric_limits< double >::max();
	aabb[ 2 ] = std::numeric_limits< double >::max();

	aabb[ 3 ] = std::numeric_limits< double >::min();
	aabb[ 4 ] = std::numeric_limits< double >::min();
	aabb[ 5 ] = std::numeric_limits< double >::min();

	for( int pp = 0; pp < numThresholdPoint ; pp++ )
	{
	    aabb[ 0 ] = std::min( points[ sortedPoints[ pp ].first ].x() , aabb[ 0 ] );	
		aabb[ 1 ] = std::min( points[ sortedPoints[ pp ].first ].y() , aabb[ 1 ] );
		aabb[ 2 ] = std::min( points[ sortedPoints[ pp ].first ].z() , aabb[ 2 ] );

		aabb[ 3 ] = std::max( points[ sortedPoints[ pp ].first ].x() , aabb[ 3 ] );
		aabb[ 4 ] = std::max( points[ sortedPoints[ pp ].first ].y() , aabb[ 4 ] );
		aabb[ 5 ] = std::max( points[ sortedPoints[ pp ].first ].z() , aabb[ 5 ] );
	}

	
}

void addTriangleInternal( vtkSmartPointer< vtkIdList > triangle , vtkSmartPointer< vtkPolyData > mesh , int id1 , int id2 , int id3 )
{
    triangle->Reset();

    triangle->InsertNextId( id1 );
    triangle->InsertNextId( id2 );
    triangle->InsertNextId( id3 );

    mesh->InsertNextCell( VTK_TRIANGLE , triangle );
}


void Math3D::generateAACube3D( tr::Point3d center , double length ,  vtkSmartPointer< vtkPolyData > cube  )
{
    vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();

    points->Allocate( 8 );

    points->InsertNextPoint( 0 + center.x() , 0 + center.y() , 0 + center.z() );
    points->InsertNextPoint( length + center.x() , 0 + center.y() , 0 + center.z() );
    points->InsertNextPoint( length + center.x() , length + center.y() , 0 + center.z() );
    points->InsertNextPoint( 0 + center.x() , length + center.y() , 0 + center.z() );

    points->InsertNextPoint( 0 + center.x() , 0 + center.y() , length + center.z() );
    points->InsertNextPoint( length + center.x() , 0 + center.y() , length  + center.z());
    points->InsertNextPoint( length + center.x() , length + center.y() , length + center.z() );
    points->InsertNextPoint( 0 + center.x() , length + center.y() , length  + center.z());

    cube->Allocate( 12 );

    vtkSmartPointer< vtkIdList > triangle = vtkSmartPointer< vtkIdList >::New();

    addTriangleInternal( triangle , cube , 0 , 5 , 4 );
    addTriangleInternal( triangle , cube , 0 , 1 , 5 );

    addTriangleInternal( triangle , cube , 1 , 6 , 5 );
    addTriangleInternal( triangle , cube , 1 , 2 , 6 );

    addTriangleInternal( triangle , cube , 3 , 7 , 6 );
    addTriangleInternal( triangle , cube , 3 , 6 , 2 );

    addTriangleInternal( triangle , cube , 0 , 2 , 7 );
    addTriangleInternal( triangle , cube , 0 , 7 , 3 );

    addTriangleInternal( triangle , cube , 0 , 2 , 1 );
    addTriangleInternal( triangle , cube , 0 , 3 , 2 );

    addTriangleInternal( triangle , cube , 1 , 6 , 7 );
    addTriangleInternal( triangle , cube , 1 , 5 , 6 );

    cube->SetPoints( points );
}



Math3D::~Math3D()
{

}

}
