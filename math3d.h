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


#ifndef MATH3D_H
#define MATH3D_H

#include "math3d.h"
#include "opencv/cv.h"
#include "eigenincludes.h"
#include "vtkPolyData.h"
#include "generaltypes.h"
#include "vtkSmartPointer.h"

namespace tr{

// typedef Eigen::Vector3d Point3d;

// class Point2d
// {
// 
//   double m_x , m_y;
// 
// public:
// 
//   Point2d(){ m_x = 0; m_y = 0; }
// 
//   Point2d( double x , double y  ): m_x( x ) , m_y( y ){}
// 
//   Eigen::Vector2d operator - ( Point2d point ){ return Eigen::Vector2d( m_x - point.m_x , m_y - point.m_y ); }
// 
//   Point2d operator + ( Eigen::Vector2d vec ){ return Point2d( m_x + vec.x() , m_y + vec.y() ); }
// 
//   Point2d operator + ( Point2d point ){ return Point2d( m_x + point.x() , m_y + point.y() ); }
// 
//   double& x(){ return m_x; }
//   double& y(){ return m_y; }
// 
// 
//   void set_x( double x ){ m_x = x; }
//   void set_y( double y ){ m_y = y; }
// 
//   void normalize()
//   {
//     double _norm = sqrt( m_x * m_x + m_y * m_y );
// 
//     m_x /= _norm;
//     m_y /= _norm;
//   }
// 
// };

class Math3D
{

public:
    Math3D();
    
    static cv::Point3f ray_plane_intersection( cv::Vec3f &normal , float origin_distance , cv::Point3f &point , cv::Vec3f &direction );

    static void copy_poly_data( vtkPolyData *src , vtkPolyData *dst );
    
    static double lineIntersection( Point3d &point1, Vector3d &direction1, Point3d &point2, Vector3d &direction2  );

	static void denseCenter( std::vector< tr::Point3d > &points , tr::Point3d &center , double denseFraction = 0.7 );

    static void denseCenter( std::vector< tr::Point3d > &points , tr::Point3d &center , double aabb[ 6 ] , double denseFraction = 0.7  ); 

    static void generateAACube3D( tr::Point3d center  , double length, vtkSmartPointer< vtkPolyData > cube );
    
    virtual ~Math3D();
};


typedef Eigen::Vector3d Point3d;

// class Point3d 
// {
//   
//   double m_x , m_y , m_z;
//   
//   Eigen::Vector3d m_vec;
//   
// public:
//   
//   Point3d(){ //m_x = 0; m_y = 0; m_z = 0; 
//              m_vec( 0 ) = 0;
// 	     m_vec( 1 ) = 0;
// 	     m_vec( 2 ) = 0;
//            }
//   
//   Point3d( double x , double y , double z )
//          { 
// 	   m_vec( 0 ) = x;
// 	   m_vec( 1 ) = y;
// 	   m_vec( 2 ) = z; 
// 	   
// 	 } 
//   
//   Point3d( const Point3d &pt ){ m_vec( 0 ) = pt.x() ; m_vec( 1 ) = pt.y(); m_vec( 2 ) = pt.z(); }
//   
//   Eigen::Vector3d operator - ( Point3d point ){ return ( m_vec - point.vec() ); }
//   
//   double& x(){ return m_vec( 0 ); }
//   double& y(){ return m_vec( 1 ); }
//   double& z(){ return m_vec( 2 ); }
//   
//   void set_x( double x ){ m_vec( 0 ) = x; }
//   void set_y( double y ){ m_vec( 1 ) = y; }
//   void set_z( double z ){ m_vec( 2 ) = z; }
//   
//   Eigen::Vector3d& vec(){ return m_vec; };
// };  
//  

}

#endif // MATH3D_H
