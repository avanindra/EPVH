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


#ifndef MATH2D_H
#define MATH2D_H

#include "opencv/cv.h"
#include "generaltypes.h"



namespace tr{  
  
  
inline  float distance( cv::Point2f point1, cv::Point2f point2 )
{
  
  float dist = sqrt( ( point1.x - point2.x ) * ( point1.x - point2.x ) + 
                ( point1.y - point2.y ) * ( point1.y - point2.y ) );
		    
   return dist;		    
	    
}  


inline  float distance( Point2 &point1, Point2 &point2 )
{
  
  float dist = ( point1 - point2 ).norm();
		    
   return dist;		    
	    
}  


//inline  float distance( cv::Point2f point1, cv::Point2f point2 )
//{

//  float dist = sqrt( ( point1.x - point2.x ) * ( point1.x - point2.x ) +
//                ( point1.y - point2.y ) * ( point1.y - point2.y ) );

//   return dist;

//}

inline  float distance( cv::Point2f& point1, Point2& point2 )
{
  
  float dist = sqrt( ( point1.x - point2.x() ) * ( point1.x - point2.x() ) + 
                ( point1.y - point2.y() ) * ( point1.y - point2.y() ) );
		    
   return dist;		    
	    
} 

void lineEquation( cv::Point2f& point1, cv::Point2f& point2 , float &m , float &c );
void lineEquation( cv::Point& point1, cv::Point& point2 , float &m , float &c );

double triangleArea( Point2 &pt1 , Point2 &pt2 , Point2 &pt3 );

double triangleArea( Eigen::Vector2f &pt1 , Eigen::Vector2f &pt2 , Eigen::Vector2f &pt3 );

double triangleArea( cv::Point2f pt1 , cv::Point2f pt2 , cv::Point2f pt3 );


void closestPointOnLine( tr::Point2& point , tr::Point2 &origin , tr::Vector2 &direction , tr::Point2 &closestPoint 
);
float lineToLineIntersection2D( tr::Point2 &point1, tr::Point2 &point2, tr::Point2 &point3, tr::Point2 &point4, tr::Point2& intersection );

float lineToLineIntersection2D( cv::Point2f point1, cv::Point2f point2, cv::Point2f point3, cv::Point2f point4, cv::Point2f& intersection );

float lineToLineIntersection2D( tr::Point2f point1, tr::Point2f point2, tr::Point2f point3, tr::Point2f point4, std::pair< float , float > &coefficients );

double lineToLineIntersection2D(Eigen::Vector2d point1, Eigen::Vector2d point2, Eigen::Vector2d point3, Eigen::Vector2d point4, std::pair< double, double > &coefficients);

void lineToLineIntersection2D( Point2 &point1, Point2 &point2, Point2 &point3, Point2 &point4, Point2& intersection , std::pair< float , float > &coefficients );

void lineToLineIntersection2D( cv::Point2f point1, cv::Point2f point2, cv::Point2f point3, cv::Point2f point4, cv::Point2f& intersection , std::pair< float , float > &coefficients );

void lineToLineIntersection2D(cv::Point2f point1, cv::Point2f point2, cv::Point2f point3, cv::Point2f point4, cv::Point2f& intersection, std::pair< double , double > &coefficients);

void lineToLineIntersection2D( Eigen::Vector2d point1, Eigen::Vector2d point2, Eigen::Vector2d point3, 
	                           Eigen::Vector2d point4, Eigen::Vector2d& intersection, std::pair< double, double > &coefficients);


class Math2D
{

public:
    Math2D();
    
    static double distance( double *p1 , double *p2 );
    
    static float distance(cv::Point2f point1, cv::Point2f point2);
    
    static double distance( cv::Point2d point1 , cv::Point2d point2 );
    
    static void create_cross_product_matrix( cv::Mat &vec , cv::Mat &product_matrix );
    
    static void unitize( cv::Vec3f &vec );
    
    static float modulus( cv::Vec3f &vec );
    
    static double modulus( cv::Vec3d &vec );
    
    static cv::Point2f parametric_point_wrt_triangle( cv::Point &triangle_point1 , cv::Point &triangle_point2 ,
                                         cv::Point &triangle_point3 , cv::Point2f &input_point 
                                       );
    static cv::Vec4f planeEquation(cv::Point3f &pt1, cv::Point3f &pt2, cv::Point3f &pt3);
    static cv::Point2f weightWrtThreePoints( cv::Point3f &pt1, cv::Point3f &pt2, cv::Point3f &pt3,
                                             cv::Point3f &plane_pt);
    virtual ~Math2D();
};

}
#endif // MATH2D_H
