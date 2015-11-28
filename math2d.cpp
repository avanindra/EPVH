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


#include "math2d.h"
#include "math.h"
#include "Eigen/Dense"

namespace tr{  
  
  
double triangleArea(Point2 &pt1, Point2 &pt2, Point2 &pt3)
{
  double c = ( pt1 - pt2 ).norm();
  double a = ( pt2 - pt3 ).norm();
  double b = ( pt3 - pt1 ).norm();
  
  double s = ( a + b + c ) / 2;
  
  return sqrt( s * ( s - a ) * ( s - b ) * ( s - c ) );
}

double triangleArea(Eigen::Vector2f& pt1, Eigen::Vector2f& pt2, Eigen::Vector2f& pt3)
{
  double c = ( pt1 - pt2 ).norm();
  double a = ( pt2 - pt3 ).norm();
  double b = ( pt3 - pt1 ).norm();
  
  double s = ( a + b + c ) / 2;
  
  return sqrt( s * ( s - a ) * ( s - b ) * ( s - c ) );
}


double triangleArea( cv::Point2f pt1 , cv::Point2f pt2 , cv::Point2f pt3 )
{
	double c = distance( pt1 , pt2 );//( pt1 - pt2 ).norm();
    double a = distance( pt2 , pt3 );//( pt2 - pt3 ).norm();
    double b = distance( pt3 , pt1 );//( pt3 - pt1 ).norm();
  
  double s = ( a + b + c ) / 2;
  
  return sqrt( s * ( s - a ) * ( s - b ) * ( s - c ) );
}


void lineEquation(cv::Point2f &point1, cv::Point2f &point2, float &m, float &c)
{
    m = ( point2.y - point1.y ) / ( point2.x - point1.x ) ;

    c = point1.y - m * point1.x;
}


void lineEquation( cv::Point &point1, cv::Point &point2, float &m, float &c )
{
    m = ( point2.y - point1.y ) / ( point2.x - point1.x ) ;

    c = point1.y - m * point1.x;
}
  

void closestPointOnLine( Point2 &point , Point2 &origin , Vector2 &direction , Point2 &closestPoint )
{
  closestPoint = origin + direction * ( point - origin ).dot( direction );
}  


float lineToLineIntersection2D( Point2 &point1 , Point2 &point2 , Point2 &point3 , Point2 &point4 , Point2 &intersection)
{
  float delX12 = point1.x() - point2.x();
  float delX34 = point3.x() - point4.x();
  float delX42 = point4.x() - point2.x();
  float delY12 = point1.y() - point2.y();
  float delY34 = point3.y() - point4.y();
  float delY42 = point4.y() - point2.y();
  
  float v = ( delX12 * delY42 - delX42 * delY12 ) / ( delX34 * delY12 - delY34 * delX12 ) ;
  
  intersection = v * point3 + ( 1 - v ) * point4;
  
  return v;
}


float lineToLineIntersection2D( cv::Point2f point1 , cv::Point2f point2 , cv::Point2f point3 , cv::Point2f point4 , cv::Point2f& intersection )
{
  float delX12 = point1.x - point2.x;
  float delX34 = point3.x - point4.x;
  float delX42 = point4.x - point2.x;
  float delY12 = point1.y - point2.y;
  float delY34 = point3.y - point4.y;
  float delY42 = point4.y - point2.y;
  
  float v = ( delX12 * delY42 - delX42 * delY12 ) / ( delX34 * delY12 - delY34 * delX12 ) ;
  
  intersection.x = v * point3.x + ( 1 - v ) * point4.x;
  intersection.y = v * point3.y + ( 1 - v ) * point4.y;
  
  return v;
}


float lineToLineIntersection2D( tr::Point2f point1, tr::Point2f point2, tr::Point2f point3, tr::Point2f point4, std::pair< float, float >& coefficients)
{
  float delX12 = point1.x() - point2.x();
  float delX34 = point3.x() - point4.x();
  float delX42 = point4.x() - point2.x();
  float delY12 = point1.y() - point2.y();
  float delY34 = point3.y() - point4.y();
  float delY42 = point4.y() - point2.y(); 
  
  coefficients.first = ( delX42 * delY34 - delX34 * delY42 ) / ( delX12 * delY34 - delY12 * delX34 );
  coefficients.second = ( delX12 * delY42 - delX42 * delY12 ) / ( delX34 * delY12 - delY34 * delX12 ) ;

  return coefficients.first;
}

double lineToLineIntersection2D(Eigen::Vector2d point1, Eigen::Vector2d point2, Eigen::Vector2d point3, Eigen::Vector2d point4, std::pair< double, double > &coefficients)
{
	float delX12 = point1.x() - point2.x();
	float delX34 = point3.x() - point4.x();
	float delX42 = point4.x() - point2.x();
	float delY12 = point1.y() - point2.y();
	float delY34 = point3.y() - point4.y();
	float delY42 = point4.y() - point2.y();

	coefficients.first = (delX42 * delY34 - delX34 * delY42) / (delX12 * delY34 - delY12 * delX34);
	coefficients.second = (delX12 * delY42 - delX42 * delY12) / (delX34 * delY12 - delY34 * delX12);

	return coefficients.first;
}

void lineToLineIntersection2D( Point2 &point1 , Point2 &point2 , Point2 &point3 , Point2 &point4 , Point2& intersection , std::pair< float , float >& coefficients)
{
  float delX12 = point1.x() - point2.x();
  float delX34 = point3.x() - point4.x();
  float delX42 = point4.x() - point2.x();
  float delY12 = point1.y() - point2.y();
  float delY34 = point3.y() - point4.y();
  float delY42 = point4.y() - point2.y(); 
  
  coefficients.first = ( delX42 * delY34 - delX34 * delY42 ) / ( delX12 * delY34 - delY12 * delX34 );
  coefficients.second = ( delX12 * delY42 - delX42 * delY12 ) / ( delX34 * delY12 - delY34 * delX12 ) ;
  
  intersection = coefficients.second * point3 + ( 1 - coefficients.second )  * point4;
  
}


void lineToLineIntersection2D( cv::Point2f point1 , cv::Point2f point2, cv::Point2f point3, cv::Point2f point4, cv::Point2f& intersection, std::pair< float, float >& coefficients)
{
  float delX12 = point1.x - point2.x;
  float delX34 = point3.x - point4.x;
  float delX42 = point4.x - point2.x;
  float delY12 = point1.y - point2.y;
  float delY34 = point3.y - point4.y;
  float delY42 = point4.y - point2.y; 
  
  coefficients.first = ( delX42 * delY34 - delX34 * delY42 ) / ( delX12 * delY34 - delY12 * delX34 );
  coefficients.second = ( delX12 * delY42 - delX42 * delY12 ) / ( delX34 * delY12 - delY34 * delX12 ) ;
  
  intersection.x = coefficients.second * point3.x + ( 1 - coefficients.second ) * point4.x;
  intersection.y = coefficients.second * point3.y + ( 1 - coefficients.second ) * point4.y;
  
}




Math2D::Math2D()
{

}

double Math2D::distance(double* p1, double* p2)
{
  double _distance = ( p1[ 0 ] - p2[ 0 ] ) * ( p1[ 0 ] - p2[ 0 ] ) + 
                     ( p1[ 1 ] - p2[ 1 ] ) * ( p1[ 1 ] - p2[ 1 ] );
  
  _distance = sqrt( _distance );
  
  return _distance;
}




float Math2D::distance(cv::Point2f point1, cv::Point2f point2)
{
  float _distance = ( point1.x - point2.x ) * ( point1.x - point2.x ) + 
                    ( point1.y - point2.y ) * ( point1.y - point2.y );
		    
  _distance = sqrt( _distance );		    
		    
  return _distance;		    
}

double Math2D::distance(cv::Point2d point1, cv::Point2d point2)
{
  double _distance = ( point1.x - point2.x ) * ( point1.x - point2.x ) + 
                     ( point1.y - point2.y ) * ( point1.y - point2.y );
		    
  _distance = sqrt( _distance );		    
		    
  return _distance;
}



void Math2D::create_cross_product_matrix(cv::Mat& vec, cv::Mat& product_matrix)
{  	
  assert( vec.rows * vec.cols == 3 );
  
  product_matrix.create( 3 , 3 , CV_64FC1 ); 
  
  double *v = ( double * )vec.data;
  
  double *p = ( double * )product_matrix.data;

  const double a1 = v[ 0 ], a2 = v[ 1 ], a3 = v[ 2 ];

  p[ 0 ] = 0; 
  p[ 1 ] = -a3; 
  p[ 2 ] = a2; 

  p[ 3 ] = a3; 
  p[ 4 ] = 0; 
  p[ 5 ] = -a1; 

  p[ 6 ] = -a2; 
  p[ 7 ] = a1; 
  p[ 8 ] = 0; 
  
}

void Math2D::unitize(cv::Vec3f& vec)
{
  float _mod = vec[ 0 ] * vec[ 0 ] + vec[ 1 ] * vec[ 1 ] + vec[ 2 ] * vec[ 2 ] ; 
  
  _mod = sqrt( _mod );
  
  float _inv_mod = 1.0 / _mod;
  
  vec[ 0 ] *= _inv_mod;
  vec[ 1 ] *= _inv_mod;
  vec[ 2 ] *= _inv_mod;
}


float Math2D::modulus(cv::Vec3f& vec)
{
   float _mod = vec[ 0 ] * vec[ 0 ] + vec[ 1 ] * vec[ 1 ] + vec[ 2 ] * vec[ 2 ] ; 
  
   _mod = sqrt( _mod );
   
   return _mod;
} 

double Math2D::modulus(cv::Vec3d& vec)
{
   double _mod = vec[ 0 ] * vec[ 0 ] + vec[ 1 ] * vec[ 1 ] + vec[ 2 ] * vec[ 2 ] ; 
  
   _mod = sqrt( _mod );
   
   return _mod;
}


cv::Vec4f Math2D::planeEquation(cv::Point3f &pt1, cv::Point3f &pt2, cv::Point3f &pt3)
{
    Eigen::Matrix<float,3,3> matA;
    Eigen::Matrix<float,3,3> matB;
    Eigen::Matrix<float,3,3> matC;
    Eigen::Matrix<float,3,3> matD;

    matA<<1,pt1.y,pt1.z,1,pt2.y,pt2.z,1,pt3.y,pt3.z;
    matB<<pt1.x,1,pt1.z,pt2.x,1,pt2.z,pt3.x,1,pt3.z;
    matC<<pt1.x,pt1.y,1,pt2.x,pt2.y,1,pt3.x,pt3.y,1;
    matD<<pt1.x,pt1.y,pt1.z,pt2.x,pt2.y,pt2.z,pt3.x,pt3.y,pt3.z;
   // std::cout<<"matrix "<<matA<<std::endl;
    cv::Vec4f plane;
    plane[0]= matA.determinant();
    plane[1]= matB.determinant();
    plane[2]= matC.determinant();
    plane[3]= -1*matD.determinant();

   // std::cout<< " plane eq:: "<<plane[0]<<", "<<plane[1]<<", "<<plane[2]<<", "<<plane[3]<<std::endl;
    return plane;
}


cv::Point2f Math2D::weightWrtThreePoints(cv::Point3f &pt1, cv::Point3f &pt2, cv::Point3f &pt3, cv::Point3f &plane_pt)
{
    Eigen::Vector3d v0,v1,v2,p;
    v0(0)=pt1.x;
    v0(1)=pt1.y;
    v0(2)=pt1.z;

    v1(0)=pt2.x;
    v1(1)=pt2.y;
    v1(2)=pt2.z;

    v2(0)=pt3.x;
    v2(1)=pt3.y;
    v2(2)=pt3.z;

    p(0)=plane_pt.x;
    p(1)=plane_pt.y;
    p(2)=plane_pt.z;

    cv::Point2f weight;

//    std::cout<<" (p.dot(v0.cross(v1))) "<<(p.dot(v0.cross(v1)))<<std::endl;
//    std::cout<<" (v2.dot(v0.cross(v1))) "<<(v2.dot(v0.cross(v1)))<<std::endl;

//    std::cout<<"(p.dot(v1.cross(v2))) "<<(p.dot(v1.cross(v2)))<<std::endl;
//    std::cout<<" (v0.dot(v1.cross(v2))) "<<(v0.dot(v1.cross(v2)))<<std::endl;

    //TODO:: if check for denominator =0 is required
    weight.x = (p.dot(v0.cross(v1)))/(v2.dot(v0.cross(v1)));
    weight.y = (p.dot(v1.cross(v2)))/(v0.dot(v1.cross(v2)));

    return weight;
}

cv::Point2f Math2D::parametric_point_wrt_triangle(cv::Point &triangle_point1, cv::Point &triangle_point2, cv::Point &triangle_point3, cv::Point2f &input_point)
{
//  float a = triangle_point1.x - triangle_point3.x;
//  float b = triangle_point2.x - triangle_point3.x;
//  float c = triangle_point3.x - input_point.x;
  
//  float d = triangle_point1.y - triangle_point3.y;
//  float e = triangle_point2.y - triangle_point3.y;
//  float f = triangle_point3.y - input_point.y;
  
  Eigen::Vector2d v0 , v1 , v2;
  
  v0( 0 ) = triangle_point3.x - triangle_point1.x;
  v0( 1 ) = triangle_point3.y - triangle_point1.y;
  
  v1( 0 ) = triangle_point2.x - triangle_point1.x;
  v1( 1 ) = triangle_point2.y - triangle_point1.y;
  
  v2( 0 ) = input_point.x - triangle_point1.x;
  v2( 1 ) = input_point.y - triangle_point1.y;
  
  double dot00 = v0.dot( v0 );
  double dot01 = v0.dot( v1 );
  double dot02 = v0.dot( v2 );
  double dot11 = v1.dot( v1 );
  double dot12 = v1.dot( v2 );
  
  // Compute barycentric coordinates
  double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
  double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
  double v = (dot00 * dot12 - dot01 * dot02) * invDenom;
  
  return cv::Point2f( v , u );
//   cv::Point2f parametric_point;
//   
//   float denom = b * f - c * e;
//   
//   if( a == 0 )
//   {
//     parametric_point.y = - c / b;
//     
//     parametric_point.x = ( -f + parametric_point.y * e ) / d;
//     
//     return parametric_point;
//   }
//   
//   if( b == 0 )
//   {
//     parametric_point.x = - c / a;
//     
//     parametric_point.y = ( -f + parametric_point.x * d ) / e;
//     
//     return parametric_point;
//   }
//   
//   if( d == 0 )
//   {
//     parametric_point.y = - f / e;
//     
//     parametric_point.x = ( -c + parametric_point.y * b ) / a;
//     
//     return parametric_point;
//   }
//   
//   if( e == 0 )
//   {
//     parametric_point.x = - f / d;
//     
//     parametric_point.y = ( -c + parametric_point.x * a ) / b;
//     
//     return parametric_point;
//   }
//   
//   parametric_point.x = ( a * e - d * b ) / denom;
//   parametric_point.y = -( parametric_point.x * a + c ) / b;
  
//   return parametric_point;
}



Math2D::~Math2D()
{

}

}

