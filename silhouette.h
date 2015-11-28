/*
    Copyright (c) <year>, <copyright holder>
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

    THIS SOFTWARE IS PROVIDED BY <copyright holder> ''AS IS'' AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL <copyright holder> BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef SILHOUETTE_H
#define SILHOUETTE_H

#include "opencvincludes.h"
#include "string"
#include "vector"
#include <generaltypes.h>


namespace tr{
  
class Silhouette
{
  int height , width;
  
  std::vector< cv::Point2f > contour_points;
  std::vector< cv::Point2f > boundary_points;
  std::vector< std::vector< cv::Point2f > > m_boundary_contours;
  std::vector< bool > mIsForeground;
  std::vector< std::vector< cv::Point2f > > mHierarchyContours;
  std::vector< cv::Vec4i > mHierarchy;
  
  cv::Mat bounding_image;
  
  bool fast_access_enabled;
  
  cv::Point2f bb_left_corner;

  std::string image_path;
  
  int pixel_margin;
  
  bool m_is_computed , m_use_outer_boundary;
  
  cv::Rect m_grab_cut_bb;
  
  std::vector< cv::Point > m_foreground_points;
  std::vector< cv::Point > m_background_points;
  
  float compression_ratio;
  
protected:
  
  void compute_bounding_box();
  
  void build_boundary();
  void build_boundary2();
  
protected:
  
  void preprocess_silhouette_image( cv::Mat &silhoutte_image  , cv::Vec3b foreground_color ,
                                   cv::Vec3b back_ground_color );
  
  void preprocess_silhouette_image( cv::Mat &silhoutte_image  , uchar foreground_color ,
                                   uchar back_ground_color );
  
  void compute_boundary_with_margin( std::vector< cv::Point2f > &boundary , int margin = 2 );
  
  void display_background_foreground( std::vector< cv::Point > &foreground , std::vector< cv::Point2f > &background );
  
public:
  
  Silhouette();
  
  bool is_computed(){ return m_is_computed; }
  
  bool is_inside( cv::Point2f point);
  
  bool isInside( Point2 &point );
  
  void set_image_path( std::string path );
  
  void set_foreground_background_points( std::vector< cv::Point2f > &foreground , std::vector< cv::Point2f > &background );
  
  void set_bounding_box( cv::Rect &bounding_box );
  
  void get_bounding_box( cv::Rect &bounding_box );
  
  void getBoundingBoxImage( cv::Mat &boundingImage );
  
  cv::Mat& boundingImage();
  
  std::string get_image_path();
  
  void build_silhouette_from_grab_cut( bool use_outer_boundary = false );
  
  void save_contour2(std::string file_path);
  
  void load_contour( std::string file_path );
  
  void load_contour2(std::vector< cv::Point2f >& point_sequence);
  
    bool load_contour2(std::string file_path);
  
  bool load_from_silhouette_image( cv::Mat &image , std::string correspoding_image_path, cv::Vec3b foreground_color ,
                                  cv::Vec3b back_ground_color  );
  
  bool load_from_silhouette_image( cv::Mat &image , std::string correspoding_image_path, uchar foreground_color ,
                                  uchar back_ground_color );
  
  
  
  std::vector< cv::Point2f > get_boundary_points();
  std::vector< std::vector< cv::Point2f > > get_boundary_contours();
  std::vector< bool > getForegroundContourInfo();
  
  void getHierarchyContours( std::vector< std::vector< cv::Point2f > > &hierarchyContours , std::vector< cv::Vec4i > &hierarchy );
  
  void display_bounding_box_image();
  
  void intersection_with_line( cv::Point2f point1 , cv::Point2f point2 , 
			       std::vector< cv::Point2f > &intersection_points );
  
  void intersection_with_line(  cv::Vec3f line, std::vector< cv::Point2f >& intersection_points);
  
  void intersectWithLine(  cv::Vec3f line, std::vector< cv::Point2f >& intersection_points );
  
  bool doesLineIntersect( cv::Vec3f line );
  
  void set_scale( float scale );
  
  float get_scale(); 
  
  void clear();
  
  ~Silhouette();
};

}

#endif // SILHOUTTE_H
