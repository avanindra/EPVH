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

#include "silhouette.h"
#include "iostream"
#include "fstream"

namespace tr
{
  
 Silhouette::Silhouette()
 {
  pixel_margin = 8;
  
  bb_left_corner.x = 0;
  bb_left_corner.y = 0;
  
  m_is_computed = false;
  
  height = 0;
  width = 0;
  
  compression_ratio = 1.0;

 }
 
 
void Silhouette::load_contour(std::string file_path)
{

}


void Silhouette::load_contour2(std::vector< cv::Point2f >& point_sequence)
{
  contour_points = point_sequence;
  
  compute_bounding_box();
  
  std::vector< cv::Point > _contour_points;
 
  int _num_pts = point_sequence.size();

  cv::Point *_points2 = new cv::Point[ _num_pts ];
  for( int pp = 0; pp < _num_pts; pp++ )
  {
   cv::Point2f _pt1 = point_sequence[ pp ];
   
   cv::Point2f _pt2 = point_sequence[ ( pp + 1 ) % _num_pts ];
   
   cv::Point _pt;
   
   _pt.x = point_sequence[ pp ].x - bb_left_corner.x;
   _pt.y = point_sequence[ pp ].y - bb_left_corner.y;
   
   _contour_points.push_back( _pt );

   _points2[ pp ] = _pt;
 }

  const cv::Point *_points = _points2;
  
  bounding_image.create( height , width , CV_8UC1 );
  
  bounding_image.setTo( cv::Scalar(0) );
  
  cv::fillPoly( bounding_image ,&_points, &_num_pts , 1 , cv::Scalar( 255 ) );
  
  for( int pp = 0; pp < _num_pts; pp++ )
  {
   cv::Point2f _pt1 = point_sequence[ pp ];
   
   cv::Point2f _pt2 = point_sequence[ ( pp + 1 ) % _num_pts ];
   
   _pt1.x = _pt1.x - bb_left_corner.x;
   _pt1.y = _pt1.y - bb_left_corner.y;

   _pt2.x = _pt2.x - bb_left_corner.x;
   _pt2.y = _pt2.y - bb_left_corner.y;
   
   cv::line(  bounding_image ,_pt1, _pt2, cv::Scalar( 100 )  );
   
   contour_points[ pp ] = _contour_points[ pp ];
  }
 
 
 
 build_boundary2();
 
}

void Silhouette::build_boundary()
{
  cv::Point2f _start = contour_points[ 0 ];
  
  int x = _start.x;
  int y = _start.y;
  
  uchar *_data = bounding_image.data;
  
  cv::Mat _copy_image;
  
  bounding_image.copyTo( _copy_image );
  
  
  uchar *_copy_data = _copy_image.data;
  
  
  if( *( _data + y * width + x ) != 100 )
  {
    if( *( _data + y * width + x + 1 ) == 100 )
    {
      x++;
      boundary_points.push_back( cv::Point2f( x + bb_left_corner.x , y + bb_left_corner.y ) );
      
      *( _copy_data + y * width + x  ) = 120;
    }
    else if( *( _data + ( y + 1 ) * width + x  ) == 100  )
    {
      y++;
      
      boundary_points.push_back( cv::Point2f( x + bb_left_corner.x , y + bb_left_corner.y ) );
      
      *( _copy_data + y * width + x  ) = 120;
    }
    else if(  *( _data + ( y + 1 ) * width + (x + 1 )  ) == 100 )
    {
      x++;
      y++;
      
      boundary_points.push_back( cv::Point2f( x + bb_left_corner.x , y + bb_left_corner.y ) );
      
      *( _copy_data + y * width + x  ) = 120;
    }
    else if( *( _data + y * width + x - 1 ) == 100 )
    {
      x--;
      
      boundary_points.push_back( cv::Point2f( x + bb_left_corner.x , y + bb_left_corner.y ) );
      
      *( _copy_data + y * width + x  ) = 120;
      
    }
    else if( *( _data + ( y - 1 ) * width + x  ) == 100  )
    {
      y--;
      
      boundary_points.push_back( cv::Point2f(  x + bb_left_corner.x , y + bb_left_corner.y  ) );
      
      *( _copy_data + y * width + x  ) = 120;
      
    }
    else if(  *( _data + ( y - 1 ) * width + (x - 1 )  ) == 100 )
    {
      x--;
      y--;
      
      boundary_points.push_back( cv::Point2f( x + bb_left_corner.x , y + bb_left_corner.y ) );
      
      *( _copy_data + y * width + x  ) = 120;
    }
  }
  else
  {
    boundary_points.push_back( cv::Point2f( x + bb_left_corner.x , y + bb_left_corner.y ) );  
    
    *( _copy_data + y * width + x  ) = 120;
  }
  
  
  while( 1 )
  {
    bool exit = true;
    
    if( *( _data + y * width + x + 1 ) == 100 && *( _copy_data  + y * width + x + 1 ) != 120 )
    {
      x++;
      boundary_points.push_back( cv::Point2f( x + bb_left_corner.x , y + bb_left_corner.y ) );
      
      *( _copy_data + y * width + x  ) = 120;
      
      exit = false;
    }
    else if( *( _data + ( y + 1 ) * width + x  ) == 100 && *( _copy_data + ( y + 1 ) * width + x   ) != 120 )
    {
      y++;
      
      boundary_points.push_back( cv::Point2f( x + bb_left_corner.x , y + bb_left_corner.y ) );
      
      *( _copy_data + y * width + x  ) = 120;
      
      exit = false;
    }
    else if(  *( _data + ( y + 1 ) * width + (x + 1 )  ) == 100 && *( _copy_data + ( y + 1 ) * width + (x + 1  ) ) != 120 )
    {
      x++;
      y++;
      
      boundary_points.push_back( cv::Point2f( x + bb_left_corner.x , y + bb_left_corner.y ) );
      
      *( _copy_data + y * width + x  ) = 120;
      
      exit = false;
    }
    else if( *( _data + y * width + x - 1 ) == 100 && *( _copy_data + y * width + x - 1 ) != 120 )
    {
      x--;
      
      boundary_points.push_back( cv::Point2f( x + bb_left_corner.x , y + bb_left_corner.y ));
      
      *( _copy_data + y * width + x  ) = 120;
      
      exit = false;
    }
    else if( *( _data + ( y - 1 ) * width + x  ) == 100 && *( _copy_data + ( y - 1 ) * width + x  ) != 120 )
    {
      y--;
      
      boundary_points.push_back( cv::Point2f( x + bb_left_corner.x , y + bb_left_corner.y ) );
      
      *( _copy_data + y * width + x  ) = 120;
      
      exit = false;
    }
    else if(  *( _data + ( y - 1 ) * width + (x - 1 )  ) == 100 && *( _copy_data + ( y - 1 ) * width + (x - 1 )  ) != 120 )
    {
      x--;
      y--;
      
      boundary_points.push_back( cv::Point2f( x + bb_left_corner.x , y + bb_left_corner.y ) );
      
      *( _copy_data + y * width + x  ) = 120;
      
      exit = false;
    }
    
    
    if( exit )
      break;
    
  }
  
}

void Silhouette::build_boundary2()
{
  
  uchar *_ptr = bounding_image.data;
  
  cv::Mat clone = bounding_image.clone();
  
  std::vector< std::vector< cv::Point > > _contours , _contours2;  
  
  std::vector< cv::Vec4i > hierarchy;
  
  cv::findContours( clone , _contours , hierarchy , CV_RETR_CCOMP , CV_CHAIN_APPROX_SIMPLE );

  cv::Mat clone2 = bounding_image.clone();

  cv::findContours(  clone2 , _contours2 , mHierarchy , CV_RETR_TREE , CV_CHAIN_APPROX_SIMPLE  );

  mHierarchyContours.clear();

  mHierarchyContours.resize( _contours2.size() );

  for( int ss = 0; ss < _contours2.size() ; ss++ )
  {
	  int size = _contours2[ ss ].size();

	  mHierarchyContours[ ss ].resize( size );

	  for( int ss2 = 0; ss2 < size ; ss2++ )
	  {
		  mHierarchyContours[ ss ][ ss2 ].x = _contours2[ ss ][ ss2 ].x;
		  mHierarchyContours[ ss ][ ss2 ].y = _contours2[ ss ][ ss2 ].y;
	  }
  }
  
#if 0  


  
  for( int y = 0; y < height; y++ )
    for( int x = 0; x < width; x++ )
    {
      if( *( _ptr + y * width + x ) == 100 )
      {
	boundary_points.push_back( cv::Point2f( x + bb_left_corner.x , y + bb_left_corner.y ) );
      }
    }
#else

int num_contours = _contours.size();

 mIsForeground.resize( num_contours , true );

m_boundary_contours.clear();

boundary_points.clear();

m_boundary_contours.resize( num_contours );

float inverse_scale = 1.0 / compression_ratio;

for( int contour = 0; contour < num_contours; contour++ )
{
  int num_contour_points = _contours[ contour ].size();
  
  if( hierarchy[ contour ][ 3 ] >= 0 )
  {
    mIsForeground[ contour ] = false;
  }
  
  for( int pp = 0; pp < num_contour_points; pp++ )
  {
     cv::Point2f pt = _contours[ contour ][ pp ];
     
     pt.x += bb_left_corner.x;
     pt.y += bb_left_corner.y;
     
     pt.x *= inverse_scale;
     pt.y *= inverse_scale;
    
     boundary_points.push_back( pt ); 
     
     m_boundary_contours[ contour ].push_back( pt );
  }  
}


int numHiewrarchyContours = mHierarchyContours.size();

for( int cc = 0; cc < numHiewrarchyContours; cc++ )
{
    int num_contour_points = mHierarchyContours[ cc ].size();

    for( int pp = 0; pp < num_contour_points; pp++ )
    {
       cv::Point2f &pt = mHierarchyContours[ cc ][ pp ];

       pt.x += bb_left_corner.x;
       pt.y += bb_left_corner.y;

       pt.x *= inverse_scale;
       pt.y *= inverse_scale;
    }
}



#endif

}

void Silhouette::preprocess_silhouette_image(cv::Mat& Silhouette_image, cv::Vec3b foreground_color, cv::Vec3b back_ground_color)
{
  int h = Silhouette_image.rows;
  
  int w = Silhouette_image.cols;
  
  uchar *_data = Silhouette_image.data;
  
  for( int x = 2; x < w - 2; x++ )
    for( int y = 2; y < h - 2; y++ )
    {
      int _r = *( _data + y * w * 3 + 3 * x );
      int _b = *( _data + y * w * 3 + 3 * x + 1);
      int _g = *( _data + y * w * 3 + 3 * x + 2);
      
      if( _r > 10 && _b > 10 && _g > 10 )
      {
	*( _data + y * w * 3 + 3 * x ) = foreground_color[ 0 ];
	*( _data + y * w * 3 + 3 * x + 1 ) = foreground_color[ 1 ];
	*( _data + y * w * 3 + 3 * x + 2) =  foreground_color[ 2 ];
	
      }
      else
      {
        *( _data + y * w * 3 + 3 * x ) = back_ground_color[ 0 ];
	*( _data + y * w * 3 + 3 * x + 1 ) = back_ground_color[ 1 ];
	*( _data + y * w * 3 + 3 * x + 2) =  back_ground_color[ 2 ];
      }
      
    }
}


void Silhouette::preprocess_silhouette_image(cv::Mat& Silhouette_image, uchar foreground_color, uchar back_ground_color)
{
  int h = Silhouette_image.rows;
  
  int w = Silhouette_image.cols;
  
  uchar *_data = Silhouette_image.data;
  
  for( int x = 2; x < w - 2; x++ )
    for( int y = 2; y < h - 2; y++ )
    {
      int _c = *( _data + y * w  + x );
      
      if( _c > 10  )
      {
	*( _data + y * w + x ) = foreground_color;
	
      }
      else
      {
        *( _data + y * w  + x ) = back_ground_color;
      }
      
    }
}


void Silhouette::compute_boundary_with_margin(std::vector< cv::Point2f >& boundary, int margin)
{
  int _num_boundary_points = boundary_points.size();
  
  cv::Mat _buffer;
  
  bounding_image.copyTo( _buffer );
  
  uchar margin_boundary_color = 150;
  
  uchar *b = _buffer.data;
  
  uchar *data = bounding_image.data;
  
  std::vector< cv::Point2f > _current_ring = boundary_points;
  
  int _ring = 0;
  
  while( 1 )
  {
    _ring++;
    
    std::vector< cv::Point2f > _next_ring;
    
    _num_boundary_points = _current_ring.size();
  
    for( int bb = 0; bb < _num_boundary_points; bb++ )
    {
      int x = _current_ring[ bb ].x - bb_left_corner.x;
      int y = _current_ring[ bb ].y - bb_left_corner.y;
    
      for( int i = -1 ; i <= 1; i++ )
       for( int j = -1; j <= 1; j++ )
       {
	 int x_b = x + i ;
	 int y_b = y + j ;
	
	 if( x_b > 0 && x_b < width && y_b > 0 && y_b < height )
	 {
	   if( data[ y_b * width + x_b ] == 0 && b[ y_b * width + x_b ] == 0 )
	   {
	     b[ y_b * width + x_b ] = margin_boundary_color;

	     _next_ring.push_back( cv::Point2f( x_b + bb_left_corner.x , y_b + bb_left_corner.y ) );
	   }
	 }
       }
  }
  
  _current_ring = _next_ring;
  
  if( _ring == margin )
    break;
  
  }
  
  boundary = _current_ring;
  
}


void Silhouette::display_background_foreground(std::vector< cv::Point >& foreground, std::vector< cv::Point2f >& background)
{
  cv::Mat _image = cv::imread( image_path );
  
  int _num_foreground_pts = foreground.size();
  
  int _num_background_pts = background.size();
  
  for( int ff = 0; ff < _num_foreground_pts; ff++ )
  {
    cv::circle( _image , foreground[ ff ] , 2 , cv::Scalar( 0 , 255 , 0 ) );  
  }
  
  for( int bb = 0; bb < _num_background_pts; bb++ )
  {
    cv::circle( _image , background[ bb ] , 2 , cv::Scalar( 0 , 0 , 255 ) );
  }
  
  cv::namedWindow( "foreground and background" , 0 );
  
  cv::imshow( "foreground and background" , _image );
  
  cv::waitKey();
  
}



void Silhouette::save_contour2(std::string file_path)
{
  
  std::fstream fs;
 
  fs.open( file_path.c_str() , std::ios::out|std::ios::binary );
  
  fs<<bb_left_corner.x<<" "<<bb_left_corner.y<<std::endl;
  
  fs<<width<<" "<<height<<std::endl;
  
  fs<<pixel_margin<<std::endl;
  
  int _size = width * height;
  
  fs.write( (char * )bounding_image.data , _size );

  fs<<image_path<<std::endl;
  
  fs.close();
  
}

bool Silhouette::load_contour2(std::string file_path)
{
  
  std::fstream fs;
 
  fs.open( file_path.c_str() , std::ios::in|std::ios::binary);
  
  char _line1[ 1028 ] , _line2[ 1028 ] , _line3[ 1028 ] , _line4[ 4096 ];
    
  fs.getline( _line1 , 1028 ,'\n' );
  
  std::string _end_quote = "Silhouette_DATA_END";
  
  std::string _temp = _line1; 
  
  int _pos = _temp.find_first_of(" ");
    
  bb_left_corner.x = atof( _temp.substr( 0 , _pos ).c_str() );
  
  _temp = _temp.substr( _pos + 1 , _temp.length() - _pos - 1 );
    
  bb_left_corner.y = atof( _temp.c_str() );
  
  fs.getline( _line2 , 1028 ,'\n' );
  
  _temp = _line2; 
  
  _pos = _temp.find_first_of(" ");
    
  width = atoi( _temp.substr( 0 , _pos ).c_str() );
  
  _temp = _temp.substr( _pos + 1 , _temp.length() - _pos - 1 );
    
  height = atoi( _temp.c_str() );
  
  if( width < 5 && height < 5 )
  {
    clear();
    
    return false;
  }
  
  
  fs.getline( _line3 , 1028 ,'\n' );
  
  _temp = _line3;
    
  pixel_margin = atoi( _temp.c_str() );
  
  bounding_image.create( height , width , CV_8UC1 );
  
  char *_data = ( char * )bounding_image.data;
  
  fs.getline( _data , width * height ,'\n' );
  
  fs.getline( _line4 , 4096 , '\n' );
  
  image_path = _line4;
  
  m_grab_cut_bb.width = width;
  m_grab_cut_bb.height = height;
  
  m_grab_cut_bb.x = bb_left_corner.x;
  m_grab_cut_bb.y = bb_left_corner.y;
  
  
  build_boundary2();  
  
  m_is_computed = true;
  
  return true;
}


bool Silhouette::load_from_silhouette_image(   cv::Mat &image, std::string correspoding_image_path, cv::Vec3b foreground_color,
  			                    cv::Vec3b back_ground_color)
{
  
   preprocess_silhouette_image( image , foreground_color , back_ground_color );
  
   image_path = correspoding_image_path;
   
   m_is_computed = true;
   
   int _height = image.rows;
   int _width = image.cols;
   
   uchar *_data = image.data;
   
   int _fr = foreground_color[ 2 ];
   int _fg = foreground_color[ 1 ];
   int _fb = foreground_color[ 0 ];
   
   int _br = back_ground_color[ 2 ];
   int _bg = back_ground_color[ 1 ];
   int _bb = back_ground_color[ 0 ];
   
   cv::Point2f _left( _width , 0 ) , _right( 0 , 0 ) , _top( 0 , _height  ) , _bottom( 0 , 0 );
   
   for( int x = 0; x < _width; x++ )
     for( int y = 0; y < _height; y++ )
     {
       int _b = *( _data + y * _width * 3 + 3 * x );
       int _g = *( _data + y * _width * 3 + 3 * x + 1);
       int _r = *( _data + y * _width * 3 + 3 * x + 2);
       
       if( _r == _fr && _b == _fb && _g == _fb )
       {
	 if( x < _left.x )
	 {
	   _left = cv::Point2f( x , y );  
	 }
	 
	 if( x > _right.x )
	 {
	   _right = cv::Point2f( x , y );  
	 }
	 
	 if( y < _top.y )
	 {
	   _top = cv::Point2f( x , y );  
	 }
	 
	 if( y >_bottom.y )
	 {
	   _bottom = cv::Point2f( x , y );  
	 }
       }
     }
     
   bb_left_corner.x = _left.x - pixel_margin;
   bb_left_corner.y = _top.y - pixel_margin;
   
   width = ceil( _right.x - _left.x ) + 2 * pixel_margin;
   
   height = ceil( _bottom.y - _top.y ) + 2 * pixel_margin;

   if( height <= 4 || width <= 4 )
   {
     clear();

     return false;
   }
   
   bounding_image.create( height , width , CV_8UC1 );
   
   bounding_image.setTo( cv::Scalar( 0 ) );
   
   uchar *_data2 = bounding_image.data;
   
   for( int x = 0; x < width; x++ )
     for( int y = 0; y < height; y++ )
     {
       int _x = x + bb_left_corner.x;
       int _y = y + bb_left_corner.y;
       
       if( _x < 1 || _x >= _width - 1 || _y < 1 || _y >= _height - 1 )
	 continue;
       
       int _r = *( _data + _y * _width * 3 + 3 * _x );
       int _g = *( _data + _y * _width * 3 + 3 * _x + 1 );
       int _b = *( _data + _y * _width * 3 + 3 * _x + 2 );
       
       if( _r == _fr && _b == _fb && _g == _fb )
       {
	 bool is_boundary = false;
	 
	 for( int x2 = _x - 1; x2 <= _x + 1; x2++ )
	   for( int y2 = _y - 1; y2 <= _y + 1; y2++ )
	   {
	       int _r2 = *( _data + y2 * _width * 3 + 3 * x2 );
               int _g2 = *( _data + y2 * _width * 3 + 3 * x2 + 1 );
               int _b2 = *( _data + y2 * _width * 3 + 3 * x2 + 2 );	       
	       
	       if(  !(_r2 == _fr && _b2 == _fb && _g2 == _fb ) )
	       {
		 is_boundary = true; 
	       }
	   }
	 
	 if( is_boundary )
	 {
	   *( _data2 + y * width + x ) = 100;  
	 }
	 else
	 {
	   *( _data2 + y * width + x ) = 255;
	 }
       }
       
       
     }    
     
     build_boundary2();    
     
     m_grab_cut_bb.x = bb_left_corner.x;
     m_grab_cut_bb.y = bb_left_corner.y;
     
     m_grab_cut_bb.width = bounding_image.cols;
     m_grab_cut_bb.height = bounding_image.rows;

     return true;
     
}

bool Silhouette::load_from_silhouette_image(cv::Mat &image, std::string correspoding_image_path, uchar foreground_color, uchar back_ground_color)
{
 
   preprocess_silhouette_image( image , foreground_color , back_ground_color );
  
   image_path = correspoding_image_path;
   
   m_is_computed = true;
   
   int _height = image.rows;
   int _width = image.cols;
   
   uchar *_data = image.data;
   
   int _fc = foreground_color;
   
   int _bc = back_ground_color;
   
   std::cout<<"channels : "<<image.channels()<<std::endl;
   
   cv::Point2f _left( _width , 0 ) , _right( 0 , 0 ) , _top( 0 , _height  ) , _bottom( 0 , 0 );
   
   for( int x = 0; x < _width; x++ )
     for( int y = 0; y < _height; y++ )
     {
       uchar c = *( _data + y * _width  + x );
       
       if( c == _fc )
       {
	 if( x < _left.x )
	 {
	   _left = cv::Point2f( x , y );  
	 }
	 
	 if( x > _right.x )
	 {
	   _right = cv::Point2f( x , y );  
	 }
	 
	 if( y < _top.y )
	 {
	   _top = cv::Point2f( x , y );  
	 }
	 
	 if( y >_bottom.y )
	 {
	   _bottom = cv::Point2f( x , y );  
	 }
       }
     }
     
   bb_left_corner.x = _left.x - pixel_margin;
   bb_left_corner.y = _top.y - pixel_margin;
   
   width = ceil( _right.x - _left.x ) + 2 * pixel_margin;
   
   height = ceil( _bottom.y - _top.y ) + 2 * pixel_margin;
   
   bounding_image.create( height , width , CV_8UC1 );
   
   bounding_image.setTo( cv::Scalar( 0 ) );
   
   uchar *_data2 = bounding_image.data;
   
   for( int x = 0; x < width; x++ )
     for( int y = 0; y < height; y++ )
     {
       int _x = x + bb_left_corner.x;
       int _y = y + bb_left_corner.y;

       if( _x < 0 || _y < 0 || _x >= _width || _y >= _height )
           continue;
       
       uchar c = *( _data + _y * _width  +  _x );
       
       if( c == _fc  )
       {
	 bool is_boundary = false;
	 
	 for( int x2 = _x - 1; x2 <= _x + 1; x2++ )
	   for( int y2 = _y - 1; y2 <= _y + 1; y2++ )
	   {
           if( x2 < 0 || y2 < 0 || x2 >= _width || y2 >= _height )
               continue;

	       int _c2 = *( _data + y2 * _width + x2 );       
	       
	       if(  _c2 != _fc )
	       {
             is_boundary = true;
	       }
	   }
	 
	 if( is_boundary )
	 {
	   *( _data2 + y * width + x ) = 100;  
	 }
	 else
	 {
	   *( _data2 + y * width + x ) = 255;
	 }
       }
       
       
     }    
     
     build_boundary2();
     
     m_grab_cut_bb.x = bb_left_corner.x;
     m_grab_cut_bb.y = bb_left_corner.y;
     
     m_grab_cut_bb.width = bounding_image.cols;
     m_grab_cut_bb.height = bounding_image.rows;

     return true;
}


bool Silhouette::is_inside(cv::Point2f point)
{
  int x = point.x * compression_ratio - bb_left_corner.x;
  int y = point.y * compression_ratio - bb_left_corner.y;
  
  if( x < 0 || x >= width || y < 0 || y >= height )
    return false;
  
  return bounding_image.at< uchar >( y , x ) > 0;
  
}


bool Silhouette::isInside(Point2 &point)
{
  int x = point.x() * compression_ratio - bb_left_corner.x;
  int y = point.y() * compression_ratio - bb_left_corner.y;
  
  if( x < 0 || x >= width || y < 0 || y >= height )
    return false;
  
  return bounding_image.at< uchar >( y , x ) > 0;
}


void Silhouette::set_image_path(std::string path)
{
  image_path = path;
}


void Silhouette::compute_bounding_box()
{
  
  int _num_pts = contour_points.size();
  
  cv::Point2f _left_corner , _left , _right , _top , _bottom;
  
  float _height = 0, _width = 0; 
  
  for( int pp = 0; pp < _num_pts; pp++ )
  {
     if( pp == 0 )
     {
       _left = contour_points[ pp ];
       _right = contour_points[ pp ];
       
       _top = contour_points[ pp ];
       _bottom = contour_points[ pp ];  
       
     }
     else
     {
       if( _left.x > contour_points[ pp ].x )
       {
	  _left = contour_points[ pp ]; 
       }
       
       if( _right.x < contour_points[ pp ].x  )
       {
	 _right = contour_points[ pp ]; 
       }
       
       if( _top.y > contour_points[ pp ].y )
       {
	 _top = contour_points[ pp ]; 
       }
       
       if( _bottom.y < contour_points[ pp ].y )
       {
	 _bottom = contour_points[ pp ]; 
       }     
       
     }
  }
  
  _height = ceil( _bottom.y - _top.y );
  _width = ceil( _right.x - _left.x );
  
  //leave 4 pixel margin at boundries  
  _height += 2 * pixel_margin;
  
  _width += 2 * pixel_margin;
  
  height = _height;
  
  width = _width;
  
  bb_left_corner.x = _left.x - pixel_margin;
  bb_left_corner.y = _top.y - pixel_margin; 

}

std::vector< cv::Point2f > Silhouette::get_boundary_points()
{
  return boundary_points;
}

std::vector< std::vector< cv::Point2f > > Silhouette::get_boundary_contours()
{
  return m_boundary_contours;
}


std::vector< bool > Silhouette::getForegroundContourInfo()
{
  return mIsForeground;
}


void Silhouette::getHierarchyContours( std::vector< std::vector< cv::Point2f > >& hierarchyContours , std::vector< cv::Vec4i >& hierarchy )
{
  hierarchyContours = mHierarchyContours;
  hierarchy = mHierarchy;
}



void Silhouette::display_bounding_box_image()
{
  cv::namedWindow( "contour_image" );
  
  cv::imshow( "contour_image" , bounding_image );
  
  cv::waitKey();  
}

void Silhouette::intersection_with_line(  cv::Point2f point1, cv::Point2f point2, 
					 std::vector< cv::Point2f >& intersection_points)
{
  
  cv::Point2f _pt1 , _pt2;
  
  _pt1.x = point1.x * compression_ratio - bb_left_corner.x;
  _pt1.y = point1.y * compression_ratio - bb_left_corner.y;
  
  _pt2.x = point2.x * compression_ratio - bb_left_corner.x;
  _pt2.y = point2.y * compression_ratio - bb_left_corner.y;
  
  cv::LineIterator _line( bounding_image , _pt1 , _pt2 );
  
  int _num_pixels =  _line.count;
  
  for( int pp = 0; pp < _num_pixels; pp++ )
  {
    uchar *_color =  *_line;
    
    if( *_color == 100 )
    {
      cv::Point2f _pt = _line.pos();
      
      _pt.x += bb_left_corner.x;
      
      _pt.y += bb_left_corner.y;
      
      intersection_points.push_back( _pt );  
    }
    
    _line++;
  }
  
}

void Silhouette::intersection_with_line(  cv::Vec3f line, std::vector< cv::Point2f >& intersection_points)
{  
  cv::Point2f _end1 , _end2( 0 , 0 );
  
  cv::Point2f _pt1 , _pt2 , _pt3 , _pt4;
  
  float _x_end1 = bb_left_corner.x - pixel_margin , _x_end2 = bb_left_corner.x + width + pixel_margin;
  float _y_end1 = bb_left_corner.y - pixel_margin , _y_end2 = bb_left_corner.y + height + pixel_margin;
  
  _pt1.x = _x_end1;
  
  _pt1.y = - ( line[ 2 ] + line[ 0 ] * _x_end1 ) / line[ 1 ] ;
  
  _pt2.x = _x_end2;
  
  _pt2.y = - ( line[ 2 ] + line[ 0 ] * _x_end2 ) / line[ 1 ] ;
  
  _pt3.y = _y_end1;
  
  _pt3.x = - ( line[ 2 ] + line[ 1 ] * _y_end1 ) / line[ 0 ] ;
  
  _pt4.y = _y_end2;
  
  _pt4.x = - ( line[ 2 ] + line[ 1 ] * _y_end2 ) / line[ 0 ] ;
  
  
  if( _pt1.y > _y_end1 &&  _pt1.y < _y_end2 )
  {    
    _end1 = _pt1;
  }
  
  if( _pt2.y > _y_end1 &&  _pt2.y < _y_end2 )
  {
    _end2 = _pt2;
  }
  
  if(  _pt3.x > _x_end1 &&  _pt3.x < _x_end2  )
  {
    if( _end2.x < _pt1.x )
    {
      _end2 = _pt3;
    }
    else
    {
      _end1 = _pt3;
    }
  }
  
  
  if(  _pt4.x > _x_end1 &&  _pt4.x < _x_end2  )
  {
    if( _end2.x < _pt1.x )
    {
      _end2 = _pt4;
    }
    else
    {
      _end1 = _pt4;
    }
  }
  
  intersection_with_line( _end1 , _end2 , intersection_points );
  
}


// void Silhouette::intersectWithLine(cv::Vec3f line, std::vector< cv::Point2f >& intersection_points)
// {
//     float l1 = line( 0 );
//     float l2 = line( 1 );
// 
//     cv::Point pt1 , pt2;
// 
//     if( abs( l2 ) < abs( l1 ) )
//     {
//         pt1.y = 0 ;
//         pt2.y = h ;
// 
//         pt1.x = -( line( 1 ) * pt1.y + line( 2 )  ) / line( 0 );
//         pt2.x = -( line( 1 ) * pt2.y + line( 2 )  ) / line( 0 );
//     }
//     else
//     {
//         pt1.x = 0;
//         pt2.x = w;
// 
//         pt1.y = -( line( 0 ) * pt1.x + line( 2 )  ) / line( 1 );
//         pt2.y = -( line( 0 ) * pt2.x + line( 2 )  ) / line( 1 );
//     }
//     
//     pt1.x = ( pt1.x - offset.x ) * silhouettScale;
//     pt1.y = ( pt1.y - offset.y ) * silhouettScale;
//     
//     pt2.x = ( pt2.x - offset.x ) * silhouettScale;
//     pt2.y = ( pt2.y - offset.y ) * silhouettScale;
//     
//     
//   
//   cv::LineIterator lineIter( binImage , pt1 , pt2 );
//   
//   int numPixels = lineIter.count;
//   
//   if( numPixels == 0 )
//     return;
//   
//   int prevVal = 0;
//   
//   
//   for( int pix = 0; pix < numPixels; pix++ )
//   {
//     uchar val = *( *lineIter );
//     
//     if( ( val && !prevVal )  || ( !val && prevVal ) )
//     {
//       intersections.push_back( lineIter.pos() );
//       
//       prevVal = val;
//     }
//     
//     lineIter++;
//     
//   }
//   
//   if( prevVal )
//   {
//     intersections.push_back( lineIter.pos() );
//   }
//   
//   assert( intersections.size() % 2 == 0 );
//   
// }
// 

void Silhouette::set_scale(float scale)
{
  compression_ratio = scale;
}

float Silhouette::get_scale()
{
  return compression_ratio;
}



void Silhouette::clear()
{  
  if( m_is_computed )
  {
    contour_points.clear();
    boundary_points.clear();
    bounding_image.release();
  
    width = 0;
    height = 0;
  
    m_is_computed = false;
  }
  
}

void Silhouette::set_bounding_box(cv::Rect& bounding_box)
{
  m_grab_cut_bb = bounding_box;
  
  m_grab_cut_bb.x *= compression_ratio;
  m_grab_cut_bb.y *= compression_ratio;
  
  m_grab_cut_bb.width *= compression_ratio;
  m_grab_cut_bb.height *= compression_ratio;
}

void Silhouette::get_bounding_box(cv::Rect& bounding_box)
{
  bounding_box = m_grab_cut_bb;
  
  m_grab_cut_bb.x /= compression_ratio;
  m_grab_cut_bb.y /= compression_ratio;
  
  m_grab_cut_bb.width /= compression_ratio;
  m_grab_cut_bb.height /= compression_ratio;
}

void Silhouette::getBoundingBoxImage(cv::Mat& boundingImage)
{
//   boundingImage = bounding_image;
  
  bounding_image.copyTo( boundingImage );
}


cv::Mat& Silhouette::boundingImage()
{
  return bounding_image;
}


std::string Silhouette::get_image_path()
{
  return image_path;
}



void Silhouette::set_foreground_background_points(std::vector< cv::Point2f >& foreground, std::vector< cv::Point2f >& background)
{
  m_foreground_points.clear();
  m_background_points.clear();
  
  m_foreground_points.resize( foreground.size() );   
  m_background_points.resize( background.size() );
  
  int bg_size = background.size() , fg_size = foreground.size();
  
  for( int bb = 0; bb < bg_size; bb++ )
  {
    m_background_points[ bb ].x = background[ bb ].x * compression_ratio;  
    m_background_points[ bb ].y = background[ bb ].y * compression_ratio;
  }  
  
  for( int ff = 0; ff < fg_size; ff++ )
  {
    m_foreground_points[ ff ].x = foreground[ ff ].x * compression_ratio;
    m_foreground_points[ ff ].y = foreground[ ff ].y * compression_ratio;    
  }
  
}

bool Silhouette::doesLineIntersect( cv::Vec3f line )
{
   cv::Point2f end1 , end2( 0 , 0 );
   
   bool isSilhouetteIntersected = false;
   
   float _x_end1 = bb_left_corner.x , _x_end2 = bb_left_corner.x + width + 2 * pixel_margin;
   
   end1.x = _x_end1;
  
   end1.y = - ( line[ 2 ] + line[ 0 ] * _x_end1 ) / line[ 1 ] ;
  
   end2.x = _x_end2;
  
   end2.y = - ( line[ 2 ] + line[ 0 ] * _x_end2 ) / line[ 1 ] ;
   
   cv::Rect bb;
   
   end1.x -= bb_left_corner.x;
   end1.y -= bb_left_corner.y;
   
   end2.x -= bb_left_corner.x;
   end2.y -= bb_left_corner.y;
   
   cv::LineIterator lineIter( bounding_image , end1 , end2  );
   
   int numPixels = lineIter.count;
   
   for( int pp = 0; pp < numPixels; pp++ )
   {
      if( *( *lineIter ) )
      {
	isSilhouetteIntersected = true; 
      }  
      
      lineIter++;
   }   
   
   return isSilhouetteIntersected;
}


void Silhouette::build_silhouette_from_grab_cut( bool use_outer_boundary )
{
  cv::Mat image_orig = cv::imread( image_path );
  
  cv::Size new_size ;
  
  new_size.height = image_orig.rows * compression_ratio;
  new_size.width = image_orig.cols * compression_ratio;
  
  cv::Mat image;
  
  cv::resize( image_orig , image , new_size );
  
  cv::Mat mask( image.size() , CV_8UC1 );
  
  cv::Mat bgdModel, fgdModel;
  
  m_grab_cut_bb.x = std::max( m_grab_cut_bb.x , 2 );
  m_grab_cut_bb.x = std::min( m_grab_cut_bb.x , image.cols - 2 );
  
  m_grab_cut_bb.y = std::max( m_grab_cut_bb.y , 2 );
  m_grab_cut_bb.y = std::min( m_grab_cut_bb.y , image.rows - 2 );

  m_grab_cut_bb.width = std::min( m_grab_cut_bb.width , image.cols - m_grab_cut_bb.x  - 3 );
  m_grab_cut_bb.height = std::min( m_grab_cut_bb.height , image.rows - m_grab_cut_bb.y  - 3 );
  
  mask.setTo( cv::GC_BGD );
   
  ( mask( m_grab_cut_bb ) ).setTo( cv::Scalar( cv::GC_PR_FGD ) );
  
  int _num_fg_pts = m_foreground_points.size();
  
  int _num_bg_pts = m_background_points.size();
  
  cv::grabCut( image , mask, m_grab_cut_bb  , bgdModel, fgdModel, 1 , cv::GC_INIT_WITH_RECT);
  
  for( int ff = 0; ff < _num_fg_pts; ff++ )
  {
    if( use_outer_boundary )
    {
      if( !is_inside( m_foreground_points[ ff ] ) )
        continue;
    }
    
    cv::Point pt;
    
    pt.x  = m_foreground_points[ ff ].x;
    pt.y =  m_foreground_points[ ff ].y;
    
    cv::circle( mask, pt , 2 , cv::GC_FGD, 2 );
  }
  
  
if( use_outer_boundary )
{
  m_background_points.clear();
  
  std::vector< cv::Point2f > _background_points;
  
  compute_boundary_with_margin( _background_points , 10);
  
  _num_bg_pts = _background_points.size();
  
  for( int bb = 0; bb < _num_bg_pts; bb++ )
  {    
    cv::Point pt;
    
    pt.x  = _background_points[ bb ].x;
    pt.y =  _background_points[ bb ].y;
    
    cv::circle( mask, pt , 1 , cv::GC_BGD, 1 );
  }
  
}
else
{
  _num_bg_pts = m_background_points.size();
  
  for( int bb = 0; bb < _num_bg_pts; bb++ )
  {
    
    cv::Point pt;
    
    pt.x  = m_background_points[ bb ].x;
    pt.y =  m_background_points[ bb ].y;
    
    cv::circle( mask, pt , 1 , cv::GC_BGD, 1 );
  }
} 

  cv::grabCut( image , mask, m_grab_cut_bb  , bgdModel, fgdModel, 3 , cv::GC_INIT_WITH_MASK);
  
  cv::Mat binMask;
  
  if( mask.empty() || mask.type()!=CV_8UC1 )
      CV_Error( CV_StsBadArg, "comMask is empty or has incorrect type (not CV_8UC1)" );
    
    binMask.create( mask.size(), CV_8UC1 );
   
    binMask = mask & 1;
    
  int width = mask.cols;
  int height = mask.rows;
  
  uchar *mask_data = binMask.data;    
  
  for( int r = 0; r < height; r++ )
    for( int c = 0; c < width; c++ )
    {
      int id = r * width + c;
      
      if( mask_data[ id ] == 1 )
      {
	mask_data[ id ] = 255;
      }
    }
    

  clear();
  
  load_from_silhouette_image( binMask , image_path , 255 , 0 );
  
}

 
 
Silhouette::~Silhouette()
{
}


  
}
