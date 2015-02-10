   /*
 * Copyright (c) 2012, avanindra <avanindra.singh@gmail.com>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *     names of its contributors may be used to endorse or promote products
 *     derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY avanindra <email> ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL avanindra <email> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include "segmentclipper.h"
#include <utility-2d/math2d.h>
#include "util/opencvincludes.h"
#include "utility-3d/display3droutines.h"
#include "vtkPolyData.h"
#include "vtkIdList.h"
#include "vtkUnsignedCharArray.h"
#include "vtkCellData.h"
#include "vtkSmartPointer.h"

#define TR_VALID_INTERSECTION( p )   ( p.first >= 0 && p.first < 1 && p.second >= 0 && p.second < 1 )


bool indexValSortPredicate( const IndexValPair &p1 , const IndexValPair &p2  )
{
  return p1.second < p2.second;
}


SegmentClipper::SegmentClipper()
{
  
}


void SegmentClipper::add( int contourId , int segmentId , tr::Point2f end1 , tr::Point2f end2 )
{
  IndexValPair hPair , vPair;
  
  hPair.first.mContourId = contourId;
  hPair.first.mSegmentId = segmentId;
  
  vPair.first.mSegmentId = segmentId;
  vPair.first.mContourId = contourId;
  
  if( end1.x() < end2.x() )
  {
    hPair.first.mHorizontalPose = 1;
    
    hPair.second = end1.x();
    
    mHorizontalArray.push_back( hPair );
    
  }
  else
  {
    hPair.first.mHorizontalPose = 0;
    
    hPair.second = end2.x();
    
    mHorizontalArray.push_back( hPair );
    
  }

  
  if( end1.y() < end2.y() )
  {
    vPair.first.mVerticalPose = 1;
    
    vPair.second = end1.y();
    
    mVerticalArray.push_back( vPair );
    
  }
  else
  {
    vPair.first.mVerticalPose = 0;
    
    vPair.second = end2.y();
    
    mVerticalArray.push_back( vPair );
    
  }
  
  
}

void SegmentClipper::add( int contourId , int segmentId , cv::Point2f end1 , cv::Point2f end2 )
{
  IndexValPair hPair , vPair;
  
  hPair.first.mContourId = contourId;
  hPair.first.mSegmentId = segmentId;
  
  vPair.first.mSegmentId = segmentId;
  vPair.first.mContourId = contourId;
  
  if( end1.x < end2.x )
  {
    hPair.first.mHorizontalPose = 1;
    
    hPair.second = end1.x;
    
    mHorizontalArray.push_back( hPair );
    
  }
  else
  {
    hPair.first.mHorizontalPose = 0;
    
    hPair.second = end2.x;
    
    mHorizontalArray.push_back( hPair );
    
  }

  
  if( end1.y < end2.y )
  {
    vPair.first.mVerticalPose = 1;
    
    vPair.second = end1.y;
    
    mVerticalArray.push_back( vPair );    
  }
  else
  {
    vPair.first.mVerticalPose = 0;
    
    vPair.second = end2.y;
    
    mVerticalArray.push_back( vPair ); 
  }

}



void SegmentClipper::build()
{
  std::sort( mHorizontalArray.begin() , mHorizontalArray.end() , indexValSortPredicate );
  std::sort( mVerticalArray.begin() , mVerticalArray.end() , indexValSortPredicate ); 

  mHorizontalStack.resize( mHorizontalArray.size() );
  mVerticalStack.resize( mVerticalArray.size() );
  
  int numElems = mHorizontalArray.size();
  
  
  mUsedStack.resize( mHorizontalArray.size() );
  
  mStackTop = -1;
  
  for( int ee = 0; ee < numElems; ee++ )
  {
    int contourId1 = mHorizontalArray[ ee ].first.mContourId;
    int segmentId1 = mHorizontalArray[ ee ].first.mSegmentId;
    
    std::pair< int , int > segment , segment2;
    
    segment.first = contourId1;
    segment.second = segmentId1;
    
    mHorizontalStack[ ee ].push_back( segment );
    
    tr::Point2f &firstEnd1 = mPolygon[ contourId1 ][ segmentId1  ];
    tr::Point2f &secondEnd1 = mPolygon[ contourId1 ][ ( segmentId1 + 1 ) % mPolygon[ contourId1 ].size() ];
    
    float xMin , xMax;
    
    if( firstEnd1.x() < secondEnd1.x() )
    {
      xMax =  secondEnd1.x();
    }
    else
    {
      xMax =  firstEnd1.x();
    }
    
    int id = ee;
    
    while( 1 )
    {
      id++;
      
      if( id >= numElems  )
      {
	    break;
      }
      
      if( xMax > mHorizontalArray[ id ].second )
      {
	    mHorizontalStack[ id ].push_back( segment );
      }
      else
      {
	    break;
      }
    }
    
    int contourId2 = mVerticalArray[ ee ].first.mContourId;
    int segmentId2 = mVerticalArray[ ee ].first.mSegmentId;
    
    segment2.first = contourId2;
    segment2.second = segmentId2;
    
    mVerticalStack[ ee ].push_back( segment2 );
    
    tr::Point2f &firstEnd2 = mPolygon[ contourId2 ][ segmentId2  ];
    tr::Point2f &secondEnd2 = mPolygon[ contourId2 ][ ( segmentId2 + 1 ) % mPolygon[ contourId2 ].size() ];
    
    float yMin , yMax;
    
    if( firstEnd2.y() < secondEnd2.y() )
    {
      yMax = secondEnd2.y();
    }
    else
    {
      yMax =  firstEnd2.y();
    }
    
    id = ee;
    
    while( 1 )
    {
      id++;
      
      if( id >= numElems  )
      {
	break;
      }
      
      if( yMax > mVerticalArray[ id ].second )
      {
	mVerticalStack[ id ].push_back( segment2 );
      }
      else
      {
	break;
      }
    }
    
  }
  
}


bool SegmentClipper::search( tr::Point2f end1 , tr::Point2f end2 , int &index )
{
  float delX = std::abs( end1.x() - end2.x() );
  float delY = std::abs( end1.y() - end2.y() );
  
  
  if( delX > delY )
  {
    float yMin = std::min( end1.y() , end2.y() );
    
    index = binarySearch( yMin , mVerticalArray );    
    
    return false;
  }
  else
  {
    float xMin = std::min( end1.x() , end2.x() );
    
    index = binarySearch( xMin , mHorizontalArray );
    
    return true;
  }
  
}



int SegmentClipper::binarySearch(float val, const std::vector< IndexValPair >& searchArray)
{
   int validIndex = -1;
   
   int size = searchArray.size();
   
   int currentIndex = size / 2;  
   
   
   if( size < 2 )
     return validIndex;
   
   
   if( val <= searchArray.front().second   )
   {
      validIndex = 0; 
     
      return validIndex;
   }
   
   if( val >= searchArray.back().second )
   {
     validIndex = searchArray.size() - 1; 
     
     return validIndex;
   }

   int end1 = 0;
   
   int end2 = size - 1;
   
   validIndex = currentIndex;
   
   while( 1 )
   { 

     if( searchArray[ validIndex ].second > val && searchArray[ validIndex - 1  ].second <= val )
     {
       break;
     }
     else if( searchArray[ validIndex ].second >= val )
     {
       
       end2 = validIndex;
       
       validIndex = ( validIndex + end1 ) / 2 ;
       
       if( validIndex == end2 )
	 break;
     }
     else
     {
       end1 = validIndex;
       
       validIndex = ( validIndex + end2 ) / 2;
       
       if( validIndex == end1 )
	 break;
     }
     
   }
   
   
   return validIndex;
}



bool SegmentClipper::isInsideSegment(tr::Point2f point, std::pair< int, int > strip)
{
  bool isInside = false;
  
  if( strip.first >= mPolygon.size() )
  {
    return false;
  }
  
  int numContourPoints = mPolygon[ strip.first ].size(); 
  
  if( strip.second > numContourPoints )
  {
    return false;
  }
  
  tr::Point2f &point1 = mPolygon[ strip.first ][ strip.second ];
  tr::Point2f &point2 = mPolygon[ strip.first ][ ( strip.second + 1 ) % numContourPoints ];
  
  float val = ( point.y() - point1.y() ) * ( point2.x() - point1.x() ) - ( point2.y() - point1.y() ) * ( point.x() - point1.x() );

  return ( val < 0 );
}


void SegmentClipper::displaySegments( std::vector< std::pair< tr::Point2f , tr::Point2f >  > &strips )
{
  
  
  vtkSmartPointer< vtkPolyData > edgePolyData = vtkSmartPointer< vtkPolyData >::New();
    
  vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();    
    
  colors->SetName("colors");
  colors->SetNumberOfComponents(3);
  
  
  int numEdges = 0 ;
  
  for( int ee1 = 0; ee1 < mPolygon.size(); ee1++ )
  {
    numEdges += mPolygon[ ee1 ].size();
  }
  
  vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();
  
  int numPoints = 2 * numEdges + 2 * strips.size();
  
  points->Allocate( numPoints );
    
  int numCells = numEdges + strips.size();
  
  edgePolyData->Allocate( numCells );
  
  vtkSmartPointer< vtkIdList > cell = vtkSmartPointer< vtkIdList >::New();
  
  for( int ee1 = 0; ee1 < mPolygon.size(); ee1++ )
    for( int ee2 = 0; ee2 < mPolygon[ ee1 ].size(); ee2++ ) 
  {
    
    tr::Point2f &pt1f = mPolygon[ ee1 ][ ee2 ];
    tr::Point2f &pt2f = mPolygon[ ee1 ][ ( ee2 + 1 ) % mPolygon[ ee1 ].size() ];
    
    tr::Point3d pt1( pt1f.x() , pt1f.y() , 0 ) , pt2( pt2f.x() , pt2f.y() , 0 );
    
    int id1 = points->InsertNextPoint( pt1.data() );
    int id2 = points->InsertNextPoint( pt2.data() );
    
    cell->Reset();
    
    cell->InsertNextId( id1 );
    cell->InsertNextId( id2 );   
    
    edgePolyData->InsertNextCell( VTK_LINE , cell );
    
    const unsigned char _color[] = { 255 , 255 , 255 };
  
    colors->InsertNextTupleValue(_color);
  }
  
  
  for( int ee = 0; ee < strips.size(); ee++ )
  {
    
    tr::Point2f &pt1f = strips[ ee ].first;
    tr::Point2f &pt2f = strips[ ee ].second;
    
    tr::Point3d pt1( pt1f.x() , pt1f.y() , 0 ) , pt2( pt2f.x() , pt2f.y() , 0 );
    
    int id1 = points->InsertNextPoint( pt1.data() );
    int id2 = points->InsertNextPoint( pt2.data() );
    
    cell->Reset();
    
    cell->InsertNextId( id1 );
    cell->InsertNextId( id2 );
    
    edgePolyData->InsertNextCell( VTK_LINE , cell );
    
    const unsigned char _color[] = { 255 , 0 , 0 };
  
    colors->InsertNextTupleValue(_color);
    
  }
  
  edgePolyData->SetPoints( points );
  
  edgePolyData->GetCellData()->SetScalars( colors );
  
  tr::Display3DRoutines::displayPolyData( edgePolyData );
  
}


bool SegmentClipper::isInsideForeground( tr::Point2f point )
{    

  std::vector< uchar > isUsed( mContourHierarchy.size() , 0 );
  
  bool isInside = false;
  
  cv::Point2f pt( point.x() , point.y() );
  
  int numContours = mContourHierarchy.size();

  for( int hh = 0; hh < numContours ; hh++ )
  {
    int isInsideCurrentChain = false;
    
    if( mContourHierarchy[ hh ][ 3 ] == -1 )
    {

      int child = mContourHierarchy[ hh ][ 2 ];
      
      isInsideCurrentChain = ( cv::pointPolygonTest( mCvPolygon[ hh ] , pt , false  ) > 0 );
 
      int next = child;
      
      if( !isInsideCurrentChain  )
      {
	    continue;
      }
      
      if( child != -1 )
      {
	    isInsideCurrentChain = ( isInsideCurrentChain && !checkChildren( hh , pt ) );
      }
      
    }  
    
    if( isInsideCurrentChain )
    {
      isInside = isInsideCurrentChain;
      
      break;
    }
  }
  
  return isInside;
  
}


bool SegmentClipper::checkChildren( int id, cv::Point2f point )
{
  int child = mContourHierarchy[ id ][ 2 ];
  
  bool isInsideChildren = false;
  
  if( child != -1 )
  {
    isInsideChildren = checkSiblings( child , point );
  }
  

  return isInsideChildren;
}


bool SegmentClipper::checkSiblings( int id , cv::Point2f point )
{
  bool isInsideSiblings = false;
  
  int sibling = id;
  
  while( sibling != -1 )
  {
    isInsideSiblings = ( cv::pointPolygonTest( mCvPolygon[ id ] , point , false  ) > 0 );
    
    if( isInsideSiblings )
    { 
      int child = mContourHierarchy[ id ][ 2 ];
      
      if( child != -1  )
      {
        isInsideSiblings = ( isInsideSiblings && !checkChildren( sibling , point ) );	
      }
      
      break;
    }
    
    sibling = mContourHierarchy[ sibling ][ 1 ];
  }
  
  return isInsideSiblings;
}


  
  


void SegmentClipper::addContours( std::vector< std::vector< cv::Point > > & contours)
{
  int numContours = contours.size();
  
  mPolygon.resize( numContours );

  mUsedIndices.resize( numContours );

  int offset = 0;
  
  xMin = std::numeric_limits< float >::max();
  xMax = std::numeric_limits< float >::min();
  
  yMin = std::numeric_limits< float >::max();
  yMax = std::numeric_limits< float >::min();
  
  for( int cc = 0; cc < numContours; cc++ )
  {
    int numPoints = contours[ cc ].size();
    
    mPolygon[ cc ].resize( numPoints );
    
    mUsedIndices[ cc ].resize( numPoints );
    
    for( int pp = 0; pp < numPoints; pp++ )
    {
      mPolygon[ cc ][ pp ].x() = contours[ cc ][ pp ].x;
      mPolygon[ cc ][ pp ].y() = contours[ cc ][ pp ].y;
      
      xMin = std::min( xMin , mPolygon[ cc ][ pp ].x() );
      xMax = std::max( xMax , mPolygon[ cc ][ pp ].x() );
      
      yMin = std::min( yMin , mPolygon[ cc ][ pp ].y() );
      yMax = std::max( yMax , mPolygon[ cc ][ pp ].y() );
      
      
      add( cc , pp , contours[ cc ][ pp ] ,  contours[ cc ][ ( pp + 1 ) % numPoints ] ); 
      
      mUsedIndices[ cc ][ pp ] = true;
    }    
  
  } 
  
  mCvPolygon.resize( numContours );
  
  for( int cc = 0; cc < numContours; cc++ )
  {
    int numPoints = contours[ cc ].size();
    
    mCvPolygon[ cc ].resize( numPoints );
    
    for( int pp = 0; pp < numPoints; pp++ )
    {
      mCvPolygon[ cc ][ pp ] = contours[ cc ][ pp ];    
    }    
    
    offset += numPoints;
    
  } 

  
}


void SegmentClipper::addContours( std::vector< std::vector< cv::Point2f > >& contours )
{
  int numContours = contours.size();
  
  mPolygon.resize( numContours );
  mCvPolygon.resize( numContours );
  
  mUsedIndices.resize( numContours );

    xMin = std::numeric_limits< float >::max();
  xMax = std::numeric_limits< float >::min();
  
  yMin = std::numeric_limits< float >::max();
  yMax = std::numeric_limits< float >::min();
  
  for( int cc = 0; cc < numContours; cc++ )
  {
    int numPoints = contours[ cc ].size();
    
    mPolygon[ cc ].resize( numPoints );
    
    mUsedIndices[ cc ].resize( numPoints );

    for( int pp = 0; pp < numPoints; pp++ )
    {
      mPolygon[ cc ][ pp ].x() = contours[ cc ][ pp ].x;
      mPolygon[ cc ][ pp ].y() = contours[ cc ][ pp ].y;
      
      xMin = std::min( xMin , mPolygon[ cc ][ pp ].x() );
      xMax = std::max( xMax , mPolygon[ cc ][ pp ].x() );
      
      yMin = std::min( yMin , mPolygon[ cc ][ pp ].y() );
      yMax = std::max( yMax , mPolygon[ cc ][ pp ].y() );
      
      add( cc , pp , contours[ cc ][ pp ] , contours[ cc ][ ( pp + 1 ) % numPoints ] );
      
      mUsedIndices[ cc ][ pp ] = true;
      
    }    
   
  }
  
  for( int cc = 0; cc < numContours; cc++ )
  {
    int numPoints = contours[ cc ].size();
    
    mCvPolygon[ cc ].resize( numPoints );
    
    for( int pp = 0; pp < numPoints; pp++ )
    {
      mCvPolygon[ cc ][ pp ] = contours[ cc ][ pp ];    
    }    
    
  } 

  
}



void SegmentClipper::addContours( std::vector< std::vector< tr::Point2f > >& contours)
{
  mPolygon = contours;
 
  int numContours = contours.size();
   
  mUsedIndices.resize( numContours );
   
  xMin = std::numeric_limits< float >::max();
  xMax = std::numeric_limits< float >::min();
  
  yMin = std::numeric_limits< float >::max();
  yMax = std::numeric_limits< float >::min();
  
  for( int cc = 0; cc < numContours; cc++ )
  {
    int numPoints = contours[ cc ].size();
    
    mUsedIndices[ cc ].resize( numPoints );
    
    for( int pp = 0; pp < numPoints; pp++ )
    {      
      add( cc , pp , mPolygon[ cc ][ pp ] , mPolygon[ cc ][ ( pp + 1 ) % numPoints ] ); 
      
      xMin = std::min( xMin , mPolygon[ cc ][ pp ].x() );
      xMax = std::max( xMax , mPolygon[ cc ][ pp ].x() );
      
      yMin = std::min( yMin , mPolygon[ cc ][ pp ].y() );
      yMax = std::max( yMax , mPolygon[ cc ][ pp ].y() );
      
       mUsedIndices[ cc ][ pp ] = true;
    } 

  } 
  
  for( int cc = 0; cc < numContours; cc++ )
  {
    int numPoints = contours[ cc ].size();
    
    mCvPolygon[ cc ].resize( numPoints );
    
    for( int pp = 0; pp < numPoints; pp++ )
    {
      mCvPolygon[ cc ][ pp ].x = contours[ cc ][ pp ].x();
      mCvPolygon[ cc ][ pp ].y = contours[ cc ][ pp ].y();
    }    
 
  } 
  
}


void SegmentClipper::setContourHierarchy( std::vector< cv::Vec4i >& contourHierarchy )
{
  mContourHierarchy = contourHierarchy;

}





bool intersectionStripPredicate( const std::pair< float , std::pair< int , int >  >  &obj1 , const std::pair< float , std::pair< int , int >  >  &obj2 )
{
   return obj1.first > obj2.first;
}



bool SegmentClipper::clipSegment2(tr::Point2f end1, tr::Point2f end2, 
				  std::vector< std::pair< tr::Point2f, tr::Point2f > >& clippedSegments, 
				  std::vector< std::pair< int, int > >& stripIds)
{
  
  mSelectedSegments.clear();
  mColors.clear();
  
  while( mStackTop >= 0 )
  {
    int contourId = mUsedStack[ mStackTop ].first;
    int segmentId = mUsedStack[ mStackTop ].second;
    
    mUsedIndices[ contourId ][ segmentId ] = true;
	  
    mStackTop--;
  }
  
  int index = -1;
  
  int size = mHorizontalArray.size();
  
  float xMax , xMin , yMax , yMin;
  
  if( end1.x() > end2.x() )
  {
    xMax = end1.x();
    xMin = end2.x();
  }
  else
  {
    xMax = end2.x();
    xMin = end1.x();
  }
  
  if( end1.y() > end2.y() )
  {
    yMax = end1.y();
    yMin = end2.y();
  }
  else
  {
    yMax = end2.y();
    yMin = end1.y();
  }
  
  std::vector< std::pair< float , std::pair< int , int > >  > clippedPoints;
  
  if( search( end1 , end2 , index  ) )
  {
    
    for( int idx = index - 1; idx < size; idx++ )
    {
      if( idx < 0 )
	    idx = 0;
      
      if( mHorizontalArray[ idx ].second > xMax )
	      break;
      
      int numSegments = mHorizontalStack[ idx ].size();
      
      for( int ss = 0; ss < numSegments; ss++ )
      {
	    int contourId = mHorizontalStack[ idx ][ ss ].first;
	    int segmentId = mHorizontalStack[ idx ][ ss ].second;

	    if( mUsedIndices[ contourId ][ segmentId ] )
	    {
	      tr::Point2f &segEnd1 = mPolygon[ contourId ][ segmentId ];
	      tr::Point2f &segEnd2 = mPolygon[ contourId ][ ( segmentId + 1 ) % mPolygon[ contourId ].size() ];
	  
	      mUsedIndices[ contourId ][ segmentId ] = false;
	  
	      mStackTop++;
	  
	      mUsedStack[ mStackTop ].first = contourId;
	      mUsedStack[ mStackTop ].second = segmentId;
	  
	     float segYMin , segYMax;
	  
	     if( segEnd1.y() < segEnd2.y() )
	     {
	       segYMin = segEnd1.y();
	       segYMax = segEnd2.y();
	     }
	     else
	     {
	       segYMin = segEnd2.y();
	       segYMax = segEnd1.y();
	     }
	  
	  if( !( segYMax < yMin || segYMin > yMax  )  )
	  {
	    std::pair< int , int > segment;
	  
	    segment.first = contourId;
	    segment.second = segmentId;
	  
	    std::pair< float , float > intersection;
	  
	    //compute the intersection
	    tr::lineToLineIntersection2D( end1 , end2 , segEnd1 , segEnd2 , intersection );
	  
	    if( TR_VALID_INTERSECTION( intersection ) )
	    {
			//if( mDisplayInfo )
			//	std::cout<<" intersection : "<<intersection.first<<" "<<intersection.second<<std::endl;

	      tr::Point2f intersectPoint = end1 * intersection.first  + ( 1 - intersection.first ) * end2;
	    
	      std::pair< float , std::pair< int , int >  > intersectionPoint;
	    
	      intersectionPoint.first = intersection.first;
	    
	      intersectionPoint.second.first = contourId;
	      intersectionPoint.second.second = segmentId;
	    
	      clippedPoints.push_back( intersectionPoint );
	      
	      mSelectedSegments.push_back( mHorizontalStack[ idx ][ ss ] );
	      mColors.push_back( cv::Vec3f( 255 , 0 , 0 ) );
	     }
	  }
	  
	}
	
      }
      
    }
  }
  else
  {
    for( int idx = index - 1; idx < size; idx++ )
    {
      if( idx < 0 )
	idx = 0;
      
      if( mVerticalArray[ idx ].second > yMax )
	break;
      
      int numSegments = mVerticalStack[ idx ].size();
      
    for( int ss = 0; ss < numSegments; ss++ )
    {
	int contourId = mVerticalStack[ idx ][ ss ].first;
	int segmentId = mVerticalStack[ idx ][ ss ].second;
	

	
	if( mUsedIndices[ contourId ][ segmentId ] )
	{
	  
	  
	  tr::Point2f &segEnd1 = mPolygon[ contourId ][ segmentId ];
	  tr::Point2f &segEnd2 = mPolygon[ contourId ][ ( segmentId + 1 ) % mPolygon[ contourId ].size() ];
	  
	  mUsedIndices[ contourId ][ segmentId ] = false;
	  
	  mStackTop++;
	  
	  mUsedStack[ mStackTop ].first = contourId;
	  mUsedStack[ mStackTop ].second = segmentId;
	  
	  float segXMin , segXMax;
	  
	  if( segEnd1.x() < segEnd2.x() )
	  {
	    segXMin = segEnd1.x();
	    segXMax = segEnd2.x();
	  }
	  else
	  {
	    segXMin = segEnd2.x();
	    segXMax = segEnd1.x();
	  }

	  if( !( ( segXMax < xMin ) || ( segXMin > xMax )  )  )
	  { 
	    
	    std::pair< int , int > segment;
	  
	    segment.first = contourId;
	    segment.second = segmentId;
	  
	    std::pair< float , float > intersection;
	  
	    //compute the intersection
	    tr::lineToLineIntersection2D( end1 , end2 , segEnd1 , segEnd2 , intersection );
	  
	    if( TR_VALID_INTERSECTION( intersection ) )
	    {

			//if( mDisplayInfo )
			//	std::cout<<" intersection : "<<intersection.first<<" "<<intersection.second<<std::endl;

	      tr::Point2f intersectPoint = end1 * intersection.first  + ( 1 - intersection.first ) * end2;
	    
	      std::pair< float , std::pair< int , int >  > intersectionPoint;
	    
	      intersectionPoint.first = intersection.first;
	    
	      intersectionPoint.second.first = contourId;
	      intersectionPoint.second.second = segmentId;
	    
	      clippedPoints.push_back( intersectionPoint );
	    
	      std::pair< int , int > strip;
	    
	      strip.first = contourId;
	      strip.second = segmentId;
	      
	      mSelectedSegments.push_back( mVerticalStack[ idx ][ ss ] );
	      mColors.push_back( cv::Vec3f( 255 , 0 , 0 ) );
	     }
	  }
	  
	}
	
      }
      
    }
  }

   if( clippedPoints.size() > 0 )
  {
    std::sort( clippedPoints.begin() , clippedPoints.end() , intersectionStripPredicate  );
  }
  else
  {
    if( isInsideForeground( end1 ) )
    {
      stripIds.push_back( std::pair< int , int >( -1 , -1 ) );
      stripIds.push_back( std::pair< int , int >( -1 , -1 ) );
      
      std::pair< tr::Point2f , tr::Point2f > segment;
	
      segment.first = end1;
      segment.second = end2;      
      
      clippedSegments.push_back( segment );
      
    }
    
    //return;
    
  }
  
  int numIntersections = clippedPoints.size();

  tr::Point2f prevPoint = end1;
  
  for( int ii = 0; ii < numIntersections; ii++ )
  {
    float coeff = clippedPoints[ ii ].first;

    tr::Point2f currentPoint = end1 * coeff + ( 1 - coeff ) * end2;
    
    if( isInsideSegment( prevPoint , clippedPoints[ ii ].second ) )
    {

      if( ii == 0 )
      {
	   stripIds.push_back( std::pair< int , int >( -1 , -1 ) );
      }
      
      std::pair< tr::Point2f , tr::Point2f > segment;
	
      segment.first = prevPoint;
      segment.second = currentPoint;
	
      clippedSegments.push_back( segment );

      prevPoint = currentPoint;
      
      stripIds.push_back(  clippedPoints[ ii ].second );

	  //if( mDisplayInfo )
	  //{
		 // std::cout<< clippedPoints[ ii ].second.first << " " <<clippedPoints[ ii ].second.second <<" "<<coeff<<" "<<mPolygon[ 0 ].size() <<std::endl;
	  //}
      
    }
    else
    {  
      prevPoint = currentPoint;
      
      stripIds.push_back(  clippedPoints[ ii ].second );  

	  //if( mDisplayInfo )
	  //{
		 // std::cout<< clippedPoints[ ii ].second.first << " " <<clippedPoints[ ii ].second.second <<" "<<coeff << mPolygon[ 0 ].size() <<std::endl;
	  //}
      
      if( ii == numIntersections - 1  )
      {
	    stripIds.push_back( std::pair< int , int >( -1 , -1 ) );	
	    std::pair< tr::Point2f , tr::Point2f > segment;	
	    segment.first = currentPoint;
	    segment.second = end2;	
	    clippedSegments.push_back( segment );

      }     
      
    }
    
    
  }
  
return true;
  
}


bool SegmentClipper::clipRay(tr::Point2f end1, tr::Point2f dir, std::vector< std::pair< tr::Point2f, tr::Point2f > >& clippedSegments, std::vector< std::pair< int, int > >& stripIds)
{
  //cv::clipLine();
  
  tr::Point2f pt1 , pt2;
  
  float offset = 5.0;
  
  if( abs( dir.x() ) < abs( dir.y() ) )
  {
    pt1.x() = ( dir.x() / dir.y() ) * ( yMin - offset ) + ( end1.x() - ( dir.x() / dir.y() ) * end1.y() ); 
    pt1.y() = ( yMin - offset );
    
    pt2.x() = ( dir.x() / dir.y() ) * ( yMax + offset ) + ( end1.x() - ( dir.x() / dir.y() ) * end1.y() ); 
    pt2.y() = ( yMax + offset );
    
  }
  else
  {
    pt1.y() = ( dir.y() / dir.x() ) * ( xMin - offset ) + ( end1.y() - ( dir.y() / dir.x() ) * end1.x() ); 
    pt1.x() = ( xMin - offset );
    
    pt2.y() = ( dir.y() / dir.x() ) * ( xMax + offset ) + ( end1.y() - ( dir.y() / dir.x() ) * end1.x() ); 
    pt2.x() = ( xMax + offset );
  }
  
  if( (pt1 - end1).norm() > (pt2 - end1).norm() )
  {
    tr::Point2f temp = pt1;
    
    pt1 = pt2;
    
    pt2 = temp;
  }
  
  offset = 2.0;
  
  tr::Point2f rpt1( xMin - offset , yMin - offset ) , rpt2( xMax + offset , yMin - offset  ) , rpt3( xMax + offset , yMax + offset ) , rpt4( xMin - offset , yMax + offset );
  
  std::pair< float , float > i1 , i2 , i3 , i4;
  
  //compute the intersection
  tr::lineToLineIntersection2D( pt1 , pt2 , rpt1 , rpt2 , i1 );
  tr::lineToLineIntersection2D( pt1 , pt2 , rpt2 , rpt3 , i2 );
  tr::lineToLineIntersection2D( pt1 , pt2 , rpt3 , rpt4 , i3 );
  tr::lineToLineIntersection2D( pt1 , pt2 , rpt4 , rpt1 , i4 );
  
  std::vector< double > intersections;
  
  if( TR_VALID_INTERSECTION( i1 ) )
  {
    intersections.push_back( i1.first );
  }
  
  if( TR_VALID_INTERSECTION( i2 ) )
  {
    intersections.push_back( i2.first );
  }
  
  if( TR_VALID_INTERSECTION( i3 ) )
  {
    intersections.push_back( i3.first );
  }
  
  if( TR_VALID_INTERSECTION( i4 ) )
  {
    intersections.push_back( i4.first );
  }
  
  if( intersections.size() == 0 )
  {
    return false;
  }
  else
  {
    double d1 , d2;
    
    if( intersections[ 0 ] > intersections[ 1 ] )
    {
      d1 = intersections[ 0 ];
      d2 = intersections[ 1 ];
    }
    else
    {
      d1 = intersections[ 1 ];
      d2 = intersections[ 0 ];
    }
    
    tr::Point2f clipEnd1 = pt1 * d1 + pt2 * ( 1 - d1 );
    tr::Point2f clipEnd2 = pt1 * d2 + pt2 * ( 1 - d2 );
    
    std::vector< std::pair< tr::Point2f , tr::Point2f > > initialStrips;
    
    
    std::pair< tr::Point2f , tr::Point2f > strip;
    
    strip.first = clipEnd1;
    strip.second = clipEnd2;
    
    initialStrips.push_back( strip );
    
    if( mDisplayInfo )
    displaySegments( initialStrips );
    
    clipSegment2( clipEnd1 , clipEnd2 , clippedSegments , stripIds );
    
    if( mDisplayInfo )
    displaySegments( clippedSegments );
    
  }

  return true;
  
}



float SegmentClipper::firstClip(  tr::Point2f end1 , tr::Point2f end2 , 
				                  tr::Point2f& clippedEnd, std::pair< int, int >& strip  )
{
  
 while( mStackTop >= 0 )
  {
    int contourId = mUsedStack[ mStackTop ].first;
    int segmentId = mUsedStack[ mStackTop ].second;
    
    mUsedIndices[ contourId ][ segmentId ] = true;
	  
    mStackTop--;
  }
  
  strip.first = -1;
  strip.second = -1;
    
  clippedEnd = tr::Point2f( -1 , -1 );
  
  int index = -1;
  
  int size = mHorizontalArray.size();
  
  float xMax , xMin , yMax , yMin;
  
  if( end1.x() > end2.x() )
  {
    xMax = end1.x();
    xMin = end2.x();
  }
  else
  {
    xMax = end2.x();
    xMin = end1.x();
  }
  
  if( end1.y() > end2.y() )
  {
    yMax = end1.y();
    yMin = end2.y();
  }
  else
  {
    yMax = end2.y();
    yMin = end1.y();
  }
  
  std::vector< std::pair< float , std::pair< int , int > >  > clippedPoints;
  
  if( search( end1 , end2 , index  ) )
  {
    
    for( int idx = index - 1; idx < size; idx++ )
    {
      if( idx < 0 )
	idx = 0;
      
      if( mHorizontalArray[ idx ].second > xMax )
	break;
      
      int numSegments = mHorizontalStack[ idx ].size();
      
      for( int ss = 0; ss < numSegments; ss++ )
      {
	int contourId = mHorizontalStack[ idx ][ ss ].first;
	int segmentId = mHorizontalStack[ idx ][ ss ].second;
	
	mSelectedSegments.push_back( mHorizontalStack[ idx ][ ss ] );
	mColors.push_back( cv::Vec3f( 255 , 0 , 0 ) );
	
	if( mUsedIndices[ contourId ][ segmentId ] )
	{
	  tr::Point2f &segEnd1 = mPolygon[ contourId ][ segmentId ];
	  tr::Point2f &segEnd2 = mPolygon[ contourId ][ ( segmentId + 1 ) % mPolygon[ contourId ].size() ];
	  
	  mUsedIndices[ contourId ][ segmentId ] = false;
	  
	  mStackTop++;
	  
	  mUsedStack[ mStackTop ].first = contourId;
	  mUsedStack[ mStackTop ].second = segmentId;
	  
	  float segYMin , segYMax;
	  
	  if( segEnd1.y() < segEnd2.y() )
	  {
	    segYMin = segEnd1.y();
	    segYMax = segEnd2.y();
	  }
	  else
	  {
	    segYMin = segEnd2.y();
	    segYMax = segEnd1.y();
	  }
	  
	  if( !( segYMax < yMin || segYMin > yMax  )  )
	  {
	    std::pair< int , int > segment;
	  
	    segment.first = contourId;
	    segment.second = segmentId;
	  
	    std::pair< float , float > intersection;
	  
	    //compute the intersection
	    tr::lineToLineIntersection2D( end1 , end2 , segEnd1 , segEnd2 , intersection );
	  
	    if( TR_VALID_INTERSECTION( intersection ) )
	    {
	      tr::Point2f intersectPoint = end1 * intersection.first  + ( 1 - intersection.first ) * end2;
	    
	      std::pair< float , std::pair< int , int >  > intersectionPoint;
	    
	      intersectionPoint.first = intersection.first;
	    
	      intersectionPoint.second.first = contourId;
	      intersectionPoint.second.second = segmentId;
	    
	      clippedPoints.push_back( intersectionPoint );
	    
	      std::pair< int , int > strip;
	    
	      strip.first = contourId;
	      strip.second = segmentId;
	     }
	  }
	  
	}
	
      }
      
    }
  }
  else
  {
    for( int idx = index - 1; idx < size; idx++ )
    {
      if( idx < 0 )
	idx = 0;
      
      if( mVerticalArray[ idx ].second > yMax )
	break;
      
      int numSegments = mVerticalStack[ idx ].size();
      
      for( int ss = 0; ss < numSegments; ss++ )
      {
	int contourId = mVerticalStack[ idx ][ ss ].first;
	int segmentId = mVerticalStack[ idx ][ ss ].second;
	
	mSelectedSegments.push_back( mVerticalStack[ idx ][ ss ] );
	mColors.push_back( cv::Vec3f( 255 , 0 , 0 ) );
	
	if( mUsedIndices[ contourId ][ segmentId ] )
	{
	  
	  
	  tr::Point2f &segEnd1 = mPolygon[ contourId ][ segmentId ];
	  tr::Point2f &segEnd2 = mPolygon[ contourId ][ ( segmentId + 1 ) % mPolygon[ contourId ].size() ];
	  
	  mUsedIndices[ contourId ][ segmentId ] = false;
	  
	  mStackTop++;
	  
	  mUsedStack[ mStackTop ].first = contourId;
	  mUsedStack[ mStackTop ].second = segmentId;
	  
	  float segXMin , segXMax;
	  
	  if( segEnd1.x() < segEnd2.x() )
	  {
	    segXMin = segEnd1.x();
	    segXMax = segEnd2.x();
	  }
	  else
	  {
	    segXMin = segEnd2.x();
	    segXMax = segEnd1.x();
	  }

	  if( !( ( segXMax < xMin ) || ( segXMin > xMax )  )  )
	  { 
	    
	    std::pair< int , int > segment;
	  
	    segment.first = contourId;
	    segment.second = segmentId;
	  
	    std::pair< float , float > intersection;
	  
	    //compute the intersection
	    tr::lineToLineIntersection2D( end1 , end2 , segEnd1 , segEnd2 , intersection );
	  
	    if( TR_VALID_INTERSECTION( intersection ) )
	    {
	      tr::Point2f intersectPoint = end1 * intersection.first  + ( 1 - intersection.first ) * end2;
	    
	      std::pair< float , std::pair< int , int >  > intersectionPoint;
	    
	      intersectionPoint.first = intersection.first;
	    
	      intersectionPoint.second.first = contourId;
	      intersectionPoint.second.second = segmentId;
	    
	      clippedPoints.push_back( intersectionPoint );
	    
	      std::pair< int , int > strip;
	    
	      strip.first = contourId;
	      strip.second = segmentId;
	     }
	  }
	  
	}
	
      }
      
    }
  }
  
  tr::Point2f prevPoint = end1;
  
  //if( mDisplayInfo )
  //  std::cout<<" num clipped points : "<<clippedPoints.size()<<std::endl;

  if( clippedPoints.size() > 0 )
  {
    std::sort( clippedPoints.begin() , clippedPoints.end() , intersectionStripPredicate  );
    
    float coeff = -1;
    
    for( int ii = 0; ii < clippedPoints.size(); ii++ )
    {
      coeff = clippedPoints[ ii ].first;
      
      tr::Point2f currentPoint = end1 * coeff + ( 1 - coeff ) * end2;
      
      if( isInsideSegment( prevPoint , clippedPoints[ ii ].second ) )
      {
	    clippedEnd = currentPoint;
	    strip = clippedPoints[ ii ].second;

	    break;
      }
      else
      {
	   prevPoint = currentPoint;
      }
      
    }   
 
    return ( 1 - coeff );
  }
  else
  {
    return -1;
  }

}






void SegmentClipper::clear()
{
  mHorizontalArray.clear();
  mVerticalArray.clear();
}


