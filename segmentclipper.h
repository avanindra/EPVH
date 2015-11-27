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


#ifndef SEGMENTCLIPPER_H
#define SEGMENTCLIPPER_H

#include <vector>
#include <algorithm>
#include <iostream>
#include "opencvincludes.h"
#include "types/generaltypes.h"





struct Segment
{
  int mContourId , mSegmentId;

  uchar mHorizontalPose , mVerticalPose;
};


typedef std::pair< Segment , float > IndexValPair;


class SegmentClipper
{  
  std::vector< IndexValPair > mHorizontalArray , mVerticalArray;  
  //std::vector< IndexValPair > mHorizontalArray2 , mVerticalArray2;
  
  std::vector< std::vector< std::pair< int , int > > > mHorizontalStack , mVerticalStack;
  
  std::vector< std::vector< bool > > mUsedIndices;
  
  std::vector< std::pair< int , int > > mUsedStack;
  
  float xMin , xMax , yMin , yMax;
  
  int mStackTop;
  
  std::vector< std::vector< tr::Point2f > > mPolygon;
  
  std::vector< std::vector< cv::Point2f > > mCvPolygon;
  
  std::vector< cv::Vec4i > mContourHierarchy;
  
 // std::vector< int > mIndexOffset;
  
  
  
public:
  
  std::vector< std::pair< int , int > > mSelectedSegments;
  std::vector< cv::Vec3f > mColors;

  bool mDisplayInfo;
  
  SegmentClipper();
  
  void build();
  
  void addContours( std::vector< std::vector< cv::Point > >  &contours );
  
  void addContours( std::vector< std::vector< cv::Point2f > >  &contours );
  
  void addContours(  std::vector< std::vector< tr::Point2f > >  &contours );
  
  void setContourHierarchy( std::vector< cv::Vec4i > &contourHierarchy );
  
  void add( int contourId , int segmentId , tr::Point2f end1 , tr::Point2f end2 );
  
  void add( int contourId , int segmentId , cv::Point2f end1 , cv::Point2f end2 );
  
  bool search( tr::Point2f end1 , tr::Point2f end2 , int &index );
  
  
  bool clipSegment2( tr::Point2f end1 , tr::Point2f end2 , 
		    std::vector< std::pair< tr::Point2f , tr::Point2f > > &clippedSegments , 
		    std::vector< std::pair< int , int > > &stripIds  );
  
  bool clipRay( tr::Point2f end1 , tr::Point2f dir , std::vector< std::pair< tr::Point2f , tr::Point2f > > &clippedSegments , 
		        std::vector< std::pair< int , int > > &stripIds );
  
  float firstClip( tr::Point2f end1 , tr::Point2f end2 , 
		           tr::Point2f &clippedEnd , std::pair< int , int > &strip  );

  void displaySegments( bool horizontal = true , bool forward = true );
  
  void displaySegments( std::vector< std::pair< tr::Point2f , tr::Point2f >  > &strips );
  
  bool isInsideForeground( tr::Point2f point );
  
  void clear();
  
protected:

  int binarySearch( float val , const std::vector< IndexValPair > &searchArray );  
  
  bool isInsideSegment( tr::Point2f point , std::pair< int , int > strip );
  
  
  
  bool checkChildren( int id , cv::Point2f point );
  
  bool checkSiblings( int id , cv::Point2f point );
  
  
  
};


#endif //SEGMENTCLIPPER
