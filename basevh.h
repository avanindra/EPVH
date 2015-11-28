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

#ifndef BASEVH_H
#define BASEVH_H

//#include "calibration/calibration_info.h"
#include "queue"
#include "camerainfo.h"
#include "opencvincludes.h"
#include "list"

 

namespace tr{
  
class Generator;  

  
 struct Vertex
{
  Vertex *mLeft , *mRight , *mOpp ; 
  
  Generator *mLeftGen , *mRightGen , *mGen;
  
  Eigen::Vector3d  mCoords;
  
  Eigen::Vector3d  mLE;
  
  Eigen::Vector3d  mLeftCoords , mRightCoords;
  
  int mId;
  
  uchar mIsGenVertex;
  
  int mIgnoreIndex;
  
}; 
  
  
typedef std::vector< cv::Point2f > contour; 
  
struct Generator
{
  int generatorId;
  
  int camId;
  
  int contourId;
  
  int stripId;
  
  int leftViewEdgeId , rightViewEdgeId;
  
  Eigen::Vector3d  normal;
  
  Eigen::Vector3d  leftRay , rightRay;
  
  uchar normalComputed;
  
  std::vector< Vertex* > mVertices , mAllVertices , mTriplePointVertices ;
  
  std::vector< bool > mUsedVertex;
  
  Generator *mLeftGen , *mRightGen;  
};  


struct ViewingEdgeInfo
{
  Eigen::Vector3d  point1 , point2;
  
  int camId , camId1 , camId2;
  
  int contourId , contourId1 , contourId2;
  
  uchar orientation;
  
  uchar isOccluding;
  
  int stripId , stripId1 , stripId2;
  
  int generatorId1 , generatorId2;
  
  Generator *leftGenerator , *rightGenerator;
  
  Generator *generator1 , *generator2;

  int frontVertexId , backVertexId;
  
};

struct StripContourMap
{
  int contourStripId;
  int contourId;
};

typedef std::pair< int , double > EdgeDistancePair;

class EdgeDistancePredicate
{

public:

    bool operator ()( EdgeDistancePair& pair1 , EdgeDistancePair& pair2 )
    {
      return pair1.second > pair2.second;
    }
};

typedef std::list< EdgeDistancePair > EdgeList;


typedef std::priority_queue< EdgeDistancePair , std::vector< EdgeDistancePair > , EdgeDistancePredicate > PriorityEdgeQueue;


#define validIntersection(x) ( x.first > 0 && x.first < 1 && x.second > 0 && x.second < 1 )
  

class BaseVH
{
  
protected:
  
  CameraInfo *mCalibration;
  std::vector< int > mSilhouetteCameras;
  std::vector< int > mMostOrthogonalCamera;
  std::vector< uchar > mIsCameraUsed;
  std::vector< std::vector< StripContourMap > > mStripContourMap;  
  std::vector< std::vector< contour > > mObjectContours;
  std::vector< std::vector< cv::Vec4i >  > mContourHierarchies;
  std::vector< std::vector< bool > > mIsOccludingContour;
  std::vector< std::vector< std::vector< Generator* >  > > mGenerators;
  std::vector< cv::Mat > mBoundingBoxImages;
  std::vector< cv::Mat > mGeneratorImages;
  std::vector< cv::Point2f > mOffset;
  std::vector< double > mScale , mInvScale ;  
  
  void buildMostOrthogonalCameras();  
  void buildPrimitives();
  void buildPrimitives( int camId );
  void buildPrimitives(int camId, std::vector< Eigen::Vector2d >& contour);
  void buildPrimitivesFromImageBoundary( int camId );
  int getMostOrthogonalUnusedCamera( int camId );
  void setCameraToUsed( int camId );
  void buildGeneratorImages();
  
  void stripEdgeIntersection( int camId , int contourId , int stripId , cv::Point2f end1 , cv::Point2f end2 , 
                              cv::Point2f &intersectionPoint , std::pair< double , double > &coefficientPairs );  

  void stripEdgeIntersection( int camId, int contourId, int stripId, Eigen::Vector2d end1, Eigen::Vector2d end2,
	  Eigen::Vector2d& intersectionPoint, std::pair< double, double >& coefficientPairs);
  
  
  double distanceToTheStrip( Eigen::Vector2d epipole , Eigen::Vector2d secondPoint , int camId , int contourId , int stripId );
  
  
  inline bool isInsideToEdge( Eigen::Vector2d point1, Eigen::Vector2d point2, Eigen::Vector2d point );
  
  bool isInsideStrip( int camId, int contourId , int stripId , Eigen::Vector2d point );

  void displayContour3D( std::vector< std::vector< cv::Point2f > > contours );

  void filterContoursByEdgeAngle( std::vector< std::vector< cv::Point2f > > &contours ,double angleThreshold = 2 );
  
  
public:
    BaseVH();
   
    void loadCameraInfo( CameraInfo *cameraInfo );
    void setSilhouetteCameras( std::vector< int > &silhoutteCameras );
	void clearGenerators();
    
    virtual bool compute() = 0;
    
    
    ~BaseVH();

};

}

#endif // BASEVH_H
