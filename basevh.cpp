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

#include "basevh.h"
#include "math2d.h"


#include "vtkincludes.h"
#include "display3droutines.h"


namespace tr{
  
 

BaseVH::BaseVH()
{

}



void BaseVH::setSilhouetteCameras(std::vector< int >& silhouetteCameras)
{
  mSilhouetteCameras = silhouetteCameras;
}





void BaseVH::buildMostOrthogonalCameras()
{
    int numSilhouetteCams = mSilhouetteCameras.size();

    int numCams = mCalibration->get_cam_count();

    mMostOrthogonalCamera.clear();

    mMostOrthogonalCamera.resize( numCams , -1 );
    
    mIsCameraUsed.clear();
    
    mIsCameraUsed.resize( numCams , 0 );

    for( int cc = 0; cc < numCams; cc++ )
    {
        int camId1 = cc;

        Vector3d cameraRay1;

        double width1 = mCalibration->get_image_width( camId1 );
        double height1 = mCalibration->get_image_height( camId1 );

        Eigen::Vector2d center1( width1 / 2 , height1 / 2 );

        mCalibration->getRay( camId1 , center1 , cameraRay1 );

        double prevVal = 1;

        int mostOrthogonalCam = -1;

        for( int sc = 0; sc < numSilhouetteCams; sc++ )
        {
            int camId2 = mSilhouetteCameras[ sc ];

            double width2 = mCalibration->get_image_width( camId2 );
            double height2 = mCalibration->get_image_height( camId2 );

            Eigen::Vector2d center2( width2 / 2 , height2 / 2 );

            Vector3d cameraRay2;

            mCalibration->getRay( camId2 , center2 , cameraRay2 );

            double dot = cameraRay1.dot( cameraRay2 );


            if( abs( dot ) < prevVal )
            {
              prevVal = abs( dot );

              mostOrthogonalCam = camId2;
            }

        }

        mMostOrthogonalCamera[ camId1 ] = mostOrthogonalCam;
    }
}



/** 
 
 Need to add few changes :
 
 1) Remove bounding image , as it is useless.
 2) Remove offset for the same reason
 3) Add contour hierarchies and do not filter the existing contours.
 
 **/

void BaseVH::buildPrimitives()
{
    int numsilhouettes = mSilhouetteCameras.size();

    int numCams = mCalibration->get_cam_count();

    mScale.resize( numCams , 1.0f );
    mInvScale.resize( numCams , 1.0f );
    mBoundingBoxImages.resize( numCams );
    mGeneratorImages.resize( numCams );
    mOffset.resize( numCams );

    mIsOccludingContour.resize( numCams );
    mObjectContours.resize( numCams );
    mGenerators.resize( numCams );
    mStripContourMap.resize( numCams );
    mContourHierarchies.resize( numCams );


    buildMostOrthogonalCameras();

    for( int sc = 0; sc < numsilhouettes; sc++ )
    {
		int sCamId = mSilhouetteCameras[ sc ];

		buildPrimitives( sCamId , mCalibration->mContours[ sCamId ]  );
    }

}


void BaseVH::buildPrimitives( int camId )
{

	int numsilhouettes = mSilhouetteCameras.size();

    int numCams = mCalibration->get_cam_count();



	
	mScale[ camId ] = mCalibration->silhouette( camId ).get_scale();

    mInvScale[ camId ] = 1.0f / mScale[ camId ] ;
	
	cv::Rect boundingBox;

    mCalibration->silhouette( camId ).get_bounding_box( boundingBox );	
		
	cv::Point offset( boundingBox.x , boundingBox.y );

    mOffset[ camId ] = cv::Point2f( 0 , 0 );

    std::vector< int > neighborCams( numsilhouettes - 1 );

    int id = 0;

    for( int sc2 = 0; sc2 < numsilhouettes; sc2++ )
    {
		if( camId == mSilhouetteCameras[ sc2 ] )
          continue;

        neighborCams[ id ] = mSilhouetteCameras[ sc2 ];

        id++;
     }

     mCalibration->silhouette( camId ).getBoundingBoxImage( mBoundingBoxImages[ camId ] );

     std::vector< std::vector< cv::Point2f > > allContours;
     std::vector< cv::Vec4i > contourHierarchy;
  
     mCalibration->silhouette( camId ).getHierarchyContours( allContours , contourHierarchy );
	
     int numAllContours = allContours.size();

      filterContoursByEdgeAngle( allContours );
	
	mScale[ camId ] = 1.0;
	mInvScale[ camId ] = 1.0;
	
	for( int cc = 0; cc < numAllContours; cc++ )
	{
	  int numPts =  allContours[ cc ].size();
	  
	  if( numPts > 20 )
	  {	    
	    cv::approxPolyDP( allContours[ cc ] , allContours[ cc ] , 2.0 , true );
	  }
	}


	mContourHierarchies[ camId ] = contourHierarchy;
	mObjectContours[ camId ] =  allContours;

    int numContours = mObjectContours[ camId ].size();
    int numContourPoints = 0;

    int stripId = 0;
	int contourId = 0;

    for( int contour = 0; contour < numContours; contour++ )
    {
        int numStrips =  mObjectContours[ camId ][ contour ].size();
        std::vector< Generator* > generators;

        generators.resize( numStrips );

        for( int strip = 0; strip < numStrips; strip++ )
        {
	       Eigen::Vector3d  ray;
	  
	       mCalibration->getRay( camId , mObjectContours[ camId ][ contour ][ strip ] , ray );
	  
           generators[ strip ] = new Generator();
           generators[ strip ]->leftViewEdgeId = -1;
           generators[ strip ]->rightViewEdgeId = -1;

           generators[ strip ]->stripId = strip;
           generators[ strip ]->camId = camId;
           generators[ strip ]->contourId = contour;
           generators[ strip ]->normalComputed = false;	    

           stripId++;
         }
                
         for( int strip = 0; strip < numStrips; strip++ )
         {
	   generators[ strip ]->mRightGen = generators[ ( strip + 1 ) % numStrips ];
	   generators[ strip ]->mLeftGen = generators[ ( strip - 1 + numStrips ) % numStrips ];
	 }

           mGenerators[ camId ].push_back( generators );

           numContourPoints += numStrips;
                
           contourId++;

        }
        
        
    for( int contour = 0; contour < numContours; contour++ )
    {
        int numStrips =  mObjectContours[ camId ][ contour ].size();
	
        for( int strip = 0; strip < numStrips; strip++ )
        {
	       Eigen::Vector3d  ray;
	  
	       mCalibration->getRay( camId , mObjectContours[ camId ][ contour ][ strip ] , ray );
	  
	       mGenerators[ camId ][ contour ][ strip ]->leftRay = ray;
	   
           int idx = ( strip - 1 + numStrips ) % numStrips;

           mGenerators[ camId ][ contour ][ idx ]->rightRay = ray;
         }                

     }   

        

    stripId = 0;
	
	int numFilteredContours = mObjectContours[ camId ].size();

	for( int contour = 0; contour < numFilteredContours; contour++ )
    {
        stripId += mObjectContours[ camId ][ contour ].size();        
    }
        
    mStripContourMap[ camId ].resize( stripId );
	
	stripId = 0;
	
    for( int contour = 0; contour < numFilteredContours; contour++ )
    {
       int numStrips = mObjectContours[ camId ][ contour ].size();

       for( int strip = 0; strip < numStrips; strip++ )
       {
           mStripContourMap[ camId ][ stripId ].contourId = contour;
           mStripContourMap[ camId ][ stripId ].contourStripId = strip;

           stripId++;
       }

    }	 
}

void BaseVH::buildPrimitives( int camId , std::vector< Eigen::Vector2d >& contour )
{
	int numsilhouettes = mSilhouetteCameras.size();

	int numCams = mCalibration->get_cam_count();

	mScale[camId] =1;

	mInvScale[camId] = 1.0;

	//cv::Rect boundingBox;

	//mCalibration->silhouette(camId).get_bounding_box(boundingBox);

	mOffset[camId] = cv::Point2f(0, 0);

	std::vector< int > neighborCams(numsilhouettes - 1);

	int id = 0;

	for (int sc2 = 0; sc2 < numsilhouettes; sc2++)
	{
		if (camId == mSilhouetteCameras[sc2])
			continue;

		neighborCams[id] = mSilhouetteCameras[sc2];

		id++;
	}

	//mCalibration->silhouette(camId).getBoundingBoxImage(mBoundingBoxImages[camId]);

	std::vector< std::vector< cv::Point2f > > allContours;
	std::vector< cv::Vec4i > contourHierarchy;

	//mCalibration->silhouette(camId).getHierarchyContours(allContours, contourHierarchy);



	int npts = contour.size();

	allContours.resize(1);
	allContours[0].resize(npts);
	contourHierarchy.resize(1);

	for (int pp = 0; pp < npts; pp++)
	{
		allContours[0][pp].x = contour[pp](0);
		allContours[0][pp].y = contour[pp](1);
	}

	contourHierarchy[ 0 ][0] = -1;
	contourHierarchy[0][1] = -1;
	contourHierarchy[0][2] = -1;
	contourHierarchy[0][3] = -1;

	int numAllContours = allContours.size();

	//filterContoursByEdgeAngle(allContours);

	mScale[camId] = 1.0;
	mInvScale[camId] = 1.0;

	for (int cc = 0; cc < numAllContours; cc++)
	{
		int numPts = allContours[cc].size();

		if (numPts > 20)
		{
			cv::approxPolyDP(allContours[cc], allContours[cc], 1.0, true);
		}
	}

	mContourHierarchies[camId] = contourHierarchy;
	mObjectContours[camId] = allContours;

	int numContours = mObjectContours[camId].size();
	int numContourPoints = 0;

	int stripId = 0;
	int contourId = 0;

	for (int contour = 0; contour < numContours; contour++)
	{
		int numStrips = mObjectContours[camId][contour].size();
		std::vector< Generator* > generators;

		generators.resize(numStrips);

		for (int strip = 0; strip < numStrips; strip++)
		{
			Eigen::Vector3d  ray;

			mCalibration->getRay(camId, mObjectContours[camId][contour][strip], ray);

			generators[strip] = new Generator();
			generators[strip]->leftViewEdgeId = -1;
			generators[strip]->rightViewEdgeId = -1;

			generators[strip]->stripId = strip;
			generators[strip]->camId = camId;
			generators[strip]->contourId = contour;
			generators[strip]->normalComputed = false;

			stripId++;
		}

		for (int strip = 0; strip < numStrips; strip++)
		{
			generators[strip]->mRightGen = generators[(strip + 1) % numStrips];
			generators[strip]->mLeftGen = generators[(strip - 1 + numStrips) % numStrips];
		}

		mGenerators[camId].push_back(generators);

		numContourPoints += numStrips;

		contourId++;

	}


	for (int contour = 0; contour < numContours; contour++)
	{
		int numStrips = mObjectContours[camId][contour].size();

		for (int strip = 0; strip < numStrips; strip++)
		{
			Eigen::Vector3d  ray;

			mCalibration->getRay(camId, mObjectContours[camId][contour][strip], ray);

			mGenerators[camId][contour][strip]->leftRay = ray;

			int idx = (strip - 1 + numStrips) % numStrips;

			mGenerators[camId][contour][idx]->rightRay = ray;
		}

	}



	stripId = 0;

	int numFilteredContours = mObjectContours[camId].size();

	for (int contour = 0; contour < numFilteredContours; contour++)
	{
		stripId += mObjectContours[camId][contour].size();
	}

	mStripContourMap[camId].resize(stripId);

	stripId = 0;

	for (int contour = 0; contour < numFilteredContours; contour++)
	{
		int numStrips = mObjectContours[camId][contour].size();

		for (int strip = 0; strip < numStrips; strip++)
		{
			mStripContourMap[camId][stripId].contourId = contour;
			mStripContourMap[camId][stripId].contourStripId = strip;

			stripId++;
		}

	}

}


void BaseVH::buildPrimitivesFromImageBoundary( int camId )
{
 
	int numsilhouettes = mSilhouetteCameras.size();

    int numCams = mCalibration->get_cam_count();
	
    mOffset[ camId ] = cv::Point2f( 0 , 0 );

	 std::vector< std::vector< cv::Point2f > > allContours;
     std::vector< cv::Vec4i > contourHierarchy;

	 std::vector< cv::Point2f > contour( 4 );

	 int w , h;

	 w = mCalibration->get_image_width( camId  );
	 h = mCalibration->get_image_height( camId  );

	 contour[ 0 ].x = 0;
	 contour[ 0 ].y = 0;
	
	 contour[ 1 ].x = 0;
	 contour[ 1 ].y = h;

	 contour[ 2 ].x = w;
	 contour[ 2 ].y = h;

	 contour[ 3 ].x = w;
	 contour[ 3 ].y = 0;

	 cv::Vec4i hr;

	 hr[ 0 ] = -1;
	 hr[ 1 ] = -1;
	 hr[ 2 ] = -1;
	 hr[ 3 ] = -1;

	allContours.push_back( contour );

	contourHierarchy.push_back( hr );

	mScale[ camId ] = 1.0;
	mInvScale[ camId ] = 1.0;
	
	
	mContourHierarchies[ camId ] = contourHierarchy;
	mObjectContours[ camId ] =  allContours;

    int numContours = mObjectContours[ camId ].size();
    int numContourPoints = 0;

    int stripId = 0;
	int contourId = 0;

    for( int contour = 0; contour < numContours; contour++ )
    {
        int numStrips =  mObjectContours[ camId ][ contour ].size();
        std::vector< Generator* > generators;

        generators.resize( numStrips );

        for( int strip = 0; strip < numStrips; strip++ )
        {
                    generators[ strip ] = new Generator();
                    generators[ strip ]->leftViewEdgeId = -1;
                    generators[ strip ]->rightViewEdgeId = -1;

                    generators[ strip ]->stripId = strip;
                    generators[ strip ]->camId = camId;
                    generators[ strip ]->contourId = contour;
                    generators[ strip ]->normalComputed = false;	    

                    stripId++;
                }
                
                for( int strip = 0; strip < numStrips; strip++ )
                {
		            generators[ strip ]->mRightGen = generators[ ( strip + 1 ) % numStrips ];
		            generators[ strip ]->mLeftGen = generators[ ( strip - 1 + numStrips ) % numStrips ];
		        }

                mGenerators[ camId ].push_back( generators );

                numContourPoints += numStrips;
                
                contourId++;
        }


	for( int contour = 0; contour < numContours; contour++ )
    {
        int numStrips =  mObjectContours[ camId ][ contour ].size();
	
        for( int strip = 0; strip < numStrips; strip++ )
        {
	       Eigen::Vector3d  ray;
	  
	       mCalibration->getRay( camId , mObjectContours[ camId ][ contour ][ strip ] , ray );
	  
	       mGenerators[ camId ][ contour ][ strip ]->leftRay = ray;
	   
           int idx = ( strip - 1 + numStrips ) % numStrips;

           mGenerators[ camId ][ contour ][ idx ]->rightRay = ray;
         }                

     } 

        

    stripId = 0;
	
	int numFilteredContours = mObjectContours[ camId ].size();

	for( int contour = 0; contour < numFilteredContours; contour++ )
    {
        stripId += mObjectContours[ camId ][ contour ].size();        
    }
        
    mStripContourMap[ camId ].resize( stripId );
	
	stripId = 0;
	
    for( int contour = 0; contour < numFilteredContours; contour++ )
    {
       int numStrips = mObjectContours[ camId ][ contour ].size();

       for( int strip = 0; strip < numStrips; strip++ )
       {
           mStripContourMap[ camId ][ stripId ].contourId = contour;
           mStripContourMap[ camId ][ stripId ].contourStripId = strip;

           stripId++;
       }

    }	  
}

int BaseVH::getMostOrthogonalUnusedCamera(int camId)
{  
    Vector3d cameraRay1;

    double width1 = mCalibration->get_image_width( camId );
    double height1 = mCalibration->get_image_height( camId );

    cv::Point2f center1( width1 / 2 , height1 / 2 );

    mCalibration->getRay( camId , center1 , cameraRay1 );

    double prevVal = 1;

    int mostOrthogonalCam = -1;
    
    int numSilhouetteCams = mSilhouetteCameras.size();

    for( int sc = 0; sc < numSilhouetteCams; sc++ )
    {
        int camId2 = mSilhouetteCameras[ sc ];
	
	if( camId2 == camId )
	  continue;

        double width2 = mCalibration->get_image_width( camId2 );
        double height2 = mCalibration->get_image_height( camId2 );

        cv::Point2f center2( width2 / 2 , height2 / 2 );

        Vector3d cameraRay2;

        mCalibration->getRay( camId2 , center2 , cameraRay2 );

        double dot = cameraRay1.dot( cameraRay2 );


        if( abs( dot ) < prevVal && !mIsCameraUsed[ camId2 ] )
        {
            prevVal = dot;

            mostOrthogonalCam = camId2;
        }

    }
    
    return mostOrthogonalCam;


}

void BaseVH::setCameraToUsed(int camId)
{
  mIsCameraUsed[ camId ] = 1;
}


void BaseVH::stripEdgeIntersection( int camId, int contourId, int stripId, cv::Point2f end1, cv::Point2f end2,
                                               cv::Point2f& intersectionPoint, std::pair< double, double >& coefficientPairs )
{   
    int numStrips = mObjectContours[ camId ][ contourId ].size();

    tr::lineToLineIntersection2D( mObjectContours[ camId ][ contourId ][ stripId ] ,
                                  mObjectContours[ camId ][ contourId ][ ( stripId + 1 ) % numStrips ] ,
                                  end1 , end2 , intersectionPoint , coefficientPairs );

}

void BaseVH::stripEdgeIntersection( int camId, int contourId, int stripId, Eigen::Vector2d end1, Eigen::Vector2d end2,
                                               Eigen::Vector2d& intersectionPoint, std::pair< double, double >& coefficientPairs )
{
    int numStrips = mObjectContours[ camId ][ contourId ].size();

    Eigen::Vector2d point1 , point2;

    point1.x() = ( mObjectContours[ camId ][ contourId ][ stripId ].x + mOffset[ camId ].x ) * mInvScale[ camId ];
    point1.y() = ( mObjectContours[ camId ][ contourId ][ stripId ].y + mOffset[ camId ].y ) * mInvScale[ camId ];

    int stripId2 = ( stripId + 1 ) % numStrips;

    point2.x() = ( mObjectContours[ camId ][ contourId ][ stripId2 ].x + mOffset[ camId ].x ) * mInvScale[ camId ];
    point2.y() = ( mObjectContours[ camId ][ contourId ][ stripId2 ].y + mOffset[ camId ].y ) * mInvScale[ camId ];


    tr::lineToLineIntersection2D( point1 , point2 , end1 , end2 , intersectionPoint , coefficientPairs );

}


double BaseVH::distanceToTheStrip( Eigen::Vector2d epipole , Eigen::Vector2d secondPoint , int camId , int contourId , int stripId )
{
  Eigen::Vector2d end1 , end2;
  
  end1.x() = ( mObjectContours[ camId ][ contourId ][ stripId ].x + mOffset[ camId ].x ) * mInvScale[ camId ];
  end1.y() = ( mObjectContours[ camId ][ contourId ][ stripId ].y + mOffset[ camId ].y ) * mInvScale[ camId ];

  int numStrips = mObjectContours[ camId ][ contourId ].size();
  
  end2.x() = ( mObjectContours[ camId ][ contourId ][ ( stripId + 1 ) % numStrips ].x + mOffset[ camId ].x ) * mInvScale[ camId ];
  end2.y() = ( mObjectContours[ camId ][ contourId ][ ( stripId + 1 ) % numStrips ].y + mOffset[ camId ].y ) * mInvScale[ camId ];
  
  Eigen::Vector2d intersectionPoint;
  
  tr::lineToLineIntersection2D( epipole , secondPoint , end1 , end2 , intersectionPoint );
  
  double dist = ( epipole - intersectionPoint ).norm();
  
  return dist;
}


bool BaseVH::isInsideToEdge( Eigen::Vector2d point1 , Eigen::Vector2d point2 , Eigen::Vector2d point )
{
  double val = ( point.y() - point1.y() ) * ( point2.x() - point1.x() ) - ( point2.y() - point1.y() ) * ( point.x() - point1.x() );

  return ( val < 0 );
}



bool BaseVH::isInsideStrip(  int camId, int contourId , int stripId, Eigen::Vector2d point )
{
    
	int numStrips = mObjectContours[ camId ][ contourId ].size();
  
    Eigen::Vector2d point1 , point2;

    point1.x() = ( mObjectContours[ camId ][ contourId ][ stripId ].x + mOffset[ camId ].x ) * mInvScale[ camId ];
    point1.y() = ( mObjectContours[ camId ][ contourId ][ stripId ].y + mOffset[ camId ].y ) * mInvScale[ camId ];

    int stripId2 = ( stripId + 1 ) % numStrips;

    point2.x() = ( mObjectContours[ camId ][ contourId ][ stripId2 ].x + mOffset[ camId ].x ) * mInvScale[ camId ];
    point2.y() = ( mObjectContours[ camId ][ contourId ][ stripId2 ].y + mOffset[ camId ].y ) * mInvScale[ camId ];
  
  
   return isInsideToEdge(  point1 , point2 , point  );

}


void BaseVH::displayContour3D( std::vector< std::vector< cv::Point2f > > contours )
{

  vtkSmartPointer< vtkPolyData > edgePolyData =  vtkSmartPointer< vtkPolyData >::New();
  vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();    
    
  colors->SetName("colors");
  colors->SetNumberOfComponents(3);
  
  int numEdges = 0 ;
  
  for( int ee1 = 0; ee1 < contours.size(); ee1++ )
  {
    numEdges += contours[ ee1 ].size();
  }
  
  vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();
  
  int numPoints = 2 * numEdges;
  
  points->Allocate( numPoints );
    
  int numCells = numEdges ;
  
  edgePolyData->Allocate( numCells );
  
  vtkSmartPointer< vtkIdList > cell = vtkSmartPointer< vtkIdList >::New();
  
  std::vector< std::vector< cv::Vec3f > > colorVecs;
  
  colorVecs.resize( contours.size() );

  bool alternate = false;
  
  for( int ee1 = 0; ee1 < contours.size(); ee1++ )
    for( int ee2 = 0; ee2 < contours[ ee1 ].size() ; ee2++ ) 
  {
	  if( ee2 == 150  )
	  {
		  colorVecs[ ee1 ].push_back( cv::Vec3f( 0 , 255 , 255 ) );
		  continue;
	  }
// 
// 	  if( ee2 == 60 || ee2 == 61  )
// 	  {
// 		  colorVecs[ ee1 ].push_back( cv::Vec3f( 255 , 255 , 255 ) );
// 		  continue;
// 	  }

	  if( alternate )
	  {
        colorVecs[ ee1 ].push_back( cv::Vec3f( 0 , 0 , 255 ) );

		alternate = false;
	  }
	  else
	  {
		colorVecs[ ee1 ].push_back( cv::Vec3f( 255 , 0 , 0 ) );

		alternate = true;
	  }
  }
  
  
  for( int ee1 = 0; ee1 < contours.size(); ee1++ )
    for( int ee2 = 0; ee2 < contours[ ee1 ].size() ; ee2++ ) 
  {
    
    cv::Point2f &pt1f = contours[ ee1 ][ ee2 ];
    cv::Point2f &pt2f = contours[ ee1 ][ ( ee2 + 1 ) % contours[ ee1 ].size() ];
    
    Eigen::Vector3d  pt1( pt1f.x , pt1f.y , 0 ) , pt2( pt2f.x , pt2f.y , 0 );

    //if( ee2 == 150 )
    // std::cout<< pt1.transpose() <<" -- "<< pt2.transpose() << std::endl;
    
    int id1 = points->InsertNextPoint( pt1.data() );
    int id2 = points->InsertNextPoint( pt2.data() );
    
    cell->Reset();
    
    cell->InsertNextId( id1 );
    cell->InsertNextId( id2 );
    
    edgePolyData->InsertNextCell( VTK_LINE , cell );
    
    const unsigned char _color[] = { colorVecs[ ee1 ][ ee2 ][ 0 ] , colorVecs[ ee1 ][ ee2 ][ 1 ] , colorVecs[ ee1 ][ ee2 ][ 2 ]};
  
    colors->InsertNextTupleValue(_color);
  }
  
  
  edgePolyData->SetPoints( points );
  
  edgePolyData->GetCellData()->SetScalars( colors );

  tr::Display3DRoutines::displayPolyData( edgePolyData );

}



void BaseVH::filterContoursByEdgeAngle( std::vector< std::vector< cv::Point2f > > &contours , double angleThreshold )
{
  int numContours = contours.size();

  double thresholdVal = - std::cos( ( angleThreshold / 180 ) * CV_PI );

  std::vector< std::vector< cv::Point2f > > filteredContours;

  double averageContourLength = 0;

  int totalNumEdges = 0;

  for( int cc = 0; cc < numContours; cc++ )
  {
	  int numEdges = contours[ cc ].size();

	  std::vector< cv::Point2f > filteredContour;

	  int id1 = -1 , id2 = -1 , id3 = -1;

	  for( int ee = 0; ee < numEdges; ee++ )
	  {
		  cv::Point2f &pt1 = contours[ cc ][ ( ee - 1 + numEdges ) % numEdges ];
		  cv::Point2f &pt2 = contours[ cc ][ ee ];

		  tr::Vector2f vec1( pt2.x - pt1.x , pt2.y - pt1.y );

		  averageContourLength += vec1.norm();

	  }

	  totalNumEdges += numEdges;
  }

  if( totalNumEdges < 4 )
	  return;

  averageContourLength /= totalNumEdges;

  double lengthThreshold = 1e-4 * averageContourLength;

  for( int cc = 0; cc < numContours; cc++ )
  {
	  int numEdges = contours[ cc ].size();

	  std::vector< cv::Point2f > filteredContour;

	  int id1 = -1 , id2 = -1 , id3 = -1;

	  cv::Point2f prevPoint( 0 , 0 );

	  for( int ee = 0; ee < numEdges; ee++ )
	  {
		  cv::Point2f &pt1 = contours[ cc ][ ( ee - 1 + numEdges ) % numEdges ];
		  cv::Point2f &pt2 = contours[ cc ][ ee ];
		  cv::Point2f &pt3 = contours[ cc ][ ( ee + 1 ) % numEdges ];

		  tr::Vector2f vec1( pt2.x - pt1.x , pt2.y - pt1.y ) , vec2( pt3.x - pt2.x , pt3.y - pt2.y );

		  vec1.normalize();
		  vec2.normalize();  

		  double val = vec1.dot( vec2 );

		  if( val < thresholdVal )
		  {
		    continue;
		  }

		  if( filteredContour.size() > 0 )
		  {
			  tr::Vector2f vec( prevPoint.x - pt2.x , prevPoint.y - pt2.y );

			  if( vec.norm() < lengthThreshold )
			  {
			    continue;
			  }
		  }

		  prevPoint = pt2;

		  filteredContour.push_back( pt2 );
	  }

	  filteredContours.push_back( filteredContour );
  }

  contours = filteredContours;

}

void BaseVH::loadCameraInfo(CameraInfo *cameraInfo)
{
	mCalibration = cameraInfo;
}

void BaseVH::clearGenerators()
{
  int size1 = mGenerators.size();

  for( int ss = 0; ss < size1 ; ss++ )
  {
      int size2 = mGenerators[ ss ].size();

      for( int ss2 = 0; ss2 < size2; ss2++ )
	  {
	     int size3 = mGenerators[ ss ][ ss2 ].size();

		 for( int ss3 = 0; ss3 < size3 ; ss3++ )
		 {
		   if( mGenerators[ ss ][ ss2 ][ ss3 ] )
		   {
		     delete mGenerators[ ss ][ ss2 ][ ss3 ];

			 mGenerators[ ss ][ ss2 ][ ss3 ] = 0;
		   }
		 }
	  }
  }

  mGenerators.clear();
}


BaseVH::~BaseVH()
{
	clearGenerators();
}

}
